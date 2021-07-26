!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!>  @brief    Module for field writing routines
!>  @details  Holds all routines for writing LFRic fields
!>
module lfric_xios_write_mod

  use clock_mod,            only: clock_type
  use constants_mod,        only: i_def, dp_xios, str_def, &
                                  str_max_filename, xios_max_int
  use field_mod,            only: field_type, field_proxy_type
  use field_parent_mod,     only: field_parent_proxy_type
  use field_collection_mod, only: field_collection_type, &
                                  field_collection_iterator_type
  use field_parent_mod,     only: field_parent_type
  use fs_continuity_mod,    only: W3
  use io_mod,               only: ts_fname
  use integer_field_mod,    only: integer_field_type, integer_field_proxy_type
  use r_solver_field_mod,   only: r_solver_field_type, r_solver_field_proxy_type
  use log_mod,              only: log_event,         &
                                  log_scratch_space, &
                                  LOG_LEVEL_INFO,    &
                                  LOG_LEVEL_WARNING, &
                                  LOG_LEVEL_ERROR
#ifdef UNIT_TEST
  use lfric_xios_mock_mod,  only: xios_send_field,      &
                                  xios_get_domain_attr, &
                                  xios_get_axis_attr
#else
  use xios,                 only: xios_send_field,      &
                                  xios_get_domain_attr, &
                                  xios_get_axis_attr
#endif

  implicit none

  private
  public :: checkpoint_write_xios,    &
            write_field_node,         &
            write_field_single_face,  &
            write_field_face,         &
            write_field_edge,         &
            write_state,              &
            write_checkpoint

contains

!>  @brief    I/O handler for writing an XIOS netcdf checkpoint
!>  @details  Note this routine accepts a filename but doesn't use it - this is
!>            to keep the interface the same for all methods
!>
!>  @param[in]      xios_field_name  XIOS identifier for the field
!>  @param[in]      file_name        Name of the file to write into
!>  @param[in,out]  field_proxy      A field proxy to be written
!>
subroutine checkpoint_write_xios(xios_field_name, file_name, field_proxy)

  implicit none

  character(len=*),               intent(in) :: xios_field_name
  character(len=*),               intent(in) :: file_name
  class(field_parent_proxy_type), intent(in) :: field_proxy

  integer(i_def)             :: undf
  real(dp_xios), allocatable :: send_field(:)

  undf = field_proxy%vspace%get_last_dof_owned()
  allocate(send_field(undf))

  ! Different field kinds are selected to access data
  select type(field_proxy)

    type is (field_proxy_type)
    send_field = field_proxy%data(1:undf)

    type is (integer_field_proxy_type)
    if ( any( abs(field_proxy%data) > xios_max_int) ) then
      call log_event( 'Data for integer field "'// trim(adjustl(xios_field_name)) // &
                      '" contains values too large for 16-bit precision', LOG_LEVEL_WARNING )
    end if
    send_field = real( field_proxy%data(1:undf), dp_xios )

    type is (r_solver_field_proxy_type)
    send_field = field_proxy%data(1:undf)

    class default
    call log_event( "Invalid type for input field proxy", LOG_LEVEL_ERROR )

  end select

  call xios_send_field(xios_field_name, reshape (send_field, (/1, undf/)))

end subroutine checkpoint_write_xios

!>  @brief  Output W0 field data to the node UGRID via XIOS
!>
!>  @param[in]     xios_field_name  XIOS identifier for the field
!>  @param[inout]  field_proxy      A field proxy to be written
!>
subroutine write_field_node(xios_field_name, field_proxy)

  implicit none

  character(len=*),               intent(in) :: xios_field_name
  class(field_parent_proxy_type), intent(in) :: field_proxy

  integer(i_def) :: i, undf
  integer(i_def) :: domain_size, axis_size
  real(dp_xios), allocatable :: send_field(:)

  undf = field_proxy%vspace%get_last_dof_owned()

  ! Get the expected horizontal domain size for the rank
  call xios_get_domain_attr('node', ni=domain_size)
  ! Get the expected vertical axis size
  call xios_get_axis_attr("vert_axis_full_levels", n_glo=axis_size)

  ! Size the arrays to be what is expected
  allocate(send_field(domain_size*axis_size))

  ! All data are scalar fields

  ! We need to reshape the raw field data to get the correct data layout for UGRID
  ! At the moment field array data is 1D with levels ordered sequentially
  ! This is only true for current scalar fields on lowest order fs and may change

  ! Different field kinds are selected to access data, which is arranged into blocks
  ! based on model level
  select type(field_proxy)

    type is (field_proxy_type)
    do i = 0, axis_size-1
      send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
                 field_proxy%data(i+1:undf:axis_size)
    end do

    type is (integer_field_proxy_type)
    if ( any( abs(field_proxy%data) > xios_max_int) ) then
      call log_event( 'Data for integer field "'// trim(adjustl(xios_field_name)) // &
                      '" contains values too large for 16-bit precision', LOG_LEVEL_WARNING )
    end if
    do i = 0, axis_size-1
      send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
                 real( field_proxy%data(i+1:undf:axis_size), dp_xios )
    end do

    type is (r_solver_field_proxy_type)
    do i = 0, axis_size-1
      send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
                 field_proxy%data(i+1:undf:axis_size)
    end do

    class default
      call log_event( "Invalid type for input field proxy", LOG_LEVEL_ERROR )

  end select

  ! Reshape into 2D horizontal + vertical levels for output
  call xios_send_field(xios_field_name, &
                      reshape (send_field, (/domain_size, axis_size/) ))

  deallocate(send_field)

end subroutine write_field_node

!>  @brief  Write W2H field data to the half-level edge UGRID via XIOS
!>
!>  @param[in]     xios_field_name  XIOS identifier for the field
!>  @param[inout]  field_proxy      A field proxy to be written
!>
subroutine write_field_edge(xios_field_name, field_proxy)

  implicit none

  character(len=*),               intent(in) :: xios_field_name
  class(field_parent_proxy_type), intent(in) :: field_proxy

  integer(i_def) :: i, undf
  integer(i_def) :: domain_size, axis_size
  real(dp_xios), allocatable :: send_field(:)

  undf = field_proxy%vspace%get_last_dof_owned()

  ! Get the expected horizontal and vertical axis size
  call xios_get_domain_attr('edge', ni=domain_size)
  call xios_get_axis_attr("vert_axis_half_levels", n_glo=axis_size)

  ! Size the arrays to be what is expected
  allocate(send_field(domain_size*axis_size))

  ! All data are scalar fields

  ! We need to reshape the raw field data to get the correct data layout for UGRID
  ! At the moment field array data is 1D with levels ordered sequentially
  ! This is only true for current scalar fields on lowest order fs and may change

  ! Different field kinds are selected to access data, which is arranged into blocks
  ! based on model level

  select type(field_proxy)

    type is (field_proxy_type)
    do i = 0, axis_size-1
      send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
                 field_proxy%data(i+1:undf:axis_size)
    end do

    type is (integer_field_proxy_type)
    if ( any( abs(field_proxy%data) > xios_max_int) ) then
      call log_event( 'Data for integer field "'// trim(adjustl(xios_field_name)) // &
                      '" contains values too large for 16-bit precision', LOG_LEVEL_WARNING )
    end if
    do i = 0, axis_size-1
      send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
                 real( field_proxy%data(i+1:undf:axis_size), dp_xios )
    end do

    type is (r_solver_field_proxy_type)
      do i = 0, axis_size-1
        send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
                   field_proxy%data(i+1:undf:axis_size)
      end do

    class default
      call log_event( "Invalid type for input field proxy", LOG_LEVEL_ERROR )

  end select

  ! Reshape into 2D horizontal + vertical levels for output
  call xios_send_field( xios_field_name, &
                        reshape( send_field, (/domain_size, axis_size/) ) )

  deallocate(send_field)

end subroutine write_field_edge

!>  @brief  Write single-level W3 field data to the face UGRID domain via XIOS
!>
!>  @param[in]     xios_field_name  XIOS identifier for the field
!>  @param[inout]  field_proxy      A field proxy to be written
!>
subroutine write_field_single_face(xios_field_name, field_proxy)

  implicit none

  character(len=*),               intent(in) :: xios_field_name
  class(field_parent_proxy_type), intent(in) :: field_proxy

  integer(i_def) :: i, undf, ndata
  integer(i_def) :: domain_size
  real(dp_xios), allocatable :: send_field(:)

  undf = field_proxy%vspace%get_last_dof_owned()
  ndata = field_proxy%vspace%get_ndata()

  ! Get the expected horizontal size
  ! all 2D fields are nominally in W3, hence half levels
  call xios_get_domain_attr('face', ni=domain_size)

  ! Size the array to be what is expected
  allocate(send_field(domain_size*ndata))

  ! If the fields do not have the same size, then exit with error
  if ( size(send_field) /= undf ) then
    call log_event( "Global size of model field /= size of field from file", &
                    LOG_LEVEL_ERROR )
  end if

  ! Different field kinds are selected to access data - data is re-ordered into
  ! slabs according to ndata
  select type(field_proxy)

    type is (field_proxy_type)
    do i = 0, ndata-1
      send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
                              field_proxy%data(i+1:(ndata*domain_size)+i:ndata)
    end do

    type is (integer_field_proxy_type)
    if ( any( abs(field_proxy%data) > xios_max_int) ) then
      call log_event( 'Data for integer field "'// trim(adjustl(xios_field_name)) // &
                      '" contains values too large for 16-bit precision', LOG_LEVEL_WARNING )
    end if
    do i = 0, ndata-1
      send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
                              field_proxy%data(i+1:(ndata*domain_size)+i:ndata)
    end do

    type is (r_solver_field_proxy_type)
    do i = 0, ndata-1
      send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
                              field_proxy%data(i+1:(ndata*domain_size)+i:ndata)
    end do

    class default
      call log_event( "Invalid type for input field proxy", LOG_LEVEL_ERROR )

  end select

  call xios_send_field(xios_field_name, reshape(send_field, (/1, domain_size*ndata/)))

  deallocate(send_field)

end subroutine write_field_single_face

!>  @brief  Write W3/WTheta field data to the full/half-level face UGRIDs via
!>          XIOS
!>
!>  @param[in]     xios_field_name  XIOS identifier for the field
!>  @param[inout]  field_proxy      A field proxy to be written
!>
subroutine write_field_face(xios_field_name, field_proxy)

  implicit none

  character(len=*),               intent(in) :: xios_field_name
  class(field_parent_proxy_type), intent(in) :: field_proxy

  integer(i_def) :: i, undf
  integer(i_def) :: fs_id
  integer(i_def) :: domain_size, axis_size
  real(dp_xios), allocatable :: send_field(:)

  ! Field must be cast to kind to get function space ID
  select type(field_proxy)
    type is (field_proxy_type)
    fs_id = field_proxy%vspace%which()

    type is (integer_field_proxy_type)
    fs_id = field_proxy%vspace%which()

    type is (r_solver_field_proxy_type)
    fs_id = field_proxy%vspace%which()

    class default
    call log_event( "Invalid type for input field proxy", LOG_LEVEL_ERROR )

  end select

  ! Size the arrays to be what is expected
  undf = field_proxy%vspace%get_last_dof_owned()

  ! Get the expected horizontal and vertical axis size
  if ( fs_id == W3 ) then
    call xios_get_domain_attr('face', ni=domain_size)
    call xios_get_axis_attr("vert_axis_half_levels", n_glo=axis_size)
  else
    call xios_get_domain_attr('face', ni=domain_size)
    call xios_get_axis_attr("vert_axis_full_levels", n_glo=axis_size)
  end if

  allocate(send_field(domain_size*axis_size))

  ! We need to reshape the raw field data to get the correct data layout for UGRID
  ! At the moment field array data is 1D with levels ordered sequentially
  ! This is only true for current scalar fields on lowest order fs and may change

  ! Different field kinds are selected to access data, which is arranged into blocks
  ! based on model level
  select type(field_proxy)

    type is (field_proxy_type)
      do i = 0, axis_size-1
        send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
                   field_proxy%data(i+1:undf:axis_size)
      end do

    type is (integer_field_proxy_type)
      if ( any( abs(field_proxy%data) > xios_max_int) ) then
        call log_event( 'Data for integer field "'// trim(adjustl(xios_field_name)) // &
                        '" contains values too large for 16-bit precision', LOG_LEVEL_WARNING )
      end if
      do i = 0, axis_size-1
        send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
                   real( field_proxy%data(i+1:undf:axis_size), dp_xios )
      end do

    type is (r_solver_field_proxy_type)
      do i = 0, axis_size-1
        send_field(i*(domain_size)+1:(i*(domain_size)) + domain_size) = &
                   field_proxy%data(i+1:undf:axis_size)
      end do

    class default
      call log_event( "Invalid type for input field proxy", LOG_LEVEL_ERROR )

  end select

  ! Reshape into 2D horizontal + vertical levels for output
  call xios_send_field(xios_field_name, &
                       reshape (send_field, (/domain_size, axis_size/) ))

  deallocate(send_field)

end subroutine write_field_face

!>  @brief    Write a collection of fields
!>  @details  Iterate over a field collection and write each field if it is
!>            enabled for writing
!>
!>  @param[in]           state   A collection of fields
!>  @param[in,optional]  prefix  A prefix to be added to the field name to
!>                               create the XIOS field ID
!>  @param[in,optional]  suffix  A suffix to be added to the field name to
!>                               create the XIOS field ID
!>
subroutine write_state(state, prefix, suffix)

  implicit none

  type(field_collection_type), intent(inout) :: state
  character(len=*), optional,  intent(in)    :: prefix
  character(len=*), optional,  intent(in)    :: suffix

  type(field_collection_iterator_type) :: iter
  character(str_def)                   :: xios_field_id

  class(field_parent_type), pointer :: fld => null()

  iter = state%get_iterator()
  do
    if ( .not.iter%has_next() ) exit
    fld => iter%next()
    select type(fld)
      type is (field_type)
        if ( fld%can_write() ) then
          write(log_scratch_space,'(3A,I6)') &
              "Writing ", trim(adjustl(fld%get_name()))
          call log_event(log_scratch_space,LOG_LEVEL_INFO)

          ! Construct the XIOS field ID from the LFRic field name and optional arguments
          xios_field_id = trim(adjustl(fld%get_name()))
          if ( present(prefix) ) xios_field_id = trim(adjustl(prefix)) // trim(adjustl(xios_field_id))
          if ( present(suffix) ) xios_field_id = trim(adjustl(xios_field_id)) // trim(adjustl(suffix))

          call fld%write_field(xios_field_id)
        else

          call log_event( 'Write method for '// trim(adjustl(fld%get_name())) // &
                      ' not set up', LOG_LEVEL_INFO )

        end if
      type is (integer_field_type)
        if ( fld%can_write() ) then
          write(log_scratch_space,'(3A,I6)') &
              "Writing ", trim(adjustl(fld%get_name()))
          call log_event(log_scratch_space,LOG_LEVEL_INFO)

          ! Construct the XIOS field ID from the LFRic field name and optional arguments
          xios_field_id = trim(adjustl(fld%get_name()))
          if ( present(prefix) ) xios_field_id = trim(adjustl(prefix)) // trim(adjustl(xios_field_id))
          if ( present(suffix) ) xios_field_id = trim(adjustl(xios_field_id)) // trim(adjustl(suffix))

          call fld%write_field(xios_field_id)
        else

          call log_event( 'Write method for '// trim(adjustl(fld%get_name())) // &
                      ' not set up', LOG_LEVEL_INFO )

        end if

      type is (r_solver_field_type)
        if ( fld%can_write() ) then
          write(log_scratch_space,'(3A,I6)') &
              "Writing ", trim(adjustl(fld%get_name()))
          call log_event(log_scratch_space,LOG_LEVEL_INFO)

          ! Construct the XIOS field ID from the LFRic field name and optional
          ! arguments
          xios_field_id = trim(adjustl(fld%get_name()))
          if ( present(prefix) ) xios_field_id = trim(adjustl(prefix)) // trim(adjustl(xios_field_id))
          if ( present(suffix) ) xios_field_id = trim(adjustl(xios_field_id)) // trim(adjustl(suffix))

          call fld%write_field(xios_field_id)
        else

          call log_event( 'Write method for '// trim(adjustl(fld%get_name())) // &
                      ' not set up', LOG_LEVEL_INFO )

        end if

    end select
  end do

  nullify(fld)

end subroutine write_state

!>  @brief    Write a checkpoint from a collection of fields
!>  @details  Iterate over a field collection and checkpoint each field
!>            if it is enabled for checkpointing
!>
!>  @param[in]  state  Fields to checkpoint.
!>  @param[in]  clock  Model time
!>  @param[in]  checkpoint_stem_name  The checkpoint file stem name
!>
subroutine write_checkpoint( state, clock, checkpoint_stem_name )

  implicit none

  type(field_collection_type), intent(inout) :: state
  class(clock_type),           intent(in)    :: clock
  character(len=*),            intent(in)    :: checkpoint_stem_name

  type(field_collection_iterator_type) :: iter

  class(field_parent_type), pointer :: fld => null()

  iter = state%get_iterator()
  do
     if ( .not.iter%has_next() ) exit
     fld => iter%next()
     select type(fld)
     type is (field_type)
        if ( fld%can_checkpoint() ) then
           write(log_scratch_space,'(2A)') &
                "Checkpointing ", trim(adjustl(fld%get_name()))
           call log_event(log_scratch_space, LOG_LEVEL_INFO)
           call fld%write_checkpoint( "checkpoint_" //                    &
                                      trim(adjustl(fld%get_name())),      &
                                      trim(ts_fname(checkpoint_stem_name, &
                                      "",                                 &
                                      trim(adjustl(fld%get_name())),      &
                                      clock%get_step(),                   &
                                      "")) )
        else if ( fld%can_write() ) then
           write(log_scratch_space,'(2A)') &
                "Writing checkpoint for ", trim(adjustl(fld%get_name()))
           call log_event(log_scratch_space, LOG_LEVEL_INFO)
           call fld%write_field( "checkpoint_" // trim(adjustl(fld%get_name())) )
        else
           call log_event( 'Writing not set up for '// trim(adjustl(fld%get_name())), &
                          LOG_LEVEL_INFO )
        end if
     type is (integer_field_type)
        if ( fld%can_checkpoint() ) then
           write(log_scratch_space,'(2A)') &
                "Checkpointing ", trim(adjustl(fld%get_name()))
           call log_event(log_scratch_space, LOG_LEVEL_INFO)
           call fld%write_checkpoint( "checkpoint_"                      &
                                     // trim(adjustl(fld%get_name())),   &
                                     trim(ts_fname(checkpoint_stem_name, &
                                     "",                                 &
                                     trim(adjustl(fld%get_name())),      &
                                     clock%get_step(),                   &
                                     "")) )
        else if ( fld%can_write() ) then
           write(log_scratch_space,'(2A)') &
                "Writing checkpoint for ", trim(adjustl(fld%get_name()))
           call log_event(log_scratch_space, LOG_LEVEL_INFO)
           call fld%write_field( "checkpoint_" // trim(adjustl(fld%get_name())) )
        else
           call log_event( 'Writing not set up for '// trim(adjustl(fld%get_name())), &
                LOG_LEVEL_INFO )
        end if
     class default
        call log_event('write_checkpoint:Invalid type of field, not supported supported',LOG_LEVEL_ERROR)
     end select
  end do

  nullify(fld)

end subroutine write_checkpoint

end module lfric_xios_write_mod
