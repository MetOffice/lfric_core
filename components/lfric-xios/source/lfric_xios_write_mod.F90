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
  use constants_mod,        only: i_def, l_def, str_def, str_max_filename
  use lfric_xios_constants_mod, &
                            only: dp_xios, xios_max_int
  use field_real32_mod,     only: field_real32_type, field_real32_proxy_type
  use field_real64_mod,     only: field_real64_type, field_real64_proxy_type
  use field_parent_mod,     only: field_parent_proxy_type
  use field_collection_iterator_mod, &
                            only: field_collection_iterator_type
  use field_collection_mod, only: field_collection_type
  use field_parent_mod,     only: field_parent_type
  use fs_continuity_mod,    only: W3
  use io_mod,               only: ts_fname
  use integer_field_mod,    only: integer_field_type, integer_field_proxy_type
  use lfric_xios_utils_mod, only: prime_io_mesh_is
  use lfric_xios_format_mod, &
                            only: format_field
  use mesh_mod,             only: mesh_type
  use log_mod,              only: log_event,         &
                                  log_scratch_space, &
                                  LOG_LEVEL_INFO,    &
                                  LOG_LEVEL_WARNING, &
                                  LOG_LEVEL_ERROR
#ifdef UNIT_TEST
  use lfric_xios_mock_mod,  only: xios_send_field,      &
                                  xios_get_domain_attr, &
                                  xios_get_axis_attr,   &
                                  lfric_xios_mock_pull_in
#else
  use lfric_xios_mock_mod,  only: lfric_xios_mock_pull_in
  use xios,                 only: xios_send_field,      &
                                  xios_get_domain_attr, &
                                  xios_get_axis_attr
#endif

  implicit none

  private
  public :: checkpoint_write_xios,    &
            write_field_generic,      &
            write_state,              &
            write_checkpoint

contains

!>  @brief  Write field data to UGRIDs via XIOS
!>
!>  @param[in]     field_name       Field name (for error reporting only)
!>  @param[in]     field_proxy      A field proxy to be written
!>
subroutine write_field_generic(field_name, field_proxy)
  use lfric_xios_diag_mod,        only:  get_field_domain_ref
  implicit none

  character(len=*), optional,     intent(in) :: field_name
  class(field_parent_proxy_type), intent(in) :: field_proxy

  integer(i_def) :: undf
  integer(i_def) :: hdim          ! horizontal dimension, domain size
  integer(i_def) :: vdim          ! vertical dimension
  real(dp_xios), allocatable :: xios_data(:)
  logical(l_def) :: legacy

  undf = field_proxy%vspace%get_last_dof_owned() ! total dimension

  vdim = field_proxy%vspace%get_ndata() * size(field_proxy%vspace%get_levels())

  hdim = undf/vdim

  ! detect field with legacy checkpointing domain
  legacy = (index(get_field_domain_ref(field_name), 'checkpoint_') == 1)

  ! sanity check
  if (.not. legacy .and. .not. (hdim*vdim == undf)) then
    call log_event('assertion failed for field ' // field_name                &
      // ': hdim*vdim == undf', log_level_error)
  end if

  allocate(xios_data(undf))

  call format_field(xios_data, field_name, field_proxy, vdim, hdim, legacy)

  if (legacy) then
    call xios_send_field( field_name, reshape(xios_data, (/ 1, undf /)) )
  else
    call xios_send_field( field_name, reshape(xios_data, (/vdim, hdim/)) )
    ! The shape is only necessary for the mock implementation, and
    ! the only thing that matters is the product of the dimensions.
  end if

  deallocate(xios_data)

end subroutine write_field_generic

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

    type is (field_real32_proxy_type)
    send_field = field_proxy%data(1:undf)

    type is (field_real64_proxy_type)
    send_field = field_proxy%data(1:undf)

    type is (integer_field_proxy_type)
    if ( any( abs(field_proxy%data(1:undf)) > xios_max_int) ) then
      call log_event( 'Data for integer field "'// trim(adjustl(xios_field_name)) // &
                      '" contains values too large for 16-bit precision', LOG_LEVEL_WARNING )
    end if
    send_field = real( field_proxy%data(1:undf), dp_xios )

    class default
    call log_event( "Invalid type for input field proxy", LOG_LEVEL_ERROR )

  end select

  call xios_send_field("checkpoint_"//trim(xios_field_name), reshape (send_field, (/1, undf/)))

end subroutine checkpoint_write_xios

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

  ! Create the iter iterator on the state collection
  call iter%initialise(state)
  do
    if ( .not.iter%has_next() ) exit
    fld => iter%next()
    select type(fld)
      type is (field_real32_type)
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
      type is (field_real64_type)
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
!>  @param[in,optional]  prefix  A prefix to be added to the field name to
!>                               create the XIOS field ID
!>  @param[in,optional]  suffix  A suffix to be added to the field name to
!>                               create the XIOS field ID
!>
subroutine write_checkpoint( state, clock, checkpoint_stem_name, prefix, suffix )

  implicit none

  type(field_collection_type), intent(inout) :: state
  class(clock_type),           intent(in)    :: clock
  character(len=*),            intent(in)    :: checkpoint_stem_name
  character(len=*), optional,  intent(in)    :: prefix
  character(len=*), optional,  intent(in)    :: suffix

  type(field_collection_iterator_type) :: iter

  class(field_parent_type), pointer    :: fld => null()

  character(str_def)                   :: xios_field_id

  ! Create the iter iterator on the state collection
  call iter%initialise(state)
  do
     if ( .not.iter%has_next() ) exit
     fld => iter%next()
     ! Construct the XIOS field ID from the LFRic field name and optional arguments
     xios_field_id = trim(adjustl(fld%get_name()))
     if ( present(prefix) ) xios_field_id = trim(adjustl(prefix)) // trim(adjustl(xios_field_id))
     if ( present(suffix) ) xios_field_id = trim(adjustl(xios_field_id)) // trim(adjustl(suffix))
     select type(fld)
     type is (field_real32_type)
        if ( fld%can_checkpoint() ) then
           write(log_scratch_space,'(2A)') &
                "Checkpointing ", xios_field_id
           call log_event(log_scratch_space, LOG_LEVEL_INFO)
           call fld%write_checkpoint( xios_field_id,      &
                                      trim(ts_fname(checkpoint_stem_name, &
                                      "",                                 &
                                      xios_field_id,      &
                                      clock%get_step(),                   &
                                      "")) )
        else if ( fld%can_write() ) then
           write(log_scratch_space,'(2A)') &
                "Writing checkpoint for ", xios_field_id
           call log_event(log_scratch_space, LOG_LEVEL_INFO)
           call fld%write_field( "checkpoint_" // xios_field_id )
        else
           call log_event( 'Writing not set up for '// xios_field_id, &
                          LOG_LEVEL_INFO )
        end if
     type is (field_real64_type)
        if ( fld%can_checkpoint() ) then
           write(log_scratch_space,'(2A)') &
                "Checkpointing ", xios_field_id
           call log_event(log_scratch_space, LOG_LEVEL_INFO)
           call fld%write_checkpoint( xios_field_id,      &
                                      trim(ts_fname(checkpoint_stem_name, &
                                      "",                                 &
                                      xios_field_id,      &
                                      clock%get_step(),                   &
                                      "")) )
        else if ( fld%can_write() ) then
           write(log_scratch_space,'(2A)') &
                "Writing checkpoint for ", xios_field_id
           call log_event(log_scratch_space, LOG_LEVEL_INFO)
           call fld%write_field( "checkpoint_" // xios_field_id )
        else
           call log_event( 'Writing not set up for '// xios_field_id, &
                          LOG_LEVEL_INFO )
        end if
     type is (integer_field_type)
        if ( fld%can_checkpoint() ) then
           write(log_scratch_space,'(2A)') &
                "Checkpointing ", xios_field_id
           call log_event(log_scratch_space, LOG_LEVEL_INFO)
           call fld%write_checkpoint( trim(adjustl(fld%get_name()) ),     &
                                      trim(ts_fname(checkpoint_stem_name, &
                                      "",                                 &
                                      xios_field_id,      &
                                      clock%get_step(),                   &
                                      "")) )
        else if ( fld%can_write() ) then
           write(log_scratch_space,'(2A)') &
                "Writing checkpoint for ", xios_field_id
           call log_event(log_scratch_space, LOG_LEVEL_INFO)
           call fld%write_field( "checkpoint_" // xios_field_id )
        else
           call log_event( 'Writing not set up for '// xios_field_id, &
                LOG_LEVEL_INFO )
        end if
     class default
        call log_event('write_checkpoint:Invalid type of field, not supported supported',LOG_LEVEL_ERROR)
     end select
  end do

  nullify(fld)

end subroutine write_checkpoint

end module lfric_xios_write_mod
