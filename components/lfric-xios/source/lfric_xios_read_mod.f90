!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!>  @brief    Module for field reading routines
!>  @details  Holds all routines for reading LFRic fields. All routines are set
!>            up to read data with dimension ordering according to the
!>            recommendations in the NetCDF CF standard.
!>
module lfric_xios_read_mod

  use constants_mod,                 only: i_def, dp_xios, str_def, r_def, &
                                           r_solver
  use field_mod,                     only: field_type, field_proxy_type
  use field_collection_mod,          only: field_collection_type, &
                                           field_collection_iterator_type
  use field_parent_mod,              only: field_parent_type, &
                                           field_parent_proxy_type
  use files_config_mod,              only: checkpoint_stem_name
  use fs_continuity_mod,             only: W3, WTheta
  use integer_field_mod,             only: integer_field_type, &
                                           integer_field_proxy_type
  use r_solver_field_mod,            only: r_solver_field_type, &
                                           r_solver_field_proxy_type
  use io_mod,                        only: ts_fname
  use log_mod,                       only: log_event,         &
                                           log_scratch_space, &
                                           LOG_LEVEL_INFO,    &
                                           LOG_LEVEL_ERROR
  use xios,                          only: xios_field,           &
                                           xios_recv_field,      &
                                           xios_get_handle,      &
                                           xios_get_attr,        &
                                           xios_get_domain_attr, &
                                           xios_get_axis_attr

  implicit none

  private
  public :: checkpoint_read_xios,    &
            read_field_node,         &
            read_field_edge,         &
            read_field_face,         &
            read_field_single_face,  &
            read_state,              &
            read_checkpoint,         &
            read_field_time_var,     &
            read_time_data

contains

!>  @brief    I/O handler for reading an XIOS netcdf checkpoint
!>  @details  Note this routine accepts a filename but doesn't use it - this is
!>           to keep the interface the same for all methods
!>
!>  @param[in]      xios_field_name  XIOS identifier for the field
!>  @param[in]      file_name        Name of the file to read
!>  @param[in,out]  field_proxy      A field proxy to read data into
!>
subroutine checkpoint_read_xios(xios_field_name, file_name, field_proxy)

  implicit none

  character(len=*),               intent(in)    :: xios_field_name
  character(len=*),               intent(in)    :: file_name
  class(field_parent_proxy_type), intent(inout) :: field_proxy

  integer(i_def) :: undf

  ! We only read in up to undf for the partition
  undf = field_proxy%vspace%get_last_dof_owned()

  select type(field_proxy)

    type is (field_proxy_type)
    call xios_recv_field(xios_field_name, field_proxy%data(1:undf))

    class default
    call log_event( "Invalid type for input field proxy", LOG_LEVEL_ERROR )

  end select

end subroutine checkpoint_read_xios

!>  @brief  Read node UGRID data on to a W0 field via XIOS
!>
!>  @param[in]     xios_field_name  XIOS identifier for the field
!>  @param[inout]  field_proxy      A field proxy to read data into
!>
subroutine read_field_node(xios_field_name, field_proxy)

  implicit none

  character(len=*),               intent(in)    :: xios_field_name
  class(field_parent_proxy_type), intent(inout) :: field_proxy

  integer(i_def) :: i, undf
  integer(i_def) :: domain_size, axis_size
  real(dp_xios), allocatable :: recv_field(:)

  undf = field_proxy%vspace%get_last_dof_owned()

  ! Get the expected horizontal and vertical axis size
  call xios_get_domain_attr('node', ni=domain_size)
  call xios_get_axis_attr("vert_axis_full_levels", n_glo=axis_size)

  ! Size the arrays to be what is expected
  allocate(recv_field(domain_size*axis_size))

  ! Read the data from XIOS to a temporary 1D array
  call xios_recv_field(xios_field_name, recv_field)

  ! Reshape the data to what we require for the LFRic field
  select type(field_proxy)

    type is (field_proxy_type)
    do i = 0, axis_size-1
      field_proxy%data( i+1 : undf : axis_size )  = &
                  real(recv_field( i*(domain_size)+1 : (i+1)*domain_size ), r_def)
    end do

    type is (integer_field_proxy_type)
    do i = 0, axis_size-1
      field_proxy%data( i+1 : undf : axis_size )  = &
                  int(recv_field( i*(domain_size)+1 : (i+1)*domain_size ), i_def)
    end do

    class default
      call log_event( "Invalid type for input field proxy", LOG_LEVEL_ERROR )

  end select

  deallocate(recv_field)

end subroutine read_field_node

!>  @brief  Read half-level edge UGRID data on to a W2H field via XIOS
!>
!>  @param[in]     xios_field_name  XIOS identifier for the field
!>  @param[inout]  field_proxy      A field proxy to read data into
!>
subroutine read_field_edge(xios_field_name, field_proxy)

  implicit none

  character(len=*),               intent(in)    :: xios_field_name
  class(field_parent_proxy_type), intent(inout) :: field_proxy

  integer(i_def) :: i, undf
  integer(i_def) :: domain_size, axis_size
  real(dp_xios), allocatable :: recv_field(:)

  undf = field_proxy%vspace%get_last_dof_owned()

  ! Get the expected horizontal and vertical axis size
  call xios_get_domain_attr('edge', ni=domain_size)
  call xios_get_axis_attr("vert_axis_half_levels", n_glo=axis_size)

  ! Size the arrays to be what is expected
  allocate(recv_field(domain_size*axis_size))

  ! Read the data from XIOS to a temporary 1D array
  call xios_recv_field(xios_field_name, recv_field)

  ! Reshape the data to what we require for the LFRic field
  ! Note the conversion from dp_xios to r_def or i_def
  select type(field_proxy)

    type is (field_proxy_type)
    do i = 0, axis_size-1
      field_proxy%data( i+1 : undf : axis_size )  = &
                  real(recv_field( i*(domain_size)+1 : (i+1)*domain_size ), r_def)
    end do

    type is (integer_field_proxy_type)
    do i = 0, axis_size-1
      field_proxy%data( i+1 : undf : axis_size )  = &
                  int( recv_field( i*(domain_size)+1 : (i+1)*domain_size ), i_def)
    end do

    type is (r_solver_field_proxy_type)
    do i = 0, axis_size-1
      field_proxy%data( i+1 : undf : axis_size )  = &
               real(recv_field( i*(domain_size)+1 : (i+1)*domain_size ), r_solver)
    end do

    class default
      call log_event( "Invalid type for input field proxy", LOG_LEVEL_ERROR )

  end select

  deallocate(recv_field)

end subroutine read_field_edge

!>  @brief  Read full/half-level face UGRID data on to a W3/WTheta field via XIOS
!>
!>  @param[in]  xios_field_name  XIOS identifier for the field
!>  @param[in]  field_proxy      A field proxy to read data into
!>
subroutine read_field_face(xios_field_name, field_proxy)

  implicit none

  character(len=*),               intent(in)    :: xios_field_name
  class(field_parent_proxy_type), intent(inout) :: field_proxy

  integer(i_def) :: i, undf
  integer(i_def) :: fs_id
  integer(i_def) :: domain_size, axis_size
  real(dp_xios), allocatable :: recv_field(:)

  ! Get the size of undf as we only read in up to last owned
  undf = field_proxy%vspace%get_last_dof_owned()
  fs_id = field_proxy%vspace%which()

  ! get the horizontal / vertical domain sizes
  if ( fs_id == W3 ) then
    call xios_get_domain_attr('face', ni=domain_size)
    call xios_get_axis_attr("vert_axis_half_levels", n_glo=axis_size)
  else
    call xios_get_domain_attr('face', ni=domain_size)
    call xios_get_axis_attr("vert_axis_full_levels", n_glo=axis_size)
  end if

  ! Size the array to be what is expected
  allocate(recv_field(domain_size*axis_size))

  ! Read the data into a temporary array - this should be in the correct order
  ! as long as we set up the horizontal domain using the global index
  call xios_recv_field(xios_field_name, recv_field)

  ! Different field kinds are selected to access data, which is arranged to get the
  ! correct data layout for the LFRic field - the reverse of what is done for writing
  ! Note the conversion from dp_xios to r_def, i_def or r_solver
  select type(field_proxy)

    type is (field_proxy_type)
    do i = 0, axis_size-1
      field_proxy%data(i+1:undf:axis_size) = &
             real(recv_field(i*(domain_size)+1:(i*(domain_size)) + domain_size), r_def)
    end do

    type is (integer_field_proxy_type)
    do i = 0, axis_size-1
      field_proxy%data(i+1:undf:axis_size) = &
             int( recv_field(i*(domain_size)+1:(i*(domain_size)) + domain_size), i_def)
    end do

    type is (r_solver_field_proxy_type)
    do i = 0, axis_size-1
      field_proxy%data(i+1:undf:axis_size) = &
             real( recv_field(i*(domain_size)+1:(i*(domain_size)) + domain_size), r_solver)
    end do

    class default
      call log_event( "Invalid type for input field proxy", LOG_LEVEL_ERROR )

  end select

  deallocate(recv_field)

end subroutine read_field_face

!>  @brief  Read face domain UGRID data on to a single-level W3 field via XIOS
!>
!>  @param[in]  xios_field_name  XIOS identifier for the field
!>  @param[in]  field_proxy      A field proxy to read data into
!>
subroutine read_field_single_face(xios_field_name, field_proxy)

  implicit none

  character(len=*),               intent(in)    :: xios_field_name
  class(field_parent_proxy_type), intent(inout) :: field_proxy

  integer(i_def) :: i, undf, ndata
  integer(i_def) :: domain_size
  real(dp_xios), allocatable :: recv_field(:)

  ! Get the size of undf as we only read in up to last owned
  undf = field_proxy%vspace%get_last_dof_owned()
  ndata = field_proxy%vspace%get_ndata()

  ! Get the expected horizontal size
  ! all 2D fields are nominally in W3, hence half levels
  call xios_get_domain_attr('face', ni=domain_size)

  ! Size the array to be what is expected
  allocate(recv_field(domain_size*ndata))

  ! If the fields do not have the same size, then exit with error
  if ( size(recv_field) /= undf ) then
    call log_event( "Global size of model field /= size of field from file", &
                    LOG_LEVEL_ERROR )
  end if

  ! Read the data into a temporary array - this should be in the correct order
  ! as long as we set up the horizontal domain using the global index
  call xios_recv_field(xios_field_name, recv_field)

  ! We need to reshape the incoming data to get the correct data layout for the LFRic
  ! multi-data field
  ! We need to cast to the field kind to access the data
  select type(field_proxy)

    ! Pass the correct data from the recieved field to the field proxy data
    ! Note the conversion from dp_xios to r_def, i_def or r_solver
    type is (field_proxy_type)
    do i = 0, ndata-1
      field_proxy%data(i+1:(ndata*domain_size)+i:ndata) = &
              real(recv_field(i*(domain_size)+1:(i*(domain_size)) + domain_size), r_def)
    end do

    type is (integer_field_proxy_type)
    do i = 0, ndata-1
      field_proxy%data(i+1:(ndata*domain_size)+i:ndata) = &
              int( recv_field(i*(domain_size)+1:(i*(domain_size)) + domain_size), i_def)
    end do

    type is (r_solver_field_proxy_type)
    do i = 0, ndata-1
      field_proxy%data(i+1:(ndata*domain_size)+i:ndata) = &
              real(recv_field(i*(domain_size)+1:(i*(domain_size)) + domain_size), r_solver)
    end do

    class default
      call log_event( "Invalid type for input field proxy", LOG_LEVEL_ERROR )

  end select

  deallocate(recv_field)

end subroutine read_field_single_face

!>  @brief  Read a time-varying field in UGRID format using XIOS
!>
!>  @param[in]     xios_field_name  XIOS identifier for the field
!>  @param[inout]  field_proxy      A field proxy to be read into
!>  @param[in]     time_index       The indices of the time 'columns' to be
!>                                  read in
!>
subroutine read_field_time_var(xios_field_name, field_proxy, time_indices)

  implicit none

  character(len=*),       intent(in)    :: xios_field_name
  type(field_proxy_type), intent(inout) :: field_proxy
  integer(i_def),         intent(in)    :: time_indices(:)

  integer(i_def) :: undf, fs_id, i, j, k, nlayers, ndata, time_index, vert_levels
  integer(i_def) :: domain_size, vert_axis_size, time_axis_size, start_index
  real(dp_xios), allocatable :: recv_field(:)
  real(r_def),   allocatable :: ndata_slice(:)
  real(r_def),   allocatable :: time_slice(:)
  real(r_def),   allocatable :: field_data(:)

  ! Get the number of layers to distiniguish between 2D and 3D fields
  nlayers = field_proxy%vspace%get_nlayers()
  ! Get the size of the multi-data field ndata axis (ndata is multi-data
  ! multiplied by length of time window so we divide by that)

  ndata = field_proxy%vspace%get_ndata() / size(time_indices)

  ! Get the size of undf as we only read in up to last owned
  undf = field_proxy%vspace%get_last_dof_owned()
  fs_id = field_proxy%vspace%which()

  ! get the horizontal / vertical / time domain sizes
  if ( fs_id == W3 ) then
    call xios_get_domain_attr( 'face', ni=domain_size )
    call xios_get_axis_attr( 'vert_axis_half_levels', n_glo=vert_axis_size )
    call xios_get_axis_attr( 'monthly_axis', n_glo=time_axis_size )
  else if ( fs_id == WTheta ) then
    call xios_get_domain_attr( 'face', ni=domain_size )
    call xios_get_axis_attr( 'vert_axis_full_levels', n_glo=vert_axis_size )
    call xios_get_axis_attr( "monthly_axis", n_glo=time_axis_size )
  else
    call log_event( 'Time varying fields only readable for W3 or WTheta function spaces', &
                     LOG_LEVEL_ERROR )
  end if

  ! Define vertical levels based on whether we are on a 2D mesh
  if ( nlayers == 1 ) then
    vert_levels = 1
  else
    vert_levels = vert_axis_size
  end if

  ! Size the various array slices
  allocate( recv_field( domain_size * vert_levels * time_axis_size * ndata ) )
  allocate( ndata_slice( domain_size * vert_levels * time_axis_size ) )
  allocate( time_slice( domain_size * vert_levels ) )
  allocate( field_data( undf ) )

  ! Read the data into a temporary array
  call xios_recv_field( trim(xios_field_name)//'_data', recv_field )

  ! Incoming data is shaped with multi-data axis first, then time axis, so set
  ! up an array for each multi-data level
  do i = 0, ndata - 1

    !Get first ndata slice - note the conversion from double precision to r_def
    ndata_slice = real( recv_field( i * ( domain_size * vert_levels * time_axis_size ) + 1 :  &
                                  ( i + 1 ) * ( domain_size * vert_levels * time_axis_size ) ), &
                                  kind=r_def )

    ! Reshape data into a single array for each time entry in the time window
    do j = 0, size(time_indices)-1
      time_index = time_indices(j+1)

      ! Get correct time-entry from current multi-data level
      time_slice = ndata_slice( ( time_index - 1 ) * ( domain_size * vert_levels ) + 1 :  &
                                 ( time_index ) * ( domain_size * vert_levels ) )

      ! We require multi-data fields with vertical levels to be multi-data first
      if ( ndata /= 1 .and. vert_levels /= 1 .and. &
           .not. field_proxy%is_ndata_first() ) then
        write( log_scratch_space,'(A,A)' ) "Only ndata_first ordering supported for read_field_time_var: "// &
                                      trim( xios_field_name )
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )

      else
        do k = 0, vert_levels-1
          start_index = (k * ndata * size(time_indices)) + (i * size(time_indices)) + j + 1

          field_data( start_index : undf : ndata*size(time_indices)*vert_levels ) &
                      = time_slice(k*(domain_size)+1:(k+1)*domain_size)

        end do

      end if

    end do

  end do

  ! Pass reshaped data array to field object via proxy
  field_proxy%data = field_data

  ! Set halos dirty here as for parallel read we only read in data for owned
  ! dofs and the halos will not be set
  call field_proxy%set_dirty()

  deallocate( recv_field )
  deallocate( ndata_slice )
  deallocate( time_slice )
  deallocate( field_data )

end subroutine read_field_time_var

!>  @brief  Read time data using XIOS
!>
!>  @param[in]   time_id    The XIOS ID for the time data
!>  @param[out]  time_data  The array of time data
!>
subroutine read_time_data(time_id, time_data)

  implicit none

  character(len=*),           intent(in)  :: time_id
  real(dp_xios), allocatable, intent(out) :: time_data(:)

  ! Local variables for XIOS interface
  type(xios_field)   :: time_hdl
  integer(i_def)     :: time_axis_size
  character(str_def) :: axis_id, time_units

  ! Set up axis size and units from XIOS configuration
  call xios_get_handle( time_id, time_hdl )
  call xios_get_attr( time_hdl, unit=time_units, axis_ref=axis_id )

  call xios_get_axis_attr( axis_id, n_glo=time_axis_size )
  allocate( time_data( time_axis_size ) )

  ! Read the time data from the ancil file
  call xios_recv_field( time_id, time_data )

  ! Set time units to seconds across a single year - the calendar type needs to
  ! be integrated into this eventually
  if ( time_units == 'days' ) then
    time_data = mod( time_data, 365.0_dp_xios )
    time_data = time_data * 3600.0_dp_xios * 24.0_dp_xios
  else if ( time_units == 'hours' ) then
    time_data = mod( time_data, 365.0_dp_xios * 24.0_dp_xios )
    time_data = time_data * 3600.0_dp_xios
  else if ( time_units == 'months' ) then
    time_data = (real(time_data-1, dp_xios)/12.0_dp_xios) * 365.0_dp_xios * 24.0_dp_xios * 3600.0_dp_xios
  else
    write( log_scratch_space,'(A,A)' ) "Invalid units for time axis: "// &
                                      trim( time_units )
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

end subroutine read_time_data

!>  @brief    Read into a collection of fields
!>  @details  Iterate over a field collection and read each field
!>            into a collection, if it is enabled for read
!>
!>  @param[in,out]       state   The collection of fields to populate
!>  @param[in,optional]  prefix  A prefix to be added to the field name to
!>                               create the XIOS field ID
!>  @param[in,optional]  suffix  A suffix to be added to the field name to
!>                               create the XIOS field ID
!>
subroutine read_state(state, prefix, suffix)

  implicit none

  type( field_collection_type ), intent(inout) :: state
  character( len=* ), optional,  intent(in)    :: prefix
  character( len=* ), optional,  intent(in)    :: suffix

  type( field_collection_iterator_type) :: iter
  character( str_def )                  :: xios_field_id

  class( field_parent_type ), pointer :: fld => null()

  iter = state%get_iterator()

  do
    if ( .not.iter%has_next() ) exit
    fld => iter%next()
    select type(fld)
      type is (field_type)
        if ( fld%can_read() ) then
          call log_event( 'Reading '//trim(adjustl(fld%get_name())), &
                          LOG_LEVEL_INFO )

          ! Construct the XIOS field ID from the LFRic field name and optional arguments
          xios_field_id = trim(adjustl(fld%get_name()))
          if ( present(prefix) ) xios_field_id = trim(adjustl(prefix)) // trim(adjustl(xios_field_id))
          if ( present(suffix) ) xios_field_id = trim(adjustl(xios_field_id)) // trim(adjustl(suffix))

          call fld%read_field(xios_field_id)
        else
          call log_event('Read method for  '//trim(adjustl(fld%get_name()))// &
                         ' not set up', LOG_LEVEL_INFO )
        end if

      type is (integer_field_type)
        if ( fld%can_read() ) then
          call log_event( &
            'Reading '//trim(adjustl(fld%get_name())), &
            LOG_LEVEL_INFO)

          ! Construct the XIOS field ID from the LFRic field name and optional arguments
          xios_field_id = trim(adjustl(fld%get_name()))
          if ( present(prefix) ) xios_field_id = trim(adjustl(prefix)) // trim(adjustl(xios_field_id))
          if ( present(suffix) ) xios_field_id = trim(adjustl(xios_field_id)) // trim(adjustl(suffix))

          call fld%read_field(xios_field_id)
        else
          call log_event( 'Read method for  '// trim(adjustl(fld%get_name())) // &
                          ' not set up', LOG_LEVEL_INFO )
        end if

      type is (r_solver_field_type)
        if ( fld%can_read() ) then
          call log_event( 'Reading '//trim(adjustl(fld%get_name())), &
            LOG_LEVEL_INFO)

          ! Construct the XIOS field ID from the LFRic field name and optional
          ! arguments
          xios_field_id = trim(adjustl(fld%get_name()))
          if ( present(prefix) ) xios_field_id = trim(adjustl(prefix)) // trim(adjustl(xios_field_id))
          if ( present(suffix) ) xios_field_id = trim(adjustl(xios_field_id)) // trim(adjustl(suffix))

          call fld%read_field(xios_field_id)
        else
          call log_event( 'Read method for  '// trim(adjustl(fld%get_name())) // &
                          ' not set up', LOG_LEVEL_INFO )
        end if

    end select
  end do

  nullify(fld)

end subroutine read_state

!>  @brief    Read from a checkpoint into a collection of fields
!>  @details  Iterate over a field collection and read each field
!>            into a collection, if it is enabled for checkpointing
!>
!>  @param[in]  state     The collection of fields to populate
!>  @param[in]  timestep  The current timestep
!>
subroutine read_checkpoint(state, timestep)

  implicit none

  type( field_collection_type ), intent(inout) :: state
  integer(i_def),                intent(in)    :: timestep

  type( field_collection_iterator_type) :: iter

  class( field_parent_type ), pointer :: fld => null()

  iter = state%get_iterator()
  do
    if ( .not.iter%has_next() ) exit
    fld => iter%next()
    select type(fld)
      type is (field_type)
        if ( fld%can_checkpoint() ) then
          call log_event( 'Reading checkpoint file to restart '// &
                           trim(adjustl(fld%get_name())), LOG_LEVEL_INFO)
          call fld%read_checkpoint( "restart_"//trim(adjustl(fld%get_name())), &
                                    trim(ts_fname(checkpoint_stem_name, "",    &
                                    trim(adjustl(fld%get_name())),timestep,"")) )
        else
          call log_event( 'Checkpointing for  '// trim(adjustl(fld%get_name())) // &
                          ' not set up', LOG_LEVEL_INFO )
        end if
      type is (integer_field_type)
        if ( fld%can_checkpoint() ) then
          call log_event( 'Reading checkpoint file to restart '// &
                           trim(adjustl(fld%get_name())), LOG_LEVEL_INFO)
          call fld%read_checkpoint( "restart_"//trim(adjustl(fld%get_name())), &
                                    trim(ts_fname(checkpoint_stem_name, "",    &
                                    trim(adjustl(fld%get_name())),timestep,"")) )
        else
          call log_event( 'Checkpointing for  '// trim(adjustl(fld%get_name())) // &
                          ' not set up', LOG_LEVEL_INFO )
        end if

      type is (r_solver_field_type)
        if ( fld%can_checkpoint() ) then
          call log_event( 'Reading checkpoint file to restart '// &
                           trim(adjustl(fld%get_name())), LOG_LEVEL_INFO)
          call fld%read_checkpoint( "restart_"//trim(adjustl(fld%get_name())), &
                                    trim(ts_fname(checkpoint_stem_name, "",    &
                                    trim(adjustl(fld%get_name())),timestep,"")))
        else
          call log_event( 'Checkpointing for  '// trim(adjustl(fld%get_name())) // &
                          ' not set up', LOG_LEVEL_INFO )
        end if

    end select
  end do

  nullify(fld)

end subroutine read_checkpoint

end module lfric_xios_read_mod
