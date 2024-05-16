!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!>  @brief    Module for field reading routines.
!>  @details  Holds all routines for reading LFRic fields. All routines are set
!>            up to read data with dimension ordering according to the
!>            recommendations in the NetCDF CF standard.
!>
module lfric_xios_read_mod

  use, intrinsic :: iso_fortran_env, only : real32, real64

  use constants_mod,            only: i_def, l_def, str_def, r_def, rmdi, &
                                      LARGE_REAL_NEGATIVE
  use lfric_xios_constants_mod, only: dp_xios
  use field_mod,                only: field_type, field_proxy_type
  use field_real32_mod,         only: field_real32_type, field_real32_proxy_type
  use field_real64_mod,         only: field_real64_type, field_real64_proxy_type
  use field_collection_iterator_mod, &
                                only: field_collection_iterator_type
  use field_collection_mod,     only: field_collection_type
  use field_parent_mod,         only: field_parent_type, &
                                      field_parent_proxy_type
  use fs_continuity_mod,        only: W3, WTheta, W2H, W2, &
                                      is_fs_horizontally_continuous
  use integer_field_mod,        only: integer_field_type, &
                                      integer_field_proxy_type
  use io_mod,                   only: ts_fname
  use lfric_xios_utils_mod,     only: prime_io_mesh_is
  use lfric_xios_format_mod,    only: inverse_format_field

  use mesh_mod,                 only: mesh_type
  use log_mod,                  only: log_event,         &
                                      log_scratch_space, &
                                      LOG_LEVEL_INFO,    &
                                      LOG_LEVEL_ERROR
#ifdef UNIT_TEST
  use lfric_xios_mock_mod,      only: xios_recv_field,      &
                                      xios_get_domain_attr, &
                                      xios_get_axis_attr,   &
                                      xios_get_field_attr,  &
                                      lfric_xios_mock_pull_in
#else
  use lfric_xios_mock_mod,      only: lfric_xios_mock_pull_in
  use xios,                     only: xios_recv_field,      &
                                      xios_get_domain_attr, &
                                      xios_get_axis_attr,   &
                                      xios_get_field_attr
#endif

  implicit none

  private
  public :: checkpoint_read_xios,    &
            read_field_generic,      &
            read_state,              &
            read_checkpoint,         &
            read_field_time_var

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
  integer(i_def) :: fs_id

  ! We only read in up to undf for the partition
  undf = field_proxy%vspace%get_last_dof_owned()

  select type(field_proxy)

    type is (field_real64_proxy_type)
      call xios_recv_field("restart_"//trim(xios_field_name), field_proxy%data(1:undf))
      call field_proxy%set_dirty()
      ! Ensure annexed dofs for continuous fields are initialised
      fs_id = field_proxy%vspace%which()
      if (is_fs_horizontally_continuous(fs_id)) then
        call field_proxy%halo_exchange(depth=1)
      end if

    class default
      call log_event( "Invalid type for input field proxy", LOG_LEVEL_ERROR )

  end select

end subroutine checkpoint_read_xios

!>  @brief    Post-processing after reading field data
!>  @details  Performs a halo swap if necessary
!>
!>  @param[in, out] field_proxy        A field proxy to be written
!>
subroutine post_read(field_proxy)

  implicit none

  class(field_parent_proxy_type), intent(inout) :: field_proxy

  call field_proxy%set_dirty()

  if (is_fs_horizontally_continuous(field_proxy%vspace%which())) then
    ! set the annexed dofs on continuous fields
    call field_proxy%halo_exchange(depth=1)
  end if

end subroutine post_read

!>  @brief   Read field data from UGRIDs via XIOS
!>  @param[in]     field_name       Field name (for error reporting only)
!>  @param[in]     field_proxy      A field proxy to be written
!>
subroutine read_field_generic(xios_field_name, field_proxy)
  use lfric_xios_diag_mod,        only: get_field_domain_ref
  implicit none

  character(len=*),               intent(in)    :: xios_field_name
  class(field_parent_proxy_type), intent(inout) :: field_proxy

  integer(i_def) :: undf
  integer(i_def) :: hdim          ! horizontal dimension, domain size
  integer(i_def) :: vdim          ! vertical dimension
  real(dp_xios), allocatable :: xios_data(:)
  logical(l_def) :: legacy

  undf = field_proxy%vspace%get_last_dof_owned() ! total dimension

  vdim = field_proxy%vspace%get_ndata() * size(field_proxy%vspace%get_levels())

  hdim = undf/vdim

  ! detect field with legacy checkpointing domain
  legacy = (index(get_field_domain_ref(xios_field_name), 'checkpoint_') == 1)

  ! sanity check
  if (.not. legacy .and. .not. (hdim*vdim == undf)) then
    call log_event('assertion failed for field ' // xios_field_name                &
      // ': hdim*vdim must equal undf', log_level_error)
  end if

  allocate(xios_data(undf))

  ! receive the field data from XIOS
  call xios_recv_field(xios_field_name, xios_data)
  ! inverse of xios formatting
  call inverse_format_field(xios_data, xios_field_name, field_proxy, vdim, hdim, legacy)
  ! deal with halo data
  call post_read(field_proxy)

  deallocate(xios_data)

end subroutine read_field_generic

!>  @brief  Read a time-varying field, with given time dimension, in UGRID format using XIOS
!>
!>  @param[in]     xios_field_name  XIOS identifier for the field
!>  @param[inout]  field_proxy      A field proxy to be read into
!>  @param[in]     time_index       The indices of the time 'columns' to be
!>                                  read in
!>  @param[in]     time_axis_size   Placeholder time axis size needed before #3265
!>
subroutine read_field_time_var(xios_field_name, field_proxy, time_indices, time_axis_size)

  implicit none

  character(len=*),       intent(in)    :: xios_field_name
  type(field_proxy_type), intent(inout) :: field_proxy
  integer(i_def),         intent(in)    :: time_indices(:)
  integer(i_def),         intent(in)    :: time_axis_size

  integer(i_def) :: undf, fs_id, i, j, k, nlayers, ndata, time_index, vert_levels
  integer(i_def) :: domain_size, vert_axis_size, start_index
  real(dp_xios), allocatable :: recv_field(:)
  real(r_def),   allocatable :: ndata_slice(:)
  real(r_def),   allocatable :: time_slice(:)
  real(r_def),   allocatable :: field_data(:)

  type(mesh_type), pointer   :: mesh => null()

  ! Call error if field not on prime mesh
  mesh => field_proxy%vspace%get_mesh()

  fs_id = field_proxy%vspace%which()
  if ( fs_id /= W3 .and. fs_id /= WTheta .and. fs_id /= W2H ) then
    call log_event( 'Time varying fields only readable for W3, WTheta or W2H function spaces', &
                     LOG_LEVEL_ERROR )
  end if

  ! Get the number of layers to distiniguish between 2D and 3D fields
  nlayers = field_proxy%vspace%get_nlayers()
  ! Get the size of the multi-data field ndata axis (ndata is multi-data
  ! multiplied by length of time window so we divide by that)

  ndata = field_proxy%vspace%get_ndata() / size(time_indices)

  ! Get the size of undf as we only read in up to last owned
  undf = field_proxy%vspace%get_last_dof_owned()
  fs_id = field_proxy%vspace%which()

  ! get the horizontal / vertical / time domain sizes
  if (prime_io_mesh_is(mesh)) then
    if ( fs_id == W3 ) then
      call xios_get_domain_attr( 'face', ni=domain_size )
      call xios_get_axis_attr( 'vert_axis_half_levels', n_glo=vert_axis_size )
    else if ( fs_id == WTheta ) then
      call xios_get_domain_attr( 'face', ni=domain_size )
      call xios_get_axis_attr( 'vert_axis_full_levels', n_glo=vert_axis_size )
    else if ( fs_id == W2H ) then
      call xios_get_domain_attr( 'edge', ni=domain_size )
      call xios_get_axis_attr( 'vert_axis_half_levels', n_glo=vert_axis_size )
    else
      call log_event( 'Time varying fields only readable for W3, WTheta or W2H function spaces', &
                      LOG_LEVEL_ERROR )
    end if
  else
    if ( fs_id == W3 ) then
      call xios_get_domain_attr( trim(adjustl(mesh%get_mesh_name()))//"_face", ni=domain_size )
      call xios_get_axis_attr( 'vert_axis_half_levels', n_glo=vert_axis_size )
    else if ( fs_id == WTheta ) then
      call xios_get_domain_attr( trim(adjustl(mesh%get_mesh_name()))//"_face", ni=domain_size )
      call xios_get_axis_attr( 'vert_axis_full_levels', n_glo=vert_axis_size )
    else if ( fs_id == W2H ) then
      call xios_get_domain_attr( trim(adjustl(mesh%get_mesh_name()))//"_edge", ni=domain_size )
      call xios_get_axis_attr( 'vert_axis_half_levels', n_glo=vert_axis_size )
    else
      call log_event( 'Time varying fields only readable for W3, WTheta or W2H function spaces', &
                      LOG_LEVEL_ERROR )
    end if
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

      ! Temporary fix for out-of-bounds issues with time_indices (in the
      ! linear tests). The linear input data only contains 2 values, but is
      ! read into a 3-element array. The last value ends up undefined. Catch
      ! the undefined value by looking for a value <1 or >100000 and resetting
      ! it to 1. It's a hack, but it should be fixed properly by ticket #4348.
      time_index = max(time_index,1)
      if( time_index > 10000) time_index=1

      ! Get correct time-entry from current multi-data level
      time_slice = ndata_slice( ( time_index - 1 ) * ( domain_size * vert_levels ) + 1 :  &
                                 ( time_index ) * ( domain_size * vert_levels ) )

      ! We require multi-data fields with vertical levels to be multi-data first
      if ( ndata /= 1 .and. vert_levels /= 1 .and. &
           .not. field_proxy%vspace%is_ndata_first() ) then
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

  ! Use our own mdi indicator
  where (field_data(1:undf) == LARGE_REAL_NEGATIVE) field_data = rmdi

  ! Pass reshaped data array to field object via proxy
  field_proxy%data( 1 : undf ) = field_data( 1 : undf )
  call field_proxy%set_dirty()

  ! Halo exchange necessary to ensure annexed dofs contain safe initial data
  ! This is only needed for horizontally continuous fields
  if (fs_id == W2H) then
    call field_proxy%halo_exchange(depth=1)
  end if

  deallocate( recv_field )
  deallocate( ndata_slice )
  deallocate( time_slice )
  deallocate( field_data )

end subroutine read_field_time_var

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

  ! Create the iter iterator on the state collection
  call iter%initialise(state)

  do
    if ( .not.iter%has_next() ) exit
    fld => iter%next()
    select type(fld)
      type is (field_real32_type)
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

      type is (field_real64_type)
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

    end select
  end do

  nullify(fld)

end subroutine read_state

!>  @brief    Read from a checkpoint into a collection of fields
!>  @details  Iterate over a field collection and read each field
!>            into a collection, if it is enabled for checkpointing
!>
!>  @param[in]  state                 The collection of fields to populate
!>  @param[in]  timestep              The current timestep
!>  @param[in]  checkpoint_stem_name  The checkpoint file stem name
!>  @param[in,optional]  prefix  A prefix to be added to the field name to
!>                               create the XIOS field ID
!>  @param[in,optional]  suffix  A suffix to be added to the field name to
!>                               create the XIOS field ID
!>
subroutine read_checkpoint(state, timestep, checkpoint_stem_name, prefix, suffix)

  implicit none

  type( field_collection_type ), intent(inout) :: state
  integer(i_def),                intent(in)    :: timestep
  character(len=*),              intent(in)    :: checkpoint_stem_name
  character(len=*), optional,    intent(in)    :: prefix
  character(len=*), optional,    intent(in)    :: suffix

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

          call log_event( 'Reading checkpoint file to restart '// &
               xios_field_id, LOG_LEVEL_INFO )
          call fld%read_checkpoint( xios_field_id, &
               trim(ts_fname(checkpoint_stem_name, "",    &
               xios_field_id, timestep,"")) )
       else if ( fld%can_read() ) then
          write(log_scratch_space,'(2A)') &
               "Reading UGRID checkpoint for ", xios_field_id
          call log_event(log_scratch_space, LOG_LEVEL_INFO)
          call fld%read_field( "restart_" // xios_field_id )
       else
          call log_event( 'Reading not set up for  '// xios_field_id, &
               LOG_LEVEL_INFO )
       end if
    type is (field_real64_type)
       if ( fld%can_checkpoint() ) then

          call log_event( 'Reading checkpoint file to restart '// &
              xios_field_id, LOG_LEVEL_INFO )
          call fld%read_checkpoint( xios_field_id, &
               trim(ts_fname(checkpoint_stem_name, "",    &
               xios_field_id, timestep,"")) )
       else if ( fld%can_read() ) then
          write(log_scratch_space,'(2A)') &
               "Reading UGRID checkpoint for ", xios_field_id
          call log_event(log_scratch_space, LOG_LEVEL_INFO)
          call fld%read_field( "restart_" // xios_field_id )
       else
          call log_event( 'Reading not set up for  '// xios_field_id, &
               LOG_LEVEL_INFO )
       end if
    type is (integer_field_type)
       if ( fld%can_checkpoint() ) then
          call log_event( 'Reading checkpoint file to restart '// &
               xios_field_id, LOG_LEVEL_INFO )
          call fld%read_checkpoint( xios_field_id, &
               trim(ts_fname(checkpoint_stem_name, "",    &
               xios_field_id,timestep,"")) )
       else if ( fld%can_read() ) then
          write(log_scratch_space,'(2A)') &
               "Reading UGRID checkpoint for ", xios_field_id
          call log_event(log_scratch_space, LOG_LEVEL_INFO)
          call fld%read_field( "restart_" // xios_field_id )
       else
          call log_event( 'Reading not set up for  '// xios_field_id, &
               LOG_LEVEL_INFO )
       end if
    class default
       call log_event('read_checkpoint:Invalid type of field, not supported',LOG_LEVEL_ERROR)
    end select
  end do

  nullify(fld)

end subroutine read_checkpoint

end module lfric_xios_read_mod
