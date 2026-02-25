!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>  @brief Module containing a file type for XIOS interface
!>
module lfric_xios_file_mod

  use, intrinsic :: iso_fortran_env, only: real64
  use constants_mod,                 only: i_def, l_def, str_def,        &
                                           str_max_filename, RMDI, IMDI
  use field_mod,                     only: field_type
  use field_parent_mod,              only: field_parent_type
  use field_collection_mod,          only: field_collection_type
  use field_collection_iterator_mod, only: field_collection_iterator_type
  use file_mod,                      only: file_type,      &
                                           file_mode_read, &
                                           file_mode_write
  use lfric_xios_process_output_mod, only: process_output_file
  use lfric_xios_field_mod,          only: lfric_xios_field_type
  use lfric_xios_diag_mod,           only: file_is_tagged
  use log_mod,                       only: log_event, log_level_error, &
                                           log_level_trace, log_level_debug
  use mesh_mod,                      only: mesh_type
  use mod_wait,                      only: init_wait
  use lfric_xios_diag_mod,           only: get_file_name
  use lfric_xios_temporal_mod,       only: temporal_type
  use xios,                          only: xios_file, xios_is_valid_file,    &
                                           xios_add_child, xios_get_handle,  &
                                           xios_set_attr, xios_filegroup,    &
                                           xios_get_file_attr,               &
                                           xios_set_file_attr,               &
                                           xios_fieldgroup, xios_duration,   &
                                           xios_is_valid_fieldgroup,         &
                                           xios_is_defined_fieldgroup_attr,  &
                                           xios_get_fieldgroup_attr,         &
                                           xios_set_fieldgroup_attr,         &
                                           xios_date, xios_get_current_date, &
                                           xios_get_start_date,              &
                                           xios_get_timestep, operator(+),   &
                                           operator(*), operator(<=)

  implicit none

private

integer, parameter :: undef_freq = -999 ! uninitialised frequency value
integer, parameter :: no_freq = 0       ! do not override the frequency
integer, parameter :: undef_file_convention = -999
character(20), parameter :: undef_group = "unset"
character(20), parameter :: diag_main_file = "lfric_diag"
character(20), parameter :: diag_field_group = "diagnostic_fields"
character(20), parameter :: diag_int_field_group = "diag_int_fields"

!> @brief Container for file properties need by XIOS
!>
type, public, extends(file_type) :: lfric_xios_file_type
  private

  !> Due to the pattern of initialising an object within LFRic and then later
  !! using the information in that object to initialise the corresponding object
  !! within the XIOS context we need to keep two versions of the same
  !! information present.
  !>
  !> Model (pre-XIOS) representations
  !> The ID of the XIOS file
  character(str_def)          :: xios_id
  !> The path to the file
  character(str_max_filename) :: path
  !> The I/O mode of open file (read or write)
  integer(i_def)              :: io_mode
  !> How the contents of the file evolve with time
  integer(i_def)              :: operation
  !> The file convention
  integer(i_def)              :: file_convention = undef_file_convention
  !> The file frequency in timesteps
  integer(i_def)              :: freq_ts = undef_freq
  !> @todo field_group is slated for removal, it is a leftover placeholder
  !!       needed to make checkpointing work for lfric_atm/gungho, but once
  !!       they are upgraded to use the "fields_in_file" API this can be removed.
  character(str_def)          :: field_group = undef_group
  !> The XIOS ID of the field group contained within the file
  character(str_def)          :: field_group_id
  !> Flag denoting if the file has been closed
  logical :: is_closed = .false.
  !> Flag denoting if it is a diagnostics file
  logical :: is_diag = .false.
  !> Flag denoting if the always-on sampling mode is selected
  logical :: diag_always_on_sampling = .true.
  !> Will the file be read again from the beginning once the end is reached
  logical :: cyclic = .false.
  !>
  logical :: context_init_read = .true.
  !> Temporal controller for file
  type(temporal_type) :: temporal

  !> XIOS representations
  !> Internal XIOS representation of the file
  type(xios_file)             :: handle
  !> The file frequency as an XIOS duration
  type(xios_duration)         :: frequency
  !> The date of the next operation to be performed on the file by XIOS.
  type(xios_date)             :: next_operation
  !> An array of lfric_xios_field_type objects attached to this file
  type(lfric_xios_field_type), allocatable :: fields(:)


contains
  procedure, public :: file_new
  procedure, public :: file_open
  procedure, public :: file_close
  procedure, public :: register_with_context
  procedure, public :: mode_is_read
  procedure, public :: mode_is_write
  procedure, public :: recv_fields
  procedure, public :: send_fields
  final             :: lfric_xios_file_final

end type lfric_xios_file_type

interface lfric_xios_file_type
  module procedure lfric_xios_file_constructor
end interface

integer(i_def), public, parameter :: OPERATION_ONCE       = 1938
integer(i_def), public, parameter :: OPERATION_TIMESERIES = 3406

integer(i_def), public, parameter :: CONVENTION_CF    =  6542
integer(i_def), public, parameter :: CONVENTION_UGRID =  3698

contains

!> @brief Registration of diagnostics file
!> @param[in] xios_id                  XIOS ID of the file
!> @param[in] freq_ts                  Output frequency for main diagnostics file in timesteps
!> @param[in] is_main                  Is this the main diagnostics file (lfric_diag)?
!> @param[in] always_on_sampling       Is the always-on sampling mode selected?
subroutine register_diagnostics_file(xios_id, freq_ts, is_main, always_on_sampling)
  implicit none
  character(len=*),  intent(in) :: xios_id
  integer(i_def),    intent(in) :: freq_ts
  logical(l_def),    intent(in) :: is_main
  logical(l_def),    intent(in) :: always_on_sampling
  type(xios_duration) :: ts
  type(xios_duration) :: freq
  logical(l_def), save :: first_time = .true.

  ! Create CF-compliant time description
  call xios_set_file_attr(xios_id, time_counter="exclusive", &
                                   time_counter_name="time")
  if (is_main) then
    ! Override XIOS file output frequency
    if (freq_ts == undef_freq) then
      call log_event("uninitialised frequency for file "//trim(xios_id), &
                     log_level_error)
    end if
    if (freq_ts <= 0) then
      ! freq_ts == 0 to be allowed and ignored in a future release
      call log_event("non-positive frequency for file "//trim(xios_id), &
                      log_level_error)
    end if
    if (freq_ts /= no_freq) then
      call xios_get_timestep(ts)
      freq = freq_ts * ts
      call xios_set_file_attr(xios_id, output_freq=freq)
    end if
  end if
  ! configure diagnostic group, but only once
  if (first_time) then
    first_time = .false.
    if (xios_is_valid_fieldgroup(diag_field_group)) then
      call xios_set_fieldgroup_attr(diag_field_group, default_value=real(RMDI, real64))
    end if
    if (xios_is_valid_fieldgroup(diag_int_field_group)) then
      call xios_set_fieldgroup_attr(diag_int_field_group, default_value=real(IMDI, real64))
    end if
    if (.not. always_on_sampling) then
      if (xios_is_valid_fieldgroup(diag_field_group)) then
        call xios_set_fieldgroup_attr(diag_field_group, &
                                      check_if_active=.true.)
      end if
      if (xios_is_valid_fieldgroup(diag_int_field_group)) then
        call xios_set_fieldgroup_attr(diag_int_field_group, &
                                      check_if_active=.true.)
      end if
    end if
  end if
end subroutine register_diagnostics_file

!> Constructs an LFRic-XIOS file type object
!>
!> @param[in] file_name      The name/path of the file
!> @param[in] xios_id        The XIOS ID of the file
!> @param[in] io_mode        Enum denoting whether file is input or output
!> @param[in] freq           The frequency in timesteps that the file is
!!                           read-from/written-to
!> @param[in] operation      Enum denoting the kind of I/O done by the file
!> @param[in] cyclic         Whether the file is cyclic (for time series files)
!> @param[in] field_group_id XIOS ID of the field group contained in the file
!> @param[in] fields_in_file Array of fields contained in the file
!> @param[in] is_diag        Is it a diagnostics file?
!> @param[in] diag_always_on_sampling Is the always-on sampling mode selected?
function lfric_xios_file_constructor( file_name, xios_id, io_mode, freq,      &
                                      operation, cyclic, field_group_id,      &
                                      fields_in_file, is_diag,                &
                                      diag_always_on_sampling,                &
                                      file_convention ) result(self)

  implicit none

  type(lfric_xios_file_type)                :: self
  character(len=*),              intent(in) :: file_name
  character(len=*),              intent(in) :: xios_id
  integer(i_def),                intent(in) :: io_mode
  integer(i_def),      optional, intent(in) :: freq
  integer(i_def),      optional, intent(in) :: operation
  logical,             optional, intent(in) :: cyclic
  character(len=*),    optional, intent(in) :: field_group_id
  type(field_collection_type), &
                       optional, intent(in) :: fields_in_file
  logical(l_def),      optional, intent(in) :: is_diag
  logical(l_def),      optional, intent(in) :: diag_always_on_sampling
  integer(i_def),      optional, intent(in) :: file_convention

  type(field_collection_iterator_type) :: iter
  class(field_parent_type), pointer    :: fld => null()

  integer(i_def) :: field_index

  self%path = file_name

  self%xios_id = xios_id

  self%io_mode = io_mode

  if (present(operation)) self%operation = operation
  if (present(file_convention)) self%file_convention = file_convention
  if (present(cyclic) .and. self%io_mode == FILE_MODE_READ) then
    self%cyclic = cyclic
  end if

  if (present(freq)) then
    if (freq < 0) then ! we are going to allow freq = 0 (= no_freq)
      call log_event( "XIOS files cannot have negative frequency", &
                      log_level_error )
    end if
    self%freq_ts = freq
  end if

  if (present(field_group_id)) self%field_group = field_group_id

  ! Set up XIOS fields representing attached field collection
  if (present(fields_in_file)) then
    if (.not. present(operation)) then
      call log_event( "XIOS file containing model fields must have a defined operation", &
                      log_level_error )
    end if
    allocate(self%fields(fields_in_file%get_length()))
    self%field_group_id = trim(self%xios_id)//"_fields"
    call iter%initialise(fields_in_file)
    do field_index = 1, fields_in_file%get_length()
      fld => iter%next()
      self%fields(field_index) = lfric_xios_field_type(fld, fieldgroup_id=self%field_group_id)
    end do
  end if

  if (present(is_diag)) self%is_diag = is_diag
  if (present(diag_always_on_sampling)) &
    self%diag_always_on_sampling = diag_always_on_sampling
  return

end function lfric_xios_file_constructor

!> @brief  Initialises an LFRic-XIOS file object
!> @param[in] file_name Path to the file
subroutine file_new(self, file_name)

  implicit none

  class(lfric_xios_file_type), intent(inout) :: self
  character(len=*),            intent(in)    :: file_name

  self%path = file_name

end subroutine file_new

!> @brief  Opens the LFRic definition for file reading
!>         after it has been opened by XIOS.
!> @param[in] file_name Path to the file
!> @param[in] file_mode Action on file.
subroutine file_open(self, file_name, file_mode)

  implicit none

  class(lfric_xios_file_type), intent(inout) :: self
  character(len=*),            intent(in)    :: file_name
  integer(i_def),    optional, intent(in)    :: file_mode

  ! TODO Write orography

end subroutine file_open

!> @brief Performs final routines before closure of XIOS file
subroutine file_close(self)

  implicit none

  class(lfric_xios_file_type), intent(inout) :: self

  if (self%is_closed) return

  if ( self%io_mode == FILE_MODE_WRITE ) then
    call log_event( "Waiting for XIOS to close file ["//trim(self%path)//".nc]", &
                    log_level_debug )
    call init_wait()
    call process_output_file(trim(self%path)//".nc")
  end if

  self%is_closed = .true.

end subroutine file_close

!> @brief Registers a file with an XIOS context
subroutine register_with_context(self)

  implicit none

  class(lfric_xios_file_type),  intent(inout) :: self

  type(xios_filegroup)   :: file_definition
  type(xios_fieldgroup)  :: field_group_hdl, file_fields
  type(xios_duration)    :: timestep_duration
  type(xios_date)        :: start_date

  integer(i_def) :: i, record_offset

  call log_event( "Registering XIOS file ["//trim(self%xios_id)//"]", &
                  log_level_trace )

  if (self%is_diag) then
    if (xios_is_valid_file(self%xios_id)) then
      self%path = get_file_name(self%xios_id) ! if there is a name specified, we use it
      call register_diagnostics_file(self%xios_id, self%freq_ts, &
        self%xios_id == diag_main_file, self%diag_always_on_sampling)
      return
    end if
  end if

  ! Register or get handle of file from XIOS
  if (xios_is_valid_file(trim(self%xios_id))) then
    call xios_get_handle( trim(self%xios_id), self%handle )
  else
    call xios_get_handle("file_definition", file_definition)
    call xios_add_child(file_definition, self%handle, trim(self%xios_id))
  end if

  ! Set file path
  call xios_set_attr( self%handle, name=trim(adjustl(self%path)) )

  ! Set I/O mode (no need for case defalt as we will fall back to XIOS's
  ! default behaviour)
  call xios_set_attr( self%handle, type="one_file" )
  select case(self%io_mode)
  case (FILE_MODE_READ)
    call xios_set_attr( self%handle, mode="read" )
  case (FILE_MODE_WRITE)
    call xios_set_attr( self%handle, mode="write" )
  end select

  ! Set XIOS file convention
  select case(self%file_convention)
    case (CONVENTION_CF)
      call xios_set_attr( self%handle, convention="CF" )
    case (CONVENTION_UGRID)
      call xios_set_attr( self%handle, convention="UGRID" )
  end select

  ! Create CF-compliant time description
  if (self%operation == OPERATION_TIMESERIES) then
    call xios_set_attr( self%handle, time_counter="exclusive", &
                                     time_counter_name="time" )
  else
    call xios_set_attr( self%handle, time_counter="none")
  end if

  ! Set XIOS duration object second value equal to file output frequency
  call xios_get_timestep(timestep_duration)
  if (.not. self%freq_ts == undef_freq) then
    self%frequency = self%freq_ts * timestep_duration
    call xios_set_attr(self%handle, output_freq=self%frequency)
  else
    ! If frequency is uninitialised, get it from XIOS
    call xios_get_file_attr(self%xios_id, output_freq=self%frequency)
  end if

  ! Set the date of the first operation
  call xios_get_start_date(start_date)
  self%next_operation = start_date + self%frequency

  ! Set up fields in file
  if (allocated(self%fields)) then
    call xios_add_child(self%handle, file_fields, self%field_group_id)

    ! Set the temporal operation for fields in the file
    select case(self%operation)
    case(OPERATION_ONCE)
      call xios_set_attr(file_fields, operation="once")
    case(OPERATION_TIMESERIES)
      call xios_set_attr(file_fields, operation="instant")
    end select

    ! Iterate over field collection and register fields
    do i = 1, size(self%fields)
      call self%fields(i)%register(field_read_access=self%mode_is_read())
    end do

    ! Set up time axis if needed
    if ( self%io_mode == FILE_MODE_READ .and. &
         self%operation == OPERATION_TIMESERIES ) then

      ! Initialise file temporal control. Get the record offset from
      ! the temporal object initialiser which will tell XIOS which time entry
      ! to start reading from
      call xios_set_attr(self%handle, cyclic=self%cyclic)
      call self%temporal%initialise( self%path, self%fields, self%frequency, &
                                     self%cyclic, record_offset )
      call xios_set_attr(self%handle, record_offset=record_offset)
    end if

    ! Enable field collection
    call xios_set_attr(file_fields, enabled=.true.)

  end if

  ! LEGACY
  ! If there is an associated field group, enable it
  if ( .not. trim(self%field_group) == undef_group ) then
    call xios_get_handle( trim(self%field_group), field_group_hdl )
    call xios_set_attr( field_group_hdl, enabled=.true. )
  end if

  ! Enable file
  call xios_set_attr( self%handle, enabled=.true. )

end subroutine register_with_context

!> Returns true if the file is in "read" mode
function mode_is_read(self) result(file_mode_is_read)

  implicit none

  class(lfric_xios_file_type), intent(inout) :: self
  logical :: file_mode_is_read

  file_mode_is_read = (self%io_mode == FILE_MODE_READ)

end function mode_is_read

!> Returns true if the file is in "write" mode
function mode_is_write(self) result(file_mode_is_write)

  implicit none

  class(lfric_xios_file_type), intent(inout) :: self
  logical :: file_mode_is_write

  file_mode_is_write = (self%io_mode == FILE_MODE_WRITE)

end function mode_is_write

!> @brief Reads all registered field entries in a file
subroutine recv_fields(self)

  implicit none

  class(lfric_xios_file_type),  intent(inout) :: self

  type(xios_date)   :: model_date
  integer(i_def)    :: i

  ! If file is closed or if there are no field descriptions in the file then
  ! don't do any operations
  if (.not. allocated(self%fields) .or. self%is_closed) return

  if (self%io_mode /= FILE_MODE_READ) then
    call log_event( "Can't perform read operations on file ["//     &
                    trim(self%xios_id)//"]: file not in read mode", &
                    log_level_error)
  end if

  call xios_get_current_date(model_date)
  if (self%context_init_read .or. self%next_operation <= model_date) then

    ! Iterate over field collection and read fields
    do i = 1, size(self%fields)
      call self%fields(i)%recv()
    end do

    ! If file should only be operated on once, close it, else set the time for
    ! the next operation
    if (self%operation == OPERATION_ONCE) then
      call self%file_close()
    else if (.not. self%context_init_read) then
      self%next_operation = self%next_operation + self%frequency
    end if

    ! Advance the time axis if present
    if ( (self%io_mode == FILE_MODE_READ) .and. &
         (self%operation == OPERATION_TIMESERIES) ) then
      if (.not. self%temporal%advance()) call self%file_close()
    end if

    self%context_init_read = .false.

  end if

end subroutine recv_fields

!> @brief Writes all fields registered with a file
subroutine send_fields(self)

  implicit none

  class(lfric_xios_file_type),  intent(inout) :: self

  type(xios_date)   :: model_date
  integer(i_def)    :: i

  ! If file is closed or if there are no field descriptions in the file then
  ! don't do any operations
  if (.not. allocated(self%fields) .or. self%is_closed) return

  ! Do not try to write if file is not in write mode
  if (self%io_mode /= FILE_MODE_WRITE) then
    call log_event( "Cannot perform write operations on file ["//     &
                    trim(self%xios_id)//"]: file not in write mode", &
                    log_level_error)
  end if

  call xios_get_current_date(model_date)
  if (self%next_operation <= model_date) then

    call log_event( "Writing fields to file ["//trim(self%xios_id)//"]", &
                    log_level_trace)

    ! Iterate over field collection and write fields
    do i = 1, size(self%fields)
      call self%fields(i)%send()
    end do

    ! If file should only be operated on once, close it, else set the time for
    ! the next operation
    if (self%operation == OPERATION_ONCE) then
      call self%file_close()
    else
      self%next_operation = self%next_operation + self%frequency
    end if

  end if

end subroutine send_fields

!> Finaliser for lfric_xios_file_type object
subroutine lfric_xios_file_final(self)

  implicit none

  type(lfric_xios_file_type), intent(inout) :: self

  if (allocated(self%fields)) deallocate(self%fields)

end subroutine lfric_xios_file_final

end module lfric_xios_file_mod
