!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing the JEDI State emulator class.
!>
!> @details This module holds a JEDI State emulator class that includes only
!>          the functionality required by the LFRic-JEDI model interface on the
!>          LFRic-API. This includes i) read/write to a defined set of fields,
!>          ii) ability to interoperate between LFRic fields and JEDI fields
!>          (Atlas fields here), and iii) storage of model_data instance to be
!>          used for IO and model time stepping.
module jedi_state_mod

  use, intrinsic :: iso_fortran_env, only : real64
  use atlas_field_emulator_mod,      only : atlas_field_emulator_type
  use atlas_field_interface_mod,     only : atlas_field_interface_type
  use jedi_lfric_datetime_mod,       only : jedi_datetime_type
  use jedi_lfric_duration_mod,       only : jedi_duration_type
  use jedi_geometry_mod,             only : jedi_geometry_type
  use driver_model_data_mod,         only : model_data_type
  use jedi_state_config_mod,         only : jedi_state_config_type
  use jedi_lfric_field_meta_mod,     only : jedi_lfric_field_meta_type
  use field_collection_mod,          only : field_collection_type
  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_INFO,     &
                                            LOG_LEVEL_ERROR
  use jedi_lfric_fake_nl_driver_mod, only : mesh, twod_mesh
  use constants_mod,                 only : i_def, l_def, str_def

  use jedi_lfric_fake_nl_driver_mod,   only : model_clock

  implicit none

  private

type, public :: jedi_state_type
  private

  !> These fields emulate a set of Atlas fields are and used purely for testing
  type ( atlas_field_emulator_type ), allocatable :: fields(:)

  !> An object that stores the field meta-data associated with the fields
  type( jedi_lfric_field_meta_type )              :: field_meta_data

  !> Interface field linking the Atlas emulator fields and LFRic fields in the
  !> model data (to do field copies)
  type( atlas_field_interface_type ), allocatable :: fields_to_model_data(:)

  !> Model data that stores the model_data fields to propagate
  type( model_data_type ),  public                :: model_data

  !> Field collection to perform IO
  type( field_collection_type ), public           :: io_collection

  ! date: '2018-04-14T21:00:00Z'
  ! mpas define the formatting for this as:
  ! dateTimeString = '$Y-$M-$D_$h:$m:$s'
  type( jedi_datetime_type )                      :: state_time

  !> The jedi_geometry object
  type( jedi_geometry_type ), pointer             :: geometry => null()

contains

  !> Jedi state initialiser.
  procedure :: initialise => state_initialiser_read
  procedure :: state_initialiser

  !> Initialise the model_data
  procedure, private :: initialise_model_data

  !> Setup the atlas_field_interface_type that enables copying between Atlas
  !> field emulators and the LFRic fields in the model_data
  procedure, private :: setup_interface_to_model_data

  !> Setup the atlas_field_interface_type that enables copying between Atlas
  !> field emulators and the fields in a LFRic field collection
  procedure, private :: setup_interface_to_field_collection

  !> Copy the data in the internal Atlas field emulators from the LFRic fields
  !> stored in the a field_collection
  procedure, public :: from_lfric_field_collection

  !> Copy the data in the internal Atlas field emulators to the LFRic fields
  !> stored in the a field_collection
  procedure, public :: to_lfric_field_collection

  !> Return the curent time
  procedure, public :: valid_time

  !> Read model fields from file into fields
  procedure, public :: read_file

  !> @todo Write model fields to file from the fields
  !> procedure, public :: write_file

  !> Create the model_data
  procedure, public :: create_model_data

  !> Update the curent time
  procedure, public :: update_time

  !> Print field
  procedure, public :: print_field

  !> Copy the data in the LFRic fields stored in the model_data to the internal
  !> Atlas field emulators
  procedure, public :: from_model_data

  !> Finalizer
  final             :: jedi_state_destructor

end type jedi_state_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------------
! Methods required by the JEDI (OOPS) model interface
!-------------------------------------------------------------------------------

!> @brief    Initialiser via read for jedi_state_type
!>
!> @param [in] geometry The geometry object required to construct the state
!> @param [in] config   The configuration object including the required
!>                      information to construct a state and read a file to
!>                      initialise the fields
subroutine state_initialiser_read( self, program_name, geometry, config )

  implicit none

  class( jedi_state_type ),           intent(inout) :: self
  character(len=*),                   intent(in)    :: program_name
  type( jedi_geometry_type ), target, intent(in)    :: geometry
  type( jedi_state_config_type ),     intent(inout) :: config

  call self%state_initialiser( geometry, config )

  ! Tick out of initialisation state
  if ( .not. model_clock%tick() ) then
    call log_event( 'The LFRic model_clock has stopped.', LOG_LEVEL_ERROR )
  end if

  ! Initialise the Atlas field emulators via the model_data or the
  ! io_collection
  if ( .not. config%use_pseudo_model ) then
    ! This calls the models initialise method that populates model_data and
    ! does a copy from model_data to the fields
    call self%initialise_model_data()
  else
    ! We are not running the non-linear model so read the file directly and
    ! do a copy from io_collection to the fields
    call self%read_file( config%state_time, config%read_file_prefix )
  end if

end subroutine state_initialiser_read


!> @brief    Initialiser for jedi_state_type
!>
!> @param [in] geometry The geometry object required to construct the state
!> @param [in] config   A configuration object including the required
!>                      information to construct a state
subroutine state_initialiser( self, geometry, config )

  use fs_continuity_mod,     only : W3, Wtheta

  implicit none

  class( jedi_state_type ),        intent(inout) :: self
  type( jedi_geometry_type ), target, intent(in) :: geometry
  type( jedi_state_config_type ),  intent(inout) :: config

  ! Local
  integer(i_def) :: n_horizontal
  integer(i_def) :: n_levels
  integer(i_def) :: n_layers
  integer(i_def) :: ivar
  integer(i_def) :: n_variables
  integer(i_def) :: fs_id
  logical(l_def) :: twod_field

  ! Setup
  self%field_meta_data = config%field_meta_data
  self%state_time = config%state_time
  self%geometry => geometry
  n_variables = self%field_meta_data%get_n_variables()

  allocate( self%fields(n_variables) )

  n_horizontal=geometry%get_n_horizontal()
  n_layers=geometry%get_n_layers()

  do ivar=1,n_variables

    fs_id      = self%field_meta_data%get_variable_function_space(ivar)
    twod_field = self%field_meta_data%get_variable_is_2d(ivar)

    select case (fs_id)
    case (W3)
      if (twod_field) then
        n_levels = 1
      else
        n_levels = n_layers
      end if

    case (Wtheta)
      n_levels = n_layers + 1

    case default
      log_scratch_space = 'The requested LFRic function space is not supported.'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )

    end select

    call self%fields(ivar)%initialise(        &
                                n_levels,     &
                                n_horizontal, &
                                self%field_meta_data%get_variable_name(ivar) )
  end do

  ! Setup the io_collection (empty ready for read/write operations)
  call self%io_collection%initialise(name = 'io_collection', table_len=100)

  ! If running the model, create model data and link to fields .
  if ( .not. config%use_pseudo_model ) then
    call self%create_model_data()
  end if

end subroutine state_initialiser

!> @brief    A method to initialise the model_data and copy into the Atlas
!>           field emulators
!>
subroutine initialise_model_data( self )

  use jedi_lfric_fake_nl_init_mod, only : initialise_da_model_data

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Read the state into the LFRic model data
  call initialise_da_model_data( self%model_data )

  ! Copy model_data fields to the Atlas field emulators
  call self%from_model_data()

end subroutine initialise_model_data

!> @brief    Returns the current time of the increment
!>
!> @return   time   The current time of the state
function valid_time( self ) result(time)

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  type( jedi_datetime_type )              :: time

  time = self%state_time

end function valid_time

!> @brief    A method to update the internal Atlas field emulators
!>
!> @param [in] read_time   The time to be read
!> @param [in] file_prefix Character array that specifies the file to read from
subroutine read_file( self, read_time, file_prefix )

  use jedi_lfric_tests_io_update_mod, only : update_io_field_collection
  use lfric_xios_read_mod,  only : read_state

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  type( jedi_datetime_type ),  intent(in) :: read_time
  character(len=*),            intent(in) :: file_prefix

  ! Local
  character( len=str_def ), allocatable :: variable_names(:)

  ! Set the clock to the desired read time
  call set_clock(self, read_time)

  ! Ensure the io_collection contains the variables defined in the list
  ! stored in the type
  call update_io_field_collection(self%io_collection, mesh, twod_mesh, &
                                  self%field_meta_data)

  ! Read the state into the io_collection
  call read_state( self%io_collection, prefix=file_prefix )

  ! Copy model_data fields to the Atlas field emulators
  call self%field_meta_data%get_variable_names(variable_names)

  call self%from_lfric_field_collection(variable_names, self%io_collection)

end subroutine read_file

!> Write fields stored in this state from file
!subroutine write_file()
! TBD ...
!end subroutine write_file

!------------------------------------------------------------------------------
! Local methods to support LFRic-JEDI implementation
!------------------------------------------------------------------------------

!> @brief    Create an instance of the model_data for jedi_state_type
!>
subroutine create_model_data( self )

  use jedi_lfric_fake_nl_init_mod,   only : create_da_model_data
  use jedi_lfric_fake_nl_driver_mod, only : mesh, twod_mesh

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Create model data and then link the Atlas fields to them
  call create_da_model_data( mesh, twod_mesh, self%model_data )
  call self%setup_interface_to_model_data()

end subroutine create_model_data

!> @brief    Setup fields_to_model_data variable that enables copying between
!>           Atlas field emulators and the LFRic fields in the model_data
!>
subroutine setup_interface_to_model_data( self )

  use jedi_lfric_tests_utils_mod, only: get_model_field
  use field_mod,                  only: field_type

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Local
  integer(i_def)            :: ivar
  type(field_type), pointer :: lfric_field_ptr
  real(real64),     pointer :: atlas_data_ptr(:,:)
  integer(i_def),   pointer :: horizontal_map_ptr(:)
  integer(i_def)            :: n_variables

  type( field_collection_type ), pointer :: depository

  nullify(depository)
  depository => self%model_data%get_field_collection("depository")

  n_variables = self%field_meta_data%get_n_variables()

  ! Allocate space for the interface fields
  if ( allocated( self%fields_to_model_data ) ) then
    deallocate( self%fields_to_model_data )
  endif
  allocate( self%fields_to_model_data( n_variables ) )

  ! Link the Atlas emulator fields with lfric fields
  call self%geometry%get_horizontal_map(horizontal_map_ptr)
  do ivar=1, n_variables

    ! Get the required data
    call get_model_field( self%field_meta_data%get_variable_name(ivar), &
                          depository, lfric_field_ptr )

    atlas_data_ptr => self%fields(ivar)%get_data()

    call self%fields_to_model_data(ivar)%initialise( atlas_data_ptr,     &
                                                     horizontal_map_ptr, &
                                                     lfric_field_ptr )

  end do

end subroutine setup_interface_to_model_data

!> @brief    Copy from model_data to the Atlas field emulators
!>
subroutine from_model_data( self )

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Local
  integer(i_def) :: ivar

  !> @todo Will need some sort of transform for winds and
  !>       possibly other higher order elements.
  !>
  !>       call transform_winds(model_data)

  ! copy to the Atlas emulator fields
  do ivar = 1, size(self%fields_to_model_data)
    call self%fields_to_model_data(ivar)%copy_from_lfric()
  end do

end subroutine from_model_data

!> @brief    Setup Atlas-LFRic interface_fields that enables copying
!>           between Atlas field emulators and the LFRic fields in
!>           io_collection
!>
!> @param [inout] interface_fields  The Atlas-LFRic interafce object
!> @param [in]    variable_names    The name of the fields to setup
!> @param [inout] field_collection  The field collection to link the Atlas fields to
subroutine setup_interface_to_field_collection( self, &
                                                interface_fields, &
                                                variable_names, &
                                                field_collection )

  use jedi_lfric_tests_utils_mod, only: get_model_field
  use field_mod,                  only: field_type

  implicit none

  class( jedi_state_type ),           intent(inout) :: self
  type( atlas_field_interface_type ), intent(inout) :: interface_fields(:)
  character( len=str_def ),              intent(in) :: variable_names(:)
  type( field_collection_type ),      intent(inout) :: field_collection

  ! Local
  integer(i_def)            :: ivar
  type(field_type), pointer :: lfric_field_ptr
  real(real64),     pointer :: atlas_data_ptr(:,:)
  integer(i_def),   pointer :: horizontal_map_ptr(:)
  integer(i_def)            :: n_variables
  logical                   :: all_variables_exists

  ! Check that the state contains all the required fields
  all_variables_exists = &
        self%field_meta_data%check_variables_exist( variable_names )
  if ( .not. all_variables_exists ) then
    log_scratch_space = &
        'The Atlas state does not conatin all the required fields.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

  ! Link the Atlas emulator fields with lfric fields
  call self%geometry%get_horizontal_map( horizontal_map_ptr )
  n_variables=size( variable_names )
  do ivar = 1, n_variables
    ! Get the required data and setup field interface
    call get_model_field( variable_names(ivar), &
                          field_collection, lfric_field_ptr )
    atlas_data_ptr => self%fields(ivar)%get_data()
    call interface_fields(ivar)%initialise( atlas_data_ptr,     &
                                            horizontal_map_ptr, &
                                            lfric_field_ptr )
  end do

end subroutine setup_interface_to_field_collection

!> @brief    Copy from a field_collection to the Atlas field emulators
!>
!>
!> @param [in]    variable_names    The name of the fields to setup
!> @param [inout] field_collection  The field collection to copy from
subroutine from_lfric_field_collection( self, &
                                        variable_names, &
                                        field_collection )

  implicit none

  class( jedi_state_type ),      intent(inout) :: self
  character( len=str_def ),         intent(in) :: variable_names(:)
  type( field_collection_type ), intent(inout) :: field_collection

  ! Local
  integer(i_def) :: ivar
  integer(i_def) :: n_variables
  type( atlas_field_interface_type ), allocatable :: interface_fields(:)

  ! Allocate space for the interface fields
  n_variables = size(variable_names)
  allocate( interface_fields(n_variables) )

  ! Link internal Atlas field emulators to LFRic fields
  call self%setup_interface_to_field_collection( interface_fields, &
                                                 variable_names, &
                                                 field_collection )

  ! Copy the LFRic fields in the field_collection to the Atlas field emulators
  do ivar = 1, n_variables
    call interface_fields(ivar)%copy_from_lfric()
  end do

end subroutine from_lfric_field_collection

!> @brief    Copy to a field_collection from the Atlas field emulators
!>
!> @param [in]    variable_names    The name of the fields to setup
!> @param [inout] field_collection  The field collection to copy to
subroutine to_lfric_field_collection( self, variable_names, field_collection )

  implicit none

  class( jedi_state_type ),      intent(inout) :: self
  character( len=str_def ),         intent(in) :: variable_names(:)
  type( field_collection_type ), intent(inout) :: field_collection

  ! Local
  integer(i_def) :: ivar
  integer(i_def) :: n_variables
  type( atlas_field_interface_type ), allocatable :: interface_fields(:)

  ! Allocate space for the interface fields
  n_variables = size( variable_names )
  allocate( interface_fields(n_variables) )

  ! Link internal Atlas field emulators to LFRic fields
  call self%setup_interface_to_field_collection( interface_fields, &
                                                 variable_names, &
                                                 field_collection )

  ! Copy the LFRic fields in the field_collection to the Atlas field emulators
  do ivar = 1, n_variables
    call interface_fields(ivar)%copy_to_lfric()
  end do

end subroutine to_lfric_field_collection

!> @brief    Update the state time by a single time-step
!>
!> @param [in] time_step  Update the state time by the time_step
!>                        duration
subroutine update_time( self, time_step )

  implicit none

  class( jedi_state_type ), intent(inout) :: self
  type( jedi_duration_type ),  intent(in) :: time_step

  self%state_time = self%state_time + time_step

end subroutine update_time

!> @brief    Set the LFRic clock to the time specified by the input time
!>
!> @param [in] new_time  The time to be used to update the LFRic clock
subroutine set_clock( self, new_time )

  use timestepping_config_mod,       only : dt

  implicit none

  class( jedi_state_type ),   intent(inout) :: self
  type( jedi_datetime_type ), intent(in)    :: new_time

  type( jedi_duration_type ) :: time_difference
  type( jedi_duration_type ) :: time_step

  logical(l_def) :: clock_stopped

  call time_step%init( int( dt, kind=i_def ) )

  time_difference = new_time - self%state_time

  if ( time_difference == 0_i_def ) then
    write ( log_scratch_space, '(A)' ) &
      "New time is the same as the current state time"
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  else if ( time_difference < 0_i_def ) then
    write ( log_scratch_space, '(A)' ) &
      "The xios clock can not go backwards."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! Tick the clock to the required time but first check that the clock is
  ! running in the case its not being ticked (i.e. new_time==self%state_time).
  ! clock_stopped=.not.clock%is_running() didnt work when I tried it - set to
  ! false initially for now.

  clock_stopped = .false.

  do while ( new_time%is_ahead( self%state_time ) )
    clock_stopped = .not. model_clock%tick()
    self%state_time = self%state_time + time_step
  end do

  ! Check the clock is still running
  if ( clock_stopped ) then
    write ( log_scratch_space, '(A)' ) &
      "State::set_clock::The LFRic clock has stopped."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

end subroutine set_clock

!> @brief    The jedi_state_type finalizer
!>
subroutine jedi_state_destructor( self )

  implicit none

  type( jedi_state_type ), intent(inout) :: self

  self%geometry => null()
  if ( allocated(self%fields ) ) then
    deallocate(self%fields )
  end if

end subroutine jedi_state_destructor

!! For testing

!> Print the 1st point in each field
subroutine print_field( self )

  implicit none

  class( jedi_state_type ), intent(inout) :: self

  ! Local
  real(real64), pointer :: atlas_data_ptr(:,:)
  integer(i_def)        :: ivar
  character(str_def)    :: iso_datetime

  ! Printing data
  call log_event( "State print ----", LOG_LEVEL_INFO )
  call self%state_time%to_string( iso_datetime )
  write ( log_scratch_space, '(2A)' ) 'Time: ', iso_datetime
  call log_event( log_scratch_space, LOG_LEVEL_INFO )
  do ivar = 1, self%field_meta_data%get_n_variables()
    atlas_data_ptr => self%fields(ivar)%get_data()
    write ( log_scratch_space, '(2A,F20.10)' ) &
      trim(self%field_meta_data%get_variable_name(ivar)), &
      ", atlas_data_ptr(1,1) = ", atlas_data_ptr(1,1)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  end do

end subroutine print_field

end module jedi_state_mod
