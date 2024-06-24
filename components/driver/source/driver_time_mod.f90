!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Module controlling the initialisation and finalisation of time related
!> functionality for a modelDB state object.
!-----------------------------------------------------------------------------
module driver_time_mod

  use constants_mod,      only: str_def, i_def, r_second, i_timestep
  use driver_modeldb_mod, only: modeldb_type
  use log_mod,            only: log_event, LOG_LEVEL_ERROR
  use model_clock_mod,    only: model_clock_type
  use namelist_mod,       only: namelist_type
  use step_calendar_mod,  only: step_calendar_type

  implicit none

  private
  public :: init_time, final_time

  interface init_time
    module procedure init_time
    module procedure init_time_deprecated
  end interface init_time

  interface final_time
    module procedure final_time
    module procedure final_time_deprecated
  end interface final_time

contains


  !> @brief Initialise model clock and calendar for a model state
  !>
  !> @param[out] modeldb Model state object
  !=================================================================
  subroutine init_time( modeldb )

    implicit none

    class(modeldb_type), intent(inout) :: modeldb

    ! Locals
    !--------
    integer(i_def) :: rc

    integer(i_timestep) :: first
    integer(i_timestep) :: last

    type(namelist_type), pointer :: time_nml
    type(namelist_type), pointer :: timestepping_nml

    character(str_def) :: timestep_start
    character(str_def) :: timestep_end
    character(str_def) :: calendar_origin
    character(str_def) :: calendar_start

    real(r_second) :: timestep_length
    real(r_second) :: spinup_period

    ! -------------------------------
    ! Extract namelist variables
    ! -------------------------------
    time_nml         => modeldb%configuration%get_namelist('time')
    timestepping_nml => modeldb%configuration%get_namelist('timestepping')

    call time_nml%get_value( 'timestep_start',  timestep_start )
    call time_nml%get_value( 'timestep_end',    timestep_end )
    call time_nml%get_value( 'calendar_origin', calendar_origin )
    call time_nml%get_value( 'calendar_start',  calendar_start )

    call timestepping_nml%get_value( 'dt', timestep_length )
    call timestepping_nml%get_value( 'spinup_period', spinup_period )

    nullify( time_nml, timestepping_nml )

    ! Instantiate the calendar
    !---------------------------------
    if ( allocated(modeldb%calendar) ) deallocate (modeldb%calendar)
    allocate( modeldb%calendar,                             &
              source = step_calendar_type( calendar_origin, &
                                           calendar_start ), stat=rc )

    if (rc /= 0) then
      call log_event( "Unable to allocate calendar", LOG_LEVEL_ERROR )
    end if

    ! Instantiate the model clock
    !---------------------------------
    first = modeldb%calendar%parse_instance(timestep_start)
    last  = modeldb%calendar%parse_instance(timestep_end)

    if ( allocated(modeldb%clock) ) deallocate (modeldb%clock)
    allocate( modeldb%clock,                              &
              source = model_clock_type( first, last,     &
                                         timestep_length, &
                                         spinup_period ), stat=rc )
    if (rc /= 0) then
      call log_event( "Unable to allocate model clock", LOG_LEVEL_ERROR )
    end if

  end subroutine init_time


  !> @brief Finalise the clock and calendar of a model state
  !>
  !> @param[in out] modeldb  Model state object
  !=================================================================
  subroutine final_time( modeldb )

    implicit none

    class(modeldb_type), intent(inout) :: modeldb

    if ( allocated(modeldb%clock) )    deallocate(modeldb%clock)
    if ( allocated(modeldb%calendar) ) deallocate(modeldb%calendar)

  end subroutine final_time


!====================================================================
! DEPRECATED ROUTINES
!====================================================================

  !> Initialise model clock and calendar from configuration
  !>
  !> @param[out] clock    The model clock
  !> @param[out] calendar The model calendar
  subroutine init_time_deprecated(clock, calendar)

    use calendar_mod,       only: calendar_type
    use constants_mod,      only: str_def, i_def, r_second, i_timestep
    use log_mod,            only: log_event, log_scratch_space, &
                                  LOG_LEVEL_ERROR,              &
                                  LOG_LEVEL_WARNING
    use model_clock_mod,    only: model_clock_type
    use step_calendar_mod,  only: step_calendar_type

    use time_config_mod,         only: timestep_end, timestep_start, &
                                       calendar_origin, calendar_start
    use timestepping_config_mod, only: dt, spinup_period

    implicit none

    type(model_clock_type),   allocatable, intent(out) :: clock
    class(calendar_type),     allocatable, intent(out) :: calendar

    integer(i_def) :: rc

    write(log_scratch_space,'(A)') 'Using DEPRECATED version of init_time'
    call log_event(log_scratch_space, LOG_LEVEL_WARNING)
    write(log_scratch_space,'(A)') 'Removal at milestone: DEIMOS (Sep 2024)'
    call log_event(log_scratch_space, LOG_LEVEL_WARNING)

    ! Choice of calendar here
    if (.not. allocated(calendar)) then
      allocate( calendar,                                     &
                source = step_calendar_type( calendar_origin, &
                                             calendar_start), stat=rc )
      if (rc /= 0) then
        call log_event( "Unable to allocate calendar", LOG_LEVEL_ERROR )
      end if
    end if

    ! Create the model's clock
    allocate( clock, source=model_clock_type(                                  &
                                      calendar%parse_instance(timestep_start), &
                                      calendar%parse_instance(timestep_end),   &
                                      dt, spinup_period ), stat=rc )
    if (rc /= 0) then
      call log_event( "Unable to allocate model clock", LOG_LEVEL_ERROR )
    end if

  end subroutine init_time_deprecated


  !> Finalise model clock and calendar
  !>
  !> @param[out] clock    The model clock
  !> @param[out] calendar The model calendar
  subroutine final_time_deprecated(clock, calendar)

    use calendar_mod,    only: calendar_type
    use model_clock_mod, only: model_clock_type
    use log_mod,         only: log_event, log_scratch_space, &
                               LOG_LEVEL_WARNING

    implicit none

    type(model_clock_type), allocatable, intent(inout) :: clock
    class(calendar_type),   allocatable, intent(inout) :: calendar

    write(log_scratch_space,'(A)') 'Using DEPRECATED version of final_time'
    call log_event(log_scratch_space, LOG_LEVEL_WARNING)
    write(log_scratch_space,'(A)') 'Removal at milestone: DEIMOS (Sep 2024)'
    call log_event(log_scratch_space, LOG_LEVEL_WARNING)

    if (allocated(clock))    deallocate(clock)
    if (allocated(calendar)) deallocate(calendar)

  end subroutine final_time_deprecated

end module driver_time_mod
