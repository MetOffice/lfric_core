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
                                         max(spinup_period, 0.0_r_second) ), &
                                         stat=rc )
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

end module driver_time_mod
