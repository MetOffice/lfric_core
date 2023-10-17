!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Module controlling the initialisation and finalisation of time related
!> functionality for a model using the LFRic infrastructure
module driver_time_mod

  use calendar_mod,            only: calendar_type
  use clock_mod,               only: clock_type
  use constants_mod,           only: i_def
  use log_mod,                 only: log_event, log_level_error
  use model_clock_mod,         only: model_clock_type
  use step_calendar_mod,       only: step_calendar_type
  use time_config_mod,         only: timestep_end, timestep_start, &
                                     calendar_origin, calendar_start
  use timestepping_config_mod, only: dt, spinup_period

  implicit none

  private
  public :: init_time, get_calendar

  ! Model calendar
  class(calendar_type), allocatable, target :: calendar

contains

  !> Initialise model clock and calendar from configuration
  !>
  !> @param[out] clock  The model clock
  subroutine init_time(clock)

    implicit none

    type(model_clock_type), allocatable, intent(out) :: clock

    integer(i_def) :: rc

    ! Choice of calendar here
    if (.not. allocated(calendar)) then
      allocate( calendar, source=step_calendar_type(calendar_origin,           &
                                                    calendar_start), stat=rc )
      if (rc /= 0) then
        call log_event( "Unable to allocate calendar", log_level_error )
      end if
    end if

    ! Create the model's clock
    allocate( clock, source=model_clock_type( &
                                      calendar%parse_instance(timestep_start), &
                                      calendar%parse_instance(timestep_end),   &
                                      dt, spinup_period ), stat=rc )
    if (rc /= 0) then
      call log_event( "Unable to allocate model clock", log_level_error )
    end if

  end subroutine init_time

  !> @brief  Returns the model calendar.
  function get_calendar() result(calendar_ptr)

    implicit none

    class(calendar_type), pointer :: calendar_ptr

    calendar_ptr => calendar

  end function get_calendar

end module driver_time_mod
