!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Wrap the XIOS context in an object for easier management and cleaner coding.
!>
module lfric_xios_context_mod

  use clock_mod,            only : clock_type
  use io_context_mod,       only : io_context_type, &
                                   io_context_initialiser_type
  use step_calendar_mod,    only : step_calendar_type
  use constants_mod,        only : i_native, &
                                   r_second, &
                                   l_def
  use lfric_xios_clock_mod, only : lfric_xios_clock_type
  use log_mod,              only : log_event,       &
                                   log_level_error, &
                                   log_level_info
  use xios,                 only : xios_close_context_definition, &
                                   xios_context,                  &
                                   xios_context_finalize,         &
                                   xios_context_initialize,       &
                                   xios_get_handle,               &
                                   xios_set_current_context

  implicit none

  private

  !> Manages interactions with XIOS.
  !>
  type, public, extends(io_context_type) :: lfric_xios_context_type
    private
    character(:),                 allocatable :: id
    type(xios_context)                        :: handle
    class(lfric_xios_clock_type), allocatable :: clock
  contains
    private
    procedure, public :: initialise
    procedure, public :: get_clock
    final :: finalise
  end type lfric_xios_context_type

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Set up an XIOS context object.
  !>
  !> @param [in]     id                Unique identifying string.
  !> @param [in]     communicator      MPI communicator used by context.
  !> @param [in,out] callback          Object to handle model specifics.
  !> @param [in]     start_time        Time of first step.
  !> @param [in]     finish_time       Time of last step.
  !> @param [in]     spinup_period     Number of seconds in spinup period.
  !> @param [in]     seconds_per_step  Number of seconds in a time step.
  !> @param [in]     timer_flag        Flag for use of subroutine timers.
  !>
  subroutine initialise( this, id, communicator,  &
                         callback,                &
                         start_time, finish_time, &
                         spinup_period,           &
                         seconds_per_step,        &
                         timer_flag )

    implicit none

    class(lfric_xios_context_type),     intent(inout) :: this
    character(*),                       intent(in)    :: id
    integer(i_native),                  intent(in)    :: communicator
    class(io_context_initialiser_type), intent(inout) :: callback
    character(*),                       intent(in)    :: start_time
    character(*),                       intent(in)    :: finish_time
    real(r_second),                     intent(in)    :: spinup_period
    real(r_second),                     intent(in)    :: seconds_per_step
    logical(l_def), optional,           intent(in)    :: timer_flag

    type(step_calendar_type), allocatable :: calendar
    integer(i_native)                     :: rc

    call xios_context_initialize( id, communicator )
    call xios_get_handle( id, this%handle )
    call xios_set_current_context( this%handle )

    allocate( calendar, stat=rc )
    if (rc /= 0) then
      call log_event( "Unable to allocate calendar", log_level_error )
    end if

    allocate( this%clock, stat=rc )
    if (rc /= 0) then
      call log_event( "Unable to allocate clock", log_level_error )
    end if
    call this%clock%initialise( calendar, start_time, finish_time, &
                                seconds_per_step, spinup_period,   &
                                timer_flag=timer_flag )

    !> @todo Rather than using this callback we might prefer to pass arrays
    !>       of objects which describe things to be set up. Alternatively we
    !>       could provide calls to "add axis" and "add file" with a final
    !>       "complete context".
    !>
    call callback%callback( this )

    call xios_close_context_definition()

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Clean up an XIOS context object.
  !>
  subroutine finalise( this )

    implicit none

    type(lfric_xios_context_type), intent(in) :: this

    call xios_context_finalize()

  end subroutine finalise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Gets the clock associated with this context.
  !>
  !> @return Clock object.
  !>
  function get_clock( this ) result(clock)

    implicit none

    class(lfric_xios_context_type), intent(in), target :: this
    class(clock_type), pointer :: clock

    clock => this%clock

  end function get_clock

end module lfric_xios_context_mod
