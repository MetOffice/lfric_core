!------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!
! Some of the content of this file has been produced with the assistance of
! Met Office GitHub Copilot Enterprise
!-----------------------------------------------------------------------------
!> @brief A Simple timer based upon calls to cpu time or mpi_wtime that outputs
!>        results to stdout as soon as the instance calls stop or is
!>        destroyed.
!>
!> Usage: The stopwatch type requires one instance per timing calliper.
!>
module stopwatch_mod

  use, intrinsic :: iso_fortran_env, only: real64, int64

  use log_mod, only: log_scratch_space, log_level_info, log_event

  implicit none
  private

  integer(int64), save :: crate = -1_int64

  type, public :: stopwatch_type
     private
     character(:), allocatable :: name
     real(real64)              :: start_time  = 0.0_real64
     real(real64)              :: pause_start = 0.0_real64
     real(real64)              :: paused_time = 0.0_real64
     logical                   :: running     = .false.
     logical                   :: paused      = .false.
   contains
     procedure :: start
     procedure :: stop
     procedure :: pause
     procedure :: resume
     procedure :: elapsed
     final     :: destructor
  end type stopwatch_type

contains

!=============================================================================!
!> @brief initialize a stopwatch instance and start timing
!> @param[in] name   The timing calliper's name, as will be logged when time
!>                   is output.
  subroutine start(this, name)
    class(stopwatch_type), intent(inout) :: this
    character(*),          intent(in)    :: name

    integer(int64) :: count

    this%name    = trim(name)
    this%running = .true.

    if (crate <= 0_int64) call system_clock(count_rate=crate)
    call system_clock(count=count)
    this%start_time = real(count, real64)

  end subroutine start

!=============================================================================!
!> @brief Calculates the total time taken between stopwatch start and finish
!> @result  time_taken   The time measured by the stopwatch instance
  function elapsed(this) result(time_taken)

    class(stopwatch_type), intent(inout) :: this
    real(real64)   :: time_taken
    integer(int64) :: now

    ! Close off any outstanding pause before calculating the elapsed time,
    ! otherwise the timer will not be using an accurate paused_time.
    if (this%paused) call this%resume()

    call system_clock(count=now)
    time_taken = (real(now, real64) - this%start_time) / real(crate, real64)
    time_taken = time_taken - this%paused_time

  end function elapsed

!=============================================================================!
!> @brief Instruct the stopwatch instance to start timing a new section to be
!>        subtracted from the total elapsed time when the timer is stopped.
  subroutine pause(this)

    class(stopwatch_type), intent(inout) :: this

    integer(int64) :: count

    if (.not. this%running) return
    ! If paused already, nothing to do
    if (this%paused) return
    this%paused = .true.

    if (crate <= 0_int64) call system_clock(count_rate=crate)
    call system_clock(count=count)
    this%pause_start = real(count, real64)

  end subroutine pause

!=============================================================================!
!> @brief Instruct the stopwatch instance to finish timing the paused section.
  subroutine resume(this)

    class(stopwatch_type), intent(inout) :: this
    real(real64)   :: time_taken
    integer(int64) :: now

    if (.not. this%running) return
    if (.not. this%paused) return
    this%paused = .false.

    call system_clock(count=now)
    time_taken = (real(now, real64) - this%pause_start) / real(crate, real64)

    this%paused_time = this%paused_time + time_taken

  end subroutine resume

!=============================================================================!
!> @brief Instruct the stopwatch instance to stop timing and return the total
!>        time measured.
  subroutine stop(this)

    class(stopwatch_type), intent(inout) :: this

    if (.not. this%running) return

    ! All stopwatch output is marked by the (STOPWATCH) identifier so it can
    ! be found easily
    write(log_scratch_space,'(3A,F21.4,A)') &
    '(STOPWATCH) Total time for ', this%name, ' : ', this%elapsed(), ' (s)'
    call log_event(log_scratch_space, log_level_info)

    this%paused  = .false.
    this%running = .false.

  end subroutine stop

!=============================================================================!
!> @brief Calls the stop subroutine to output the total time measured
!>        when going out of scope or being destroyed manually.
  subroutine destructor(this)
    type(stopwatch_type), intent(inout) :: this

    call stop(this)

  end subroutine destructor

end module stopwatch_mod
