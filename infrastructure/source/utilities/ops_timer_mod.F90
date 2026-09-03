!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief A Simple timer based upon calls to cpu time or mpi_wtime that outputs
!>        results to stdout as soon as the instance calls stop_timer or is
!>        destroyed.
!>
!> Usage: The ops_timer type requires one instance per timing calliper.
!>
module ops_timer_mod

  use constants_mod, only: r_double, i_def

#ifdef NO_MPI
  ! No "use mpi" in non-mpi build
#else
#ifdef LEGACY_MPI
  use mpi,     only: mpi_wtime
#else
  use mpi_f08, only: mpi_wtime
#endif
#endif

  use log_mod, only: log_scratch_space, log_level_always, log_event

  implicit none
  private

#ifdef NO_MPI
  integer(i_def), save :: crate = -1_i_def
#endif

  type, public :: ops_timer_type
     private
     character(:), allocatable :: name
     real(r_double)              :: start_time  = 0.0_r_double
     real(r_double)              :: pause_start = 0.0_r_double
     real(r_double)              :: paused_time = 0.0_r_double
     logical                   :: running     = .false.
     logical                   :: paused      = .false.
   contains
     procedure :: start_timer
     procedure :: stop_timer
     procedure :: pause_timer
     procedure :: resume_timer
     procedure :: elapsed
     final     :: destructor
  end type ops_timer_type

contains

!=============================================================================!
!> @brief initialize an ops_timer instance and start timing with either
!>        mpi_wtime or system-clock
!> @param[in] name   The timing calliper's name, as will be logged when time
!>                   is output.
  subroutine start_timer(this, name)
    class(ops_timer_type), intent(inout) :: this
    character(*),          intent(in)    :: name

    this%name    = trim(name)
    this%running = .true.

#ifdef NO_MPI
    if (crate <= 0_i_def) call system_clock(count_rate=crate)
    block
      integer(i_def) :: count
      call system_clock(count=count)
      this%start_time = real(count, r_double)
    end block
#else
    this%start_time = mpi_wtime()
#endif
  end subroutine start_timer

!=============================================================================!
!> @brief Calculates the total time taken between ops_timer start and finish
!> @result  time_taken   The time measured by the ops_timer instance
  function elapsed(this) result(time_taken)

    class(ops_timer_type), intent(inout) :: this
    real(r_double) :: time_taken
#ifdef NO_MPI
    integer(i_def) :: now
#endif

    ! Close off any outstanding pause before calculating the elapsed time,
    ! otherwise the timer will not be using an accurate paused_time.
    if (this%paused) call this%resume_timer()

#ifdef NO_MPI
    call system_clock(count=now)
    time_taken = (real(now, r_double) - this%start_time) / real(crate, r_double)
    time_taken = time_taken - this%paused_time
#else
    time_taken = (mpi_wtime() - this%start_time) - this%paused_time
#endif
  end function elapsed

!=============================================================================!
!> @brief Instruct the ops_timer instance to start timing a new section to be
!>        subtracted from the total elapsed time when the timer is stopped.
  subroutine pause_timer(this)

    class(ops_timer_type), intent(inout) :: this

    if (.not. this%running) return
    ! If paused already, nothing to do
    if (this%paused) return
    this%paused = .true.

#ifdef NO_MPI
    if (crate <= 0_i_def) call system_clock(count_rate=crate)
    block
      integer(i_def) :: count
      call system_clock(count=count)
      this%pause_start = real(count, r_double)
    end block
#else
    this%pause_start = mpi_wtime()
#endif

  end subroutine pause_timer

!=============================================================================!
!> @brief Instruct the ops_timer instance to finish timing the paused section.
  subroutine resume_timer(this)

    class(ops_timer_type), intent(inout) :: this
    real(r_double) :: time_taken

    if (.not. this%running) return
    if (.not. this%paused) return
    this%paused = .false.

#ifdef NO_MPI
    block
      integer(i_def) :: now
      call system_clock(count=now)
      time_taken = (real(now, r_double) - this%pause_start) / real(crate, r_double)
    end block
#else
    time_taken = mpi_wtime() - this%pause_start
#endif
    this%paused_time = this%paused_time + time_taken

  end subroutine resume_timer

!=============================================================================!
!> @brief Instruct the ops_timer instance to stop timing and return the total
!>        time measured.
  subroutine stop_timer(this)

    class(ops_timer_type), intent(inout) :: this

    if (.not. this%running) return

    ! All ops_timer output is marked by the (OPS TIMER) identifier so it can
    ! be found easily
    write(log_scratch_space,'(3A,F21.4,A)') &
    '(OPS TIMER) Time taken for ', this%name, ' : ', this%elapsed(), ' (s)'
    call log_event(log_scratch_space, log_level_always)

    this%running = .false.

  end subroutine stop_timer

!=============================================================================!
!> @brief Calls the stop_timer subroutine to output the total time measured
!>        when going out of scope or being destroyed manually.
  subroutine destructor(this)
    type(ops_timer_type), intent(inout) :: this

    call stop_timer(this)

  end subroutine destructor

end module ops_timer_mod
