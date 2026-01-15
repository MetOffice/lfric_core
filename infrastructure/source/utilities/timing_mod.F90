!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Provides wrapper support for profiler timings
!>
module timing_mod
  use log_mod,            only:   log_event, log_scratch_space,     &
                                  LOG_LEVEL_DEBUG, LOG_LEVEL_WARNING
  use constants_mod,      only:   i_def, imdi, cmdi, str_def

#ifdef VERNIER
  use vernier_mod,        only:   vernier_init, vernier_start,      &
                                  vernier_stop, vernier_write,      &
                                  vernier_finalize, vik

#elif defined( LEGACY_TIMER )
  use timer_mod,          only: timer, init_timer, output_timer

#endif

  implicit none

  public :: init_timing, final_timing, start_timing, stop_timing
  public :: tik

#ifdef VERNIER
  integer, parameter :: tik = vik
  integer(tik), private :: global_timing_handle

#else

  integer, parameter :: tik = i_def

#endif

#ifdef TIMING_ON
 ! LPROF enables profiler timings.
 logical, public, protected :: LPROF = .false.

#else
 ! LPROF enables profiler timings.
 ! The logical is declared as a parameter for this build
 ! so compilers can easily optimise out the profiler
 ! calliper calls from the code.
 logical, public, parameter :: LPROF = .false.

#endif

contains

!=============================================================================!
!> @brief Initialise timings and start a global calliper
!> @param[in] communicator       LFRic mpi communicator
!> @param[in] lsubroutine_timers Runtime logical controlling timer use
!> @param[in] application_name   String for the global calliper
!> @param[in] timer_output_path  Temporary string used for the legacy timer path
  subroutine init_timing( communicator, lsubroutine_timers, application_name, &
                          timer_output_path )
    use lfric_mpi_mod,      only: lfric_comm_type

    implicit none

    type( lfric_comm_type ), intent(in) :: communicator
    logical,      intent(in)            :: lsubroutine_timers
    character(*), intent(in)            :: application_name
    character(*), intent(in), optional  :: timer_output_path

#ifdef TIMING_ON
    character(str_def) :: name

    ! If timing is on, LPROF will be defined by subroutine_timers
    LPROF = lsubroutine_timers
    name = cmdi

#ifdef LEGACY_TIMER
    name = 'Timer'
    if ( LPROF ) then
      if ( present ( timer_output_path ) ) then
        call init_timer( timer_output_path )
      else
        call init_timer( 'timer.txt' )
      end if

      call timer( application_name )

    end if

#elif defined( VERNIER )
    name = 'Vernier'
    if ( LPROF ) then
      call vernier_init( communicator%get_comm_mpi_val() )
      if ( LPROF ) call vernier_start( global_timing_handle, '__' // &
                                       application_name // '__' )

    end if

#endif

  if ( LPROF ) then
    if (trim( name ) == trim( cmdi )) then
      call log_event('Subroutine timings unavailable, no profiler compiled', &
                      log_level_warning)
    else
      call log_event( trim( name ) // ' initialised', log_level_debug )
    end if
  end if

#endif

  end subroutine init_timing

!=============================================================================!
!> @brief Output and finalise timings
  subroutine final_timing( application_name )
    implicit none
    character(*), intent(in)            :: application_name

#ifdef TIMING_ON
#ifdef VERNIER
    ! If Vernier is on then it will write to a file and then finalise
    if ( LPROF ) then
      call vernier_stop( global_timing_handle )
      call vernier_write()
      write(log_scratch_space, '(A)') 'Timing Mod: Vernier has written to file'
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

      call vernier_finalize()
      write(log_scratch_space, '(A)') 'Timing Mod: Vernier finalised'
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end if

#elif defined(LEGACY_TIMER)
    if ( LPROF ) then
      call timer ( application_name )
      call output_timer()

      write(log_scratch_space, '(A)') 'Timing Mod: Legacy timing finalised'
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end if

#endif
#endif
  end subroutine final_timing

!=============================================================================!
!> @brief Start timings
!> @param[out] timing_section_handle The integer handle for timed region
!> @param[in]  timing_state_name     Name of the measured region
  subroutine start_timing( timing_section_handle, timing_section_name )

    implicit none

    integer(tik), intent(out) :: timing_section_handle
    character(*), intent(in)  :: timing_section_name

    timing_section_handle = imdi

#ifdef VERNIER
    ! If Vernier is on will start a calliper
    call vernier_start( timing_section_handle , timing_section_name )

#elif defined(LEGACY_TIMER)
    call timer( timing_section_name )

#endif

  end subroutine start_timing

!=============================================================================!
!> @brief Stop timings
!> @param[in] timing_section_handle The integer handle for timed region
!> @param[in] timing_state_name     Optional, name of the measured region
  subroutine stop_timing( timing_section_handle, timing_section_name )
    implicit none

    integer(tik), optional, intent(in) :: timing_section_handle
    character(*), optional, intent(in) :: timing_section_name

#ifdef VERNIER
    ! If Vernier is on will end a calliper
    call vernier_stop( timing_section_handle )

#elif defined(LEGACY_TIMER)
    call timer( timing_section_name )

#endif

  end subroutine stop_timing

end module timing_mod
