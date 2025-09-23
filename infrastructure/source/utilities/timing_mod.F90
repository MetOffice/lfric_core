!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Provides wrapper support for Vernier timings
!>
module timing_mod
    use log_mod,            only:   log_event, log_scratch_space,     &
                                    LOG_LEVEL_DEBUG, LOG_LEVEL_WARNING
    use constants_mod,      only:   i_def,  IMDI

#ifdef VERNIER
    !Vernier will only be loaded if the VERNIER environment variable is used
    use vernier_mod,        only:   vernier_init, vernier_start,      &
                                    vernier_stop, vernier_write,      &
                                    vernier_finalize, vik

#endif

    implicit none

    public :: init_timing, final_timing, start_timing, stop_timing
    public :: tik

#ifdef VERNIER
    !If Vernier is on then the calliper hash 'tik' is defined as Vernier's hash 'vik'
    integer, parameter :: tik = vik

#else

    integer, parameter :: tik = i_def

#endif

#ifndef TIMING_ON
    !If the timing macro is not explicity turned on then the callipers won't be called
    logical, public, parameter :: LPROF = .false.

#else
    !If the timing macro is defined/ turned on, LPROF will be defined later (by subroutine_timers)
    logical, public :: LPROF

#endif

contains

!=============================================================================!
!> @brief Initialize timings
!> @param[in] communicator  LFRic mpi communicator
!> @param[in] lsubroutine_timers Runtime logical controlling timer use
    subroutine init_timing( communicator, lsubroutine_timers )
        use lfric_mpi_mod,      only: lfric_comm_type

        implicit none

        logical, intent(in)                 :: lsubroutine_timers
        type( lfric_comm_type ), intent(in) :: communicator

#ifdef TIMING_ON
        !If timing is on, LPROF will be defined by subroutine_timers
        LPROF = lsubroutine_timers

        write(log_scratch_space, '(A)') 'Timing Mod: Runtime timing is turned on'
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

#ifdef VERNIER
        !If Timing and Vernier is on, Vernier will be initialised

        write(log_scratch_space, '(A)') 'Timing Mod: Vernier is turned on'
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

        call vernier_init( communicator%get_comm_mpi_val() )

        write(log_scratch_space, '(A)') 'Timing Mod: Vernier initialised'
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
#endif

#ifndef VERNIER
        !If Timing is on but Vernier is not on then a warning will be thrown
        write(log_scratch_space, '(A)') 'Timing Mod: Runtime timing is turned on but no profiling tool (such as Vernier) is turned on!'
        call log_event(log_scratch_space, LOG_LEVEL_WARNING)
#endif
#endif

#ifndef TIMING_ON
#ifdef VERNIER
        write(log_scratch_space, '(A)') 'Timing Mod: Vernier is on but Timing is not!'
        call log_event(log_scratch_space, LOG_LEVEL_WARNING)
#endif
#endif

    end subroutine init_timing

!=============================================================================!
!> @brief Output and finalize timings
    subroutine final_timing()

        implicit none

#ifdef TIMING_ON
#ifdef VERNIER
        !If Vernier is on then it will write to a file and then finalise
        call vernier_write()

        write(log_scratch_space, '(A)') 'Timing Mod: Vernier has written to file'
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

        call vernier_finalize()

        write(log_scratch_space, '(A)') 'Timing Mod: Vernier finalised'
        call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
#endif
#endif

    end subroutine final_timing

!=============================================================================!
!> @brief Start timings
!> @param[out] timing_section_handle The name of the section that is being timed
!> @param[in]  timing_state_name     Starting or stopping the given timing, either 'start' or 'stop'
    subroutine start_timing( timing_section_handle, timing_section_name )

        implicit none

        character(*),   intent(in)  :: timing_section_name
        integer(tik),   intent(out) :: timing_section_handle

#ifdef VERNIER
        !If Vernier is on will start a calliper
        call vernier_start( timing_section_handle , timing_section_name )
#else
        timing_section_handle = IMDI
#endif

    end subroutine start_timing

    !=============================================================================!
!> @brief Stop timings
!> @param[in] timing_section_handle The name of the section that is being timed
    subroutine stop_timing( timing_section_handle )

        implicit none

        integer(tik),  intent(in) :: timing_section_handle
        !Future callipers may require the section name as well as the handle

#ifdef VERNIER
        !If Vernier is on will end a calliper
        call vernier_stop( timing_section_handle )
#endif

    end subroutine stop_timing

end module timing_mod
