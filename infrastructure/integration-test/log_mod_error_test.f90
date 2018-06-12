!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

! A very simply program which just logs an error.
!
program log_mod_error_test

  use iso_fortran_env, only : error_unit
  use log_mod,         only : initialise_logging, finalise_logging, log_event, &
                              LOG_LEVEL_ERROR
  use mpi_mod,         only : initialise_comm, store_comm, &
                              finalise_comm, &
                              get_comm_size, get_comm_rank

  integer :: total_ranks, local_rank, comm

  ! Initialse mpi and create the default communicator: mpi_comm_world
  call initialise_comm(comm)
  ! Store the communicator for later use
  call store_comm(comm)
  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()
  call initialise_logging(local_rank, total_ranks, 'log_mod_error_test')

  ! Everything else is here purely to support the testing of this line:
  !
  call log_event( 'An error was logged.', LOG_LEVEL_ERROR )

  ! Finalise mpi and release the communicator
  call finalise_comm()

  ! Finalise the logging system
  call finalise_logging()

end program log_mod_error_test
