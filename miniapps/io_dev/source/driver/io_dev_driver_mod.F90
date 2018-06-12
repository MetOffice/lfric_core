!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the io_dev miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module io_dev_driver_mod

use constants_mod,              only: i_def
use derived_config_mod,         only: set_derived_config 
use field_mod,                  only: field_type
use finite_element_config_mod,  only: element_order
use global_mesh_collection_mod, only: global_mesh_collection, &
                                      global_mesh_collection_type
use init_io_dev_mod,            only: init_io_dev
use init_fem_mod,               only: init_fem
use init_mesh_mod,              only: init_mesh
use io_dev_mod,                 only: load_configuration
use io_mod,                     only: output_xios_nodal, &
                                      xios_domain_init
use log_mod,                    only: log_event,         &
                                      log_set_level,     &
                                      log_scratch_space, &
                                      initialise_logging, &
                                      finalise_logging, &
                                      LOG_LEVEL_ERROR,   &
                                      LOG_LEVEL_INFO,    &
                                      LOG_LEVEL_DEBUG,   &
                                      LOG_LEVEL_TRACE,   &
                                      log_scratch_space
use mod_wait
use restart_config_mod,         only: restart_filename => filename
use restart_control_mod,        only: restart_type
use timestepping_config_mod,    only: dt
use mpi_mod,                    only: initialise_comm, store_comm, &
                                      finalise_comm, &
                                      get_comm_size, get_comm_rank
use xios
use yaxt,                       only: xt_initialize, xt_finalize

implicit none

private

public initialise, run, finalise

character(*), public, parameter :: program_name = 'io_dev'
character(*), public, parameter :: xios_ctx = 'io_dev'
character(*), public, parameter :: xios_id  = 'lfric_client'

type(restart_type) :: restart

! Prognostic fields
type( field_type ) :: test_field

! Coordinate field
type(field_type), target, dimension(3) :: chi

integer(i_def) :: mesh_id
integer(i_def) :: twod_mesh_id

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise( filename )

  implicit none

  character(*), intent(in) :: filename

  integer(i_def) :: total_ranks, local_rank
  integer(i_def) :: comm = -999

  integer(i_def) :: dtime

  ! Initialse mpi and create the default communicator: mpi_comm_world
  call initialise_comm(comm)

  ! Initialise XIOS and get back the split mpi communicator
  call init_wait()
  call xios_initialize(xios_id, return_comm = comm)

  ! Save lfric's part of the split communicator for later use
  call store_comm(comm)

  ! Initialise YAXT
  call xt_initialize(comm)

  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()

  ! Store the MPI communicator for later use
  call store_comm(comm)

  call initialise_logging(local_rank, total_ranks, program_name)

  call log_event( program_name//': Running miniapp ...', LOG_LEVEL_INFO )

  call load_configuration( filename )
  call set_derived_config( .true. )

  restart = restart_type( restart_filename, local_rank, total_ranks )


  !----------------------------------------------------------------------------
  ! Mesh init
  !----------------------------------------------------------------------------
  allocate( global_mesh_collection, &
            source = global_mesh_collection_type() )

  ! Create the mesh
  call init_mesh(local_rank, total_ranks, mesh_id, twod_mesh_id)

  !----------------------------------------------------------------------------
  ! FEM init
  !----------------------------------------------------------------------------
  ! Create FEM specifics (function spaces and chi field)
  call init_fem(mesh_id, chi)


  !-----------------------------------------------------------------------------
  ! IO init
  !-----------------------------------------------------------------------------

  ! Set up XIOS domain and context
  dtime = int(dt)

  call xios_domain_init( xios_ctx, comm, dtime, &
                         restart, mesh_id, chi )

  !-----------------------------------------------------------------------------
  ! Model init
  !-----------------------------------------------------------------------------

  ! Create and initialise field
  call init_io_dev(mesh_id, chi, test_field)


  ! Full global meshes no longer required, so reclaim
  ! the memory from global_mesh_collection
  write(log_scratch_space,'(A)') &
      "Purging global mesh collection."
  call log_event( log_scratch_space, LOG_LEVEL_INFO )
  deallocate(global_mesh_collection)

  end subroutine initialise


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run()

  implicit none

  !-----------------------------------------------------------------------------
  ! Model step 
  !-----------------------------------------------------------------------------

  ! Update XIOS calendar
  call log_event( program_name//': Updating XIOS timestep', LOG_LEVEL_INFO )
  call xios_update_calendar(1)

  call log_event( program_name//': Writing XIOS output', LOG_LEVEL_INFO)
  call output_xios_nodal('test_field', test_field, mesh_id)

  end subroutine run

  !-----------------------------------------------------------------------------
  ! Driver layer finalise
  !-----------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise()

  implicit none

  integer(i_def) :: rc
  integer(i_def) :: ierr

  ! Finalise XIOS context
  call xios_context_finalize()

  ! Finalise XIOS
  call xios_finalize()

  ! Finalise YAXT
  call xt_finalize()

  ! Finalise mpi and release the communicator
  call finalise_comm()

  ! Finalise the logging system
  call finalise_logging()

  end subroutine finalise

end module io_dev_driver_mod
