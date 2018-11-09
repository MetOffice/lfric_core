!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @mainpage solver_miniapp
!> The solver API is abstract so cannot be tested unless with a particular implementation.
!! The mini-app makes an example linear operator, preconditioner and krylov subspace solver.
!! This can then be executed to test the solver api against different compilers.

program solver_miniapp

  use constants_mod,                  only : i_def
  use cli_mod,                        only : get_initial_filename
  use init_mesh_mod,                  only : init_mesh
  use init_fem_mod,                   only : init_fem  
  use init_solver_miniapp_mod,        only : init_solver_miniapp
  use yaxt,                           only : xt_initialize, xt_finalize
  use mpi_mod,                        only : initialise_comm, store_comm, &
                                             finalise_comm, &
                                             get_comm_size, get_comm_rank
  use global_mesh_collection_mod,     only : global_mesh_collection, &
                                             global_mesh_collection_type
  use field_mod,                      only : field_type
  use field_vector_mod,               only : field_vector_type
  use solver_miniapp_alg_mod,           only : solver_miniapp_alg
  use solver_miniapp_configuration_mod, only : final_configuration
  use solver_miniapp_mod,               only : load_configuration
  use log_mod,                        only : log_event,         &
                                             log_set_level,     &
                                             log_scratch_space, &
                                             initialise_logging, &
                                             finalise_logging, &
                                             LOG_LEVEL_ERROR,   &
                                             LOG_LEVEL_INFO
  use output_config_mod,              only : write_nodal_output
  use io_mod,                         only : output_nodal
  use checksum_alg_mod,               only : checksum_alg

  implicit none

  character(:), allocatable :: filename

  integer(i_def)     :: total_ranks, local_rank
  integer(i_def)     :: comm = -999

  integer(i_def)     :: mesh_id, twod_mesh_id

  ! prognostic fields
  type( field_type ), target, dimension(3) :: chi
  type( field_type ) :: field_1, field_2
  type( field_vector_type) :: fv_1
  !-----------------------------------------------------------------------------
  ! Driver layer init
  !-----------------------------------------------------------------------------

  ! Initialise MPI communicatios and get a valid communicator
  call initialise_comm(comm)

  ! Save the commmunicator for later use
  call store_comm(comm)

  ! Initialise YAXT
  call xt_initialize(comm)

  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()

  call initialise_logging(local_rank, total_ranks, 'solver_miniapp')

    ! First call to log event must occur after the call to initialise_logging
  call log_event( 'solver miniapp running...', LOG_LEVEL_INFO )

  call get_initial_filename( filename )

  call load_configuration( filename )
  deallocate( filename )

  !-----------------------------------------------------------------------------
  ! model init
  !-----------------------------------------------------------------------------

  ! Create the mesh and function space collection
  allocate( global_mesh_collection, &
            source = global_mesh_collection_type() )
  call init_mesh(local_rank, total_ranks, mesh_id, twod_mesh_id)
 ! Full global meshes no longer required, so reclaim
  ! the memory from global_mesh_collection
  write(log_scratch_space,'(A)') &
      "Purging global mesh collection."
  call log_event( log_scratch_space, LOG_LEVEL_INFO )
  deallocate(global_mesh_collection)

  call init_fem(mesh_id,chi)

  ! Create and initialise prognostic fields
  call init_solver_miniapp(mesh_id, chi, fv_1)

  ! Call an algorithm
  call solver_miniapp_alg(fv_1)

  ! Write out output file
  call log_event("solver miniapp: writing diagnostic output", LOG_LEVEL_INFO)

  ! pull the fields from the vector
  call fv_1%export_field( field_1, 1 )
  call fv_1%export_field( field_2, 2 )
  ! Original nodal output
  if ( write_nodal_output)  then
     call output_nodal('solver_field_1', 0, field_1, mesh_id)
     call output_nodal('solver_field_2', 0, field_2, mesh_id)
  end if

  !-----------------------------------------------------------------------------
  ! model finalise
  !-----------------------------------------------------------------------------

  ! Write checksums to file
  call checksum_alg('solver_miniapp', field_1, 'solver_field_1',field_2, 'solver_field_2')

  call log_event( 'solver miniapp completed', LOG_LEVEL_INFO )

  !-----------------------------------------------------------------------------
  ! Driver layer finalise
  !-----------------------------------------------------------------------------

  ! Finalise namelist configurations
  call final_configuration()

  ! Finalise YAXT
  call xt_finalize()

  ! Finalise MPI communications
  call finalise_comm()

  ! Finalise the logging system
  call finalise_logging()

end program solver_miniapp




