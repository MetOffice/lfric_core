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

  use constants_mod,                    only : i_def, PRECISION_REAL
  use convert_to_upper_mod,             only : convert_to_upper
  use cli_mod,                          only : get_initial_filename
  use create_mesh_mod,                  only : init_mesh
  use create_fem_mod,                   only : init_fem
  use init_solver_miniapp_mod,          only : init_solver_miniapp
  use yaxt,                             only : xt_initialize, xt_finalize
  use mpi_mod,                          only : initialise_comm, store_comm, &
                                               finalise_comm,               &
                                               get_comm_size, get_comm_rank
  use local_mesh_collection_mod,        only : local_mesh_collection, &
                                               local_mesh_collection_type
  use field_mod,                        only : field_type
  use field_vector_mod,                 only : field_vector_type
  use solver_miniapp_alg_mod,           only : solver_miniapp_alg
  use configuration_mod,                only : final_configuration
  use solver_miniapp_mod,               only : load_configuration, program_name
  use log_mod,                          only : log_event,            &
                                               log_set_level,        &
                                               log_scratch_space,    &
                                               initialise_logging,   &
                                               finalise_logging,     &
                                               LOG_LEVEL_ALWAYS,     &
                                               LOG_LEVEL_ERROR,      &
                                               LOG_LEVEL_WARNING,    &
                                               LOG_LEVEL_INFO,       &
                                               LOG_LEVEL_DEBUG,      &
                                               LOG_LEVEL_TRACE
  use logging_config_mod,               only : run_log_level,          &
                                               key_from_run_log_level, &
                                               RUN_LOG_LEVEL_ERROR,    &
                                               RUN_LOG_LEVEL_INFO,     &
                                               RUN_LOG_LEVEL_DEBUG,    &
                                               RUN_LOG_LEVEL_TRACE,    &
                                               RUN_LOG_LEVEL_WARNING
  use mesh_collection_mod,              only : mesh_collection, &
                                               mesh_collection_type
  use checksum_alg_mod,                 only : checksum_alg

  implicit none

  character(:), allocatable :: filename

  integer(i_def) :: total_ranks, local_rank
  integer(i_def) :: log_level
  integer(i_def) :: comm = -999

  integer(i_def) :: mesh_id

  ! prognostic fields
  type( field_type ), target, dimension(3) :: chi_xyz
  type( field_type ), target, dimension(3) :: chi_sph
  type( field_type ), target :: panel_id
  type( field_type ) :: field_1, field_2
  type( field_vector_type) :: fv_1

  !-----------------------------------------------------------------------------
  ! Driver layer init
  !-----------------------------------------------------------------------------

  ! Initialise MPI communicatios and get a valid communicator
  call initialise_comm(comm)

  ! Initialise YAXT
  call xt_initialize(comm)

  ! Save the commmunicator for later use
  call store_comm(comm)

  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()

  call initialise_logging(local_rank, total_ranks, program_name)

  call get_initial_filename( filename )
  call load_configuration( filename )
  deallocate( filename )

  select case (run_log_level)
  case( RUN_LOG_LEVEL_ERROR )
    log_level = LOG_LEVEL_ERROR
  case( RUN_LOG_LEVEL_WARNING )
    log_level = LOG_LEVEL_WARNING
  case( RUN_LOG_LEVEL_INFO )
    log_level = LOG_LEVEL_INFO
  case( RUN_LOG_LEVEL_DEBUG )
    log_level = LOG_LEVEL_DEBUG
  case( RUN_LOG_LEVEL_TRACE )
    log_level = LOG_LEVEL_TRACE
  end select

  call log_set_level( log_level )

  write(log_scratch_space,'(A)')                              &
      'Runtime message logging severity set to log level: '// &
      convert_to_upper(key_from_run_log_level(run_log_level))
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

  write(log_scratch_space,'(A)')                        &
      'Application built with '//trim(PRECISION_REAL)// &
      '-bit real numbers'
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

  !-----------------------------------------------------------------------------
  ! model init
  !-----------------------------------------------------------------------------
  call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

  ! Create the mesh and function space collection
  allocate( local_mesh_collection, &
            source = local_mesh_collection_type() )
  allocate( mesh_collection, &
            source=mesh_collection_type() )

  call init_mesh( local_rank, total_ranks, mesh_id )

  call init_fem( mesh_id, chi_xyz, chi_sph, panel_id )

  ! Create and initialise prognostic fields
  call init_solver_miniapp(mesh_id, chi_sph, panel_id, fv_1)

  ! Call an algorithm
  call solver_miniapp_alg(fv_1, chi_sph, panel_id)

  ! Write out output file
  call log_event(program_name//": writing diagnostic output", LOG_LEVEL_INFO)

  ! pull the fields from the vector
  call fv_1%export_field( field_1, 1 )
  call fv_1%export_field( field_2, 2 )

  !-----------------------------------------------------------------------------
  ! model finalise
  !-----------------------------------------------------------------------------
  call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

  ! Write checksums to file
  call checksum_alg(program_name, field_1, 'solver_field_1',field_2, 'solver_field_2')

  call log_event( program_name//': Miniapp completed', LOG_LEVEL_INFO )

  !-----------------------------------------------------------------------------
  ! Driver layer finalise
  !-----------------------------------------------------------------------------

  ! Finalise namelist configurations
  call final_configuration()

  ! Finalise YAXT
  call xt_finalize()

  ! Finalise MPI communications
  call finalise_comm()

  call log_event( program_name//' completed.', LOG_LEVEL_ALWAYS )

  ! Finalise the logging system
  call finalise_logging()

end program solver_miniapp
