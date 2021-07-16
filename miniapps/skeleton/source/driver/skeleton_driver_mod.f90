!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the skeleton miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module skeleton_driver_mod

  use base_mesh_config_mod,       only : prime_mesh_name
  use checksum_alg_mod,           only : checksum_alg
  use clock_mod,                  only : clock_type
  use configuration_mod,          only : final_configuration
  use constants_mod,              only : i_def, i_native, &
                                         PRECISION_REAL
  use convert_to_upper_mod,       only : convert_to_upper
  use create_mesh_mod,            only : init_mesh
  use create_fem_mod,             only : init_fem
  use derived_config_mod,         only : set_derived_config
  use diagnostics_io_mod,         only : write_scalar_diagnostic
  use field_mod,                  only : field_type
  use init_skeleton_mod,          only : init_skeleton
  use lfric_xios_io_mod,          only : initialise_xios
  use io_config_mod,              only : write_diag, &
                                         use_xios_io
  use io_context_mod,             only : io_context_type
  use lfric_xios_clock_mod,       only : lfric_xios_clock_type
  use lfric_xios_context_mod,     only : lfric_xios_context_type
  use local_mesh_collection_mod,  only : local_mesh_collection, &
                                         local_mesh_collection_type
  use linked_list_mod,            only : linked_list_type
  use log_mod,                    only : log_event,          &
                                         log_set_level,      &
                                         log_scratch_space,  &
                                         initialise_logging, &
                                         finalise_logging,   &
                                         LOG_LEVEL_ALWAYS,   &
                                         LOG_LEVEL_ERROR,    &
                                         LOG_LEVEL_WARNING,  &
                                         LOG_LEVEL_INFO,     &
                                         LOG_LEVEL_DEBUG,    &
                                         LOG_LEVEL_TRACE
  use mesh_collection_mod,        only : mesh_collection, &
                                         mesh_collection_type
  use mpi_mod,                    only : store_comm,    &
                                         get_comm_size, &
                                         get_comm_rank
  use simple_io_mod,              only : initialise_simple_io
  use simple_io_context_mod,      only : simple_io_context_type
  use skeleton_mod,               only : load_configuration, program_name
  use skeleton_alg_mod,           only : skeleton_alg
  use time_config_mod,            only : timestep_start, timestep_end
  use timestepping_config_mod,    only : dt, &
                                         spinup_period
  use xios,                       only : xios_context_finalize, &
                                         xios_update_calendar
  use yaxt,                       only : xt_initialize, xt_finalize

  implicit none

  private
  public initialise, run, finalise

  class(io_context_type), allocatable :: io_context

  ! Prognostic fields
  type( field_type ) :: field_1

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi_xyz
  type(field_type), target, dimension(3) :: chi_sph
  type(field_type), target               :: panel_id

  integer(i_def) :: mesh_id
  integer(i_def) :: twod_mesh_id

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise( filename, model_communicator )

    use logging_config_mod, only: run_log_level,          &
                                  key_from_run_log_level, &
                                  RUN_LOG_LEVEL_ERROR,    &
                                  RUN_LOG_LEVEL_INFO,     &
                                  RUN_LOG_LEVEL_DEBUG,    &
                                  RUN_LOG_LEVEL_TRACE,    &
                                  RUN_LOG_LEVEL_WARNING

    implicit none

    character(:),      intent(in), allocatable :: filename
    integer(i_native), intent(in)              :: model_communicator

    integer(i_def)    :: total_ranks, local_rank
    integer(i_native) :: log_level

    !Store the MPI communicator for later use
    call store_comm( model_communicator )

    ! Initialise YAXT
    call xt_initialize( model_communicator )

    ! and get the rank information from the virtual machine
    total_ranks = get_comm_size()
    local_rank  = get_comm_rank()

    call initialise_logging(local_rank, total_ranks, program_name)

    call load_configuration( filename )

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

    call set_derived_config( .true. )

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    allocate( local_mesh_collection, &
              source = local_mesh_collection_type() )

    allocate( mesh_collection, &
              source=mesh_collection_type() )

    ! Create the mesh
    call init_mesh( local_rank, total_ranks, mesh_id, &
                    twod_mesh_id = twod_mesh_id )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh_id, chi_xyz, chi_sph, panel_id )

    !-------------------------------------------------------------------------
    ! IO init
    !-------------------------------------------------------------------------

    ! If using XIOS for diagnostic output or checkpointing, then set up
    ! XIOS domain and context

    if ( use_xios_io ) then
      call initialise_xios( io_context,         &
                            program_name,       &
                            model_communicator, &
                            mesh_id,            &
                            twod_mesh_id,       &
                            chi_sph,            &
                            panel_id,           &
                            timestep_start,     &
                            timestep_end,       &
                            spinup_period,      &
                            dt )
    else
      call initialise_simple_io( io_context,     &
                                 timestep_start, &
                                 timestep_end,   &
                                 spinup_period,  &
                                 dt )
    end if

    ! Create and initialise prognostic fields
    call init_skeleton(mesh_id, twod_mesh_id, chi_xyz, chi_sph, panel_id, field_1)

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run()

    implicit none

    class(clock_type), pointer :: clock
    logical                    :: running

    clock => io_context%get_clock()
    running = clock%tick()

    ! Call an algorithm
    call skeleton_alg(field_1)

    ! Write out output file
    call log_event(program_name//": Writing diagnostic output", LOG_LEVEL_INFO)

    if (write_diag ) then
      ! Calculation and output of diagnostics
      call write_scalar_diagnostic( 'skeleton_field', field_1, &
                                    clock, mesh_id, .false. )
    end if

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise()

    implicit none

!-----------------------------------------------------------------------------
    ! Model finalise
    !-----------------------------------------------------------------------------
    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    ! Write checksums to file
    call checksum_alg(program_name, field_1, 'skeleton_field_1')

    call log_event( program_name//': Miniapp completed', LOG_LEVEL_INFO )

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

   ! Finalise XIOS context if we used it for diagnostic output or checkpointing
    if ( use_xios_io ) then
      call xios_context_finalize()
    end if

    ! Finalise namelist configurations
    call final_configuration()

    ! Finalise YAXT
    call xt_finalize()

    call log_event( program_name//' completed.', LOG_LEVEL_ALWAYS )

    ! Finalise the logging system
    call finalise_logging()

  end subroutine finalise

end module skeleton_driver_mod
