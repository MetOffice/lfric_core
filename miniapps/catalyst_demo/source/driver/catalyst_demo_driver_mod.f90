!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the catalyst_demo miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module catalyst_demo_driver_mod

  use catalyst_demo_mod,              only: load_configuration
  use checksum_alg_mod,               only: checksum_alg
  use clock_mod,                      only: clock_type
  use constants_mod,                  only: i_def, i_native
  use derived_config_mod,             only: set_derived_config
  use field_mod,                      only: field_type
  use finite_element_config_mod,      only: element_order
  use formulation_config_mod,         only: l_multigrid
  use function_space_chain_mod,       only: function_space_chain_type
  use gw_alg_mod,                     only: gravity_wave_alg_init, &
                                            gravity_wave_alg_step, &
                                            gravity_wave_alg_final
  use init_catalyst_demo_mod,         only: init_catalyst_demo
  use init_clock_mod,                 only: initialise_clock
  use lfric_xios_io_mod,              only: initialise_xios
  use io_mod,                         only: ts_fname
  use diagnostics_io_mod,             only: write_scalar_diagnostic, &
                                            write_vector_diagnostic
  use local_mesh_collection_mod,      only: local_mesh_collection, &
                                            local_mesh_collection_type
  use log_mod,                        only: log_event,          &
                                            log_set_level,      &
                                            log_scratch_space,  &
                                            initialise_logging, &
                                            finalise_logging,   &
                                            LOG_LEVEL_ERROR,    &
                                            LOG_LEVEL_INFO,     &
                                            LOG_LEVEL_DEBUG,    &
                                            LOG_LEVEL_TRACE,    &
                                            log_scratch_space
  use mesh_collection_mod,            only: mesh_collection, &
                                            mesh_collection_type
  use operator_mod,                   only: operator_type
  use io_config_mod,                  only: write_diag,           &
                                            diagnostic_frequency, &
                                            use_xios_io,          &
                                            nodal_output_on_w3,   &
                                            checkpoint_write,     &
                                            checkpoint_stem_name, &
                                            subroutine_timers
  use timer_mod,                      only: init_timer, timer, output_timer
  use timestepping_config_mod,        only: dt
  use mpi_mod,                        only: store_comm,    &
                                            get_comm_size, &
                                            get_comm_rank
  use visualisation_config_mod,       only: write_catalyst_output,   &
                                            visualisation_frequency, &
                                            visualisation_stem_name, &
                                            use_python_vis_pipeline, &
                                            python_vis_pipeline_name
  use visualisation_mod,              only: vis_field_list_type, &
                                            catalyst_initialize, &
                                            catalyst_coprocess,  &
                                            catalyst_finalize
  use xios,                           only: xios_ctx,              &
                                            xios_id,               &
                                            xios_context_finalize, &
                                            xios_domain_init,      &
                                            xios_finalize,         &
                                            xios_initialize,       &
                                            xios_update_calendar
                                            initialise_xios,       &
                                            xios_finalize,         &
                                            xios_initialize,       &
                                            xios_update_calendar

  implicit none

  private

  public initialise, run, finalise

  character(*), public, parameter   :: program_name = 'catalyst_demo'
  character(*), public, parameter   :: xios_ctx = 'catalyst_demo'

  type(restart_type) :: restart

  ! Prognostic fields
  type( field_type ) :: wind
  type( field_type ) :: buoyancy
  type( field_type ) :: pressure

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi_xyz
  type(field_type), target, dimension(3) :: chi_sph
  type(field_type), target,              :: panel_id

  integer(i_def) :: mesh_id
  integer(i_def) :: twod_mesh_id

  integer(i_def), allocatable :: multigrid_mesh_ids(:)

  ! Function space chains
  type(function_space_chain_type) :: multigrid_function_space_chain

  ! Lists of fields for visualisation
  type(vis_field_list_type) :: vis_fields

  class(clock_type), allocatable :: clock

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise( filename, model_communicator )

  implicit none

  character(*),      intent(in) :: filename
  integer(i_native), intent(in) :: model_communicator

  integer(i_def) :: total_ranks, local_rank

  ! Initialse mpi and create the default communicator: mpi_comm_world
  call initialise_comm(comm)

  ! Initialise XIOS and get back the split mpi communicator
  call init_wait()
  call xios_initialize(xios_id, return_comm = comm)

  !Store the MPI communicator for later use
  call store_comm(model_communicator)

  ! Initialise YAXT
  call xt_initialize(model_communicator)

  ! Get the rank information
  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()

  call initialise_logging(local_rank, total_ranks, program_name)

  call log_event( program_name//': Running miniapp ...', LOG_LEVEL_INFO )

  call load_configuration( filename )
  call set_derived_config( .false. )

  call initialise_clock( clock )

  restart = restart_type( restart_filename, local_rank, total_ranks )

  !----------------------------------------------------------------------------
  ! Mesh init
  !----------------------------------------------------------------------------
  if ( subroutine_timers ) then
    call init_timer()
    call timer(program_name)
  end if

  allocate( local_mesh_collection, &
       source = local_mesh_collection_type() )
  allocate( mesh_collection, &
       source=mesh_collection_type() )

  ! Create the mesh
  call init_mesh( local_rank, total_ranks, mesh_id, &
                  twod_mesh_id=twod_mesh_id,        &
                  multigrid_mesh_ids=multigrid_mesh_ids, &
                  use_multigrid=l_multigrid )

  !----------------------------------------------------------------------------
  ! FEM init
  !----------------------------------------------------------------------------
  ! Create FEM specifics (function spaces and chi field)
  call init_fem( mash_id, chi_xyz, chi_sph, panel_id )


  !----------------------------------------------------------------------------
  ! IO init
  !----------------------------------------------------------------------------

  ! If using XIOS for diagnostic output or checkpointing, then set up XIOS
  ! domain and context
  if ( use_xios_io ) then
    call initialise_xios( xios_ctx,           &
                          model_communicator, &
                          clock,              &
                          mesh_id,            &
                          twod_mesh_id,       &
                          chi_sph,            &
                          panel_id )

    ! Make sure XIOS calendar is set to timestep 1 as it starts there
    ! not timestep 0.
    call xios_update_calendar(1)

  end if

  !----------------------------------------------------------------------------
  ! Visualisation init
  !----------------------------------------------------------------------------

  if ( write_catalyst_output ) then
    call catalyst_initialize(visualisation_frequency,       &
                             trim(visualisation_stem_name), &
                             model_communicator,            &
                             use_python_vis_pipeline,       &
                             trim(python_vis_pipeline_name))
  end if

  !----------------------------------------------------------------------------
  ! Model init
  !----------------------------------------------------------------------------
  multigrid_function_space_chain = function_space_chain_type()

  ! Create function space collection and initialise prognostic fields
  call init_catalyst_demo( mesh_id, twod_mesh_id, multigrid_mesh_ids, &
                           chi_xyz, chi_sph, panel_id,                &
                           multigrid_function_space_chain,            &
                           wind, pressure, buoyancy )

  ! Output initial conditions
  ! We only want these once at the beginning of a run
  if (clock%is_initialisation()) then

    if ( use_xios_io ) then

      ! Need to ensure calendar is initialised here as XIOS has no concept of timestep 0
      call xios_update_calendar(1)

    end if

    if ( write_diag ) then
      ! Calculation and output of diagnostics
      call write_vector_diagnostic('wind', wind, &
                                   clock, mesh_id, nodal_output_on_w3)
      call write_scalar_diagnostic('pressure', pressure, &
                                   clock, mesh_id, nodal_output_on_w3)
      call write_scalar_diagnostic('buoyancy', buoyancy, &
                                   clock, mesh_id, nodal_output_on_w3)
    end if

    ! Catalyst visualisation
    if (write_catalyst_output) then
      vis_fields = vis_field_list_type()
      call vis_fields%add_vis_field(wind, 'wind')
      call vis_fields%add_vis_field(pressure, 'pressure')
      call vis_fields%add_vis_field(buoyancy, 'buoyancy')
      call catalyst_coprocess(ts_init, ts_init*dt, vis_fields, mesh_id)
    end if

  end if

  call gravity_wave_alg_init(mesh_id, wind, pressure, buoyancy)

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run()

  implicit none

  !--------------------------------------------------------------------------
  ! Model step
  !--------------------------------------------------------------------------
  do while(clock%tick())

    ! Update XIOS calendar if we are using it for diagnostic output or checkpoint
    if ( use_xios_io ) then
      call log_event( program_name//': Updating XIOS timestep', LOG_LEVEL_INFO )
      call xios_update_calendar( clock%get_step() )
    end if

    write( log_scratch_space, '("/", A, "\ ")' ) repeat( '*', 76 )
    call log_event( log_scratch_space, LOG_LEVEL_TRACE )
    write( log_scratch_space, '(A,I0)' ) 'Start of timestep ', clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    call gravity_wave_alg_step(wind, pressure, buoyancy)
    write( log_scratch_space, '(A,I0)' ) 'End of timestep ', clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '("\", A, "/")' ) repeat( '*', 76 )
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    if ( (mod(clock%get_step(), diagnostic_frequency) == 0) \
         .and. (write_diag) ) then

      call log_event("Catalyst demo: writing diagnostic output", &
                     LOG_LEVEL_INFO)

      call write_vector_diagnostic('wind', wind, &
                                   clock, mesh_id, nodal_output_on_w3)
      call write_scalar_diagnostic('pressure', pressure, &
                                   clock, mesh_id, nodal_output_on_w3)
      call write_scalar_diagnostic('buoyancy', buoyancy, &
                                   clock, mesh_id, nodal_output_on_w3)

    end if

    ! Catalyst uses its own mechanism to determine if output should be written
    if (write_catalyst_output) then
      call catalyst_coprocess(timestep, timestep*dt, vis_fields, mesh_id)
    end if

  end do

  ! Write checkpoint/restart files if required
  if( checkpoint_write ) then
    write(log_scratch_space,'(A,I6)') &
        "Checkpointing pressure at timestep ", timestep_end
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call pressure%write_checkpoint("checkpoint_pressure", trim(ts_fname(checkpoint_stem_name,&
                                  "", "pressure", timestep_end,"")))

    write(log_scratch_space,'(A,I6)') &
         "Checkpointing wind at timestep ", timestep_end
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call wind%write_checkpoint("checkpoint_wind", trim(ts_fname(checkpoint_stem_name,&
                                  "", "wind", timestep_end,"")))

    write(log_scratch_space,'(A,I6)') &
         "Checkpointing buoyancy at timestep ", timestep_end
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call buoyancy%write_checkpoint("checkpoint_buoyancy", trim(ts_fname(checkpoint_stem_name,&
                                  "", "buoyancy", timestep_end,"")))
  end if

  end subroutine run


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise()

  implicit none

  integer(i_def) :: rc
  integer(i_def) :: ierr

  !----------------------------------------------------------------------------
  ! Model finalise
  !----------------------------------------------------------------------------

  ! Write checksums to file
  call checksum_alg( program_name, wind, 'wind', buoyancy, 'buoyancy', pressure, 'pressure')

  call log_event( program_name//': Miniapp run complete', LOG_LEVEL_INFO )
  if ( subroutine_timers ) then
    call timer(program_name)
    call output_timer()
  end if

  call gravity_wave_alg_final()

  !----------------------------------------------------------------------------
  ! Driver layer finalise
  !----------------------------------------------------------------------------

  if ( write_catalyst_output ) then
    call catalyst_finalize()
  end if

  ! Finalise XIOS context if we used it for diagnostic output or checkpointing
  if ( use_xios_io ) then
    call xios_context_finalize()
  end if

  ! Finalise YAXT
  call xt_finalize()

  ! Finalise the logging system
  call finalise_logging()

end subroutine finalise

end module catalyst_demo_driver_mod
