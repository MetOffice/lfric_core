!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the gravity_wave miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module gravity_wave_driver_mod

  use checksum_alg_mod,               only: checksum_alg
  use constants_mod,                  only: i_def
  use derived_config_mod,             only: set_derived_config
  use yaxt,                           only: xt_initialize, xt_finalize
  use field_mod,                      only: field_type
  use finite_element_config_mod,      only: element_order
  use function_space_chain_mod,       only: function_space_chain_type
  use global_mesh_collection_mod,     only: global_mesh_collection, &
                                            global_mesh_collection_type
  use gravity_wave_mod,               only: load_configuration
  use gw_alg_mod,                     only: gravity_wave_alg_init, &
                                            gravity_wave_alg_step, &
                                            gravity_wave_alg_final
  use init_fem_mod,                   only: init_fem
  use init_gravity_wave_mod,          only: init_gravity_wave
  use init_mesh_mod,                  only: init_mesh
  use io_mod,                         only: output_nodal,      &
                                            output_xios_nodal, &
                                            xios_domain_init
  use log_mod,                        only: log_event,         &
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
  use operator_mod,                   only: operator_type
  use output_config_mod,              only: diagnostic_frequency, &
                                            subroutine_timers,    &
                                            write_nodal_output,   &
                                            write_xios_output
  use restart_config_mod,             only: restart_filename => filename
  use restart_control_mod,            only: restart_type
  use timer_mod,                      only: timer, output_timer
  use timestepping_config_mod,        only: dt
  use mpi_mod,                        only: initialise_comm, store_comm, &
                                            finalise_comm, &
                                            get_comm_size, get_comm_rank
  use xios

  implicit none

  private

  public initialise, run, finalise

  character(*), public, parameter   :: program_name = 'gravity_wave'
  character(*), public, parameter   :: xios_ctx = 'gravity_wave'
  character(*), public, parameter   :: xios_id  = 'lfric_client'

  type(restart_type) :: restart

  ! Prognostic fields
  type( field_type ) :: wind
  type( field_type ) :: buoyancy
  type( field_type ) :: pressure

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi

  integer(i_def) :: mesh_id
  integer(i_def) :: twod_mesh_id

  ! Function space chains
  type(function_space_chain_type) :: multigrid_function_space_chain

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise( filename )

  implicit none

  character(*), intent(in) :: filename

  integer(i_def) :: total_ranks, local_rank
  integer(i_def) :: comm = -999

  integer(i_def) :: ts_init
  integer(i_def) :: dtime

  ! Initialse mpi and create the default communicator: mpi_comm_world
  call initialise_comm(comm)

  ! Initialise XIOS and get back the split mpi communicator
  call init_wait()
  call xios_initialize(xios_id, return_comm = comm)

  !Store the MPI communicator for later use
  call store_comm(comm)

  ! Initialise YAXT
  call xt_initialize(comm)

  ! Get the rank information
  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()

  call initialise_logging(local_rank, total_ranks, program_name)

  call log_event( program_name//': Running miniapp ...', LOG_LEVEL_INFO )

  call load_configuration( filename )
  call set_derived_config( .false. )

  restart = restart_type( restart_filename, local_rank, total_ranks )

  !----------------------------------------------------------------------------
  ! Mesh init
  !----------------------------------------------------------------------------
  if ( subroutine_timers ) call timer(program_name)

  allocate( global_mesh_collection, &
       source = global_mesh_collection_type() )

  ! Create the mesh
  call init_mesh(local_rank, total_ranks, mesh_id, twod_mesh_id )

  !----------------------------------------------------------------------------
  ! FEM init
  !----------------------------------------------------------------------------
  ! Create FEM specifics (function spaces and chi field)
  call init_fem(mesh_id, chi)

  !----------------------------------------------------------------------------
  ! IO init
  !----------------------------------------------------------------------------

  ! If using XIOS for diagnostic output or checkpointing, then set up XIOS
  ! domain and context
  if ( (write_xios_output) .or. (restart%use_xios()) ) then

    dtime = int(dt)

    call xios_domain_init( xios_ctx,   &
                           comm,       &
                           dtime,      &
                           restart,    &
                           mesh_id,    &
                           chi)

    ! Make sure XIOS calendar is set to timestep 1 as it starts there
    ! not timestep 0.
    call xios_update_calendar(1)

  end if

  !----------------------------------------------------------------------------
  ! Model init
  !----------------------------------------------------------------------------
  multigrid_function_space_chain = function_space_chain_type()

  ! Create function space collection and initialise prognostic fields
  call init_gravity_wave( mesh_id, chi, multigrid_function_space_chain, &
                          wind, pressure, buoyancy, restart )

  ! Full global meshes no longer required, so reclaim
  ! the memory from global_mesh_collection
  write(log_scratch_space,'(A)') &
      "Purging global mesh collection."
  call log_event( log_scratch_space, LOG_LEVEL_INFO )
  deallocate(global_mesh_collection)

  ! Output initial conditions
  ! We only want these once at the beginning of a run
  ts_init = max( (restart%ts_start() - 1), 0 ) ! 0 or t previous.

  if (ts_init == 0) then

    if (write_nodal_output) then
      call output_nodal('wind',     ts_init, wind,     mesh_id)
      call output_nodal('pressure', ts_init, pressure, mesh_id)
      call output_nodal('buoyancy', ts_init, buoyancy, mesh_id)
    end if

    ! XIOS output
    if (write_xios_output) then

      ! Make sure XIOS calendar is set to first timestep to be computed
      call xios_update_calendar(ts_init+1)

      ! Output scalar fields
      call output_xios_nodal("init_wind",     wind,     mesh_id)
      call output_xios_nodal("init_pressure", pressure, mesh_id)
      call output_xios_nodal("init_buoyancy", buoyancy, mesh_id)

    end if

  end if

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run()

  implicit none

  integer(i_def) :: timestep

  !--------------------------------------------------------------------------
  ! Model step
  !--------------------------------------------------------------------------
  do timestep = restart%ts_start(),restart%ts_end()

    ! Update XIOS calendar if we are using it for diagnostic output or checkpoint
    if ( (write_xios_output) .or. (restart%use_xios()) ) then
      call log_event( program_name//': Updating XIOS timestep', LOG_LEVEL_INFO )
      call xios_update_calendar(timestep)
    end if

    call log_event( &
    "/****************************************************************************\ ", &
    LOG_LEVEL_TRACE )
    write( log_scratch_space, '(A,I0)' ) 'Start of timestep ', timestep
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    if (timestep == restart%ts_start()) then
      call gravity_wave_alg_init(mesh_id, wind, pressure, buoyancy)
    end if

    call gravity_wave_alg_step(wind, pressure, buoyancy)
    write( log_scratch_space, '(A,I0)' ) 'End of timestep ', timestep
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    call log_event( &
    '\****************************************************************************/ ', &
    LOG_LEVEL_INFO )
    if ( mod(timestep, diagnostic_frequency) == 0 ) then

      if (write_nodal_output) then
        call log_event(program_name//': Writing nodal diag output', LOG_LEVEL_INFO)
        call output_nodal('wind',     timestep, wind,     mesh_id)
        call output_nodal('pressure', timestep, pressure, mesh_id)
        call output_nodal('buoyancy', timestep, buoyancy, mesh_id)
      end if

      ! XIOS output
      if (write_xios_output) then

        call log_event(program_name//': Writing xios output', LOG_LEVEL_INFO)
        call output_xios_nodal('wind',     wind,     mesh_id)
        call output_xios_nodal('pressure', pressure, mesh_id)
        call output_xios_nodal('buoyancy', buoyancy, mesh_id)
      end if

    end if

  end do

  ! Write checkpoint/restart files if required
  if( restart%checkpoint() ) then
    write(log_scratch_space,'(A,I6)') &
        "Checkpointing pressure at timestep ", restart%ts_end()
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call pressure%write_checkpoint( 'checkpoint_pressure', &
                                    trim(restart%endfname("pressure")) )

    write(log_scratch_space,'(A,I6)') &
         "Checkpointing wind at timestep ", restart%ts_end()
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call wind%write_checkpoint( 'checkpoint_wind', &
                                trim(restart%endfname("wind")) )

    write(log_scratch_space,'(A,I6)') &
         "Checkpointing buoyancy at timestep ", restart%ts_end()
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call buoyancy%write_checkpoint( 'checkpoint_buoyancy', &
                                    trim(restart%endfname("buoyancy")) )
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

  ! Finalise XIOS context if we used it for diagnostic output or checkpointing
  if ( (write_xios_output) .or. (restart%use_xios()) ) then
    call xios_context_finalize()
  end if

  ! Finalise XIOS
  call xios_finalize()

  ! Finalise YAXT
  call xt_finalize()

  ! Finalise mpi and release the communicator
  call finalise_comm()

  ! Finalise the logging system
  call finalise_logging()

end subroutine finalise

end module gravity_wave_driver_mod
