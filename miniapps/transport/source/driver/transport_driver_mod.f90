!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>
!> @brief Drives the execution of the transport miniapp.
!>
module transport_driver_mod

  use checksum_alg_mod,               only: checksum_alg
  use constants_mod,                  only: i_def
  use derived_config_mod,             only: set_derived_config
  use yaxt,                           only: xt_initialize, xt_finalize
  use field_mod,                      only: field_type
  use global_mesh_collection_mod,     only: global_mesh_collection,           &
                                            global_mesh_collection_type
  use transport_configuration_mod,    only: final_configuration
  use transport_mod,                  only: transport_load_configuration
  use init_fem_mod,                   only: init_fem
  use init_transport_mod,             only: init_transport
  use init_mesh_mod,                  only: init_mesh
  use init_shifted_fields_mod,        only: init_shifted_fields
  use io_mod,                         only: xios_domain_init
  use diagnostics_io_mod,             only: write_scalar_diagnostic,          &
                                            write_vector_diagnostic      
  use diagnostics_mod,                only: write_density_diagnostic
  use log_mod,                        only: log_event,                        &
                                            log_scratch_space,                &
                                            initialise_logging,               &
                                            finalise_logging,                 &
                                            LOG_LEVEL_ERROR,                  &
                                            LOG_LEVEL_INFO,                   &
                                            LOG_LEVEL_TRACE
  use output_config_mod,              only: diagnostic_frequency,             &
                                            nodal_output_on_w3,               &
                                            write_xios_output,                &
                                            subroutine_timers
  use restart_config_mod,             only: restart_filename => filename
  use restart_control_mod,            only: restart_type
  use timer_mod,                      only: timer, output_timer
  use timestepping_config_mod,        only: dt
  use mpi_mod,                        only: initialise_comm, store_comm,      &
                                            finalise_comm,                    &
                                            get_comm_size, get_comm_rank
  use mod_wait,                       only: init_wait
  use transport_config_mod,           only: scheme,                           &
                                            transport_scheme_yz_bip_cosmic,   &
                                            transport_scheme_cosmic_3D,       &
                                            transport_scheme_horz_cosmic
  use mass_conservation_alg_mod,      only: mass_conservation
  use yz_bip_cosmic_alg_mod,          only: yz_bip_cosmic_step
  use cusph_cosmic_transport_alg_mod, only: set_winds,      &
                                            cusph_cosmic_transport_step
  use cosmic_threed_alg_mod,          only: cosmic_threed_transport_step
  use calc_dep_pts_alg_mod,           only: calc_dep_pts
  use density_inc_update_alg_mod,     only: density_inc_update_alg
  use xios

  implicit none

  private

  public :: initialise_transport, run_transport, finalise_transport

  character(*), public, parameter :: program_name = 'transport'

  type(restart_type) :: restart

  ! Prognostic fields
  type(field_type) :: wind
  type(field_type) :: density
  type(field_type) :: dep_pts_x
  type(field_type) :: dep_pts_y
  type(field_type) :: dep_pts_z

  ! Shifted Prognostic fields
  type(field_type) :: wind_shifted
  type(field_type) :: density_shifted
  type(field_type) :: increment

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi
  type(field_type), target, dimension(3) :: shifted_chi

  integer(i_def) :: mesh_id
  integer(i_def) :: twod_mesh_id
  integer(i_def) :: shifted_mesh_id

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  !! @param[in] filename Configuration namelist file
  subroutine initialise_transport( filename )

    implicit none

    character(*), intent(in) :: filename

    character(len=*), parameter :: xios_id   = "lfric_client"
    character(len=*), parameter :: xios_ctx  = "transport"

    integer(i_def) :: total_ranks, local_rank
    integer(i_def) :: comm = -999
    integer(i_def) :: timestep, ts_init, dtime

    ! Initialse mpi and create the default communicator: mpi_comm_world
    call initialise_comm( comm )

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

    call initialise_logging( local_rank, total_ranks, program_name )

    call log_event( program_name//': Running transport miniapp ...', LOG_LEVEL_INFO )

    call transport_load_configuration( filename )
    call set_derived_config( .true. )

    restart = restart_type( restart_filename, local_rank, total_ranks )

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    if ( subroutine_timers ) call timer( program_name )

    ! Mesh initialisation
    allocate( global_mesh_collection, &
              source = global_mesh_collection_type() )

    ! Create the mesh
    call init_mesh( local_rank, total_ranks, mesh_id, twod_mesh_id, shifted_mesh_id )

    ! FEM  initialisation
    call init_fem( mesh_id, chi, shifted_mesh_id, shifted_chi )

    ! Transport initialisation
    call init_transport( mesh_id, chi, wind, density, dep_pts_x, dep_pts_y, dep_pts_z, increment )
    call init_shifted_fields( shifted_mesh_id, wind_shifted, density_shifted )


    if ( (write_xios_output) ) then
      dtime = int(dt)
      call xios_domain_init( xios_ctx,     &
                             comm,         &
                             dtime,        &
                             restart,      &
                             mesh_id,      &
                             twod_mesh_id, &
                             chi)
    end if

    ! Output initial conditions
    ts_init = max( (restart%ts_start() - 1), 0 )
    if ( ts_init == 0 ) then

      if (write_xios_output) then

        ! Need to ensure calendar is initialised here as XIOS has no concept of timestep 0
        call xios_update_calendar(1)

      end if

      call write_vector_diagnostic('wind', wind, ts_init, mesh_id, nodal_output_on_w3)
      call write_scalar_diagnostic('density', density, ts_init, mesh_id, nodal_output_on_w3)

    end if

  end subroutine initialise_transport

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Performs time stepping.
  !>
  subroutine run_transport()

    implicit none

    integer(i_def) :: timestep

    !--------------------------------------------------------------------------
    ! Model step
    !--------------------------------------------------------------------------
    do timestep = restart%ts_start(), restart%ts_end()

      call log_event( &
      "/****************************************************************************\ ", &
      LOG_LEVEL_TRACE )
      write( log_scratch_space, '(A,I0)' ) 'Start of timestep ', timestep
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      ! Update XIOS calendar if we are using it for diagnostic output or checkpoint
      if ( write_xios_output ) then
        call log_event( "Transport: Updating XIOS timestep", LOG_LEVEL_INFO )
        call xios_update_calendar( timestep )
      end if

      ! Update the wind each timestep.
      call set_winds( wind, mesh_id, timestep )
      ! Calculate departure points.
      call calc_dep_pts( dep_pts_x, dep_pts_y, dep_pts_z, wind, chi )

      select case( scheme )
        case ( transport_scheme_yz_bip_cosmic )
          call yz_bip_cosmic_step( increment, density, dep_pts_y, dep_pts_z )
        case ( transport_scheme_horz_cosmic )
          call cusph_cosmic_transport_step( increment, density, dep_pts_x, dep_pts_y )
        case ( transport_scheme_cosmic_3D )
          call cosmic_threed_transport_step( increment, density, dep_pts_x, dep_pts_y, dep_pts_z )
        case default
          call log_event( "Transport mini-app: incorrect transport option chosen, "// &
                          "stopping program! ",LOG_LEVEL_ERROR )
          stop
      end select

      ! Add the increment to density
      call density_inc_update_alg(density, increment)
      call write_density_diagnostic( density, timestep )

      call mass_conservation( timestep, density )

      write( log_scratch_space, '(A,I0)' ) 'End of timestep ', timestep
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      call log_event( &
      '\****************************************************************************/ ', &
      LOG_LEVEL_INFO )

      ! Output wind and density values.
      if ( mod( timestep, diagnostic_frequency ) == 0 ) then

        call write_vector_diagnostic('wind', wind, timestep, mesh_id, nodal_output_on_w3)
        call write_scalar_diagnostic('density', density, timestep, mesh_id, nodal_output_on_w3)

      end if

    end do

  end subroutine run_transport

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Tidies up after a run.
  !>
  subroutine finalise_transport()

    implicit none

    !----------------------------------------------------------------------------
    ! Model finalise
    !----------------------------------------------------------------------------
    ! Write checksums to file
    call checksum_alg( program_name, density, 'density', wind, 'wind' )

    call log_event( program_name//': Miniapp run complete', LOG_LEVEL_INFO )
    if ( subroutine_timers ) then
      call timer( program_name )
      call output_timer()
    end if

    if ( write_xios_output ) then
      call xios_context_finalize()
    end if

    ! Finalise XIOS
    call xios_finalize()

    ! Finalise namelist configurations
    call final_configuration()

    !----------------------------------------------------------------------------
    ! Driver layer finalise
    !----------------------------------------------------------------------------
    ! Finalise YAXT
    call xt_finalize()

    ! Finalise mpi and release the communicator
    call finalise_comm()

    ! Finalise the logging system
    call finalise_logging()

  end subroutine finalise_transport

end module transport_driver_mod
