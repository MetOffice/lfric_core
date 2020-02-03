!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Handles initialisation and finalisation of infrastructure, constants
! and the gungho model simulations
module gungho_model_mod


  use constants_mod,              only : i_def, i_native, &
                                         PRECISION_REAL
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

  use io_config_mod,              only : subroutine_timers, &
                                         subroutine_counters, &
                                         use_xios_io, &
                                         write_diag, &
                                         write_dump, &
                                         write_minmax_tseries

  use timer_mod,                  only : timer, output_timer, init_timer

  use assign_orography_field_mod, only : assign_orography_field
  use field_mod,                  only : field_type, &
                                         write_interface
  use count_mod,                  only : count_type, halo_calls
  use runtime_constants_mod,      only : create_runtime_constants, &
                                         final_runtime_constants
#ifdef UM_PHYSICS
  use planet_constants_mod,       only : set_planet_constants
#endif
  use convert_to_upper_mod,       only : convert_to_upper
  use mpi_mod,                    only : store_comm, &
                                         get_comm_size, get_comm_rank
  use gungho_mod,                 only : load_configuration
  use configuration_mod,          only : final_configuration
  use derived_config_mod,         only : set_derived_config
  use yaxt,                       only : xt_initialize, xt_finalize
  use mod_wait,                   only : init_wait
  use global_mesh_collection_mod, only : global_mesh_collection, &
                                         global_mesh_collection_type
  use create_mesh_mod,            only : init_mesh, final_mesh
  use create_fem_mod,             only : init_fem, final_fem
  use timestepping_config_mod,        only : method, &
                                             method_semi_implicit, &
                                             method_rk, &
                                             dt
  use time_config_mod,            only : timestep_start

  use io_mod,                     only : xios_domain_init, &
                                         write_state, &
                                         xios_write_field_single_face

  use xios,                       only : xios_initialize,       &
                                         xios_finalize,         &
                                         xios_context_finalize, &
                                         xios_update_calendar
  use checksum_alg_mod,           only : checksum_alg
  use conservation_algorithm_mod, only : conservation_algorithm
  use field_collection_mod,       only : field_collection_type, &
                                         field_collection_iterator_type
  use gungho_model_data_mod,      only : model_data_type
  use formulation_config_mod,     only : transport_only, &
                                         use_moisture,   &
                                         use_physics
  use iter_timestep_alg_mod,      only : iter_alg_init, &
                                         iter_alg_final
  use minmax_tseries_mod,         only : minmax_tseries,      &
                                         minmax_tseries_init, &
                                         minmax_tseries_final
  use mr_indices_mod,             only : nummr
  use rk_alg_timestep_mod,        only : rk_alg_init, &
                                         rk_alg_final
  use rk_transport_mod,           only : rk_transport_init, &
                                         rk_transport_final
  use runge_kutta_init_mod,       only : runge_kutta_init, &
                                         runge_kutta_final
  use transport_config_mod,       only : scheme, &
                                         scheme_method_of_lines
  use init_altitude_mod,          only : init_altitude
  implicit none

  private
  public initialise_infrastructure, initialise_model, finalise_infrastructure, finalise_model

  contains


  !> @brief Initialises the infrastructure and sets up constants used by the model
  !> @param [inout] comm The MPI communicator for use within the model
  !>                     (not XIOS' communicator)
  !> @param [in]    filename The name of the configuration namelist file
  !> @param [in]    program_name An identifier given to the model begin run
  !> @param [inout] mesh_id The identifier given to the current 3d mesh
  !> @param [inout] twod_mesh_id The identifier given to the current 2d mesh
  !> @param [inout] chi A size 3 array of fields holding the coordinates of the mesh
  subroutine initialise_infrastructure(comm,         &
                                       filename,     &
                                       program_name, &
                                       mesh_id,      &
                                       twod_mesh_id, &
                                       chi)

    use logging_config_mod, only: run_log_level,          &
                                  key_from_run_log_level, &
                                  RUN_LOG_LEVEL_ERROR,    &
                                  RUN_LOG_LEVEL_INFO,     &
                                  RUN_LOG_LEVEL_DEBUG,    &
                                  RUN_LOG_LEVEL_TRACE,    &
                                  RUN_LOG_LEVEL_WARNING



    implicit none

    integer(i_def),     intent(inout) :: comm
    character(*),       intent(in)    :: filename
    character(*),       intent(in)    :: program_name

    integer(i_def),     intent(inout) :: mesh_id, twod_mesh_id
    type(field_type),   intent(inout) :: chi(3)

    character(len=*), parameter :: xios_id   = "lfric_client"
    character(len=*), parameter :: xios_ctx  = "gungho_atm"

    integer(i_def)    :: total_ranks, local_rank
    integer(i_def)    :: dtime
    integer(i_def)    :: ts_init
    integer(i_native) :: log_level

    type(field_type)  :: surface_altitude

    !-------------------------------------------------------------------------
    ! Initialise aspects of the infrastructure
    !-------------------------------------------------------------------------

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
    ! Initialise timers and counters
    !-------------------------------------------------------------------------
    if ( subroutine_timers ) then
      call init_timer()
      call timer('gungho')
    end if

    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    if ( subroutine_counters ) then
      allocate(halo_calls, source=count_type('halo_calls'))
      call halo_calls%counter('gungho')
    end if

    !-------------------------------------------------------------------------
    ! Initialise aspects of the grid
    !-------------------------------------------------------------------------

    allocate( global_mesh_collection, &
              source = global_mesh_collection_type() )

    ! Get the rank information from the virtual machine
    total_ranks = get_comm_size()
    local_rank  = get_comm_rank()

    ! Create the mesh
    call init_mesh(local_rank, total_ranks, mesh_id, twod_mesh_id)

    ! Create FEM specifics (function spaces and chi field)
    call init_fem(mesh_id, chi)

    ! Full global meshes no longer required, so reclaim
    ! the memory from global_mesh_collection
    deallocate(global_mesh_collection)

    !-------------------------------------------------------------------------
    ! Initialise aspects of output
    !-------------------------------------------------------------------------

    ! If using XIOS for diagnostic output or checkpointing, then set up XIOS
    ! domain and context
    if ( use_xios_io ) then

      dtime = int(dt)

      call xios_domain_init( xios_ctx,     &
                             comm,         &
                             dtime,        &
                             mesh_id,      &
                             twod_mesh_id, &
                             chi)

      ts_init = max( (timestep_start - 1), 0 ) ! 0 or t previous.

      if (ts_init == 0) then
        ! Make sure XIOS calendar is set to timestep 1 as it starts there
        ! not timestep 0.
        call xios_update_calendar(1)
      end if

    end if

    ! Set up surface altitude field - this will be used to generate orography
    ! for models with global land mass included (i.e GA/GL)
    call init_altitude( twod_mesh_id, surface_altitude )

    ! Set up orography after XIOS to enable read-in of surface altitude from
    ! ancillary file - assignment of orography from altitude to be done by #2023
    call assign_orography_field(chi, mesh_id)


    !-------------------------------------------------------------------------
    ! Setup constants
    !-------------------------------------------------------------------------

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_runtime_constants(mesh_id, twod_mesh_id, chi)

#ifdef UM_PHYSICS
    ! Set derived planet constants and presets
    call set_planet_constants()
#endif

  end subroutine initialise_infrastructure

  !---------------------------------------------------------------------------
  !> @brief Initialises the gungho application
  !> @param[in] mesh_id The identifier of the primary mesh
  !> @param[inout] model_data The working data set for the model run
  !> @param[in] timestep_start number of timestep at which this run started
  subroutine initialise_model( mesh_id, &
                               model_data, &
                               timestep_start )
    implicit none

    integer(i_def),                  intent(in)    :: mesh_id
    type( model_data_type ), target, intent(inout) :: model_data
    integer(i_def),                  intent(in)    :: timestep_start

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_collection_type ), pointer :: diagnostic_fields => null()
    type( field_type ),            pointer :: mr(:) => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()

    ! Get pointers to field collections for use downstream
    prognostic_fields => model_data%prognostic_fields
    diagnostic_fields => model_data%diagnostic_fields
    mr => model_data%mr

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')

    if (write_minmax_tseries) then
      call minmax_tseries_init('u', mesh_id)
      call minmax_tseries(u, 'u', mesh_id)
    end if

    if ( transport_only ) then

      select case( scheme )
        case ( scheme_method_of_lines )
          ! Initialise and output initial conditions for first timestep
          call runge_kutta_init()
          call rk_transport_init( mesh_id, u, rho, theta)
        case default
          call log_event("Gungho: Incorrect transport option chosen, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      end select
    else
      select case( method )
        case( method_semi_implicit )  ! Semi-Implicit
          ! Initialise and output initial conditions for first timestep
          call runge_kutta_init()
          call iter_alg_init(mesh_id, u, rho, theta, exner, mr)
          if ( write_diag ) &
           call conservation_algorithm(timestep_start, rho, u, theta, exner)
        case( method_rk )             ! RK
          ! Initialise and output initial conditions for first timestep
          call runge_kutta_init()
          call rk_alg_init(mesh_id, u, rho, theta, exner)
          if ( write_diag ) &
           call conservation_algorithm(timestep_start, rho, u, theta, exner)
        case default
          call log_event("Gungho: Incorrect time stepping option chosen, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      end select

    end if

  end subroutine initialise_model

  !> @brief Finalises infrastructure and constants used by the model
  subroutine finalise_infrastructure(program_name)

    implicit none

    character(*), intent(in) :: program_name

    !-------------------------------------------------------------------------
    ! Finalise constants
    !-------------------------------------------------------------------------

    call final_runtime_constants()

    !-------------------------------------------------------------------------
    ! Finalise timers and counters
    !-------------------------------------------------------------------------

    if ( subroutine_timers ) then
      call timer('gungho')
      call output_timer()
    end if

    if ( subroutine_counters ) then
      call halo_calls%counter('gungho')
      call halo_calls%output_counters()
    end if

    !-------------------------------------------------------------------------
    ! Finalise aspects of output
    !-------------------------------------------------------------------------

    ! Finalise XIOS context if we used it for diagnostic output or checkpointing
    if ( use_xios_io ) then
      call xios_context_finalize()
    end if

    !-------------------------------------------------------------------------
    ! Finalise aspects of the grid
    !-------------------------------------------------------------------------

    call final_mesh()
    call final_fem()

    !-------------------------------------------------------------------------
    ! Final logging before infrastructure is destroyed
    !-------------------------------------------------------------------------

    call log_event( program_name//' completed.', LOG_LEVEL_ALWAYS )

    !-------------------------------------------------------------------------
    ! Finalise infrastructure
    !-------------------------------------------------------------------------

    ! Finalise XIOS
    call xios_finalize()

    ! Finalise namelist configurations
    call final_configuration()

    ! Finalise YAXT
    call xt_finalize()

    ! Finalise the logging system
    call finalise_logging()

  end subroutine finalise_infrastructure

  !---------------------------------------------------------------------------
  !> @brief Finalise the gungho application
  !> @param[in] mesh_id The identifier of the primary mesh
  !> @param[inout] model_data The working data set for the model run
  !> @param[in] program_name An identifier given to the model begin run
  subroutine finalise_model( mesh_id, &
                             model_data, &
                             program_name )

    implicit none

    integer(i_def),                  intent(in)    :: mesh_id
    type( model_data_type ), target, intent(inout) :: model_data
    character(*),                    intent(in)    :: program_name

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_collection_type ), pointer :: diagnostic_fields => null()
    type( field_type ),            pointer :: mr(:) => null()
    type( field_collection_type ), pointer :: fd_fields => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()

    ! Pointer for tstar_2d to allow write to dump
    type( field_type ), pointer   :: tstar_2d => null()
    ! Pointer for setting I/O handlers on fields
    procedure(write_interface), pointer :: tmp_write_ptr => null()

    ! Get pointers to field collections for use downstream
    prognostic_fields => model_data%prognostic_fields
    diagnostic_fields => model_data%diagnostic_fields
    mr => model_data%mr
    fd_fields => model_data%fd_fields

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')

    ! Log fields
    call rho%log_field(   LOG_LEVEL_DEBUG, 'rho' )
    call theta%log_field( LOG_LEVEL_DEBUG, 'theta' )
    call exner%log_field( LOG_LEVEL_DEBUG, 'exner' )
    call u%log_field(     LOG_LEVEL_DEBUG, 'u' )

    ! Write checksums to file
    if (use_moisture) then
      call checksum_alg(program_name, rho, 'rho', theta, 'theta', u, 'u', &
         field_bundle=mr, bundle_name='mr')
    else
      call checksum_alg(program_name, rho, 'rho', theta, 'theta', u, 'u')
    end if

    ! Call timestep finalizers
    if ( transport_only .and. scheme == scheme_method_of_lines) then
      call rk_transport_final( rho, theta)
    end if

    if(write_minmax_tseries) call minmax_tseries_final(mesh_id)

    call runge_kutta_final()

    if ( .not. transport_only ) then
      if ( method == method_semi_implicit ) call iter_alg_final()
      if ( method == method_rk )            call rk_alg_final()
    end if

  end subroutine finalise_model

end module gungho_model_mod
