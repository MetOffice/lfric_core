!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Handles initialisation and finalisation of infrastructure, constants
! and the gungho model simulations
module gungho_model_mod

  use assign_orography_field_mod, only : assign_orography_field
  use base_mesh_config_mod,       only : prime_mesh_name
  use checksum_alg_mod,           only : checksum_alg
  use clock_mod,                  only : clock_type
  use create_fem_mod,             only : init_fem, final_fem
  use create_mesh_mod,            only : init_mesh, final_mesh
  use configuration_mod,          only : final_configuration
  use conservation_algorithm_mod, only : conservation_algorithm
  use constants_mod,              only : i_def, i_native, &
                                         PRECISION_REAL
  use convert_to_upper_mod,       only : convert_to_upper
  use count_mod,                  only : count_type, halo_calls
  use derived_config_mod,         only : set_derived_config
  use field_mod,                  only : field_type
  use field_parent_mod,           only : write_interface
  use field_collection_mod,       only : field_collection_type, &
                                         field_collection_iterator_type
  use formulation_config_mod,     only : transport_only, &
                                         use_moisture,   &
                                         use_physics
  use global_mesh_collection_mod, only : global_mesh_collection, &
                                         global_mesh_collection_type
  use gungho_mod,                 only : load_configuration
  use gungho_model_data_mod,      only : model_data_type
  use init_altitude_mod,          only : init_altitude
  use init_clock_mod,             only : initialise_clock
  use init_altitude_mod,          only : init_altitude
  use io_mod,                     only : initialise_xios
  use io_config_mod,              only : subroutine_timers,       &
                                         subroutine_counters,     &
                                         use_xios_io,             &
                                         write_conservation_diag, &
                                         write_dump,              &
                                         write_minmax_tseries
  use semi_implicit_timestep_alg_mod, &
                                  only : semi_implicit_alg_init, &
                                         semi_implicit_alg_final
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
  use moisture_conservation_alg_mod, &
                                  only : moisture_conservation_alg
  use minmax_tseries_mod,         only : minmax_tseries,      &
                                         minmax_tseries_init, &
                                         minmax_tseries_final
  use mpi_mod,                    only : store_comm,    &
                                         get_comm_size, &
                                         get_comm_rank
  use mr_indices_mod,             only : nummr
  use rk_alg_timestep_mod,        only : rk_alg_init, &
                                         rk_alg_final
  use rk_transport_mod,           only : rk_transport_init, &
                                         rk_transport_final
  use runge_kutta_init_mod,       only : runge_kutta_init, &
                                         runge_kutta_final
  use runtime_constants_mod,      only : create_runtime_constants, &
                                         final_runtime_constants
  use section_choice_config_mod,  only : radiation,         &
                                         radiation_socrates,&
                                         surface, surface_jules
  use time_config_mod,            only : timestep_start
  use timer_mod,                  only : timer, output_timer, init_timer
  use timestepping_config_mod,    only : method,               &
                                         method_semi_implicit, &
                                         method_rk
  use transport_config_mod,       only : scheme, &
                                         scheme_method_of_lines
  use write_methods_mod,          only : write_state, &
                                         write_field_single_face
  use xios,                       only : xios_context_finalize, &
                                         xios_update_calendar
  use yaxt,                       only : xt_initialize, xt_finalize

  use linked_list_mod,            only : linked_list_type
  use gungho_setup_io_mod,        only : init_gungho_files

#ifdef UM_PHYSICS
  use jules_control_init_mod,     only : jules_control_init
  use jules_physics_init_mod,     only : jules_physics_init
  use planet_constants_mod,       only : set_planet_constants
  use socrates_init_mod,          only : socrates_init
  use um_control_init_mod,        only : um_control_init
  use um_physics_init_mod,        only : um_physics_init
#endif

  implicit none

  private
  public initialise_infrastructure, &
         initialise_model,          &
         finalise_infrastructure,   &
         finalise_model

  contains


  !> @brief Initialises the infrastructure and sets up constants used by the model
  !> @param [in,out] communicator The MPI communicator for use within the model
  !>                     (not XIOS' communicator)
  !> @param [in]    filename The name of the configuration namelist file
  !> @param [in]    program_name An identifier given to the model begin run
  !> @param[out]     clock Model time
  !> @param [in,out] mesh_id The identifier given to the current 3d mesh
  !> @param [in,out] twod_mesh_id The identifier given to the current 2d mesh
  !> @param [in,out] shifted_mesh_id The identifier given to the vertically shifted 3d mesh
  !> @param [in,out] double_level_mesh_id The identifier given to the double-level 3d mesh
  !> the double-level mesh
  subroutine initialise_infrastructure(communicator,         &
                                       filename,             &
                                       program_name,         &
                                       clock,                &
                                       mesh_id,              &
                                       twod_mesh_id,         &
                                       shifted_mesh_id,      &
                                       double_level_mesh_id)

    use logging_config_mod, only: run_log_level,          &
                                  key_from_run_log_level, &
                                  RUN_LOG_LEVEL_ERROR,    &
                                  RUN_LOG_LEVEL_INFO,     &
                                  RUN_LOG_LEVEL_DEBUG,    &
                                  RUN_LOG_LEVEL_TRACE,    &
                                  RUN_LOG_LEVEL_WARNING

    implicit none

    integer(i_native),  intent(in)    :: communicator
    character(*),       intent(in)    :: filename
    character(*),       intent(in)    :: program_name
    class(clock_type), intent(out), allocatable :: clock

    integer(i_def),     intent(inout) :: mesh_id, twod_mesh_id
    integer(i_def),     intent(inout) :: double_level_mesh_id, shifted_mesh_id

    character(len=*), parameter :: xios_ctx  = "gungho_atm"

    integer(i_def)    :: total_ranks, local_rank
    integer(i_native) :: log_level

    type(field_type)  :: surface_altitude

    type(field_type), target :: chi_xyz(3)
    type(field_type), target :: chi_sph(3)
    type(field_type), target :: panel_id
    type(field_type), target :: shifted_chi_xyz(3)
    type(field_type), target :: shifted_chi_sph(3)
    type(field_type), target :: double_level_chi_xyz(3)
    type(field_type), target :: double_level_chi_sph(3)

    type(linked_list_type) :: files_list

    integer(i_def), allocatable :: multigrid_mesh_ids(:)
    integer(i_def), allocatable :: multigrid_2d_mesh_ids(:)

    !-------------------------------------------------------------------------
    ! Initialise aspects of the infrastructure
    !-------------------------------------------------------------------------

    !Store the MPI communicator for later use
    call store_comm( communicator )

    ! Initialise YAXT
    call xt_initialize( communicator )

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

    call initialise_clock( clock )

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
      call timer(program_name)
    end if

    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    if ( subroutine_counters ) then
      allocate(halo_calls, source=count_type('halo_calls'))
      call halo_calls%counter(program_name)
    end if

    !-------------------------------------------------------------------------
    ! Initialise aspects of the grid
    !-------------------------------------------------------------------------

    allocate( global_mesh_collection, &
              source = global_mesh_collection_type() )

    ! Create the mesh
    call init_mesh( local_rank, total_ranks, mesh_id,             &
                    twod_mesh_id          = twod_mesh_id,         &
                    shifted_mesh_id       = shifted_mesh_id,      &
                    double_level_mesh_id  = double_level_mesh_id, &
                    multigrid_mesh_ids    = multigrid_mesh_ids,   &
                    multigrid_2D_mesh_ids = multigrid_2D_mesh_ids )

    call init_fem( mesh_id, chi_xyz, chi_sph, panel_id,          &
                   shifted_mesh_id       = shifted_mesh_id,      &
                   shifted_chi_xyz       = shifted_chi_xyz,      &
                   shifted_chi_sph       = shifted_chi_sph,      &
                   double_level_mesh_id  = double_level_mesh_id, &
                   double_level_chi_xyz  = double_level_chi_xyz, &
                   double_level_chi_sph  = double_level_chi_sph, &
                   multigrid_mesh_ids    = multigrid_mesh_ids,   &
                   multigrid_2D_mesh_ids = multigrid_2D_mesh_ids )


    ! Full global meshes no longer required, so reclaim
    ! the memory from global_mesh_collection
    if (allocated(global_mesh_collection)) deallocate(global_mesh_collection)

    !-------------------------------------------------------------------------
    ! Initialise aspects of output
    !-------------------------------------------------------------------------

    ! If using XIOS for diagnostic output or checkpointing, then set up XIOS
    ! domain and context
    if ( use_xios_io ) then

      call init_gungho_files(files_list, clock)

      call initialise_xios( xios_ctx,     &
                            communicator, &
                            clock,        &
                            mesh_id,      &
                            twod_mesh_id, &
                            chi_xyz,      &
                            files_list )

      if (clock%is_initialisation()) then
        ! Make sure XIOS calendar is set to timestep 1 as it starts there
        ! not timestep 0.
        call xios_update_calendar(1)
      end if

    end if

    ! Set up surface altitude field - this will be used to generate orography
    ! for models with global land mass included (i.e GAL)
    call init_altitude( twod_mesh_id, surface_altitude )

    ! Assignment of orography from surface_altitude
    call assign_orography_field(chi_xyz, panel_id,                      &
                                mesh_id, surface_altitude,              &
                                spherical_coords=.false.)
    call assign_orography_field(shifted_chi_xyz, panel_id,              &
                                shifted_mesh_id, surface_altitude,      &
                                spherical_coords=.false.)
    call assign_orography_field(double_level_chi_xyz, panel_id,         &
                                double_level_mesh_id, surface_altitude, &
                                spherical_coords=.false.)
    call assign_orography_field(chi_sph, panel_id,                      &
                                mesh_id, surface_altitude,              &
                                spherical_coords=.true.)
    call assign_orography_field(shifted_chi_sph, panel_id,              &
                                shifted_mesh_id, surface_altitude,      &
                                spherical_coords=.true.)
    call assign_orography_field(double_level_chi_sph, panel_id,         &
                                double_level_mesh_id, surface_altitude, &
                                spherical_coords=.true.)

    !-------------------------------------------------------------------------
    ! Setup constants
    !-------------------------------------------------------------------------

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_runtime_constants(mesh_id, twod_mesh_id, chi_xyz, chi_sph,    &
                                  panel_id, shifted_mesh_id, shifted_chi_xyz, &
                                  shifted_chi_sph, double_level_mesh_id,      &
                                  double_level_chi_xyz, double_level_chi_sph, &
                                  surface_altitude)

#ifdef UM_PHYSICS
    ! Set derived planet constants and presets
    call set_planet_constants()

    if ( use_physics ) then
      if (radiation == radiation_socrates) then
        ! Initialisation for the Socrates radiation scheme
        call socrates_init()
      end if
      ! Initialisation of UM high-level variables
      call um_control_init(mesh_id)
      ! Initialisation of UM physics variables
      call um_physics_init()
      ! Initialisation of Jules high-level variables
      call jules_control_init()
      if (surface == surface_jules) then
        ! Initialisation of Jules physics variables
        call jules_physics_init()
      end if
    end if
#endif

  end subroutine initialise_infrastructure

  !---------------------------------------------------------------------------
  !> @brief Initialises the gungho application
  !> @param[in] clock Model time
  !> @param[in] mesh_id The identifier of the primary mesh
  !> @param[inout] model_data The working data set for the model run
  subroutine initialise_model( clock,   &
                               mesh_id, &
                               model_data )
    implicit none

    class(clock_type),               intent(in)    :: clock
    integer(i_def),                  intent(in)    :: mesh_id
    type( model_data_type ), target, intent(inout) :: model_data

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
          if ( use_moisture ) then
            call rk_transport_init( rho, theta, mr )
          else
            call rk_transport_init( rho, theta )
          end if
        case default
          call log_event("Gungho: Incorrect transport option chosen, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      end select
    else
      select case( method )
        case( method_semi_implicit )  ! Semi-Implicit
          ! Initialise and output initial conditions for first timestep
          call runge_kutta_init()
          call semi_implicit_alg_init(mesh_id, u, rho, theta, exner, mr)
          if ( write_conservation_diag ) then
           call conservation_algorithm( clock%get_step(), &
                                        rho,              &
                                        u,                &
                                        theta,            &
                                        exner )
           if ( use_moisture ) &
             call moisture_conservation_alg( clock%get_step(), &
                                             rho,              &
                                             mr,               &
                                             'Before timestep' )
          end if
        case( method_rk )             ! RK
          ! Initialise and output initial conditions for first timestep
          call runge_kutta_init()
          call rk_alg_init(mesh_id, u, rho, theta, exner)
          if ( write_conservation_diag ) then
           call conservation_algorithm( clock%get_step(), &
                                        rho,              &
                                        u,                &
                                        theta,            &
                                        exner )
           if ( use_moisture ) &
             call moisture_conservation_alg( clock%get_step(), &
                                             rho,              &
                                             mr,               &
                                             'Before timestep' )
          end if
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
      call timer(program_name)
      call output_timer()
    end if

    if ( subroutine_counters ) then
      call halo_calls%counter(program_name)
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
      if ( use_moisture ) then
        call rk_transport_final( rho, theta, mr )
      else
        call rk_transport_final( rho, theta )
      end if
    end if

    if(write_minmax_tseries) call minmax_tseries_final(mesh_id)

    call runge_kutta_final()

    if ( .not. transport_only ) then
      if ( method == method_semi_implicit ) call semi_implicit_alg_final()
      if ( method == method_rk )            call rk_alg_final()
    end if

  end subroutine finalise_model

end module gungho_model_mod
