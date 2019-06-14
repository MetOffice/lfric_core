!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>@brief Drives the execution of the GungHo model.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module gungho_driver_mod

  use checksum_alg_mod,           only : checksum_alg
  use section_choice_config_mod,  only : cloud, cloud_um, cloud_none
  use conservation_algorithm_mod, only : conservation_algorithm
  use constants_mod,              only : i_def, imdi, str_def, str_short
  use derived_config_mod,         only : set_derived_config
  use time_config_mod,            only : timestep_start, &
                                         timestep_end
  use diagnostics_io_mod,         only : write_scalar_diagnostic,     &
                                         write_vector_diagnostic
  use diagnostics_mod,            only : write_divergence_diagnostic, &
                                         write_density_diagnostic,    &
                                         write_hydbal_diagnostic
  use yaxt,                       only : xt_initialize, xt_finalize
  use field_mod,                  only : field_type,     &
                                         write_interface
  use formulation_config_mod,     only : transport_only, &
                                         use_moisture,   &
                                         use_physics
  use function_space_collection_mod, &
                                  only : function_space_collection
  use fs_continuity_mod,          only : W3, Wtheta
  use field_collection_mod,       only : field_collection_type, &
                                         field_collection_iterator_type
  use global_mesh_collection_mod, only : global_mesh_collection, &
                                         global_mesh_collection_type
  use configuration_mod,          only : final_configuration
  use gungho_mod,                 only : load_configuration
  use init_fem_mod,               only : init_fem
  use create_gungho_prognostics_mod, &
                                  only : create_gungho_prognostics
  use create_physics_prognostics_mod, &
                                  only : create_physics_prognostics
  use init_gungho_prognostics_alg_mod, &
                                  only : init_gungho_prognostics_alg
  use moist_dyn_factors_alg_mod,  only : moist_dyn_factors_alg
  use initial_cloud_alg_mod,      only : initial_cloud_alg
  use init_physics_incs_alg_mod,  only : init_physics_incs_alg
  use init_physics_prognostics_alg_mod, &
                                  only : init_physics_prognostics_alg
  use init_physics_incs_alg_mod,  only : init_physics_incs_alg
  use init_mesh_mod,              only : init_mesh
  use runtime_constants_mod,      only : create_runtime_constants
  use io_mod,                     only : xios_domain_init, &
                                         ts_fname,         &
                                         write_checkpoint, &
                                         read_checkpoint,  &
                                         write_state,      &
                                         dump_write_xios
  use io_config_mod,              only : write_diag,           &
                                         diagnostic_frequency, &
                                         use_xios_io,          &
                                         nodal_output_on_w3,   &
                                         checkpoint_write,     &
                                         checkpoint_read,      &
                                         write_dump,           &
                                         write_minmax_tseries, &
                                         subroutine_timers,    &
                                         subroutine_counters
  use files_config_mod,           only:  checkpoint_stem_name
  use initialization_config_mod,  only : init_option,                &
                                         init_option_analytic,       & 
                                         init_option_fd_start_dump,  &
                                         init_option_checkpoint_dump,&
                                         init_option_fe_start_dump,  &
                                         ancil_option,               &
                                         ancil_option_none,          &
                                         ancil_option_analytic,      & 
                                         ancil_option_aquaplanet
  use create_fd_prognostics_mod,  &
                                  only : create_fd_prognostics
  use init_fd_prognostics_mod,    &
                                  only : init_fd_prognostics_dump
  use init_ancils_mod,            only : init_analytic_ancils, &
                                         init_aquaplanet_ancils
  use map_fd_to_prognostics_mod,  only : map_fd_to_prognostics, &
                                         hydrostatic_balance
  use iter_timestep_alg_mod,      only : iter_alg_init, &
                                         iter_alg_step, &
                                         iter_alg_final
  use log_mod,                    only : log_event,          &
                                         log_set_level,      &
                                         log_scratch_space,  &
                                         initialise_logging, &
                                         finalise_logging,   &
                                         LOG_LEVEL_ERROR,    &
                                         LOG_LEVEL_INFO,     &
                                         LOG_LEVEL_DEBUG,    &
                                         LOG_LEVEL_TRACE
  use mesh_collection_mod,        only : mesh_collection
  use mod_wait
  use minmax_tseries_mod,         only : minmax_tseries, &
                                         minmax_tseries_init, &
                                         minmax_tseries_final
  use mr_indices_mod,             only : nummr, mr_names
  use moist_dyn_mod,              only : num_moist_factors, gas_law, &
                                         total_mass, water
  use rk_alg_timestep_mod,        only : rk_alg_init, &
                                         rk_alg_step, &
                                         rk_alg_final
  use rk_transport_mod,           only : rk_transport_init, &
                                         rk_transport_step, &
                                         rk_transport_final
  use runge_kutta_init_mod,       only : runge_kutta_init, runge_kutta_final
  use runtime_constants_mod,      only : final_runtime_constants
  use timer_mod,                  only : timer, output_timer
  use timestepping_config_mod,    only : method, dt, &
                                         method_semi_implicit, &
                                         method_rk
  use transport_config_mod,       only : scheme, &
                                         scheme_method_of_lines
  use xios
  use count_mod,                  only : count_type, halo_calls
  use mpi_mod,                    only : initialise_comm, store_comm, &
                                         finalise_comm, &
                                         get_comm_size, get_comm_rank


  implicit none

  private
  public initialise, run, finalise

  ! Field collections
  type( field_collection_type ) :: prognostic_fields
  type( field_collection_type ) :: derived_fields
  type( field_collection_type ) :: cloud_fields
  type( field_collection_type ) :: twod_fields
  type( field_collection_type ) :: physics_incs
  type( field_collection_type ) :: fd_fields

  ! Prognostic fields (pointers to fields in a collection)
  type( field_type ), pointer   :: u => null()
  type( field_type ), pointer   :: rho => null()
  type( field_type ), pointer   :: theta => null()
  type( field_type ), pointer   :: exner => null()

  ! Pointer for tstar_2d to allow write to dump
  type( field_type ), pointer   :: tstar_2d => null()

  ! Auxiliary prognostic fields
  ! Moisture mixing ratios
  type( field_type ) :: mr(nummr)

  ! Diagnostic fields
  type( field_type ) :: xi
  ! Moist dynamics
  type( field_type ) :: moist_dyn(num_moist_factors)

  ! Iterator for field collection
  type(field_collection_iterator_type)  :: iterator

  ! A pointer used for retrieving fields from collections
  ! when iterating over them
  type( field_type ), pointer :: field_ptr  => null()

  ! Pointer for setting I/O handlers on fields
  procedure(write_interface), pointer    :: tmp_write_ptr => null()


  ! Coordinate field
  type(field_type), target :: chi(3)

  integer(i_def) :: mesh_id      = imdi
  integer(i_def) :: twod_mesh_id = imdi

  character(str_def) :: name

  integer(i_def) :: i, fs

  integer(i_def) :: prognostic_init_choice, ancil_choice

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Sets up required state in preparation for run.
  !>@param[in] filename Name of the file containing the desired configuration 
  subroutine initialise( filename )

    implicit none

    character(:), intent(in), allocatable :: filename

    character(len=*), parameter :: xios_id   = "lfric_client"
    character(len=*), parameter :: xios_ctx  = "gungho_atm"

    integer(i_def) :: total_ranks, local_rank
    integer(i_def) :: comm = -999
    integer(i_def) :: ts_init, dtime

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

    call initialise_logging(local_rank, total_ranks, 'gungho')

    ! First call to log event must occur after the call to initialise_logging
    call log_event( 'gungho running...', LOG_LEVEL_INFO )

    call load_configuration( filename )
    call set_derived_config( .true. )

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    if ( subroutine_timers ) call timer('gungho')

    if ( subroutine_counters ) then
      allocate(halo_calls, source=count_type('halo_calls'))
      call halo_calls%counter('gungho')
    end if

    allocate( global_mesh_collection, &
              source = global_mesh_collection_type() )

    ! Create the mesh
    call init_mesh(local_rank, total_ranks, mesh_id, twod_mesh_id)

    ! Create FEM specifics (function spaces and chi field)
    call init_fem(mesh_id, chi)

    ! Full global meshes no longer required, so reclaim
    ! the memory from global_mesh_collection
    write(log_scratch_space,'(A)') &
        "Purging global mesh collection."
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    call global_mesh_collection%clear()
    deallocate(global_mesh_collection)

    !-------------------------------------------------------------------------
    ! IO init
    !-------------------------------------------------------------------------

    ! If using XIOS for diagnostic output, checkpointing or dumping, then set up
    ! XIOS domain and context

    if ( use_xios_io ) then

      dtime = int(dt)

      call xios_domain_init( xios_ctx,     &
                             comm,         &
                             dtime,        &
                             mesh_id,      &
                             twod_mesh_id, &
                             chi)

    end if

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_runtime_constants(mesh_id, twod_mesh_id, chi)

    ! Create gungho prognostics and auxilliary (diagnostic) fields
    call create_gungho_prognostics( mesh_id, prognostic_fields, mr, moist_dyn, xi )

    ! Create prognostics used by physics
    if (use_physics) then
      call create_physics_prognostics( mesh_id, twod_mesh_id, &
                                       prognostic_fields, &
                                       derived_fields, cloud_fields, &
                                       twod_fields, physics_incs )
    end if

    !-------------------------------------------------------------------------
    ! Select how to initialize model prognostic fields
    !-------------------------------------------------------------------------

    ! This way of setting up the initialisation options is not ideal, but
    ! pragmatic for now and avoids extra namelist changes. It should be
    ! reviewed in the next round of driver layer refactoring
 
    ! Get the specified namelist options for prognostic initialisation
    prognostic_init_choice = init_option
    ancil_choice = ancil_option

    ! If checkpoint reading has been specified then override these options
    if (checkpoint_read) then
      prognostic_init_choice = init_option_checkpoint_dump
      ancil_choice = ancil_option_none
    end if

    ! Initialise prognostic fields appropriately
    select case ( prognostic_init_choice )

      case ( init_option_analytic )

        ! Initialise prognostics analytically according to
        ! namelist options

        call init_gungho_prognostics_alg(prognostic_fields, mr, moist_dyn, xi)

        if (use_physics) then
          call init_physics_prognostics_alg(derived_fields, &
                                            cloud_fields, &
                                            twod_fields, &
                                            physics_incs)
        end if

      case ( init_option_checkpoint_dump )
 
        ! Initialize prognostics using a checkpoint file
        ! from a previous run

        call read_checkpoint(prognostic_fields, timestep_start-1)

        ! Update factors for moist dynamics
        call moist_dyn_factors_alg(moist_dyn, mr)

        ! if no cloud scheme, reset cloud variables
        if (use_physics) then
          if ( cloud == cloud_none ) then  
            call initial_cloud_alg(cloud_fields)
          end if
          ! Set the increments to 0 initially
          call init_physics_incs_alg(physics_incs)
        end if

      case ( init_option_fd_start_dump )

        if (use_physics) then

          ! Initialise FD prognostic fields from a UM2LFRic dump     

          ! Create FD prognostic fields
          call create_fd_prognostics(mesh_id, fd_fields)                          

          ! Read in from a UM2LFRic dump file
          call init_fd_prognostics_dump(fd_fields)

          ! Set physics increments to 0
          call init_physics_incs_alg(physics_incs)

          ! Populate prognostics from input finite difference fields
          call map_fd_to_prognostics(prognostic_fields, xi, &
                                     mr, moist_dyn,         &
                                     cloud_fields,          &
                                     twod_fields,           &
                                     fd_fields)

        else
          call log_event("Gungho: Prognostic initialisation from an FD dump not valid "// &
                          "if use_physics is .false., stopping program! ",LOG_LEVEL_ERROR)

        end if   

      case ( init_option_fe_start_dump )
        ! Initialise FE prognostic fields from an FE dump
        ! Not yet supported
        call log_event("Gungho: Prognostic initialisation from an FE dump not yet supported, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      case default
        ! No valid initialisation option selected
        call log_event("Gungho: No valid prognostic initialisation option selected, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)

    end select
    
    ! Assuming this is only relevant for physics runs at the moment
    if (use_physics) then

      ! Initialise ancillary fields
      select case ( ancil_choice )
        case ( ancil_option_none )
          call log_event( "Gungho: No ancillaries to be read for this run.", LOG_LEVEL_INFO )
        case ( ancil_option_aquaplanet )
          call log_event( "Gungho: Reading ancillaries from aquaplanet dump ", LOG_LEVEL_INFO )
          call init_aquaplanet_ancils(twod_fields)
        case ( ancil_option_analytic )
          call log_event( "Gungho: Setting ancillaries from analytic representation ", LOG_LEVEL_INFO )
          call init_analytic_ancils(twod_fields)
        case default
          ! No valid ancil option selected
          call log_event("Gungho: No valid ancillary initialisation option selected, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      end select

    end if

    ! Get pointers to fields in the prognostic fields collection
    ! for use downstream
    theta => prognostic_fields%get_field('theta')
    u => prognostic_fields%get_field('u')
    rho => prognostic_fields%get_field('rho')
    exner => prognostic_fields%get_field('exner')


    ! Initial output
    ! We only want these once at the beginning of a run
    ts_init = max( (timestep_start - 1), 0 ) ! 0 or t previous.

    if (ts_init == 0) then

      if ( use_xios_io ) then

        ! Need to ensure calendar is initialised here as XIOS has no concept of timestep 0
        call xios_update_calendar(ts_init + 1)

      end if

      if ( write_diag ) then

        ! Calculation and output of diagnostics

        ! Prognostic Scalar fields
        call write_scalar_diagnostic('rho', rho, ts_init, mesh_id, nodal_output_on_w3)
        call write_scalar_diagnostic('theta', theta, ts_init, mesh_id, nodal_output_on_w3)
        call write_scalar_diagnostic('exner', exner, ts_init, mesh_id, nodal_output_on_w3)

        ! Prognostic Vector fields
        call write_vector_diagnostic('u', u, ts_init, mesh_id, nodal_output_on_w3)
        call write_vector_diagnostic('xi', xi, ts_init, mesh_id, nodal_output_on_w3)

        ! Moisture fields
        if (use_moisture) then
          do i=1,nummr
            call write_scalar_diagnostic( trim(mr_names(i)), mr(i), &
                                          ts_init, mesh_id, nodal_output_on_w3 )
          end do
        end if

        ! Cloud fields
        if (use_physics .and. cloud == cloud_um) then

          iterator = cloud_fields%get_iterator()
          do
            if ( .not.iterator%has_next() ) exit
              field_ptr => iterator%next()
              name = trim(adjustl( field_ptr%get_name() ))
              call write_scalar_diagnostic( trim(name), field_ptr, &
                                            ts_init, mesh_id, nodal_output_on_w3 )
          end do
          field_ptr => null()

        end if

        ! Other derived diagnostics with special pre-processing
        call write_divergence_diagnostic(u, ts_init, mesh_id)
        call write_hydbal_diagnostic(theta, moist_dyn, exner, mesh_id)

      end if

      if (write_minmax_tseries) then
        call minmax_tseries_init('u', mesh_id)
        call minmax_tseries(u, 'u', mesh_id)
      end if

    end if

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Timesteps the model, calling the desired timestepping algorithm based
  !upon the configuration
  subroutine run()

    implicit none

    integer(i_def) :: timestep

    do timestep = timestep_start, timestep_end

      ! Update XIOS calendar if we are using it for diagnostic output or checkpoint
      if ( use_xios_io ) then
        call log_event( "Gungho: Updating XIOS timestep", LOG_LEVEL_INFO )
        call xios_update_calendar(timestep)
      end if

      write( log_scratch_space, '("/", A, "\ ")' ) repeat( "*", 76 )
      call log_event( log_scratch_space, LOG_LEVEL_TRACE )
      write( log_scratch_space, '(A,I0)' ) 'Start of timestep ', timestep
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      if ( transport_only ) then

        select case( scheme )
          case ( scheme_method_of_lines )
            if (timestep == timestep_start) then
              ! Initialise and output initial conditions on first timestep
              call runge_kutta_init()
              call rk_transport_init( mesh_id, u, rho, theta)
            end if
            call rk_transport_step( u, rho, theta)
          case default
          call log_event("Gungho: Incorrect transport option chosen, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
          stop
        end select
        call write_density_diagnostic(rho, timestep)
        if ( write_diag ) call conservation_algorithm(timestep, rho, u, theta, exner, xi)
      else
        select case( method )
          case( method_semi_implicit )  ! Semi-Implicit
            ! Initialise and output initial conditions on first timestep
            if (timestep == timestep_start) then
              call runge_kutta_init()
              call iter_alg_init(mesh_id, u, rho, theta, exner, mr, &
                                 twod_fields)
              if ( write_diag ) call conservation_algorithm(timestep, rho, u, theta, exner, xi)
            end if
            call iter_alg_step(u, rho, theta, exner, mr, moist_dyn, xi,      &
                               derived_fields, cloud_fields, twod_fields,    &
                               physics_incs, timestep)

          case( method_rk )             ! RK
            ! Initialise and output initial conditions on first timestep
            if (timestep == timestep_start) then
              call runge_kutta_init()
              call rk_alg_init(mesh_id, u, rho, theta, exner)
              if ( write_diag ) call conservation_algorithm(timestep, rho, u, theta, exner, xi)
            end if
            call rk_alg_step(u, rho, theta, moist_dyn, exner, xi)
          case default
            call log_event("Gungho: Incorrect time stepping option chosen, "// &
                            "stopping program! ",LOG_LEVEL_ERROR)
            stop
        end select

        if ( write_diag ) call conservation_algorithm(timestep, rho, u, theta, exner, xi)

        if(write_minmax_tseries) call minmax_tseries(u, 'u', mesh_id)

       call u%log_minmax(LOG_LEVEL_INFO, ' u')

      end if

      write( log_scratch_space, '(A,I0)' ) 'End of timestep ', timestep
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      write( log_scratch_space, '("\", A, "/ ")' ) repeat( "*", 76 )
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      ! Use diagnostic output frequency to determine whether to write
      ! diagnostics on this timestep

      if ( ( mod(timestep, diagnostic_frequency) == 0 ) .and. ( write_diag ) ) then

        call log_event("Gungho: writing diagnostic output", LOG_LEVEL_INFO)

        ! Calculation and output of diagnostics

        ! Prognostic Scalar fields
        call write_scalar_diagnostic('rho', rho, timestep, mesh_id, nodal_output_on_w3)
        call write_scalar_diagnostic('theta', theta, timestep, mesh_id, nodal_output_on_w3)
        call write_scalar_diagnostic('exner', exner, timestep, mesh_id, nodal_output_on_w3)

        ! Prognostic Vector fields
        call write_vector_diagnostic('u', u, timestep, mesh_id, nodal_output_on_w3)
        call write_vector_diagnostic('xi', u, timestep, mesh_id, nodal_output_on_w3)

        ! Moisture fields
        if (use_moisture) then
          do i=1,nummr
            call write_scalar_diagnostic( trim(mr_names(i)), mr(i), &
                                          timestep, mesh_id, nodal_output_on_w3 )
          end do
        end if

        ! Cloud fields
        if ( use_physics .and. cloud == cloud_um ) then

          iterator = cloud_fields%get_iterator()
          do
            if ( .not.iterator%has_next() ) exit
              field_ptr => iterator%next()
              name = trim(adjustl( field_ptr%get_name() ))
              call write_scalar_diagnostic( trim(name), field_ptr, &
                                            timestep, mesh_id, nodal_output_on_w3 )
          end do
          field_ptr => null()

        end if

        ! Derived physics fields (only those on W3 or Wtheta)
        if (use_physics) then

          iterator = derived_fields%get_iterator()
          do
            if ( .not.iterator%has_next() ) exit
            field_ptr => iterator%next()
            fs = field_ptr%which_function_space()
            if ( fs == W3 .or. fs == Wtheta ) then
              name = trim(adjustl( field_ptr%get_name() ))
              call write_scalar_diagnostic( trim(name), field_ptr, &
                                            timestep, mesh_id, nodal_output_on_w3 )
            end if
          end do
          field_ptr => null()

        end if

        ! Other derived diagnostics with special pre-processing
        call write_divergence_diagnostic(u, timestep, mesh_id)
        call write_hydbal_diagnostic(theta, moist_dyn, exner, mesh_id)

      end if

    end do ! end ts loop

    call runge_kutta_final()

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tidies up after a run.
  subroutine finalise()

    implicit none


    ! Log fields
    call rho%log_field(   LOG_LEVEL_DEBUG, 'rho' )
    call theta%log_field( LOG_LEVEL_DEBUG, 'theta' )
    call exner%log_field( LOG_LEVEL_DEBUG, 'exner' )
    call u%log_field(     LOG_LEVEL_DEBUG, 'u' )

    ! Write checksums to file
    if (use_moisture) then
      call checksum_alg('gungho', rho, 'rho', theta, 'theta', u, 'u', &
         field_bundle=mr, bundle_name='mr')
    else
      call checksum_alg('gungho', rho, 'rho', theta, 'theta', u, 'u')
    end if

    ! Write checkpoint files if required
    if( checkpoint_write ) then

       call write_checkpoint(prognostic_fields, timestep_end)

    end if

    !===================== Write fields to dump ======================!

    if( write_dump ) then

      ! Current dump writing is only relevant for physics runs at the moment
      if (use_physics) then

        ! For the purposes of dumping from one collection, we add a pointer
        ! to tstar to the fd_prognostics collection

        tmp_write_ptr => dump_write_xios
        tstar_2d => twod_fields%get_field('tstar')
        call tstar_2d%set_write_behaviour(tmp_write_ptr)

        call fd_fields%add_field(tstar_2d)

        call log_event("Gungho: writing FD fields to dump", LOG_LEVEL_INFO)

        ! Write prognostic fields to dump
        call write_state(fd_fields)
      
        nullify(tmp_write_ptr, tstar_2d)

      end if

    end if

    ! Call timestep finalizers
    if ( transport_only .and. scheme == scheme_method_of_lines) then
      call rk_transport_final( rho, theta)
    end if

    call theta%field_final()
    call rho%field_final()
    call exner%field_final()
    call u%field_final()
    call xi%field_final()
    do i=1, nummr
      call mr(i)%field_final()
    end do
    do i = 1, num_moist_factors
      call moist_dyn(i)%field_final()
    end do

    if (use_moisture .and. use_physics) then

      iterator = cloud_fields%get_iterator()
      do
        if ( .not.iterator%has_next() ) exit
        field_ptr => iterator%next()
        call field_ptr%field_final()
      end do
      field_ptr => null()

    endif

    if (use_physics) then

      iterator = twod_fields%get_iterator()
      do
        if ( .not.iterator%has_next() ) exit
        field_ptr => iterator%next()
        call field_ptr%field_final()
      end do
      field_ptr => null()

    endif

    ! Finalise the fd fields if we used them
    if ( prognostic_init_choice == init_option_fd_start_dump) then

      iterator = fd_fields%get_iterator()
      do
        if ( .not.iterator%has_next() ) exit
        field_ptr => iterator%next()
        call field_ptr%field_final()
      end do 
      field_ptr => null()


    end if

    if(write_minmax_tseries) call minmax_tseries_final(mesh_id)

    if ( .not. transport_only ) then
      if ( method == method_semi_implicit ) call iter_alg_final()
      if ( method == method_rk )            call rk_alg_final()
    end if
    call final_runtime_constants()

    if ( subroutine_timers ) then
      call timer('gungho')
      call output_timer()
    end if

    if ( subroutine_counters ) then
      call halo_calls%counter('gungho')
      call halo_calls%output_counters()
    end if

    nullify( theta, rho, u, exner, field_ptr )

    call log_event( 'gungho completed', LOG_LEVEL_INFO )

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

    ! Finalise XIOS context if we used it

    if ( use_xios_io ) then
      call xios_context_finalize()
    end if

    if (allocated(mesh_collection)) then
      call mesh_collection%clear()
      deallocate(mesh_collection)
    end if

    if (allocated(function_space_collection)) then
      call function_space_collection%clear()
      deallocate(function_space_collection)
    end if

    ! Finalise XIOS
    call xios_finalize()

    ! Finalise namelist configurations
    call final_configuration()

    ! Finalise YAXT
    call xt_finalize()

    ! Finalise mpi and release the communicator
    call finalise_comm()

    ! Finalise the logging system
    call finalise_logging()

  end subroutine finalise

end module gungho_driver_mod
