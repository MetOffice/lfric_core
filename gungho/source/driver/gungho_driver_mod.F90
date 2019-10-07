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

  use section_choice_config_mod,  only : cloud, cloud_none
  use constants_mod,              only : i_def, imdi, str_def, str_short, l_def
  use time_config_mod,            only : timestep_start, &
                                         timestep_end
  use field_mod,                  only : field_type
  use gungho_model_data_mod,      only : model_data_type, &
                                         finalise_model_data
  use formulation_config_mod,     only : use_physics
  use function_space_collection_mod, &
                                  only : function_space_collection
  use field_collection_mod,       only : field_collection_type
  use gungho_diagnostics_driver_mod, &
                                  only : gungho_diagnostics_driver
  use gungho_grid_mod,            only : initialise_grid, finalise_grid
  use gungho_infrastructure_mod,  only : initialise_infrastructure, &
                                         finalise_infrastructure
  use gungho_io_mod,              only : initialise_io, &
                                         finalise_io
  use create_gungho_prognostics_mod, &
                                  only : create_gungho_prognostics
  use create_physics_prognostics_mod, &
                                  only : create_physics_prognostics

  use create_fd_prognostics_mod,  &
                                  only : create_fd_prognostics
  use init_gungho_prognostics_alg_mod, &
                                  only : init_gungho_prognostics_alg
  use map_fd_to_prognostics_alg_mod,  only : map_fd_to_prognostics
  use moist_dyn_factors_alg_mod,  only : moist_dyn_factors_alg
  use initial_cloud_alg_mod,      only : initial_cloud_alg
  use init_jules_alg_mod,         only : init_jules_alg
  use init_physics_incs_alg_mod,  only : init_physics_incs_alg
  use init_physics_prognostics_alg_mod, &
                                  only : init_physics_prognostics_alg
  use update_tstar_alg_mod,       only : update_tstar_alg
  use runtime_constants_mod,      only : create_runtime_constants
  use gungho_model_mod,           only : initialise_model, & 
                                         finalise_model
  use io_mod,                     only : write_checkpoint, &
                                         read_checkpoint
  use io_config_mod,              only : write_diag,           &
                                         diagnostic_frequency, &
                                         use_xios_io,          &
                                         nodal_output_on_w3,   &
                                         checkpoint_write,     &
                                         checkpoint_read,      &
                                         subroutine_timers,    &
                                         subroutine_counters
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
  use log_mod,                    only : log_event,          &
                                         log_set_level,      &
                                         log_scratch_space,  &
                                         LOG_LEVEL_ALWAYS,   &
                                         LOG_LEVEL_ERROR,    &
                                         LOG_LEVEL_INFO,     &
                                         LOG_LEVEL_TRACE
  use mesh_collection_mod,        only : mesh_collection
#ifdef UM_PHYSICS
  use planet_constants_mod,       only : set_planet_constants
#endif
  use runtime_constants_mod,      only : final_runtime_constants
  use gungho_step_mod,            only : step
  use timer_mod,                  only : timer, output_timer
  use xios
  use count_mod,                  only : count_type, halo_calls
  use mpi_mod,                    only : initialise_comm, finalise_comm
  use gungho_mod,                 only : program_name

  implicit none

  private
  public initialise, run, finalise

  character(len=*), parameter :: xios_id   = "lfric_client"
  character(len=*), parameter :: xios_ctx  = "gungho_atm"

  
  ! Model run working data set
  type (model_data_type) :: model_data

  ! Coordinate field
  type(field_type), target :: chi(3)

  integer(i_def) :: mesh_id      = imdi
  integer(i_def) :: twod_mesh_id = imdi

  integer(i_def) :: prognostic_init_choice, ancil_choice

  logical(l_def) :: put_field

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Sets up required state in preparation for run.
  !>@param[in] filename Name of the file containing the desired configuration 
  subroutine initialise( filename )

    implicit none

    character(*), intent(in) :: filename

    integer(i_def) :: comm = -999
    integer(i_def) :: ts_init

    ! Initialse mpi and create the default communicator: mpi_comm_world
    call initialise_comm(comm)

    ! Initialise aspects of the infrastructure
    call initialise_infrastructure(comm, filename, program_name, xios_id)

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    if ( subroutine_timers ) call timer('gungho')

    call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    if ( subroutine_counters ) then
      allocate(halo_calls, source=count_type('halo_calls'))
      call halo_calls%counter('gungho')
    end if

    ! Initialise aspects of the grid
    call initialise_grid(mesh_id, twod_mesh_id, chi)

    !Initialise aspects of output
    call initialise_io(comm, mesh_id, twod_mesh_id, chi, xios_ctx)

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_runtime_constants(mesh_id, twod_mesh_id, chi)

#ifdef UM_PHYSICS
    ! Set derived planet constants and presets
    call set_planet_constants()
#endif

    !-------------------------------------------------------------------------
    ! Select how to initialise model prognostic fields
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

    ! Create gungho prognostics and auxilliary (diagnostic) fields
    call create_gungho_prognostics( mesh_id,                        & 
                                    model_data%depository,          &
                                    model_data%prognostic_fields,   &
                                    model_data%diagnostic_fields,   &
                                    model_data%mr,                  &
                                    model_data%moist_dyn )

    ! Create prognostics used by physics
    if (use_physics) then
      call create_physics_prognostics( mesh_id, twod_mesh_id,        &
                                       model_data%depository,        &
                                       model_data%prognostic_fields, &
                                       model_data%derived_fields,    &
                                       model_data%cloud_fields,      &
                                       model_data%twod_fields,       &
                                       model_data%radstep_fields,    &
                                       model_data%physics_incs,      &
                                       model_data%jules_ancils,      &
                                       model_data%jules_prognostics )
    end if

    ! Create FD prognostic fields
    select case ( prognostic_init_choice )
      case ( init_option_fd_start_dump )
        if (use_physics) call create_fd_prognostics(mesh_id, model_data%fd_fields)
    end select

    ! Initialise prognostic fields appropriately
    select case ( prognostic_init_choice )

      case ( init_option_analytic )

        ! Initialise prognostics analytically according to
        ! namelist options

        call init_gungho_prognostics_alg( model_data%prognostic_fields, &
                                          model_data%diagnostic_fields, &
                                          model_data%mr,                & 
                                          model_data%moist_dyn )

        if (use_physics) then
          call init_physics_prognostics_alg( model_data%derived_fields, &
                                             model_data%cloud_fields,   &
                                             model_data%twod_fields,    &
                                             model_data%physics_incs,   &
                                             model_data%jules_ancils,   &
                                             model_data%jules_prognostics )
        end if

      case ( init_option_checkpoint_dump )
 
        ! Initialize prognostics using a checkpoint file
        ! from a previous run

        call read_checkpoint( model_data%prognostic_fields, timestep_start-1 )

        ! Update factors for moist dynamics
        call moist_dyn_factors_alg( model_data%moist_dyn, model_data%mr )

        if (use_physics) then
          ! if no cloud scheme, reset cloud variables
          if ( cloud == cloud_none ) then  
            call initial_cloud_alg( model_data%cloud_fields )
          end if
          ! re-initialise jules fields
          call init_jules_alg( model_data%jules_ancils, model_data%jules_prognostics )
          ! Set the increments to 0 initially
          call init_physics_incs_alg( model_data%physics_incs )
        end if

      case ( init_option_fd_start_dump )

        if (use_physics) then

          ! Initialise FD prognostic fields from a UM2LFRic dump     

          ! Read in from a UM2LFRic dump file
          call init_fd_prognostics_dump( model_data%fd_fields )

          ! Initialise jules fields
          call init_jules_alg( model_data%jules_ancils, model_data%jules_prognostics )

          ! Set physics increments to 0
          call init_physics_incs_alg( model_data%physics_incs )

          ! Populate prognostics from input finite difference fields
          call map_fd_to_prognostics( model_data%prognostic_fields,          &
                                      model_data%diagnostic_fields,          &
                                      model_data%mr,                         &
                                      model_data%moist_dyn,                  &
                                      model_data%cloud_fields,               &
                                      model_data%twod_fields,                &
                                      model_data%fd_fields )

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
          call init_aquaplanet_ancils( model_data%twod_fields )
        case ( ancil_option_analytic )
          call log_event( "Gungho: Setting ancillaries from analytic representation ", LOG_LEVEL_INFO )
          call init_analytic_ancils( model_data%twod_fields )
        case default
          ! No valid ancil option selected
          call log_event("Gungho: No valid ancillary initialisation option selected, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      end select

      ! All reading has been done, map the SST into the correct
      ! location of the multi-dimensional field
      put_field = .true.
      call update_tstar_alg( model_data%twod_fields, model_data%jules_prognostics, put_field )

    end if

    ! Initial output
    ! We only want these once at the beginning of a run
    ts_init = max( (timestep_start - 1), 0 ) ! 0 or t previous.

    if (ts_init == 0) then

      if ( write_diag ) then

        ! Calculation and output of diagnostics
        call gungho_diagnostics_driver( mesh_id,                      &
                                        model_data%prognostic_fields, &
                                        model_data%diagnostic_fields, &
                                        model_data%mr,                &
                                        model_data%moist_dyn,         &
                                        model_data%cloud_fields,      &
                                        model_data%derived_fields,    &
                                        ts_init,                      &
                                        nodal_output_on_w3 )

      end if

    end if

    call initialise_model( mesh_id,                      &
                           model_data%prognostic_fields, &
                           model_data%diagnostic_fields, &
                           model_data%mr,                &
                           model_data%twod_fields,       &
                           timestep_start )

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Timesteps the model, calling the desired timestepping algorithm based
  !upon the configuration
  subroutine run()

    implicit none

    integer(i_def) :: timestep

    call log_event( 'Running '//program_name//' ...', LOG_LEVEL_ALWAYS )

    do timestep = timestep_start, timestep_end

      ! Update XIOS calendar if we are using it for diagnostic output or checkpoint
      if ( use_xios_io ) then
        call log_event( "Gungho: Updating XIOS timestep", LOG_LEVEL_INFO )
        call xios_update_calendar(timestep - timestep_start + 1)
      end if

      write( log_scratch_space, '("/", A, "\ ")' ) repeat( "*", 76 )
      call log_event( log_scratch_space, LOG_LEVEL_TRACE )
      write( log_scratch_space, '(A,I0)' ) 'Start of timestep ', timestep
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      ! Perform a timestep
      call step( mesh_id,                      &
                 twod_mesh_id,                 &
                 model_data%prognostic_fields, &
                 model_data%diagnostic_fields, &
                 model_data%mr,                &
                 model_data%moist_dyn,         &
                 model_data%derived_fields,    &
                 model_data%cloud_fields,      &
                 model_data%twod_fields,       &
                 model_data%radstep_fields,    &
                 model_data%physics_incs,      &
                 model_data%jules_ancils,      &
                 model_data%jules_prognostics, &
                 timestep )

      write( log_scratch_space, '(A,I0)' ) 'End of timestep ', timestep
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      write( log_scratch_space, '("\", A, "/ ")' ) repeat( "*", 76 )
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      ! Use diagnostic output frequency to determine whether to write
      ! diagnostics on this timestep

      if ( ( mod(timestep, diagnostic_frequency) == 0 ) .and. ( write_diag ) ) then

        call log_event("Gungho: writing diagnostic output", LOG_LEVEL_INFO)

        ! Calculation and output diagnostics
        call gungho_diagnostics_driver( mesh_id,                      &
                                        model_data%prognostic_fields, &
                                        model_data%diagnostic_fields, &
                                        model_data%mr,                &
                                        model_data%moist_dyn,         &
                                        model_data%cloud_fields,      &
                                        model_data%derived_fields,    &
                                        timestep,                     &
                                        nodal_output_on_w3 )
      end if

    end do ! end ts loop

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tidies up after a run.
  subroutine finalise()

    implicit none

    call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

    if (use_physics) then
      ! All running has been done, map the SST back into the dumped field
      put_field = .false.
      call update_tstar_alg( model_data%twod_fields, model_data%jules_prognostics, put_field )
    end if

    ! Write checkpoint files if required
    if( checkpoint_write ) then
       call write_checkpoint( model_data%prognostic_fields, timestep_end )
    end if

    call finalise_model( mesh_id,                      &
                         model_data%prognostic_fields, &
                         model_data%diagnostic_fields, &
                         model_data%mr,                &
                         model_data%twod_fields,       &
                         model_data%fd_fields,         &
                         program_name )

    call finalise_model_data( model_data )

    call final_runtime_constants()

    if ( subroutine_timers ) then
      call timer('gungho')
      call output_timer()
    end if

    if ( subroutine_counters ) then
      call halo_calls%counter('gungho')
      call halo_calls%output_counters()
    end if

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

    call finalise_io()

    call finalise_grid()

    call log_event( program_name//' completed.', LOG_LEVEL_ALWAYS )

    ! Finalise infrastructure
    call finalise_infrastructure()

    ! Finalise mpi and release the communicator
    call finalise_comm()

  end subroutine finalise

end module gungho_driver_mod
