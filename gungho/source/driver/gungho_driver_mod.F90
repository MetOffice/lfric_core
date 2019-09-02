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
  use formulation_config_mod,     only : use_physics
  use function_space_collection_mod, &
                                  only : function_space_collection
  use field_collection_mod,       only : field_collection_type
  use final_gungho_mod,           only : final_gungho
  use gungho_diagnostics_driver_mod, &
                                  only : gungho_diagnostics_driver
  use gungho_grid_mod,            only : initialise_grid
  use gungho_infrastructure_mod,  only : initialise_infrastructure, &
                                         finalise_infrastructure
  use gungho_io_mod,              only : initialise_io, &
                                         finalise_io
  use create_gungho_prognostics_mod, &
                                  only : create_gungho_prognostics
  use create_physics_prognostics_mod, &
                                  only : create_physics_prognostics
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
  use init_gungho_mod,            only : init_gungho
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
  use mr_indices_mod,             only : nummr
  use moist_dyn_mod,              only : num_moist_factors
#ifdef UM_PHYSICS
  use planet_constants_mod,       only : set_planet_constants
#endif
  use runtime_constants_mod,      only : final_runtime_constants
  use step_gungho_mod,            only : step_gungho
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

  ! Depository of shared fields
  type( field_collection_type ), target :: depository

  ! Surface fields
  type( field_collection_type ) :: jules_ancils
  type( field_collection_type ) :: jules_prognostics

  ! Field collections
  type( field_collection_type ) :: prognostic_fields
  type( field_collection_type ) :: diagnostic_fields
  type( field_collection_type ) :: derived_fields
  type( field_collection_type ) :: cloud_fields
  type( field_collection_type ) :: twod_fields
  type( field_collection_type ) :: radstep_fields
  type( field_collection_type ) :: physics_incs
  type( field_collection_type ) :: fd_fields


  ! Auxiliary prognostic fields
  ! Moisture mixing ratios
  type( field_type ) :: mr(nummr)

  ! Moist dynamics
  type( field_type ) :: moist_dyn(num_moist_factors)

  ! A pointer used for retrieving fields from collections
  ! when iterating over them
  type( field_type ), pointer :: field_ptr  => null()

  ! Coordinate field
  type(field_type), target :: chi(3)

  integer(i_def) :: mesh_id      = imdi
  integer(i_def) :: twod_mesh_id = imdi

  integer(i_def) :: i 
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

    ! Create the depository, prognostics and diagnostics field collections.
    depository = field_collection_type(name='depository')
    prognostic_fields = field_collection_type(name="prognostics")
    diagnostic_fields = field_collection_type(name="diagnostics")

    ! Create gungho prognostics and auxilliary (diagnostic) fields
    call create_gungho_prognostics( mesh_id, depository, &
                                    prognostic_fields, diagnostic_fields, &
                                    mr, moist_dyn )

    ! Create prognostics used by physics
    if (use_physics) then
      call create_physics_prognostics( mesh_id, twod_mesh_id, &
                                       depository, &
                                       prognostic_fields, &
                                       derived_fields, cloud_fields, &
                                       twod_fields, radstep_fields, &
                                       physics_incs, &
                                       jules_ancils, jules_prognostics )
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

        call init_gungho_prognostics_alg(prognostic_fields, diagnostic_fields, &
                                         mr, moist_dyn)

        if (use_physics) then
          call init_physics_prognostics_alg(derived_fields, &
                                            cloud_fields, &
                                            twod_fields, &
                                            physics_incs, &
                                            jules_ancils, &
                                            jules_prognostics)
        end if

      case ( init_option_checkpoint_dump )
 
        ! Initialize prognostics using a checkpoint file
        ! from a previous run

        call read_checkpoint(prognostic_fields, timestep_start-1)

        ! Update factors for moist dynamics
        call moist_dyn_factors_alg(moist_dyn, mr)

        if (use_physics) then
          ! if no cloud scheme, reset cloud variables
          if ( cloud == cloud_none ) then  
            call initial_cloud_alg(cloud_fields)
          end if
          ! re-initialise jules fields
          call init_jules_alg(jules_ancils, jules_prognostics)
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

          ! Initialise jules fields
          call init_jules_alg(jules_ancils, jules_prognostics)

          ! Set physics increments to 0
          call init_physics_incs_alg(physics_incs)

          ! Populate prognostics from input finite difference fields
          call map_fd_to_prognostics(prognostic_fields, diagnostic_fields, &
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

      ! All reading has been done, map the SST into the correct
      ! location of the multi-dimensional field
      put_field = .true.
      call update_tstar_alg(twod_fields, jules_prognostics, put_field)

    end if

    ! Initial output
    ! We only want these once at the beginning of a run
    ts_init = max( (timestep_start - 1), 0 ) ! 0 or t previous.

    if (ts_init == 0) then

      if ( write_diag ) then

        ! Calculation and output of diagnostics
        call gungho_diagnostics_driver( mesh_id, &
                                        prognostic_fields, &
                                        diagnostic_fields, &
                                        mr, &
                                        moist_dyn, &
                                        cloud_fields, &
                                        derived_fields, &
                                        ts_init, &
                                        nodal_output_on_w3 )

      end if

    end if

#ifdef UM_PHYSICS
    ! Set derived planet constants and presets
    call set_planet_constants()
#endif

    call init_gungho( mesh_id, &
                      prognostic_fields,diagnostic_fields, &
                      mr, twod_fields, &
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
        call xios_update_calendar(timestep)
      end if

      write( log_scratch_space, '("/", A, "\ ")' ) repeat( "*", 76 )
      call log_event( log_scratch_space, LOG_LEVEL_TRACE )
      write( log_scratch_space, '(A,I0)' ) 'Start of timestep ', timestep
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      ! Perform a timestep
      call step_gungho( mesh_id,           &
                        twod_mesh_id,      &
                        prognostic_fields, &
                        diagnostic_fields, &
                        mr,                &
                        moist_dyn,         &
                        derived_fields,    &
                        cloud_fields,      &
                        twod_fields,       &
                        radstep_fields,    &
                        physics_incs,      &
                        jules_ancils,      &
                        jules_prognostics, &
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
        call gungho_diagnostics_driver( mesh_id, &
                                        prognostic_fields, &
                                        diagnostic_fields, &
                                        mr, &
                                        moist_dyn, &
                                        cloud_fields, &
                                        derived_fields, &
                                        timestep, &
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
      call update_tstar_alg(twod_fields, jules_prognostics, put_field)
    end if

    ! Write checkpoint files if required
    if( checkpoint_write ) then
       call write_checkpoint(prognostic_fields, timestep_end)
    end if

    call final_gungho( mesh_id,           &
                       prognostic_fields, &
                       diagnostic_fields, &
                       mr,                &
                       twod_fields,       &
                       fd_fields,         &
                       program_name )

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

    if (allocated(mesh_collection)) then
      call mesh_collection%clear()
      deallocate(mesh_collection)
    end if

    if (allocated(function_space_collection)) then
      call function_space_collection%clear()
      deallocate(function_space_collection)
    end if

    call log_event( program_name//' completed.', LOG_LEVEL_ALWAYS )

    ! Finalise infrastructure
    call finalise_infrastructure()

    ! Finalise mpi and release the communicator
    call finalise_comm()

  end subroutine finalise

end module gungho_driver_mod
