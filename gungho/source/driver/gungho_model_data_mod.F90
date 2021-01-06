!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Container for Gung Ho model run working data set inclduing methods
!> to initialise, copy and finalise the data set
!>
!> This module provides a type to hold all the model fields and methods to
!> initialise (create and read), copy and finalise (write and destroy) the
!> data contained within the type.
!>
module gungho_model_data_mod

  use clock_mod,                        only : clock_type
  use mr_indices_mod,                   only : nummr
  use moist_dyn_mod,                    only : num_moist_factors
  use field_mod,                        only : field_type
  use field_parent_mod,                 only : write_interface
  use field_collection_mod,             only : field_collection_type
  use constants_mod,                    only : i_def, l_def
  use log_mod,                          only : log_event,      &
                                               LOG_LEVEL_INFO, &
                                               LOG_LEVEL_ERROR
  use formulation_config_mod,           only : use_physics
  use initialization_config_mod,        only : init_option,                 &
                                               init_option_analytic,        &
                                               init_option_fd_start_dump,   &
                                               init_option_checkpoint_dump, &
                                               init_option_fe_start_dump,   &
                                               ancil_option,                &
                                               ancil_option_none,           &
                                               ancil_option_aquaplanet,     &
                                               ancil_option_basic_gal,      &
                                               ancil_option_prototype_gal
  use io_config_mod,                    only : checkpoint_read,  &
                                               checkpoint_write, &
                                               write_dump
  use read_methods_mod,                 only : read_checkpoint,  &
                                               read_state
  use write_methods_mod,                only : write_checkpoint, &
                                               write_state
  use boundaries_config_mod,            only : limited_area
  use create_lbcs_mod,                  only : create_lbc_fields
  use create_gungho_prognostics_mod,    only : create_gungho_prognostics
  use create_physics_prognostics_mod,   only : create_physics_prognostics
  use section_choice_config_mod,        only : cloud, cloud_none
  use map_fd_to_prognostics_alg_mod,    only : map_fd_to_prognostics
  use init_gungho_prognostics_alg_mod,  only : init_gungho_prognostics_alg
  use init_physics_prognostics_alg_mod, only : init_physics_prognostics_alg
  use update_tstar_alg_mod,             only : update_tstar_alg
  use moist_dyn_factors_alg_mod,        only : moist_dyn_factors_alg
  use init_fd_prognostics_mod,          only : init_fd_prognostics_dump
#ifdef UM_PHYSICS
  use create_fd_prognostics_mod,        only : create_fd_prognostics
  use init_ancils_mod,                  only : create_fd_ancils
  use process_inputs_alg_mod,           only : process_inputs_alg
#endif
  use linked_list_mod,                  only : linked_list_type
  use variable_fields_mod,              only : init_variable_fields

  implicit none

  private

  !> Holds the working data set for a model run and other working state.
  !>
  type :: model_data_type

    private

    !> Stores all the fields used by the model - the remaining collections store
    !> pointers to the data in the depository.
    type( field_collection_type ), public :: depository

    !> @name Fields needed to time-step the model.
    !> @{

    !> All the prognostic fields (except for field arrays: auxiliary prognostic)
    type( field_collection_type ), public   :: prognostic_fields
    !> All the diagnostic fields
    type( field_collection_type ), public   :: diagnostic_fields
    !> FD fields derived from FE fields for use in physics time-stepping schemes
    type( field_collection_type ), public   :: derived_fields
    !> LBC fields - lateral boundary conditions to run a limited area model
    type( field_collection_type ), public   :: lbc_fields
    !> Fields owned by the radiation scheme
    type( field_collection_type ), public   :: radiation_fields
    !> Fields owned by the microphysics scheme
    type( field_collection_type ), public   :: microphysics_fields
    !> Fields owned by the orographic drag schemes
    type( field_collection_type ), public   :: orography_fields
    !> Fields owned by the turbulence scheme
    type( field_collection_type ), public   :: turbulence_fields
    !> Fields owned by the convection schemes
    type( field_collection_type ), public   :: convection_fields
    !> Fields owned by the cloud schemes
    type( field_collection_type ), public   :: cloud_fields
    !> Fields owned by the surface exchange scheme
    type( field_collection_type ), public   :: surface_fields
    !> Fields owned by the soil hydrology scheme
    type( field_collection_type ), public   :: soil_fields
    !> Fields owned by the snow scheme
    type( field_collection_type ), public   :: snow_fields
    !> Fields owned by the aerosol schemes
    type( field_collection_type ), public   :: aerosol_fields
    !> Array of fields containing the moisture mixing ratios
    !>  (auxiliary prognostic)
    type( field_type ), allocatable, public :: mr(:)
    !> Array of fields containing the moist dynamics (auxiliary prognostic)
    type( field_type ), allocatable, public :: moist_dyn(:)
    !> @}

    !> FD fields used to read initial conditions from LFRic-Input files
    type( field_collection_type ), public   :: fd_fields

    !> Fields used to store data read in from ancillary files
    type( field_collection_type ), public   :: ancil_fields
    !> Linked list of time axis objects used to update time-varying ancils
    type( linked_list_type ),      public   :: ancil_times_list


  end type model_data_type

  ! Set these to select how to initialize model prognostic fields
  integer(i_def) :: prognostic_init_choice, ancil_choice

  logical(l_def) :: put_field

  public model_data_type, create_model_data, finalise_model_data, &
         initialise_model_data, output_model_data

contains

  !> @brief Create the fields contained in model_data
  !> @param[inout] model_data The working data set for a model run
  !> @param[in]    mesh_id The identifier given to the current 3d mesh
  !> @param[in]    twod_mesh_id The identifier given to the current 2d mesh
  subroutine create_model_data( model_data, &
                                mesh_id,    &
                                twod_mesh_id, &
                                clock )

    implicit none

    type( model_data_type ), intent(inout) :: model_data
    integer(i_def),          intent(in)    :: mesh_id
    integer(i_def),          intent(in)    :: twod_mesh_id
    class(clock_type),       intent(in)    :: clock

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

    !-------------------------------------------------------------------------
    ! Instantiate the fields
    !-------------------------------------------------------------------------

    ! Field bundles - allocate the fields so thay can be cleared
    allocate(model_data%moist_dyn(num_moist_factors))
    allocate(model_data%mr(nummr))

    ! Create gungho prognostics and auxilliary (diagnostic) fields
    call create_gungho_prognostics( mesh_id,                        &
                                    model_data%depository,          &
                                    model_data%prognostic_fields,   &
                                    model_data%diagnostic_fields,   &
                                    model_data%mr,                  &
                                    model_data%moist_dyn )

    if (limited_area) call create_lbc_fields( mesh_id,                      &
                                              model_data%depository,        &
                                              model_data%prognostic_fields, &
                                              model_data%lbc_fields )

    ! Create prognostics used by physics
    if (use_physics) then
      call create_physics_prognostics( mesh_id, twod_mesh_id,          &
                                       clock,                          &
                                       model_data%depository,          &
                                       model_data%prognostic_fields,   &
                                       model_data%derived_fields,      &
                                       model_data%radiation_fields,    &
                                       model_data%microphysics_fields, &
                                       model_data%orography_fields,    &
                                       model_data%turbulence_fields,   &
                                       model_data%convection_fields,   &
                                       model_data%cloud_fields,        &
                                       model_data%surface_fields,      &
                                       model_data%soil_fields,         &
                                       model_data%snow_fields,         &
                                       model_data%aerosol_fields )

#ifdef UM_PHYSICS
      ! Create FD prognostic fields
      select case ( prognostic_init_choice )
        case ( init_option_fd_start_dump )
          call create_fd_prognostics(mesh_id, twod_mesh_id, &
                                     model_data%fd_fields,  &
                                     model_data%depository)
      end select

      ! Create and populate collection of fields to be read from ancillary files
      select case ( ancil_choice )
        case ( ancil_option_basic_gal, ancil_option_prototype_gal )
          call create_fd_ancils( model_data%depository,   &
                                 model_data%ancil_fields, &
                                 mesh_id, twod_mesh_id ,  &
                                 model_data%ancil_times_list )
      end select
#endif
    end if

  end subroutine create_model_data

  !-------------------------------------------------------------------------------
  !> @brief Initialises the working data set dependent of namelist configuration
  !> @param [inout] model_data The working data set for a model run
  subroutine initialise_model_data( model_data, clock )

    implicit none

    type( model_data_type ), intent(inout) :: model_data
    class(clock_type),       intent(in)    :: clock

    ! Initialise all the physics fields here. We'll then re initialise
    ! them below if need be
    if (use_physics) then
          call init_physics_prognostics_alg( model_data%radiation_fields,    &
                                             model_data%microphysics_fields, &
                                             model_data%orography_fields,    &
                                             model_data%turbulence_fields,   &
                                             model_data%convection_fields,   &
                                             model_data%cloud_fields,        &
                                             model_data%surface_fields,      &
                                             model_data%soil_fields,         &
                                             model_data%snow_fields,         &
                                             model_data%aerosol_fields )

    end if

    ! Initialise prognostic fields appropriately
    select case ( prognostic_init_choice )

      case ( init_option_analytic )

        ! Initialise prognostics analytically according to
        ! namelist options

        call init_gungho_prognostics_alg( model_data%prognostic_fields, &
                                          model_data%diagnostic_fields, &
                                          model_data%mr,                &
                                          model_data%moist_dyn )

      case ( init_option_checkpoint_dump )

        ! Initialize prognostics using a checkpoint file
        ! from a previous run

        call read_checkpoint( model_data%prognostic_fields, &
                              clock%get_first_step() - 1 )

        ! Update factors for moist dynamics
        call moist_dyn_factors_alg( model_data%moist_dyn, model_data%mr )

      case ( init_option_fd_start_dump )

        if (use_physics) then

          ! Initialise FD prognostic fields from a UM2LFRic dump

          ! Read in from a UM2LFRic dump file
          call init_fd_prognostics_dump( model_data%fd_fields )

          ! Populate prognostics from input finite difference fields
          call map_fd_to_prognostics( model_data%prognostic_fields,          &
                                      model_data%mr,                         &
                                      model_data%moist_dyn,                  &
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
          ! Update the tiled surface temperature with the calculated tstar
          put_field = .true.
          call update_tstar_alg(model_data%surface_fields, &
                                model_data%fd_fields, put_field )
        case ( ancil_option_basic_gal, ancil_option_prototype_gal )
          call log_event( "Gungho: Reading basic/proto GAL ancils ", LOG_LEVEL_INFO )
          call read_state( model_data%ancil_fields )
          call init_variable_fields( model_data%ancil_times_list, &
                                     clock, model_data%ancil_fields )
#ifdef UM_PHYSICS
          call process_inputs_alg( model_data%ancil_fields,   &
                                   model_data%fd_fields,      &
                                   model_data%surface_fields, &
                                   model_data%soil_fields,    &
                                   model_data%snow_fields)
#endif
        case default
          ! No valid ancil option selected
          call log_event("Gungho: No valid ancillary initialisation option selected, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      end select

    end if

  end subroutine initialise_model_data

  !----------------------------------------------------------------------------
  !> @brief Writes out a checkpoint and dump file dependent on namelist
  !> options
  !> @param[inout] model_data The working data set for the model run
  !> @param[in] clock Model time.
  !>
  subroutine output_model_data( model_data, &
                                clock )

    implicit none

    type( model_data_type ), intent(inout), target :: model_data
    class(clock_type),       intent(in)            :: clock

    type( field_collection_type ), pointer :: fd_fields => null()
    type( field_collection_type ), pointer :: prognostic_fields => null()

    ! Get pointers to field collections for use downstream
    fd_fields => model_data%fd_fields
    prognostic_fields => model_data%prognostic_fields

    !===================== Write fields to dump ======================!
    if( use_physics ) then

      ! Current dump writing is only relevant for physics runs at the moment
      if (write_dump) then

        call log_event("Gungho: writing FD fields to dump is not yet supported", LOG_LEVEL_ERROR)

        ! Write prognostic fields to dump
        call write_state(fd_fields)

      end if

    end if

    !=================== Write fields to checkpoint files ====================!
    if( checkpoint_write ) then
       call write_checkpoint( prognostic_fields, clock )
    end if

  end subroutine output_model_data

  !>@brief Routine to destroy all the field collections in the working data set
  !> @param[inout] model_data The working data set for a model run
  subroutine finalise_model_data( model_data )

    implicit none

      type(model_data_type), intent(inout) :: model_data

      ! Clear all the fields in each field collection
      call model_data%depository%clear()
      call model_data%prognostic_fields%clear()
      call model_data%diagnostic_fields%clear()
      call model_data%derived_fields%clear()
      call model_data%radiation_fields%clear()
      call model_data%microphysics_fields%clear()
      call model_data%orography_fields%clear()
      call model_data%turbulence_fields%clear()
      call model_data%convection_fields%clear()
      call model_data%cloud_fields%clear()
      call model_data%surface_fields%clear()
      call model_data%soil_fields%clear()
      call model_data%snow_fields%clear()
      call model_data%aerosol_fields%clear()
      call model_data%fd_fields%clear()
      call model_data%lbc_fields%clear()
      if (allocated(model_data%mr)) deallocate(model_data%mr)
      if (allocated(model_data%moist_dyn)) deallocate(model_data%moist_dyn)

      call log_event( 'finalise_model_data: all fields have been cleared', &
                       LOG_LEVEL_INFO )

  end subroutine finalise_model_data

end module gungho_model_data_mod
