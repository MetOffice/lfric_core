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

  use mr_indices_mod,                   only : nummr
  use moist_dyn_mod,                    only : num_moist_factors
  use field_mod,                        only : field_type, &
                                               write_interface
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
                                               ancil_option_analytic,       &
                                               ancil_option_aquaplanet,     &
                                               ancil_option_basic_gagl
  use io_config_mod,                    only : checkpoint_read,  &
                                               checkpoint_write, &
                                               write_dump
  use io_mod,                           only : read_checkpoint,  &
                                               write_checkpoint, &
                                               read_state,       &
                                               write_state,      &
                                               xios_write_field_single_face
  use create_gungho_prognostics_mod,    only : create_gungho_prognostics
  use create_physics_prognostics_mod,   only : create_physics_prognostics
  use create_fd_prognostics_mod,        only : create_fd_prognostics
  use section_choice_config_mod,        only : cloud, cloud_none
  use time_config_mod,                  only : timestep_start
  use map_fd_to_prognostics_alg_mod,    only : map_fd_to_prognostics
  use init_gungho_prognostics_alg_mod,  only : init_gungho_prognostics_alg
  use init_physics_prognostics_alg_mod, only : init_physics_prognostics_alg
  use init_aerosol_alg_mod,             only : init_aerosol_alg
  use update_tstar_alg_mod,             only : update_tstar_alg
  use moist_dyn_factors_alg_mod,        only : moist_dyn_factors_alg
  use initial_cloud_alg_mod,            only : initial_cloud_alg
  use init_jules_alg_mod,               only : init_jules_alg
  use init_orography_fields_alg_mod,    only : init_orography_fields_alg
  use init_physics_incs_alg_mod,        only : init_physics_incs_alg
  use init_fd_prognostics_mod,          only : init_fd_prognostics_dump
  use init_ancils_mod,                  only : init_analytic_ancils,    &
                                               init_aquaplanet_ancils,  &
                                               create_fd_ancils

  implicit none

  private

  !> Holds the working data set for a model run.
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

  end type model_data_type

  ! Set these to select how to initialize model prognostic fields
  integer(i_def) :: prognostic_init_choice, ancil_choice

  public model_data_type, create_model_data, finalise_model_data, &
         initialise_model_data, output_model_data

contains

  !> @brief Create the fields contained in model_data
  !> @param[inout] model_data The working data set for a model run
  !> @param[in]    mesh_id The identifier given to the current 3d mesh
  !> @param[in]    twod_mesh_id The identifier given to the current 2d mesh
  subroutine create_model_data( model_data, &
                                mesh_id,    &
                                twod_mesh_id )

    implicit none

    type( model_data_type ), intent(inout) :: model_data
    integer(i_def), intent(in)             :: mesh_id
    integer(i_def), intent(in)             :: twod_mesh_id

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

    ! Create prognostics used by physics
    if (use_physics) then
      call create_physics_prognostics( mesh_id, twod_mesh_id,          &
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
    end if

    ! Create FD prognostic fields
    select case ( prognostic_init_choice )
      case ( init_option_fd_start_dump )
        if (use_physics) call create_fd_prognostics(mesh_id, model_data%fd_fields)
    end select

    select case ( ancil_choice )
      case ( ancil_option_basic_gagl )
        if (use_physics) call create_fd_ancils( model_data%depository,   &
                                                model_data%ancil_fields, &
                                                mesh_id, twod_mesh_id )
    end select

  end subroutine create_model_data

  !-------------------------------------------------------------------------------
  !> @brief Initialises the working data set dependent of namelist configuration
  !> @param [inout] model_data The working data set for a model run
  subroutine initialise_model_data( model_data )

    implicit none

    type( model_data_type ), intent(inout) :: model_data

    logical(l_def) :: put_field

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

        call read_checkpoint( model_data%prognostic_fields, timestep_start-1 )

        ! Update factors for moist dynamics
        call moist_dyn_factors_alg( model_data%moist_dyn, model_data%mr )

        if (use_physics) then
          ! if no cloud scheme, reset cloud variables
          if ( cloud == cloud_none ) then
            call initial_cloud_alg( model_data%convection_fields, &
                                    model_data%cloud_fields )
          end if

        end if

      case ( init_option_fd_start_dump )

        if (use_physics) then

          ! Initialise FD prognostic fields from a UM2LFRic dump

          ! Read in from a UM2LFRic dump file
          call init_fd_prognostics_dump( model_data%fd_fields )

          ! Populate prognostics from input finite difference fields
          call map_fd_to_prognostics( model_data%prognostic_fields,          &
                                      model_data%diagnostic_fields,          &
                                      model_data%mr,                         &
                                      model_data%moist_dyn,                  &
                                      model_data%turbulence_fields,          &
                                      model_data%convection_fields,          &
                                      model_data%cloud_fields,               &
                                      model_data%surface_fields,             &
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
          call init_aquaplanet_ancils( model_data%surface_fields )
        case ( ancil_option_analytic )
          call log_event( "Gungho: Setting ancillaries from analytic representation ", LOG_LEVEL_INFO )
          call init_analytic_ancils( model_data%surface_fields )
        case ( ancil_option_basic_gagl )
          call log_event( "Gungho: Reading ancillaries from disk ", LOG_LEVEL_INFO )
          call read_state( model_data%ancil_fields )
        case default
          ! No valid ancil option selected
          call log_event("Gungho: No valid ancillary initialisation option selected, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
      end select

      ! All reading has been done, map the SST into the correct
      ! location of the multi-dimensional field
      put_field = .true.
      call update_tstar_alg( model_data%surface_fields, put_field )

    end if

  end subroutine initialise_model_data

  !-------------------------------------------------------------------------------
  !> @brief Writes out a checkpoint and dump file dependent on namelist
  !> options
  !> @param[inout] model_data The working data set for the model run
  !> @param[in] timestep number of current timestep
  subroutine output_model_data( model_data, &
                                timestep )

    implicit none

    type( model_data_type ), target, intent(inout) :: model_data
    integer(i_def),          intent(in)    :: timestep

    type( field_collection_type ), pointer :: surface_fields => null()
    type( field_collection_type ), pointer :: fd_fields => null()
    type( field_collection_type ), pointer :: prognostic_fields => null()

    type( field_type ), pointer            :: tstar_2d => null()
    procedure(write_interface), pointer    :: tmp_write_ptr => null()
    logical(l_def)                         :: put_field

    ! Get pointers to field collections for use downstream
    surface_fields => model_data%surface_fields
    fd_fields => model_data%fd_fields
    prognostic_fields => model_data%prognostic_fields

    !===================== Write fields to dump ======================!
    if( use_physics ) then

      ! All running has been done, map the SST back into the dumped field
      put_field = .false.
      call update_tstar_alg(surface_fields, put_field )

      ! Current dump writing is only relevant for physics runs at the moment
      if (write_dump) then

        ! For the purposes of dumping from one collection, we add a pointer
        ! to tstar to the fd_prognostics collection

        tmp_write_ptr => xios_write_field_single_face
        tstar_2d => surface_fields%get_field('tstar')
        call tstar_2d%set_write_behaviour(tmp_write_ptr)

        call fd_fields%add_field(tstar_2d)

        call log_event("Gungho: writing FD fields to dump", LOG_LEVEL_INFO)

        ! Write prognostic fields to dump
        call write_state(fd_fields)

        nullify(tmp_write_ptr, tstar_2d)

      end if

    end if

    !===================== Write fields to checkpoint files ======================!
    if( checkpoint_write ) then
       call write_checkpoint( prognostic_fields, timestep )
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
      if (allocated(model_data%mr)) deallocate(model_data%mr)
      if (allocated(model_data%moist_dyn)) deallocate(model_data%moist_dyn)

      call log_event( 'finalise_model_data: all fields have been cleared',     &
                       LOG_LEVEL_INFO )

  end subroutine finalise_model_data

end module gungho_model_data_mod
