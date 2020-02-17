!-------------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief create physics prognostics
!> @details Creates the physics prognostic fields
module create_physics_prognostics_mod

  use constants_mod,                  only : i_def, l_def
  use field_mod,                      only : field_type,                       &
                                             write_interface,                  &
                                             checkpoint_write_interface,       &
                                             checkpoint_read_interface
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use field_collection_mod,           only : field_collection_type
  use fs_continuity_mod,              only : W2, W3, Wtheta
  use function_space_mod,             only : function_space_type
  use io_mod,                         only : tile_order, pft_order, &
                                             sice_order, soil_order, snow_order
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO,         &
                                             LOG_LEVEL_ERROR
  use radiation_config_mod,           only : n_radstep
  use aerosol_config_mod,             only : c_aerosol,                        &
                                             c_aerosol_glomap_mode_climatology
  use section_choice_config_mod,      only : cloud, cloud_um,                  &
                                             aerosol, aerosol_um,              &
                                             radiation, radiation_socrates
  use cloud_config_mod,               only : scheme, &
                                             scheme_pc2
  use time_config_mod,                only : timestep_start, timestep_end
  implicit none

  private
  public :: create_physics_prognostics

contains
  !>@brief Routine to initialise the field objects required by the physics
  !> @param[in]    mesh_id The identifier given to the current 3d mesh
  !> @param[in]    twod_mesh_id The identifier given to the current 2d mesh
  !> @param[inout] depository Main collection of all fields in memory
  !> @param[inout] prognostic_fields A collection of the fields that make up the
  !>                                 prognostic variables in the model
  !> @param[out]   derived_fields Collection of FD fields derived from FE fields
  !> @param[out]   radition_fields Collection of fields for radiation scheme
  !> @param[out]   microphysics_fields Collection of fields for microphys scheme
  !> @param[out]   orography_fields Collection of fields for orogr drag scheme
  !> @param[out]   turbulence_fields Collection of fields for turbulence scheme
  !> @param[out]   convection_fields Collection of fields for convection scheme
  !> @param[out]   cloud_fields Collection of fields for cloud scheme
  !> @param[out]   surface_fields Collection of fields for surface scheme
  !> @param[out]   soil_fields Collection of fields for soil hydrology scheme
  !> @param[out]   snow_fields Collection of fields for snow scheme
  !> @param[out]   aerosol_fields Collection of fields for aerosol scheme
  subroutine create_physics_prognostics( mesh_id,                              &
                                         twod_mesh_id,                         &
                                         depository,                           &
                                         prognostic_fields,                    &
                                         derived_fields,                       &
                                         radiation_fields,                     &
                                         microphysics_fields,                  &
                                         orography_fields,                     &
                                         turbulence_fields,                    &
                                         convection_fields,                    &
                                         cloud_fields,                         &
                                         surface_fields,                       &
                                         soil_fields,                          &
                                         snow_fields,                          &
                                         aerosol_fields )

    implicit none

    integer(i_def), intent(in) :: mesh_id
    integer(i_def), intent(in) :: twod_mesh_id

    ! Collections of fields
    type(field_collection_type), intent(inout) :: depository
    type(field_collection_type), intent(inout) :: prognostic_fields
    type(field_collection_type), intent(out) :: derived_fields
    type(field_collection_type), intent(out) :: radiation_fields
    type(field_collection_type), intent(out) :: microphysics_fields
    type(field_collection_type), intent(out) :: orography_fields
    type(field_collection_type), intent(out) :: turbulence_fields
    type(field_collection_type), intent(out) :: convection_fields
    type(field_collection_type), intent(out) :: cloud_fields
    type(field_collection_type), intent(out) :: surface_fields
    type(field_collection_type), intent(out) :: soil_fields
    type(field_collection_type), intent(out) :: snow_fields
    type(field_collection_type), intent(out) :: aerosol_fields

    ! pointers to vector spaces
    type(function_space_type), pointer :: vector_space => null()
    type(function_space_type), pointer :: wtheta_space => null()
    type(function_space_type), pointer :: w3_space => null()
    type(function_space_type), pointer :: w2_space => null()
    type(function_space_type), pointer :: twod_space => null()
    type(function_space_type), pointer :: tile_space => null()
    type(function_space_type), pointer :: pft_space => null()
    type(function_space_type), pointer :: soil_space => null()
    type(function_space_type), pointer :: sice_space => null()
    type(function_space_type), pointer :: snow_space => null()

    type( field_type ), pointer :: theta => null()

    integer(i_def) :: theta_space
    logical(l_def) :: checkpoint_restart_flag
    logical(l_def) :: advection_flag

    call log_event( 'Create physics prognostics', LOG_LEVEL_INFO )

    theta => prognostic_fields%get_field('theta')
    theta_space = theta%which_function_space()

    if (theta_space /= Wtheta)then
      call log_event( 'Physics: requires theta variable to be in Wtheta',      &
                      LOG_LEVEL_ERROR )
    end if

    if (element_order > 0)then
      call log_event( 'Physics: requires lowest order elements',               &
                      LOG_LEVEL_ERROR )
    end if

    ! Create the vector spaces once here for re-use later
    wtheta_space => function_space_collection%get_fs(mesh_id, 0, Wtheta)
    w3_space => function_space_collection%get_fs(mesh_id, 0, W3)
    w2_space => function_space_collection%get_fs(mesh_id, 0, W2)
    twod_space => function_space_collection%get_fs(twod_mesh_id, 0, W3)
    tile_space => function_space_collection%get_fs(twod_mesh_id, tile_order, W3)
    pft_space => function_space_collection%get_fs(twod_mesh_id, pft_order, W3)
    soil_space => function_space_collection%get_fs(twod_mesh_id, soil_order, W3)
    sice_space => function_space_collection%get_fs(twod_mesh_id, sice_order, W3)
    snow_space => function_space_collection%get_fs(twod_mesh_id, snow_order, W3)
    !========================================================================
    ! Fields derived from the FE dynamical fields for use in physics
    !========================================================================
    derived_fields  =  field_collection_type(name='derived_fields')
    checkpoint_restart_flag = .false.

    ! Wtheta fields
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'w_physics',      wtheta_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'rho_in_wth',     wtheta_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'wetrho_in_wth',  wtheta_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'exner_in_wth',   wtheta_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'w_physics_star', wtheta_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'theta_star',     wtheta_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'shear',          wtheta_space, checkpoint_restart_flag )

    ! W3 fields
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u1_in_w3',      w3_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u2_in_w3',      w3_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u3_in_w3',      w3_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'theta_in_w3',   w3_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'wetrho_in_w3',  w3_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u1_in_w3_star', w3_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u2_in_w3_star', w3_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u3_in_w3_star', w3_space, checkpoint_restart_flag )

    ! W2 fields
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u_physics',      w2_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u_physics_star', w2_space, checkpoint_restart_flag )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u_star',         w2_space, checkpoint_restart_flag )

    !========================================================================
    ! Fields owned by the radiation scheme
    !========================================================================
    radiation_fields = field_collection_type(name='radiation_fields')

    ! 2D fields, need checkpointing
    checkpoint_restart_flag = .true.
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'albedo_obs_sw', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'albedo_obs_vis', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'albedo_obs_nir', twod_space, checkpoint_restart_flag, twod=.true. )

    ! 2D fields, don't need checkpointing
    checkpoint_restart_flag = .false.
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'cos_zenith_angle',   twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lit_fraction',       twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lw_down_surf', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_down_surf', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_direct_surf', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_down_blue_surf', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_direct_blue_surf', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lw_up_tile', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_up_tile', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_up_blue_tile', tile_space, checkpoint_restart_flag, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_heating_rate', wtheta_space, checkpoint_restart_flag )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lw_heating_rate', wtheta_space, checkpoint_restart_flag )

    ! Fields which need checkpointing for radiation timestepping
    if (radiation == radiation_socrates) then
      ! Checkpoint unless both the first timestep of this run and the
      ! first timestep of the next run are radiation timesteps
      checkpoint_restart_flag = &
        mod(timestep_start-1, n_radstep) /= 0 .or. &
        mod(timestep_end,     n_radstep) /= 0
    else
      checkpoint_restart_flag = .false.
    end if

    ! 2D fields
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lw_down_surf_rts', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_down_surf_rts', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_direct_surf_rts', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_down_blue_surf_rts', twod_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_direct_blue_surf_rts', twod_space, checkpoint_restart_flag,          &
      twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'cos_zenith_angle_rts', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lit_fraction_rts', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'stellar_irradiance_rts', twod_space, checkpoint_restart_flag,           &
      twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sin_stellar_declination_rts', twod_space, checkpoint_restart_flag,      &
      twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'stellar_eqn_of_time_rts', twod_space, checkpoint_restart_flag,          &
      twod=.true. )

    ! 3D fields
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_heating_rate_rts', wtheta_space, checkpoint_restart_flag )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lw_heating_rate_rts', wtheta_space, checkpoint_restart_flag )

    ! Fields on surface tiles
    checkpoint_restart_flag = .false. ! #1920 needed for this to work
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lw_up_tile_rts', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_up_tile_rts', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_up_blue_tile_rts', tile_space, checkpoint_restart_flag, twod=.true. )

    !========================================================================
    ! Fields owned by the microphysics scheme
    !========================================================================
    microphysics_fields = field_collection_type(name='microphysics_fields')

    ! 2D fields, don't need checkpointing
    checkpoint_restart_flag = .false.
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      'ls_rain',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      'ls_snow',  twod_space, checkpoint_restart_flag, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      'dtl_mphys', wtheta_space, checkpoint_restart_flag )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      'dmt_mphys', wtheta_space, checkpoint_restart_flag )

    !========================================================================
    ! Fields owned by the orographic drag schemes
    !========================================================================
    orography_fields = field_collection_type(name='orography_fields')

    ! 2D fields, need checkpointing
    checkpoint_restart_flag = .true.
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      'sd_orog', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      'grad_xx_orog', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      'grad_xy_orog', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      'grad_yy_orog', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      'peak_to_trough_orog', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      'silhouette_area_orog', twod_space, checkpoint_restart_flag, twod=.true. )

    !========================================================================
    ! Fields owned by the turbulence scheme
    !========================================================================
    turbulence_fields = field_collection_type(name='turbulence_fields')

    ! 2D fields, need checkpointing
    checkpoint_restart_flag = .true.
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'zh',      twod_space, checkpoint_restart_flag, twod=.true. )

    ! 2D fields, don't need checkpointing
    checkpoint_restart_flag = .false.
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'ntml',    twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'cumulus', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'zlcl_mixed',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'blend_height_tq',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'blend_height_uv',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'zh_nonloc',  twod_space, checkpoint_restart_flag, twod=.true. )

    ! a couple of BL arrays which are helpful to group together as higher
    ! order fields
    vector_space => function_space_collection%get_fs(twod_mesh_id, 1, W3)
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'rhokm_surf',  vector_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'bl_types',  vector_space, checkpoint_restart_flag, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'dt_bl', wtheta_space, checkpoint_restart_flag )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'dmv_bl', wtheta_space, checkpoint_restart_flag )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'rhokm_bl', wtheta_space, checkpoint_restart_flag )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'ngstress_bl', wtheta_space, checkpoint_restart_flag )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'bq_bl', wtheta_space, checkpoint_restart_flag )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'bt_bl', wtheta_space, checkpoint_restart_flag )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'dtrdz_tq_bl', wtheta_space, checkpoint_restart_flag )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'rdz_uv_bl', wtheta_space, checkpoint_restart_flag )

    ! 3D fields on W3 (rho) levels
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'rhokh_bl',     w3_space, checkpoint_restart_flag )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'moist_flux_bl',     w3_space, checkpoint_restart_flag )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'heat_flux_bl',     w3_space, checkpoint_restart_flag )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'dtrdz_uv_bl',     w3_space, checkpoint_restart_flag )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'rdz_tq_bl',     w3_space, checkpoint_restart_flag )

    !========================================================================
    ! Fields owned by the convection scheme
    !========================================================================
    convection_fields = field_collection_type(name='convection_fields')

    ! 2D fields, need checkpointing
    checkpoint_restart_flag = .true.
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'conv_rain',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'conv_snow',  twod_space, checkpoint_restart_flag, twod=.true. )

    ! 3D fields, need checkpointing
    call add_physics_field(convection_fields, depository, prognostic_fields, &
      'cca', wtheta_space, checkpoint_restart_flag)
    call add_physics_field(convection_fields, depository, prognostic_fields, &
      'ccw', wtheta_space, checkpoint_restart_flag)

    ! 2D fields, don't need checkpointing
    checkpoint_restart_flag = .false.
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'shallow_flag',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'uw0_flux',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'vw0_flux',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'lcl_height',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'parcel_top',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'level_parcel_top',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'wstar',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'thv_flux',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'parcel_buoyancy',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'qsat_at_lcl',  twod_space, checkpoint_restart_flag, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'dt_conv', wtheta_space, checkpoint_restart_flag )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'dmv_conv', wtheta_space, checkpoint_restart_flag )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'dmcl_conv', wtheta_space, checkpoint_restart_flag )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'dmcf_conv', wtheta_space, checkpoint_restart_flag )
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      'dcfl_conv', wtheta_space, checkpoint_restart_flag)
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      'dcff_conv', wtheta_space, checkpoint_restart_flag)
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      'dbcf_conv', wtheta_space, checkpoint_restart_flag)

    ! 3D fields on W3 (rho) levels
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'du_conv', w3_space, checkpoint_restart_flag )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'dv_conv', w3_space, checkpoint_restart_flag )

    !========================================================================
    ! Fields owned by the cloud scheme
    !========================================================================
    cloud_fields = field_collection_type(name='cloud_fields')

    ! 3D fields, might need checkpointing
    if (cloud == cloud_um) then
      checkpoint_restart_flag = .true.
    else
      checkpoint_restart_flag = .false.
    end if
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'area_fraction',   wtheta_space, checkpoint_restart_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'rh_crit',     wtheta_space, checkpoint_restart_flag)

    ! 3D fields, might need advecting
    if ( scheme == scheme_pc2 ) then
      advection_flag=.true.
    else
      advection_flag=.false.
    endif
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'liquid_fraction', wtheta_space, checkpoint_restart_flag,         &
      advection_flag=advection_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'ice_fraction',    wtheta_space, checkpoint_restart_flag,         &
      advection_flag=advection_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'bulk_fraction',   wtheta_space, checkpoint_restart_flag,         &
      advection_flag=advection_flag)

    !========================================================================
    ! Fields owned by the surface exchange scheme
    !========================================================================
    surface_fields = field_collection_type(name='surface_fields')

    ! 2D fields, need checkpointing
    checkpoint_restart_flag = .true.
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'tstar',   twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'z0msea',  twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'surface_conductance', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'chloro_sea', twod_space, checkpoint_restart_flag, twod=.true. )

    ! Fields on surface tiles, need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'tile_fraction', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'tile_temperature', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'canopy_water', tile_space, checkpoint_restart_flag, twod=.true. )

    ! Fields on plant functional types, need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'leaf_area_index', pft_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'canopy_height', pft_space, checkpoint_restart_flag, twod=.true. )

    ! Sea-ice category fields, need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'sea_ice_thickness', sice_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'sea_ice_pond_frac', sice_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'sea_ice_pond_depth', sice_space, checkpoint_restart_flag, twod=.true. )

    ! Fields on surface tiles, don't need checkpointing
    checkpoint_restart_flag = .false.
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'tile_heat_flux', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'tile_moisture_flux', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'alpha1_tile', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'ashtf_prime_tile', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'dtstar_tile', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'fraca_tile', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'z0h_tile', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'z0m_tile', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'rhokh_tile', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'chr1p5m_tile', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'resfs_tile', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'canhc_tile', tile_space, checkpoint_restart_flag, twod=.true. )

    ! 2D fields
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'ustar', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field(surface_fields, depository, prognostic_fields,      &
      'net_prim_prod', twod_space, checkpoint_restart_flag, twod=.true.)

    ! This needs to be order 3 as it is land tile (9) x soil levels (4) = 36
    vector_space => function_space_collection%get_fs(twod_mesh_id, 3, W3)
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'tile_water_extract', vector_space, checkpoint_restart_flag, twod=.true. )

    !========================================================================
    ! Fields owned by the soil hydrology scheme
    !========================================================================
    soil_fields = field_collection_type(name='soil_fields')

    ! 2D fields, need checkpointing
    checkpoint_restart_flag = .true.
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_albedo', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_roughness', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_thermal_cond', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_carbon_content', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'decrease_sath_cond', twod_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'mean_topog_index', twod_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'a_sat_frac', twod_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'c_sat_frac', twod_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'a_wet_frac', twod_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'c_wet_frac', twod_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'soil_sat_frac', twod_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'soil_wet_frac', twod_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'water_table', twod_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'wetness_under_soil', twod_space, checkpoint_restart_flag, twod=.true.)

    ! Fields on soil levels
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_moist_wilt', soil_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_moist_crit', soil_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_moist_sat', soil_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_cond_sat', soil_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_thermal_cap', soil_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_suction_sat', soil_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'clapp_horn_b', soil_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_temperature', soil_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_moisture', soil_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'unfrozen_soil_moisture', soil_space, checkpoint_restart_flag,          &
      twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'frozen_soil_moisture', soil_space, checkpoint_restart_flag, twod=.true. )

    ! 2D fields, don't need checkpointing
    checkpoint_restart_flag = .false.
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_moist_avail', twod_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'soil_respiration', twod_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'thermal_cond_wet_soil', twod_space, checkpoint_restart_flag, twod=.true.)

    !========================================================================
    ! Fields owned by the snow scheme
    !========================================================================
    snow_fields = field_collection_type(name='snow_fields')

    ! Fields on surface tiles, need checkpointing
    checkpoint_restart_flag = .true.
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      'tile_snow_mass', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      'tile_snow_rgrain', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      'n_snow_layers', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      'snow_depth', tile_space, checkpoint_restart_flag, twod=.true. )
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      'snowpack_density', tile_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      'snow_under_canopy', tile_space, checkpoint_restart_flag, twod=.true.)

    ! Fields on snow layers
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      'snow_layer_thickness', snow_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      'snow_layer_ice_mass', snow_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      'snow_layer_liq_mass', snow_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      'snow_layer_temp', snow_space, checkpoint_restart_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      'snow_layer_rgrain', snow_space, checkpoint_restart_flag, twod=.true.)

    ! 2D fields
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      'snow_soot', twod_space, checkpoint_restart_flag, twod=.true. )

    !========================================================================
    ! Fields owned by the aerosol scheme
    !========================================================================
    aerosol_fields = field_collection_type(name='aerosol_fields')

    ! 3D fields, might need checkpointing
    if ( ( aerosol == aerosol_um ) .and.                                       &
         ( c_aerosol == c_aerosol_glomap_mode_climatology ) ) then
      checkpoint_restart_flag = .true.
    else
      checkpoint_restart_flag = .false.
    end if
    ! Nucleation Soluble number density
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'nd_nuc_sol', wtheta_space, checkpoint_restart_flag )
    ! Nucleation Soluble sulphate aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'nuc_sol_su', wtheta_space, checkpoint_restart_flag )
    ! Nucleation Soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'nuc_sol_oc', wtheta_space, checkpoint_restart_flag )
    ! Aitken Soluble number density
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'nd_ait_sol', wtheta_space, checkpoint_restart_flag )
    ! Aitken Soluble sulphate aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'ait_sol_su', wtheta_space, checkpoint_restart_flag )
    ! Aitken Soluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'ait_sol_bc', wtheta_space, checkpoint_restart_flag )
    ! Aitken Soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'ait_sol_oc', wtheta_space, checkpoint_restart_flag )
    ! Accumulation Soluble number density
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'nd_acc_sol', wtheta_space, checkpoint_restart_flag )
    ! Accumulation Soluble sulphate aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'acc_sol_su', wtheta_space, checkpoint_restart_flag )
    ! Accumulation Soluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'acc_sol_bc', wtheta_space, checkpoint_restart_flag )
    ! Accumulation Soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'acc_sol_oc', wtheta_space, checkpoint_restart_flag )
    ! Accumulation Soluble sea salt aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'acc_sol_ss', wtheta_space, checkpoint_restart_flag )
    ! Coarse Soluble number density
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'nd_cor_sol', wtheta_space, checkpoint_restart_flag )
    ! Coarse Soluble sulphate aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'cor_sol_su', wtheta_space, checkpoint_restart_flag )
    ! Coarse Soluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'cor_sol_bc', wtheta_space, checkpoint_restart_flag )
    ! Coarse Soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'cor_sol_oc', wtheta_space, checkpoint_restart_flag )
    ! Coarse Soluble sea salt aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'cor_sol_ss', wtheta_space, checkpoint_restart_flag )
    ! Aitken Insoluble number density
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'nd_ait_ins', wtheta_space, checkpoint_restart_flag )
    ! Aitken Insoluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'ait_ins_bc', wtheta_space, checkpoint_restart_flag )
    ! Aitken Insoluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'ait_ins_oc', wtheta_space, checkpoint_restart_flag )

    !========================================================================
    ! Here we create CDNC (and later RADAER) fields
    !========================================================================
    ! No checkpoint restart for these fields
    checkpoint_restart_flag = .false.
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'cloud_drop_no_conc', wtheta_space, checkpoint_restart_flag )

  end subroutine create_physics_prognostics

  !>@brief Add field to field collection and set its write,
  !>       checkpoint and restart behaviour
  !> @param[in,out] field_collection Field collection that 'name' will be added to
  !> @param[in]     name             Name of field to be added to collection
  !> @param[in]     vector_space     Function space of field to set behaviour for
  !> @param[in]     checkpoint_restart_flag  Flag to allow checkpoint and
  !>                                        restart behaviour of field to be set
  !> @param[in]     twod             Optional flag to determine if this is a
  !>                                        2D field for diagnostic output
  !> @param[in]     advection_flag   Optional flag whether this field is to be advected

  subroutine add_physics_field(field_collection, &
                               depository, prognostic_fields, &
                               name, vector_space, &
                               checkpoint_restart_flag, twod, advection_flag)

    use io_config_mod,      only : use_xios_io, &
                                   write_diag
    use io_mod,             only : xios_write_field_face,        &
                                   xios_write_field_single_face, &
                                   checkpoint_write_xios,        &
                                   checkpoint_write_netcdf,      &
                                   checkpoint_read_netcdf,       &
                                   checkpoint_read_xios

    implicit none

    character(*), intent(in)                       :: name
    type(field_collection_type), intent(inout)     :: field_collection
    type(field_collection_type), intent(inout)     :: depository
    type(field_collection_type), intent(inout)     :: prognostic_fields
    type(function_space_type), pointer, intent(in) :: vector_space
    logical(l_def), intent(in)                     :: checkpoint_restart_flag
    logical(l_def), optional, intent(in)           :: twod
    logical(l_def), optional, intent(in)           :: advection_flag
    !Local variables
    type(field_type)                               :: new_field

    ! pointers for xios write interface
    procedure(write_interface), pointer   :: write_diag_behaviour => null()
    procedure(checkpoint_write_interface), pointer  :: checkpoint_write_behaviour => null()
    procedure(checkpoint_read_interface), pointer   :: checkpoint_read_behaviour => null()

    if ( use_xios_io) then
      checkpoint_write_behaviour => checkpoint_write_xios
      checkpoint_read_behaviour => checkpoint_read_xios
    else
      checkpoint_write_behaviour => checkpoint_write_netcdf
      checkpoint_read_behaviour => checkpoint_read_netcdf
    endif

    if (present(advection_flag)) then
      call new_field%initialise( vector_space, name=trim(name), advection_flag=advection_flag )
    else
      call new_field%initialise( vector_space, name=trim(name) )
    end if

    if (use_xios_io .and. write_diag) then
      ! All physics fields currently require output on faces...
      write_diag_behaviour => xios_write_field_face
      if (present(twod))then
        if (twod) write_diag_behaviour => xios_write_field_single_face
      end if
      call new_field%set_write_behaviour(write_diag_behaviour)
    end if
    if ( checkpoint_restart_flag ) then
      call new_field%set_checkpoint_write_behaviour(checkpoint_write_behaviour)
      call new_field%set_checkpoint_read_behaviour(checkpoint_read_behaviour)
    endif

    ! Add the field to the depository
    call depository%add_field(new_field)
    ! Put a pointer to the field in the required collection
    call field_collection%add_reference_to_field( depository%get_field(name) )
    ! If checkpointing the field, put a pointer to it in the prognostics collection
    if ( checkpoint_restart_flag ) then
      call prognostic_fields%add_reference_to_field( depository%get_field(name) )
    endif

  end subroutine add_physics_field

end module create_physics_prognostics_mod
