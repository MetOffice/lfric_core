!-------------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief create physics prognostics
!> @details Creates the physics prognostic fields
module create_physics_prognostics_mod

  use clock_mod,                      only : clock_type
  use constants_mod,                  only : i_def, l_def
  use field_mod,                      only : field_type
  use field_parent_mod,               only : write_interface, read_interface,  &
                                             checkpoint_write_interface,       &
                                             checkpoint_read_interface
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use field_collection_mod,           only : field_collection_type
  use fs_continuity_mod,              only : W2, W3, Wtheta
  use function_space_mod,             only : function_space_type
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_ERROR
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use radiation_config_mod,           only : n_radstep, cloud_representation,  &
                                             cloud_representation_combined,    &
                                             cloud_representation_conv_strat_liq_ice, &
                                             cloud_representation_split,       &
                                             l_trans_zen_correction
  use aerosol_config_mod,             only : c_aerosol,                        &
                                             c_aerosol_glomap_mode_climatology
  use section_choice_config_mod,      only : cloud, cloud_um,                  &
                                             aerosol, aerosol_um,              &
                                             radiation, radiation_socrates,    &
                                             boundary_layer,                   &
                                             boundary_layer_um,                &
                                             surface, surface_jules,           &
                                             orographic_drag,                  &
                                             orographic_drag_um,               &
                                             convection, convection_um
  use cloud_config_mod,               only : scheme, &
                                             scheme_pc2
  use surface_config_mod,             only : albedo_obs, srf_ex_cnv_gust,      &
                                             sea_alb_var_chl
  use spectral_gwd_config_mod,        only : add_cgw

  implicit none

  private
  public :: create_physics_prognostics

contains
  !>@brief Routine to initialise the field objects required by the physics
  !> @param[in]    mesh_id The identifier given to the current 3d mesh
  !> @param[in]    twod_mesh_id The identifier given to the current 2d mesh
  !> @param[in]    clock Model time.
  !> @param[in,out] depository Main collection of all fields in memory
  !> @param[in,out] prognostic_fields The prognostic variables in the model
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
  subroutine create_physics_prognostics( mesh_id,             &
                                         twod_mesh_id,        &
                                         clock,               &
                                         depository,          &
                                         prognostic_fields,   &
                                         derived_fields,      &
                                         radiation_fields,    &
                                         microphysics_fields, &
                                         orography_fields,    &
                                         turbulence_fields,   &
                                         convection_fields,   &
                                         cloud_fields,        &
                                         surface_fields,      &
                                         soil_fields,         &
                                         snow_fields,         &
                                         aerosol_fields )

#ifdef UM_PHYSICS
    use jules_control_init_mod,  only: n_surf_tile, n_sea_ice_tile, soil_lev_tile
    use jules_physics_init_mod,  only: snow_lev_tile
    use jules_surface_types_mod, only: npft
    use nlsizes_namelist_mod,    only: sm_levels
#endif

    implicit none

    integer(i_def),    intent(in) :: mesh_id
    integer(i_def),    intent(in) :: twod_mesh_id
    class(clock_type), intent(in) :: clock

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
#ifdef UM_PHYSICS
    type(function_space_type), pointer :: vector_space => null()
    type(function_space_type), pointer :: twod_space => null()
    type(function_space_type), pointer :: surft_space => null()
    type(function_space_type), pointer :: pft_space => null()
    type(function_space_type), pointer :: soil_space => null()
    type(function_space_type), pointer :: sice_space => null()
    type(function_space_type), pointer :: snow_space => null()
#endif
    type(function_space_type), pointer :: wtheta_space => null()
    type(function_space_type), pointer :: w3_space => null()
    type(function_space_type), pointer :: w2_space => null()

    type( field_type ), pointer :: theta => null()

    integer(i_def) :: theta_space
#ifdef UM_PHYSICS
    logical(l_def) :: checkpoint_flag
    logical(l_def) :: advection_flag
#endif

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
#ifdef UM_PHYSICS
    twod_space => function_space_collection%get_fs(twod_mesh_id, 0, W3)
    surft_space => function_space_collection%get_fs(twod_mesh_id, 0, W3, n_surf_tile)
    pft_space => function_space_collection%get_fs(twod_mesh_id, 0, W3, npft)
    soil_space => function_space_collection%get_fs(twod_mesh_id, 0, W3, sm_levels)
    sice_space => function_space_collection%get_fs(twod_mesh_id, 0, W3, n_sea_ice_tile)
    snow_space => function_space_collection%get_fs(twod_mesh_id, 0, W3, snow_lev_tile)
#endif
    !========================================================================
    ! Fields derived from the FE dynamical fields for use in physics
    !========================================================================
    derived_fields  =  field_collection_type(name='derived_fields')

    ! Wtheta fields
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'w_physics',      wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'rho_in_wth',     wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'wetrho_in_wth',  wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'exner_in_wth',   wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'w_physics_star', wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'theta_star',     wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'shear',          wtheta_space )

    ! W3 fields
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u1_in_w3',      w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u2_in_w3',      w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u3_in_w3',      w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'theta_in_w3',   w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'wetrho_in_w3',  w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u1_in_w3_star', w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u2_in_w3_star', w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u3_in_w3_star', w3_space )

    ! W2 fields
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u_physics',      w2_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u_physics_star', w2_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      'u_star',         w2_space )

#ifdef UM_PHYSICS
    !========================================================================
    ! Fields owned by the radiation scheme
    !========================================================================
    radiation_fields = field_collection_type(name='radiation_fields')

    ! 2D fields, might need checkpointing
    if (surface == surface_jules .and. albedo_obs) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'albedo_obs_vis', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'albedo_obs_nir', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! 3D fields, need checkpointing
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'ozone', wtheta_space, checkpoint_flag=.true. )

    ! 2D fields, don't need checkpointing
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'cos_zenith_angle',   twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lit_fraction',       twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lw_down_surf', twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_down_surf', twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_direct_surf', twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_down_blue_surf', twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_direct_blue_surf', twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lw_up_tile', surft_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_up_tile', surft_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_up_blue_tile', surft_space, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_heating_rate', wtheta_space )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lw_heating_rate', wtheta_space )

    ! Fields which need checkpointing for radiation timestepping
    !
    !> @todo There is probably a better way of handling this test which doesn't
    !>       involving passing the clock down here.
    !>
    if (radiation == radiation_socrates) then
      ! Checkpoint unless both the first timestep of this run and the
      ! first timestep of the next run are radiation timesteps
      checkpoint_flag = &
        mod(clock%get_first_step()-1, n_radstep) /= 0 .or. &
        mod(clock%get_last_step(),    n_radstep) /= 0
    else
      checkpoint_flag = .false.
    end if

    ! 2D fields
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lw_down_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_down_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_direct_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_direct_toa_rts', twod_space,                                         &
      checkpoint_flag=(checkpoint_flag .and. l_trans_zen_correction), twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_down_blue_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_direct_blue_surf_rts', twod_space, checkpoint_flag=checkpoint_flag,  &
      twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'cos_zenith_angle_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lit_fraction_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'stellar_irradiance_rts', twod_space, checkpoint_flag=checkpoint_flag,   &
      twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sin_stellar_declination_rts', twod_space,                               &
      checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'stellar_eqn_of_time_rts', twod_space, checkpoint_flag=checkpoint_flag,  &
      twod=.true. )

    ! 3D fields
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_heating_rate_rts', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lw_heating_rate_rts', wtheta_space, checkpoint_flag=checkpoint_flag )

    ! Fields on surface tiles
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'lw_up_tile_rts', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_up_tile_rts', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      'sw_up_blue_tile_rts', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    !========================================================================
    ! Fields owned by the microphysics scheme
    !========================================================================
    microphysics_fields = field_collection_type(name='microphysics_fields')

    ! 2D fields, don't need checkpointing
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      'ls_rain',  twod_space, twod=.true. )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      'ls_snow',  twod_space, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      'dtl_mphys', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      'dmt_mphys', wtheta_space )

    !========================================================================
    ! Fields owned by the orographic drag schemes
    !========================================================================
    orography_fields = field_collection_type(name='orography_fields')

    ! 2D fields, might need checkpointing
    if (boundary_layer == boundary_layer_um .or. &
         surface == surface_jules           .or. &
         orographic_drag == orographic_drag_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      'sd_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      'grad_xx_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      'grad_xy_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      'grad_yy_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      'peak_to_trough_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      'silhouette_area_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    !========================================================================
    ! Fields owned by the turbulence scheme
    !========================================================================
    turbulence_fields = field_collection_type(name='turbulence_fields')

    ! 2D fields, might need checkpointing
    if (boundary_layer == boundary_layer_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'zh',      twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! 2D fields, don't need checkpointing
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'ntml',    twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'cumulus', twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'z_lcl',  twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'inv_depth',  twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'qcl_at_inv_top',  twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'blend_height_tq',  twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'blend_height_uv',  twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'zh_nonloc',  twod_space, twod=.true. )

    ! Space for the three fields required to regrid this
    vector_space => function_space_collection%get_fs(twod_mesh_id, 0, W3, 5)
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'rhokm_surf',  vector_space, twod=.true. )
    ! Space for the 7 BL types
    vector_space => function_space_collection%get_fs(twod_mesh_id, 0, W3, 7)
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'bl_type_ind',  vector_space, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'dt_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'dmv_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'rhokm_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'ngstress_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'bq_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'bt_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'dtrdz_tq_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'rdz_uv_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'fd_taux', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'fd_tauy', wtheta_space )

    ! 3D fields on W3 (rho) levels
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'rhokh_bl',     w3_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'moist_flux_bl',     w3_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'heat_flux_bl',     w3_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'dtrdz_uv_bl',     w3_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'rdz_tq_bl',     w3_space )

    ! W2 fields, don't need checkpointing
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      'du_bl', w2_space )

    !========================================================================
    ! Fields owned by the convection scheme
    !========================================================================
    convection_fields = field_collection_type(name='convection_fields')

    ! 2D fields, might need checkpointing
    if (convection == convection_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'conv_rain',  twod_space,                                                &
      checkpoint_flag=(checkpoint_flag .and. add_cgw), twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'conv_snow',  twod_space,                                                &
      checkpoint_flag=(checkpoint_flag .and. add_cgw), twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'dd_mf_cb',  twod_space,                                                 &
      checkpoint_flag=(checkpoint_flag .and. srf_ex_cnv_gust), twod=.true. )

    ! 3D fields, might need checkpointing
    if (convection == convection_um) then
      select case (cloud_representation)
      case (cloud_representation_combined, &
            cloud_representation_conv_strat_liq_ice, &
            cloud_representation_split)
        checkpoint_flag = .true.
      case default
        checkpoint_flag = .false.
      end select
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field(convection_fields, depository, prognostic_fields, &
      'cca', wtheta_space, checkpoint_flag=checkpoint_flag)
    call add_physics_field(convection_fields, depository, prognostic_fields, &
      'ccw', wtheta_space, checkpoint_flag=checkpoint_flag)

    ! 2D fields, don't need checkpointing
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'shallow_flag',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'uw0_flux',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'vw0_flux',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'lcl_height',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'parcel_top',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'level_parcel_top',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'wstar',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'thv_flux',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'parcel_buoyancy',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'qsat_at_lcl',  twod_space, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'dt_conv', wtheta_space )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'dmv_conv', wtheta_space )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'dmcl_conv', wtheta_space )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'dmcf_conv', wtheta_space )
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      'dcfl_conv', wtheta_space )
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      'dcff_conv', wtheta_space )
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      'dbcf_conv', wtheta_space )

    ! 3D fields on W3 (rho) levels
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'du_conv', w3_space )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      'dv_conv', w3_space )

    !========================================================================
    ! Fields owned by the cloud scheme
    !========================================================================
    cloud_fields = field_collection_type(name='cloud_fields')

    ! 3D fields, might need checkpointing
    if (cloud == cloud_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'area_fraction',   wtheta_space,                                  &
      checkpoint_flag=(checkpoint_flag .and. radiation==radiation_socrates))

    ! 3D fields, might need advecting
    if ( scheme == scheme_pc2 ) then
      advection_flag=.true.
    else
      advection_flag=.false.
    endif
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'liquid_fraction', wtheta_space, checkpoint_flag=checkpoint_flag, &
      advection_flag=advection_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'ice_fraction',    wtheta_space, checkpoint_flag=checkpoint_flag, &
      advection_flag=advection_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'bulk_fraction',   wtheta_space, checkpoint_flag=checkpoint_flag, &
      advection_flag=advection_flag)

    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'rh_crit',     wtheta_space )
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'departure_exner_wth', wtheta_space, advection_flag=advection_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      'sigma_qcw',   wtheta_space )

    !========================================================================
    ! Fields owned by the surface exchange scheme
    !========================================================================
    surface_fields = field_collection_type(name='surface_fields')

    ! 2D fields, might need checkpointing
    if (surface == surface_jules) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'z0msea',  twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'surface_conductance', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'chloro_sea', twod_space,                                                &
      checkpoint_flag=(checkpoint_flag .and. sea_alb_var_chl), twod=.true. )

    ! Fields on surface tiles, might need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'tile_fraction', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'tile_temperature', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'canopy_water', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! Fields on plant functional types, might need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'leaf_area_index', pft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'canopy_height', pft_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! Sea-ice category fields, might need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'sea_ice_thickness', sice_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'sea_ice_temperature', sice_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! Fields on surface tiles, don't need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'tile_heat_flux', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'tile_moisture_flux', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'alpha1_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'ashtf_prime_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'dtstar_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'fraca_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'z0h_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'z0m_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'rhokh_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'chr1p5m_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'resfs_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'canhc_tile', surft_space, twod=.true. )

    ! 2D fields
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'ustar', twod_space, twod=.true. )
    call add_physics_field(surface_fields, depository, prognostic_fields,      &
      'net_prim_prod', twod_space, twod=.true.)

    ! Field on soil levels and land tiles
    vector_space => function_space_collection%get_fs(twod_mesh_id, 0, W3, soil_lev_tile)
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      'tile_water_extract', vector_space, twod=.true. )

    !========================================================================
    ! Fields owned by the soil hydrology scheme
    !========================================================================
    soil_fields = field_collection_type(name='soil_fields')

    ! 2D fields, might need checkpointing
    if (surface == surface_jules) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_albedo', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_roughness', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_thermal_cond', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_carbon_content', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'mean_topog_index', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'a_sat_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'c_sat_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'a_wet_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'c_wet_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'soil_sat_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'water_table', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'wetness_under_soil', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_moist_wilt', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_moist_crit', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_moist_sat', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_cond_sat', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_thermal_cap', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_suction_sat', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'clapp_horn_b', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! Fields on soil levels
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_temperature', soil_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_moisture', soil_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'unfrozen_soil_moisture', soil_space, checkpoint_flag=checkpoint_flag,  &
      twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'frozen_soil_moisture', soil_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! 2D fields, don't need checkpointing
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      'soil_moist_avail', twod_space, twod=.true. )
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'soil_respiration', twod_space, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      'thermal_cond_wet_soil', twod_space, twod=.true.)

    !========================================================================
    ! Fields owned by the snow scheme
    !========================================================================
    snow_fields = field_collection_type(name='snow_fields')

    ! Fields on surface tiles, might need checkpointing
    if (surface == surface_jules) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      'tile_snow_mass', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      'tile_snow_rgrain', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      'n_snow_layers', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      'snow_depth', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      'snowpack_density', surft_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      'snow_under_canopy', surft_space, checkpoint_flag=checkpoint_flag, twod=.true.)

    ! Fields on snow layers
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      'snow_layer_thickness', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      'snow_layer_ice_mass', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      'snow_layer_liq_mass', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      'snow_layer_temp', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      'snow_layer_rgrain', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)

    ! 2D fields
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      'snow_soot', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! Fields which don't need checkpointing
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      'snow_unload_rate', pft_space, twod=.true. )

    !========================================================================
    ! Fields owned by the aerosol scheme
    !========================================================================
    aerosol_fields = field_collection_type(name='aerosol_fields')

    ! 3D fields, might need checkpointing
    if ( ( aerosol == aerosol_um ) .and.                                       &
         ( c_aerosol == c_aerosol_glomap_mode_climatology ) ) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    ! Nucleation Soluble number density
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'nd_nuc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Nucleation Soluble sulphate aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'nuc_sol_su', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Nucleation Soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'nuc_sol_oc', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Aitken Soluble number density
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'nd_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Aitken Soluble sulphate aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'ait_sol_su', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Aitken Soluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'ait_sol_bc', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Aitken Soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'ait_sol_oc', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Accumulation Soluble number density
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'nd_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Accumulation Soluble sulphate aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'acc_sol_su', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Accumulation Soluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'acc_sol_bc', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Accumulation Soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'acc_sol_oc', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Accumulation Soluble sea salt aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'acc_sol_ss', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Coarse Soluble number density
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'nd_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Coarse Soluble sulphate aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'cor_sol_su', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Coarse Soluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'cor_sol_bc', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Coarse Soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'cor_sol_oc', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Coarse Soluble sea salt aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'cor_sol_ss', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Aitken Insoluble number density
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'nd_ait_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Aitken Insoluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'ait_ins_bc', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Aitken Insoluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'ait_ins_oc', wtheta_space, checkpoint_flag=checkpoint_flag )

    !========================================================================
    ! Aerosol fields that do not require checkpoint restart
    !========================================================================
    ! No checkpoint restart for these fields
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'cloud_drop_no_conc', wtheta_space )
    ! Sulphuric Acid aerosol MMR
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      'sulphuric', wtheta_space )
#endif

  end subroutine create_physics_prognostics

  !>@brief Add field to field collection and set its write,
  !>       checkpoint-restart and advection behaviour
  !> @param[in,out] field_collection  Field collection that 'name' will be added to
  !> @param[in,out] depository        Collection of all fields
  !> @param[in,out] prognostic_fields Collection of checkpointed fields
  !> @param[in]     name              Name of field to be added to collection
  !> @param[in]     vector_space      Function space of field to set behaviour for
  !> @param[in]     checkpoint_flag   Optional flag to allow checkpoint-
  !>                                   restart behaviour of field to be set
  !> @param[in]     twod              Optional flag to determine if this is a
  !>                                   2D field for diagnostic output
  !> @param[in]     advection_flag    Optional flag whether this field is to be advected

  subroutine add_physics_field(field_collection, &
                               depository, prognostic_fields, &
                               name, vector_space, &
                               checkpoint_flag, twod, advection_flag)

    use io_config_mod,           only : use_xios_io, &
                                        write_diag, checkpoint_write, &
                                        checkpoint_read
    use lfric_xios_read_mod,     only : read_field_face, &
                                        read_field_single_face
    use lfric_xios_write_mod,    only : write_field_face, &
                                        write_field_single_face
    use io_mod,                  only : checkpoint_write_netcdf, &
                                        checkpoint_read_netcdf

    implicit none

    character(*), intent(in)                       :: name
    type(field_collection_type), intent(inout)     :: field_collection
    type(field_collection_type), intent(inout)     :: depository
    type(field_collection_type), intent(inout)     :: prognostic_fields
    type(function_space_type), pointer, intent(in) :: vector_space
    logical(l_def), optional, intent(in)           :: checkpoint_flag
    logical(l_def), optional, intent(in)           :: twod
    logical(l_def), optional, intent(in)           :: advection_flag
    !Local variables
    type(field_type)                               :: new_field
    class(pure_abstract_field_type), pointer       :: field_ptr => null()
    logical(l_def)                                 :: twod_field, checkpointed

    ! pointers for xios write interface
    procedure(write_interface), pointer :: write_behaviour => null()
    procedure(read_interface),  pointer :: read_behaviour => null()
    procedure(checkpoint_write_interface), pointer :: checkpoint_write_behaviour => null()
    procedure(checkpoint_read_interface), pointer  :: checkpoint_read_behaviour => null()

    ! Create the new field
    if (present(advection_flag)) then
      call new_field%initialise( vector_space, name=trim(name), advection_flag=advection_flag )
    else
      call new_field%initialise( vector_space, name=trim(name) )
    end if

    ! Set checkpoint flag
    if (present(checkpoint_flag)) then
      checkpointed = checkpoint_flag
    else
      checkpointed = .false.
    end if

    ! Set read and write behaviour
    if (use_xios_io) then
      if (present(twod)) then
        twod_field = twod
      else
        twod_field = .false.
      end if
      if (twod_field) then
        write_behaviour => write_field_single_face
        read_behaviour  => read_field_single_face
      else
        write_behaviour => write_field_face
        read_behaviour  => read_field_face
      end if
      if (write_diag .or. checkpoint_write) &
        call new_field%set_write_behaviour(write_behaviour)
      if (checkpoint_read .and. checkpointed) &
        call new_field%set_read_behaviour(read_behaviour)
    else
      checkpoint_write_behaviour => checkpoint_write_netcdf
      checkpoint_read_behaviour  => checkpoint_read_netcdf
      call new_field%set_checkpoint_write_behaviour(checkpoint_write_behaviour)
      call new_field%set_checkpoint_read_behaviour(checkpoint_read_behaviour)
    endif

    ! Add the field to the depository
    call depository%add_field(new_field)
    field_ptr => depository%get_field(name)
    ! Put a pointer to the field in the required collection
    call field_collection%add_reference_to_field( field_ptr )
    ! If checkpointing the field, put a pointer to it in the prognostics collection
    if ( checkpointed ) then
      call prognostic_fields%add_reference_to_field( field_ptr )
    endif

  end subroutine add_physics_field

end module create_physics_prognostics_mod
