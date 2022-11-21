!-------------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief create physics prognostics
!> @details Creates the physics prognostic fields
module create_physics_prognostics_mod

  use checks_config_mod,              only : energy_correction,     &
                                             energy_correction_none
  use clock_mod,                      only : clock_type
  use constants_mod,                  only : i_def, l_def
  use field_mod,                      only : field_type
  use integer_field_mod,              only : integer_field_type
  use field_parent_mod,               only : write_interface, read_interface,  &
                                             checkpoint_write_interface,       &
                                             checkpoint_read_interface
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use field_collection_mod,           only : field_collection_type
  use fs_continuity_mod,              only : W2, W3, Wtheta, W2H
  use function_space_mod,             only : function_space_type
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_ERROR
  use mesh_mod,                       only : mesh_type
  use mixing_config_mod,              only : smagorinsky
  use physics_config_mod,             only : stochastic_physics_placement,     &
                                             stochastic_physics_placement_fast
  use pure_abstract_field_mod,        only : pure_abstract_field_type
  use radiation_config_mod,           only : n_radstep, cloud_representation,  &
                                             cloud_representation_combined,    &
                                             cloud_representation_conv_strat_liq_ice, &
                                             cloud_representation_split,       &
                                             topography, topography_flat,      &
                                             topography_horizon,               &
                                             n_horiz_layer, n_horiz_ang,       &
                                             l_inc_radstep
  use aerosol_config_mod,             only : glomap_mode,                      &
                                             glomap_mode_climatology,          &
                                             glomap_mode_ukca
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
  use jules_surface_config_mod,       only : srf_ex_cnv_gust
  use surface_config_mod,             only : albedo_obs, sea_alb_var_chl
  use spectral_gwd_config_mod,        only : add_cgw
  use microphysics_config_mod,        only : turb_gen_mixph
  use derived_config_mod,             only : l_esm_couple
  use esm_couple_config_mod,          only : l_esm_couple_test
#ifdef UM_PHYSICS
  use multidata_field_dimensions_mod, only :                                   &
       get_ndata_val => get_multidata_field_dimension
  use cv_run_mod,                     only:  l_conv_prog_precip,               &
                                             l_conv_prog_dtheta,               &
                                             l_conv_prog_dq,                   &
                                             adv_conv_prog_dtheta,             &
                                             adv_conv_prog_dq
#endif

  implicit none

  private
  public :: create_physics_prognostics

contains

  !>@brief Routine to initialise the field objects required by the physics
  !> @param[in]    mesh      The current 3d mesh
  !> @param[in]    twod_mesh The current 2d mesh
  !> @param[in]    clock Model time.
  !> @param[in,out] depository Main collection of all fields in memory
  !> @param[in,out] prognostic_fields The prognostic variables in the model
  !> @param[out]   adv_fields_all_outer Collection of fields that need to be advected every outer iteration
  !> @param[out]   adv_fields_last_outer Collection of fields that need to be advected at final outer iteration
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
  !> @param[out]   chemistry_fields Collection of fields for chemistry scheme
  subroutine create_physics_prognostics( mesh,                  &
                                         twod_mesh,             &
                                         clock,                 &
                                         depository,            &
                                         prognostic_fields,     &
                                         adv_fields_all_outer,  &
                                         adv_fields_last_outer, &
                                         derived_fields,        &
                                         radiation_fields,      &
                                         microphysics_fields,   &
                                         orography_fields,      &
                                         turbulence_fields,     &
                                         convection_fields,     &
                                         cloud_fields,          &
                                         surface_fields,        &
                                         soil_fields,           &
                                         snow_fields,           &
                                         chemistry_fields,      &
                                         aerosol_fields )

    implicit none

    type( mesh_type ), intent(in), pointer :: mesh
    type( mesh_type ), intent(in), pointer :: twod_mesh

    class(clock_type), intent(in) :: clock

    ! Collections of fields
    type(field_collection_type), intent(inout) :: depository
    type(field_collection_type), intent(inout) :: prognostic_fields
    type(field_collection_type), intent(inout) :: adv_fields_all_outer
    type(field_collection_type), intent(inout) :: adv_fields_last_outer
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
    type(field_collection_type), intent(out) :: chemistry_fields
    type(field_collection_type), intent(out) :: aerosol_fields

    ! pointers to vector spaces
#ifdef UM_PHYSICS
    type(function_space_type), pointer :: vector_space => null()
    type(function_space_type), pointer :: surft_space => null()
    type(function_space_type), pointer :: pft_space => null()
    type(function_space_type), pointer :: soil_space => null()
    type(function_space_type), pointer :: sice_space => null()
    type(function_space_type), pointer :: snow_space => null()
    type(function_space_type), pointer :: h_ang_space => null()
    type(function_space_type), pointer :: h_asp_space => null()
#endif
    type(function_space_type), pointer :: twod_space => null()
    type(function_space_type), pointer :: wtheta_space => null()
    type(function_space_type), pointer :: w3_space => null()
    type(function_space_type), pointer :: w2_space => null()
    type(function_space_type), pointer :: w2h_space => null()

    type( field_type ), pointer :: theta => null()

    integer(i_def) :: theta_space
#ifdef UM_PHYSICS
    logical(l_def) :: checkpoint_flag
    logical(l_def) :: checkpoint_couple
    logical(l_def) :: advection_flag
#endif

    call log_event( 'Create physics prognostics', LOG_LEVEL_INFO )

    call prognostic_fields%get_field('theta', theta)
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
    wtheta_space => function_space_collection%get_fs(mesh, 0, Wtheta)
    w3_space => function_space_collection%get_fs(mesh, 0, W3)
    w2_space => function_space_collection%get_fs(mesh, 0, W2)
    w2h_space => function_space_collection%get_fs(mesh, 0, W2H)
    twod_space  => function_space_collection%get_fs(twod_mesh, 0, W3)
#ifdef UM_PHYSICS
    surft_space => function_space_collection%get_fs(twod_mesh, 0, W3,          &
         get_ndata_val('surface_tiles'))
    pft_space   => function_space_collection%get_fs(twod_mesh, 0, W3,          &
         get_ndata_val('plant_types'))
    soil_space  => function_space_collection%get_fs(twod_mesh, 0, W3,          &
         get_ndata_val('soil_levels'))
    sice_space  => function_space_collection%get_fs(twod_mesh, 0, W3,          &
         get_ndata_val('sea_ice_tiles'))
    snow_space  => function_space_collection%get_fs(twod_mesh, 0, W3,          &
         get_ndata_val('snow_layers_and_tiles'))
#endif

    !========================================================================
    ! Fields derived from the FE dynamical fields for use in physics
    !========================================================================
    call derived_fields%initialise(name='derived_fields', table_len=100)

    ! Wtheta fields
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'velocity_w2v',      wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'w_in_wth',      wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'rho_in_wth',     wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'wetrho_in_wth',  wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'exner_in_wth',   wtheta_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'exner_wth_n',   wtheta_space )

    if ( boundary_layer == boundary_layer_um .or. &
         convection     == convection_um ) then

      call add_physics_field( derived_fields, depository, prognostic_fields,     &
        adv_fields_last_outer, &
        'theta_star',     wtheta_space )

    end if

    if ( boundary_layer == boundary_layer_um .or. &
         convection     == convection_um     .or. &
         smagorinsky ) then

      call add_physics_field( derived_fields, depository, prognostic_fields,     &
        adv_fields_last_outer, &
        'shear',          wtheta_space )

    end if

    ! W3 fields
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'u_in_w3',      w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'v_in_w3',      w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'w_in_w3',      w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'theta_in_w3',   w3_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'wetrho_in_w3',  w3_space )

    if ( boundary_layer               == boundary_layer_um .or.                  &
         convection                   == convection_um     .or.                  &
         stochastic_physics_placement == stochastic_physics_placement_fast ) then

      call add_physics_field( derived_fields, depository, prognostic_fields,     &
        adv_fields_last_outer, &
        'u_in_w3_star', w3_space )
      call add_physics_field( derived_fields, depository, prognostic_fields,     &
        adv_fields_last_outer, &
        'v_in_w3_star', w3_space )
      call add_physics_field( derived_fields, depository, prognostic_fields,     &
        adv_fields_last_outer, &
        'w_in_w3_star', w3_space )

    end if

    ! W2 fields
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'u_physics',      w2_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'u_physics_star', w2_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'u_star',         w2_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'wetrho_in_w2',   w2_space )

    ! W2H fields
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'u_in_w2h',      w2h_space )
    call add_physics_field( derived_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'v_in_w2h',      w2h_space )

    ! 2D fields
    if ( energy_correction /= energy_correction_none ) then
      call add_physics_field( derived_fields, depository,               &
                              prognostic_fields, adv_fields_last_outer, &
                              'temp_correction_field',                  &
                              twod_space, twod = .true.,                &
                              checkpoint_flag = .true. )

      call add_physics_field( derived_fields, depository,               &
                              prognostic_fields, adv_fields_last_outer, &
                              'accumulated_fluxes',                     &
                              twod_space, twod = .true. )
    end if

#ifdef UM_PHYSICS
    !========================================================================
    ! Fields owned by the radiation scheme
    !========================================================================
    call radiation_fields%initialise(name='radiation_fields', table_len=100)

    ! 2D fields, might need checkpointing
    if (surface == surface_jules .and. albedo_obs) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'albedo_obs_vis', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'albedo_obs_nir', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! 3D fields, need checkpointing
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'ozone', wtheta_space, checkpoint_flag=.true. )

    ! 2D fields, don't need checkpointing
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'lw_down_surf', twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_down_surf', twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_direct_surf', twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_down_blue_surf', twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_direct_blue_surf', twod_space, twod=.true. )

    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'cos_zenith_angle',   twod_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'lit_fraction',       twod_space, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_heating_rate', wtheta_space )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'lw_heating_rate', wtheta_space )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'dtheta_rad', wtheta_space )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'dmv_pc2_rad', wtheta_space )

    ! Fields on surface tiles, don't need checkpointing
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'lw_up_tile', surft_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_up_tile', surft_space, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_up_blue_tile', surft_space, twod=.true. )

    ! Fields that need checkpointing for the topographic correction scheme
    if (radiation == radiation_socrates .and. &
        topography /= topography_flat) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'slope_angle', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'slope_aspect', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'skyview', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    if (radiation == radiation_socrates .and. &
        topography == topography_horizon) then
      checkpoint_flag = .true.
      h_ang_space => function_space_collection%get_fs( &
                       twod_mesh, 0, W3, n_horiz_layer*n_horiz_ang)
      h_asp_space => function_space_collection%get_fs( &
                       twod_mesh, 0, W3, n_horiz_ang)
    else
      checkpoint_flag = .false.
      h_ang_space => function_space_collection%get_fs(twod_mesh, 0, W3, 1)
      h_asp_space => function_space_collection%get_fs(twod_mesh, 0, W3, 1)
    end if
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'horizon_angle', h_ang_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'horizon_aspect', h_asp_space, checkpoint_flag=checkpoint_flag, twod=.true. )

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
      adv_fields_last_outer, &
      'lw_down_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_down_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_direct_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_down_blue_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_direct_blue_surf_rts', twod_space, checkpoint_flag=checkpoint_flag,  &
      twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'lw_up_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_up_surf_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'lw_up_toa_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_up_toa_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_direct_toa_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'cos_zenith_angle_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'lit_fraction_rts', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'stellar_irradiance_rts', twod_space, checkpoint_flag=checkpoint_flag,   &
      twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sin_stellar_declination_rts', twod_space,                               &
      checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'stellar_eqn_of_time_rts', twod_space, checkpoint_flag=checkpoint_flag,  &
      twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, 'orographic_correction_rts', twod_space,          &
      checkpoint_flag=checkpoint_flag, twod=.true. )

    ! 3D fields
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_heating_rate_rts', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'lw_heating_rate_rts', wtheta_space, checkpoint_flag=checkpoint_flag )

    ! Fields on surface tiles
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'lw_up_tile_rts', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_up_tile_rts', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_up_blue_tile_rts', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )


    ! Fields which need checkpointing for radiation incremental timestepping
    checkpoint_flag = checkpoint_flag .and. l_inc_radstep

    ! 2D fields
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'lw_down_surf_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_down_surf_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_direct_surf_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_down_blue_surf_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_direct_blue_surf_rtsi', twod_space, checkpoint_flag=checkpoint_flag,  &
      twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'lw_up_surf_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_up_surf_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'lw_up_toa_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_up_toa_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_direct_toa_rtsi', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! 3D fields
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_heating_rate_rtsi', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'lw_heating_rate_rtsi', wtheta_space, checkpoint_flag=checkpoint_flag )

    ! Fields on surface tiles
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'lw_up_tile_rtsi', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_up_tile_rtsi', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( radiation_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sw_up_blue_tile_rtsi', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )


    !========================================================================
    ! Fields owned by the microphysics scheme
    !========================================================================
    call microphysics_fields%initialise(name='microphysics_fields', table_len=100)

    ! 2D fields, don't need checkpointing
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      adv_fields_last_outer, &
      'ls_rain',  twod_space, twod=.true. )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      adv_fields_last_outer, &
      'ls_snow',  twod_space, twod=.true. )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      adv_fields_last_outer, &
      'lsca_2d',  twod_space, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      adv_fields_last_outer, &
      'dtheta_mphys', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      adv_fields_last_outer, &
      'dmv_mphys', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      adv_fields_last_outer, &
      'ls_rain_3d', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      adv_fields_last_outer, &
      'ls_snow_3d', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      adv_fields_last_outer, &
      'autoconv', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      adv_fields_last_outer, &
      'accretion', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      adv_fields_last_outer, &
      'rim_cry', wtheta_space )
    call add_physics_field( microphysics_fields, depository, prognostic_fields,&
      adv_fields_last_outer, &
      'rim_agg', wtheta_space )

    !========================================================================
    ! Fields owned by the orographic drag schemes
    !========================================================================
    call orography_fields%initialise(name='orography_fields', table_len=100)

    ! 2D fields, might need checkpointing
    if (boundary_layer == boundary_layer_um .or. &
         surface == surface_jules           .or. &
         orographic_drag == orographic_drag_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'sd_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'grad_xx_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'grad_xy_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'grad_yy_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'peak_to_trough_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'silhouette_area_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! Mean X, Y Gradients for the orographic correction to radiation
    checkpoint_flag = .false.
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'grad_x_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( orography_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'grad_y_orog', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )


    !========================================================================
    ! Fields owned by the turbulence scheme
    !========================================================================
    call turbulence_fields%initialise(name='turbulence_fields', table_len=100)

    ! 2D fields, might need checkpointing

    if (boundary_layer == boundary_layer_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'zh',      twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'wvar',    wtheta_space, checkpoint_flag=turb_gen_mixph )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'gradrinr', wtheta_space )

    ! 2D fields, don't need checkpointing
    call add_integer_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'ntml',    twod_space, twod=.true. )
    call add_integer_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'cumulus', twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'z_lcl',  twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'inv_depth',  twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'qcl_at_inv_top',  twod_space, twod=.true. )
    call add_integer_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'blend_height_tq',  twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'zh_nonloc',  twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'zhsc', twod_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'bl_weight_1dbl',  twod_space, twod=.true. )
    call add_integer_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'level_ent', twod_space, twod=.true. )
    call add_integer_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'level_ent_dsc', twod_space, twod=.true. )

    ! Space for the 7 BL types
    vector_space => function_space_collection%get_fs(twod_mesh, 0, W3,         &
         get_ndata_val('boundary_layer_types'))
    call add_integer_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'bl_type_ind',  vector_space, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'bq_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'bt_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'lmix_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'dsldzm',  wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'tke_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'dtrdz_tq_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'dw_bl', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'thetal_inc_leonard', wtheta_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'mt_inc_leonard', wtheta_space )

    ! 3D fields on W3 (rho) levels
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'moist_flux_bl',     w3_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'heat_flux_bl',     w3_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'rhokh_bl', w3_space )

    ! W2 fields, don't need checkpointing
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'rhokm_w2', w2_space )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'tau_w2', w2_space )

    ! Fields on entrainment levels, don't need checkpointing
    vector_space => function_space_collection%get_fs(twod_mesh, 0, W3,         &
         get_ndata_val('entrainment_levels'))
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'ent_we_lim', vector_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'ent_t_frac', vector_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'ent_zrzi', vector_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'ent_we_lim_dsc', vector_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'ent_t_frac_dsc', vector_space, twod=.true. )
    call add_physics_field( turbulence_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'ent_zrzi_dsc', vector_space, twod=.true. )

    !========================================================================
    ! Fields owned by the convection scheme
    !========================================================================
    call convection_fields%initialise(name='convection_fields', table_len=100)

    ! 2D fields, might need checkpointing
    if (convection == convection_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'conv_rain',  twod_space,                                                &
      checkpoint_flag=(checkpoint_flag .and. add_cgw), twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'conv_snow',  twod_space,                                                &
      checkpoint_flag=(checkpoint_flag .and. add_cgw), twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
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
      adv_fields_last_outer, &
      'cca', wtheta_space, checkpoint_flag=checkpoint_flag)
    call add_physics_field(convection_fields, depository, prognostic_fields, &
      adv_fields_last_outer, &
      'ccw', wtheta_space, checkpoint_flag=checkpoint_flag)

    ! 2D fields, don't need checkpointing
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'cca_2d',  twod_space, twod=.true. )
    call add_integer_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'shallow_flag',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'uw0_flux',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'vw0_flux',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'lcl_height',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'parcel_top',  twod_space, twod=.true. )
    call add_integer_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'level_parcel_top',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'wstar',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'thv_flux',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'parcel_buoyancy',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'qsat_at_lcl',  twod_space, twod=.true. )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'cape_diluted', twod_space, twod=.true. )

    ! 3D fields, don't need checkpointing
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'dt_conv', wtheta_space )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'dmv_conv', wtheta_space )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'dmcl_conv', wtheta_space )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'dmci_conv', wtheta_space )
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'dcfl_conv', wtheta_space )
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'dcff_conv', wtheta_space )
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'dbcf_conv', wtheta_space )
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'massflux_up', wtheta_space )
    call add_physics_field(convection_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'massflux_down', wtheta_space )
    call add_physics_field( convection_fields, depository, prognostic_fields, &
      adv_fields_last_outer, &
      'conv_rain_3d', wtheta_space )
    call add_physics_field( convection_fields, depository, prognostic_fields, &
      adv_fields_last_outer, &
      'conv_snow_3d', wtheta_space )

    ! 3D fields on W3 (rho) levels
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'du_conv', w3_space )
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'dv_conv', w3_space )

    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_all_outer,  'conv_prog_dtheta', wtheta_space,                 &
      checkpoint_flag=l_conv_prog_dtheta,                                      &
      advection_flag=(l_conv_prog_dtheta .and. adv_conv_prog_dtheta))
    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_all_outer,  'conv_prog_dmv', wtheta_space,                    &
      checkpoint_flag=l_conv_prog_dq,                                          &
      advection_flag=(l_conv_prog_dq .and. adv_conv_prog_dq))

    call add_physics_field( convection_fields, depository, prognostic_fields,  &
      adv_fields_all_outer, 'conv_prog_precip', wtheta_space,                  &
      checkpoint_flag=l_conv_prog_precip,                                      &
      advection_flag=l_conv_prog_precip)

    !========================================================================
    ! Fields owned by the cloud scheme
    !========================================================================
    call cloud_fields%initialise(name='cloud_fields', table_len=100)

    ! 3D fields, might need checkpointing
    if (cloud == cloud_um) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      adv_fields_last_outer, &
      'area_fraction',   wtheta_space,                                  &
      checkpoint_flag=(checkpoint_flag .and. radiation==radiation_socrates))

    ! 3D fields, might need advecting
    if ( scheme == scheme_pc2 ) then
      advection_flag=.true.
    else
      advection_flag=.false.
    endif
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      adv_fields_all_outer, &
      'liquid_fraction', wtheta_space, checkpoint_flag=checkpoint_flag, &
      advection_flag=advection_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      adv_fields_all_outer, &
      'frozen_fraction', wtheta_space, checkpoint_flag=checkpoint_flag, &
      advection_flag=advection_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      adv_fields_all_outer, &
      'bulk_fraction',   wtheta_space, checkpoint_flag=checkpoint_flag, &
      advection_flag=advection_flag)

    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      adv_fields_last_outer, &
      'rh_crit',     wtheta_space )
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      adv_fields_last_outer, &
      'departure_exner_wth', wtheta_space, advection_flag=advection_flag)
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
       adv_fields_last_outer, &
      'sigma_mc',   wtheta_space )

    ! Fields for bimodal cloud scheme
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      adv_fields_last_outer, &
      'tau_dec_bm',  wtheta_space )
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      adv_fields_last_outer, &
      'tau_hom_bm',  wtheta_space )
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      adv_fields_last_outer, &
      'tau_mph_bm',  wtheta_space )

    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      adv_fields_last_outer, &
      'sskew_bm',     wtheta_space )
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      adv_fields_last_outer, &
      'svar_bm',     wtheta_space )
    call add_physics_field(cloud_fields, depository, prognostic_fields, &
      adv_fields_last_outer, &
      'svar_tb',     wtheta_space )

    !========================================================================
    ! Fields owned by the surface exchange scheme
    !========================================================================
    call surface_fields%initialise(name='surface_fields', table_len=100)

    ! 2D fields, might need checkpointing
    if (surface == surface_jules) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if

    ! Coupling fields might need checkpointing
    if (surface == surface_jules .and. (l_esm_couple .OR. l_esm_couple_test)) then
      checkpoint_couple = .true.
    else
      checkpoint_couple = .false.
    end if

    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'z0msea',  twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'surface_conductance', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'chloro_sea', twod_space,                                                &
      checkpoint_flag=(checkpoint_flag .and. sea_alb_var_chl), twod=.true. )

    ! Fields on surface tiles, might need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'tile_fraction', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'tile_temperature', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'screen_temperature', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'time_since_transition', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'canopy_water', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'tile_lw_grey_albedo', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    vector_space=>function_space_collection%get_fs(twod_mesh, 0, W3,           &
         get_ndata_val('land_tile_rad_band'))
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'albedo_obs_scaling', vector_space,                                      &
      checkpoint_flag=(checkpoint_flag .and. albedo_obs), twod=.true. )

    ! Fields on plant functional types, might need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'leaf_area_index', pft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'canopy_height', pft_space, checkpoint_flag=checkpoint_flag, twod=.true. )


    ! Sea-ice category fields, might need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'sea_ice_thickness', sice_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'sea_ice_temperature', sice_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'sea_ice_conductivity', sice_space, checkpoint_flag=.false., twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'melt_pond_fraction', sice_space, checkpoint_flag=.false., twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'melt_pond_depth', sice_space, checkpoint_flag=.false., twod=.true. )

    ! Fields on surface tiles, don't need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'tile_heat_flux', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'alpha1_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'ashtf_prime_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'dtstar_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'fraca_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'z0h_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'z0m_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'rhokh_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'chr1p5m_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'resfs_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'canhc_tile', surft_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'gc_tile', surft_space, twod=.true. )

   ! Fields on surface tiles used by coupler, need checkpointing in coupled models
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'tile_moisture_flux', surft_space, checkpoint_flag=checkpoint_couple, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'snowice_melt', surft_space, checkpoint_flag=checkpoint_couple, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'surf_ht_flux', surft_space, checkpoint_flag=checkpoint_couple, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'snowice_sublimation', surft_space, checkpoint_flag=checkpoint_couple, twod=.true. )

    ! 2D fields, don't need checkpointing
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'z0m_eff', twod_space, twod=.true. )
    call add_physics_field(surface_fields, depository, prognostic_fields,      &
      adv_fields_last_outer, &
      'net_prim_prod', twod_space, twod=.true.)
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'taux_ssi', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'tauy_ssi', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'z0m', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'ustar', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'wspd10m', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'urbwrr', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'urbhwr', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'urbhgt', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'urbdisp', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'urbztm', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'urbalbwl', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'urbalbrd', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'urbemisw', twod_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'urbemisr', twod_space, twod=.true. )

    ! Space for variables required for regridding to cell faces
    vector_space => function_space_collection%get_fs(twod_mesh, 0, W2,         &
         get_ndata_val('surface_regrid_vars'))
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'surf_interp_w2',  vector_space, twod=.true. )

    ! 2D fields at W2 points
    vector_space => function_space_collection%get_fs(twod_mesh, 0, W2)
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'tau_land_w2',  vector_space, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'tau_ssi_w2',  vector_space, twod=.true. )

    ! Field on soil levels and land tiles
    vector_space => function_space_collection%get_fs(twod_mesh, 0, W3,         &
         get_ndata_val('soil_levels_and_tiles'))
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'tile_water_extract', vector_space, twod=.true. )

    checkpoint_flag = .false.
    ! Sea surface currents, should not need checkpointing because they
    ! come from the coupler
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'sea_u_current', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( surface_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'sea_v_current', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )



    !========================================================================
    ! Fields owned by the soil hydrology scheme
    !========================================================================
    call soil_fields%initialise(name='soil_fields', table_len=100)

    ! 2D fields, might need checkpointing
    if (surface == surface_jules) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'soil_albedo', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'soil_roughness', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'soil_thermal_cond', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      adv_fields_last_outer, &
      'mean_topog_index', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      adv_fields_last_outer, &
      'a_sat_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      adv_fields_last_outer, &
      'c_sat_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      adv_fields_last_outer, &
      'a_wet_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      adv_fields_last_outer, &
      'c_wet_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      adv_fields_last_outer, &
      'soil_sat_frac', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      adv_fields_last_outer, &
      'water_table', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      adv_fields_last_outer, &
      'wetness_under_soil', twod_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'soil_moist_wilt', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'soil_moist_crit', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'soil_moist_sat', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'soil_cond_sat', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'soil_thermal_cap', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'soil_suction_sat', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'clapp_horn_b', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! Fields on soil levels
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'soil_temperature', soil_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'soil_moisture', soil_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'unfrozen_soil_moisture', soil_space, checkpoint_flag=checkpoint_flag,  &
      twod=.true. )
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'frozen_soil_moisture', soil_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! 2D fields, don't need checkpointing
    call add_physics_field( soil_fields, depository, prognostic_fields,       &
      adv_fields_last_outer, &
      'soil_moist_avail', twod_space, twod=.true. )
    call add_physics_field(soil_fields, depository, prognostic_fields,        &
      adv_fields_last_outer, &
      'thermal_cond_wet_soil', twod_space, twod=.true.)

    !========================================================================
    ! Fields owned by the snow scheme
    !========================================================================
    call snow_fields%initialise(name='snow_fields', table_len=100)

    ! Fields on surface tiles, might need checkpointing
    if (surface == surface_jules) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'tile_snow_mass', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'tile_snow_rgrain', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_integer_field( snow_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'n_snow_layers', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'snow_depth', surft_space, checkpoint_flag=checkpoint_flag, twod=.true. )
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'snowpack_density', surft_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'snow_under_canopy', surft_space, checkpoint_flag=checkpoint_flag, twod=.true.)

    ! Fields on snow layers
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'snow_layer_thickness', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'snow_layer_ice_mass', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'snow_layer_liq_mass', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'snow_layer_temp', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)
    call add_physics_field(snow_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'snow_layer_rgrain', snow_space, checkpoint_flag=checkpoint_flag, twod=.true.)

    ! 2D fields
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'snow_soot', twod_space, checkpoint_flag=checkpoint_flag, twod=.true. )

    ! Fields which don't need checkpointing
    call add_physics_field( snow_fields, depository, prognostic_fields,  &
      adv_fields_last_outer, &
      'snow_unload_rate', pft_space, twod=.true. )

    !========================================================================
    ! Fields owned by the chemistry scheme
    !========================================================================
    call chemistry_fields%initialise(name='chemistry_fields', table_len=100)

    ! 3D fields, might need checkpointing and advecting
    if ( aerosol == aerosol_um .and. glomap_mode == glomap_mode_ukca ) then
      checkpoint_flag = .true.
      advection_flag = .true.
    else
      checkpoint_flag = .false.
      advection_flag = .false.
    end if
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'h2o2', wtheta_space, checkpoint_flag=checkpoint_flag,                   &
      advection_flag=advection_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'dms', wtheta_space, checkpoint_flag=checkpoint_flag,                    &
      advection_flag=advection_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'so2', wtheta_space, checkpoint_flag=checkpoint_flag,                    &
      advection_flag=advection_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'h2so4', wtheta_space, checkpoint_flag=checkpoint_flag,                  &
      advection_flag=advection_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'dmso', wtheta_space, checkpoint_flag=checkpoint_flag,                   &
      advection_flag=advection_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'monoterpene', wtheta_space, checkpoint_flag=checkpoint_flag,            &
      advection_flag=advection_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'secondary_organic', wtheta_space, checkpoint_flag=checkpoint_flag,      &
      advection_flag=advection_flag )

    ! 3D fields, might need checkpointing
    ! Not advected in Offline Oxidants chemistry scheme
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'o3', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'no3', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'oh', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'ho2', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( chemistry_fields, depository, prognostic_fields,   &
      adv_fields_last_outer, &
      'h2o2_limit', wtheta_space, checkpoint_flag=checkpoint_flag )

    !========================================================================
    ! Fields owned by the aerosol scheme
    !========================================================================
    call aerosol_fields%initialise(name='aerosol_fields', table_len=100)

    ! 2D fields, might need checkpointing
    if ( aerosol == aerosol_um .and. glomap_mode == glomap_mode_ukca ) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'emiss_bc_biofuel', twod_space, twod=.true.,                             &
      checkpoint_flag=checkpoint_flag)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'emiss_bc_fossil', twod_space, twod=.true.,                              &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'emiss_dms_land', twod_space, twod=.true.,                               &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'dms_conc_ocean', twod_space, twod=.true.,                               &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'emiss_monoterp', twod_space, twod=.true.,                               &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'emiss_om_biofuel', twod_space, twod=.true.,                             &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'emiss_om_fossil', twod_space, twod=.true.,                              &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'emiss_so2_low', twod_space, twod=.true.,                                &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'emiss_so2_high', twod_space, twod=.true.,                               &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'soil_clay', twod_space, twod=.true.,                                    &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'soil_sand', twod_space, twod=.true.,                                    &
      checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'surf_wetness', twod_space, twod=.true.,                                 &
      checkpoint_flag=checkpoint_flag )

    ! 3D fields, might need checkpointing
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'emiss_bc_biomass', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'emiss_om_biomass', wtheta_space, checkpoint_flag=checkpoint_flag )
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'emiss_so2_nat', wtheta_space, checkpoint_flag=checkpoint_flag )

    ! 3D fields, might need checkpointing and/or advecting
    if ( aerosol == aerosol_um .and.                                           &
           ( glomap_mode == glomap_mode_climatology .or.                       &
             glomap_mode == glomap_mode_ukca ) ) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    if ( aerosol == aerosol_um .and. glomap_mode == glomap_mode_ukca ) then
      advection_flag = .true.
    else
      advection_flag = .false.
    end if
    ! Nucleation soluble mode number mixing ratio
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'n_nuc_sol', wtheta_space, checkpoint_flag=checkpoint_flag,              &
      advection_flag=advection_flag )
    ! Nucleation soluble H2SO4 aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'nuc_sol_su', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Nucleation soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'nuc_sol_om', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Aitken soluble mode number mixing ratio
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'n_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag,              &
      advection_flag=advection_flag )
    ! Aitken soluble H2SO4 aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'ait_sol_su', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Aitken soluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'ait_sol_bc', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Aitken soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'ait_sol_om', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Accumulation soluble mode number mixing ratio
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'n_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag,              &
      advection_flag=advection_flag )
    ! Accumulation soluble H2SO4 aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'acc_sol_su', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Accumulation soluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'acc_sol_bc', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Accumulation soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'acc_sol_om', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Accumulation soluble sea salt aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'acc_sol_ss', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Accumulation soluble dust aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'acc_sol_du', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Coarse soluble mode number mixing ratio
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'n_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag,              &
      advection_flag=advection_flag )
    ! Coarse soluble H2SO4 aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'cor_sol_su', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Coarse soluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'cor_sol_bc', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Coarse soluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'cor_sol_om', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Coarse soluble sea salt aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'cor_sol_ss', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Coarse soluble dust aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'cor_sol_du', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Aitken insoluble mode number mixing ratio
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'n_ait_ins', wtheta_space, checkpoint_flag=checkpoint_flag,              &
      advection_flag=advection_flag )
    ! Aitken insoluble black carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'ait_ins_bc', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Aitken insoluble organic carbon aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'ait_ins_om', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Accumulation insoluble mode number mixing ratio
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'n_acc_ins', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Accumulation insoluble dust aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'acc_ins_du', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Coarse insoluble mode number mixing ratio
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'n_cor_ins', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )
    ! Coarse insoluble dust aerosol mmr
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'cor_ins_du', wtheta_space, checkpoint_flag=checkpoint_flag,             &
      advection_flag=advection_flag )

    ! 3D fields, might need checkpointing
    if (aerosol == aerosol_um .and. glomap_mode == glomap_mode_ukca) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if

    ! Cloud droplet number concentration
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'cloud_drop_no_conc', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Dry diameter Aitken mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'drydp_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Dry diameter Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'drydp_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Dry diameter Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'drydp_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Dry diameter Aitken mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'drydp_ait_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Dry diameter Accumulation mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'drydp_acc_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Dry diameter Coarse mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'drydp_cor_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Wet diameter Aitken mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'wetdp_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Wet diameter Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'wetdp_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Wet diameter Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'wetdp_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Particle density Aitken mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'rhopar_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Particle density Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'rhopar_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Particle density Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'rhopar_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Particle density Aitken mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'rhopar_ait_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Particle density Accumulation mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'rhopar_acc_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Particle density Coarse mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'rhopar_cor_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume of water Aitken mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_wat_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume of water Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_wat_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume of water Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_wat_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Sulphate Aitken mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_su_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Black Carbon Aitken mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_bc_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Organic Matter Aitken mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_om_ait_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Sulphate Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_su_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Black Carbon Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_bc_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Organic Matter Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_om_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Sea Salt Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_ss_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Dust Accumulation mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_du_acc_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Sulphate Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_su_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Black Carbon Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_bc_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Organic Matter Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_om_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Sea Salt Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_ss_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Dust Coarse mode (Soluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_du_cor_sol', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Black Carbon Aitken mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_bc_ait_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Organic Matter Aitken mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_om_ait_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Dust Accumulation mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_du_acc_ins', wtheta_space, checkpoint_flag=checkpoint_flag )
    ! Partial volume component Dust Coarse mode (Insoluble)
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'pvol_du_cor_ins', wtheta_space, checkpoint_flag=checkpoint_flag )

    ! Fields on dust space, might need checkpointing
    vector_space => function_space_collection%get_fs(twod_mesh, 0, W3,         &
         get_ndata_val('dust_divisions'))
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'dust_mrel', vector_space, twod=.true. , checkpoint_flag=checkpoint_flag )

    ! Fields on dust space, don't need checkpointing
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'dust_flux', vector_space, twod=.true. )

    ! 3D fields, don't need checkpointing
    ! Sulphuric Acid aerosol MMR
    call add_physics_field( aerosol_fields, depository, prognostic_fields,     &
      adv_fields_last_outer, &
      'sulphuric', wtheta_space )
#endif

  end subroutine create_physics_prognostics

  !>@brief Add field to field collection and set its write,
  !>       checkpoint-restart and advection behaviour
  !> @param[in,out] field_collection  Field collection that 'name' will be added to
  !> @param[in,out] depository        Collection of all fields
  !> @param[in,out] prognostic_fields Collection of checkpointed fields
  !> @param[in,out] advected_fields   Collection of fields to be advected
  !> @param[in]     name              Name of field to be added to collection
  !> @param[in]     vector_space      Function space of field to set behaviour for
  !> @param[in]     checkpoint_flag   Optional flag to allow checkpoint-
  !>                                   restart behaviour of field to be set
  !> @param[in]     twod              Optional flag to determine if this is a
  !>                                   2D field for diagnostic output
  !> @param[in]     advection_flag    Optional flag whether this field is to be advected
  subroutine add_physics_field(field_collection, &
                               depository, prognostic_fields, advected_fields, &
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
    type(field_collection_type), intent(inout)     :: advected_fields
    type(function_space_type), pointer, intent(in) :: vector_space
    logical(l_def), optional, intent(in)           :: checkpoint_flag
    logical(l_def), optional, intent(in)           :: twod
    logical(l_def), optional, intent(in)           :: advection_flag
    !Local variables
    type(field_type)                               :: new_field
    type(field_type), pointer                      :: field_ptr => null()
    class(pure_abstract_field_type), pointer       :: tmp_ptr => null()
    logical(l_def)                                 :: twod_field, checkpointed
    logical(l_def)                                 :: advected

    ! pointers for xios write interface
    procedure(write_interface), pointer :: write_behaviour => null()
    procedure(read_interface),  pointer :: read_behaviour => null()
    procedure(checkpoint_write_interface), pointer :: checkpoint_write_behaviour => null()
    procedure(checkpoint_read_interface), pointer  :: checkpoint_read_behaviour => null()

    ! Create the new field
    call new_field%initialise( vector_space, name=trim(name) )

    ! Set advection flag
    if (present(advection_flag)) then
      advected = advection_flag
    else
      advected = .false.
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
    call depository%get_field(name, field_ptr)
    tmp_ptr => field_ptr
    ! Put a pointer to the field in the required collection
    call field_collection%add_reference_to_field( tmp_ptr )
    ! If checkpointing the field, put a pointer to it in the prognostics collection
    if ( checkpointed ) then
      call prognostic_fields%add_reference_to_field( tmp_ptr )
    endif
    ! If advecting the field, put a pointer to it in the advected collection
    if ( advected ) then
      call advected_fields%add_reference_to_field( tmp_ptr )
    endif

  end subroutine add_physics_field

  !>@brief Add integer field to field collection and set its write,
  !>       checkpoint-restart and advection behaviour
  !> @param[in,out] field_collection  Field collection that 'name' will be added to
  !> @param[in,out] depository        Collection of all fields
  !> @param[in,out] prognostic_fields Collection of checkpointed fields
  !> @param[in,out] advected_fields   Collection of fields to be advected
  !> @param[in]     name              Name of field to be added to collection
  !> @param[in]     vector_space      Function space of field to set behaviour for
  !> @param[in]     checkpoint_flag   Optional flag to allow checkpoint-
  !>                                   restart behaviour of field to be set
  !> @param[in]     twod              Optional flag to determine if this is a
  !>                                   2D field for diagnostic output
  !> @param[in]     advection_flag    Optional flag whether this field is to be advected
  subroutine add_integer_field(field_collection, &
                               depository, prognostic_fields, advected_fields, &
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
    type(field_collection_type), intent(inout)     :: advected_fields
    type(function_space_type), pointer, intent(in) :: vector_space
    logical(l_def), optional, intent(in)           :: checkpoint_flag
    logical(l_def), optional, intent(in)           :: twod
    logical(l_def), optional, intent(in)           :: advection_flag
    !Local variables
    type(integer_field_type)                       :: new_field
    type(integer_field_type), pointer              :: field_ptr => null()
    class(pure_abstract_field_type), pointer       :: tmp_ptr => null()
    logical(l_def)                                 :: twod_field, checkpointed
    logical(l_def)                                 :: advected

    ! pointers for xios write interface
    procedure(write_interface), pointer :: write_behaviour => null()
    procedure(read_interface),  pointer :: read_behaviour => null()
    procedure(checkpoint_write_interface), pointer :: checkpoint_write_behaviour => null()
    procedure(checkpoint_read_interface), pointer  :: checkpoint_read_behaviour => null()

    ! Create the new field
    call new_field%initialise( vector_space, name=trim(name) )

    ! Set advection flag
    if (present(advection_flag)) then
      advected = advection_flag
    else
      advected = .false.
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
    call depository%get_field(name, field_ptr)
    ! Put a pointer to the field in the required collection
    tmp_ptr => field_ptr
    call field_collection%add_reference_to_field( tmp_ptr )
    ! If checkpointing the field, put a pointer to it in the prognostics collection
    if ( checkpointed ) then
      call prognostic_fields%add_reference_to_field( tmp_ptr )
    endif
    ! If advecting the field, put a pointer to it in the advected collection
    if ( advected ) then
      call advected_fields%add_reference_to_field( tmp_ptr )
    endif

  end subroutine add_integer_field

end module create_physics_prognostics_mod
