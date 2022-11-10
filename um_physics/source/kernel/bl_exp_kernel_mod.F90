!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the explicit UM boundary layer scheme.
module bl_exp_kernel_mod

  use argument_mod,           only : arg_type,                   &
                                     GH_FIELD, GH_SCALAR,        &
                                     GH_REAL, GH_INTEGER,        &
                                     GH_READ, GH_WRITE,          &
                                     GH_READWRITE, CELL_COLUMN,  &
                                     ANY_DISCONTINUOUS_SPACE_1,  &
                                     ANY_DISCONTINUOUS_SPACE_2,  &
                                     ANY_DISCONTINUOUS_SPACE_3,  &
                                     ANY_DISCONTINUOUS_SPACE_4,  &
                                     ANY_DISCONTINUOUS_SPACE_5,  &
                                     ANY_DISCONTINUOUS_SPACE_6,  &
                                     ANY_DISCONTINUOUS_SPACE_7,  &
                                     ANY_DISCONTINUOUS_SPACE_8,  &
                                     ANY_DISCONTINUOUS_SPACE_9,  &
                                     ANY_DISCONTINUOUS_SPACE_10, &
                                     ANY_DISCONTINUOUS_SPACE_11, &
                                     ANY_DISCONTINUOUS_SPACE_12, &
                                     STENCIL, REGION
  use constants_mod,          only : i_def, i_um, r_def, r_um, rmdi
  use empty_data_mod,         only : empty_real_data
  use fs_continuity_mod,      only : W3, Wtheta
  use kernel_mod,             only : kernel_type
  use blayer_config_mod,      only : fixed_flux_e, fixed_flux_h, flux_bc_opt, &
                                     flux_bc_opt_specified_scalars, bl_mix_w
  use cloud_config_mod,       only : rh_crit_opt, rh_crit_opt_tke, scheme,    &
                                     scheme_bimodal, scheme_pc2,              &
                                     pc2ini, pc2ini_bimodal
  use radiation_config_mod,   only : topography, topography_horizon
  use convection_config_mod,  only : use_jules_flux
  use microphysics_config_mod, only: turb_gen_mixph
  use mixing_config_mod,      only : smagorinsky
  use jules_surface_config_mod, only : formdrag, formdrag_dist_drag
  use surface_config_mod,     only : albedo_obs, sea_surf_alg, &
                                     sea_surf_alg_fixed_roughness, &
                                     buddy_sea, buddy_sea_on, &
                                     emis_method_soil, emis_method_soil_fixed
  use timestepping_config_mod, only: outer_iterations
  use water_constants_mod,     only: tfs
  use jules_surface_types_mod, only: soil

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: bl_exp_kernel_type
    private
    type(arg_type) :: meta_args(152) = (/                                      &
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! theta_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3),                       &! rho_in_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! wetrho_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3),                       &! exner_in_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! exner_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3, STENCIL(REGION)),      &! u_in_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3, STENCIL(REGION)),      &! v_in_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! w_in_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! velocity_w2v
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! m_v_n
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! m_cl_n
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! m_ci_n
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      W3),                       &! height_w3
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! height_wth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! shear
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! delta
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! max_diff_smag
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! zh_2d
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! z0msea_2d
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! z0m_2d
         arg_type(GH_FIELD, GH_INTEGER,  GH_WRITE,  ANY_DISCONTINUOUS_SPACE_1),&! ntml_2d
         arg_type(GH_FIELD, GH_INTEGER,  GH_WRITE,  ANY_DISCONTINUOUS_SPACE_1),&! cumulus_2d
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2, STENCIL(REGION)),&! tile_fraction
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_3),&! leaf_area_index
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_3),&! canopy_height
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sd_orog_2d
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! peak_to_trough_orog
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! silhouette_area_orog
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_albedo
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_roughness
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_wilt
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_crit
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_sat
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_thermal_cond
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_suction_sat
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! clapp_horn_b
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! soil_respiration
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! thermal_cond_wet_soil
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_temperature
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_conductivity
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! tile_temperature
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_lw_grey_albedo
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_snow_mass
         arg_type(GH_FIELD, GH_INTEGER,  GH_READ,   ANY_DISCONTINUOUS_SPACE_2),&! n_snow_layers
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! snow_depth
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! snow_layer_thickness
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! snow_layer_ice_mass
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! snow_layer_liq_mass
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! snow_layer_temp
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! surface_conductance
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! canopy_water
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! soil_temperature
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! soil_moisture
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! unfrozen_soil_moisture
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! frozen_soil_moisture
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! tile_heat_flux
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! tile_moisture_flux
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! gross_prim_prod
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! net_prim_prod
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! cos_zen_angle
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! skyview
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! sw_up_tile
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sw_down_surf
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! lw_down_surf
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sw_down_surf_blue
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! dd_mf_cb
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! dtl_mphys
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! dmt_mphys
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! sw_heating_rate
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! lw_heating_rate
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! ozone
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! cf_bulk
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      WTHETA),                   &! cf_liquid
         arg_type(GH_FIELD, GH_REAL,  GH_READWRITE, WTHETA),                   &! rh_crit
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! dsldzm
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! wvar
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! visc_m_blend
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! visc_h_blend
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! dw_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! rhokm_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_7),&! surf_interp
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! rhokh_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! tke_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! ngstress_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! bq_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! bt_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! moist_flux_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! heat_flux_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! dtrdz_tq_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! rdz_tq_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! fd_taux
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     W3),                       &! fd_tauy
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! lmix_bl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     WTHETA),                   &! gradrinr
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! alpha1_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! ashtf_prime_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! dtstar_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! fraca_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! z0h_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! z0m_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! rhokh_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! chr1p5m_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! resfs_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! gc_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! canhc_tile
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_8),&! tile_water_extract
         arg_type(GH_FIELD, GH_INTEGER,  GH_WRITE,  ANY_DISCONTINUOUS_SPACE_1),&! blend_height_tq
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! z0m_eff
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! ustar
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_avail
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! zh_nonloc
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! zhsc_2d
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! z_lcl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! inv_depth
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! qcl_at_inv_top
         arg_type(GH_FIELD, GH_INTEGER,  GH_WRITE,  ANY_DISCONTINUOUS_SPACE_1),&! shallow_flag
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! uw0_flux
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! vw0_flux
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! lcl_height
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! parcel_top
         arg_type(GH_FIELD, GH_INTEGER,  GH_WRITE,  ANY_DISCONTINUOUS_SPACE_1),&! level_parcel_top
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! wstar_2d
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! thv_flux
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! parcel_buoyancy
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! qsat_at_lcl
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! bl_weight_1dbl
         arg_type(GH_FIELD, GH_INTEGER,  GH_WRITE,  ANY_DISCONTINUOUS_SPACE_9),&! bl_type_ind
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_3),&! snow_unload_rate
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_10),&! albedo_obs_scaling
         arg_type(GH_FIELD, GH_INTEGER,  GH_WRITE,  ANY_DISCONTINUOUS_SPACE_1), &! level_ent
         arg_type(GH_FIELD, GH_INTEGER,  GH_WRITE,  ANY_DISCONTINUOUS_SPACE_1), &! level_ent_dsc
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_11),&! ent_we_lim
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_11),&! ent_t_frac
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_11),&! ent_zrzi
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_11),&! ent_we_lim_dsc
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_11),&! ent_t_frac_dsc
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_11),&! ent_zrzi_dsc
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &! soil_clay_2d
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &! soil_sand_2d
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_12),&! dust_div_mrel
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_12),&! dust_div_flux
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ                             ), &! day_of_year
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &! urbwrr
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &! urbhwr
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &! urbhgt
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &! urbztm
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &! urbdisp
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &! urbemisw
         arg_type(GH_FIELD, GH_REAL,  GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &! urbemisr
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), &! diag__zht
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), &! diag__z0h_eff
         arg_type(GH_FIELD, GH_REAL,  GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1)  &! diag__oblen
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: bl_exp_code
  end type

  public :: bl_exp_code

contains

  !> @brief Interface to the UM BL scheme
  !> @details The UM Boundary Layer scheme does:
  !>             vertical mixing of heat, momentum and moisture,
  !>             as documented in UMDP24
  !>          NB This version uses winds in w3 space (i.e. A-grid)
  !> @param[in]     nlayers                Number of layers
  !> @param[in]     theta_in_wth           Potential temperature field
  !> @param[in]     rho_in_w3              Density field in density space
  !> @param[in]     wetrho_in_wth          Wet density field in wth space
  !> @param[in]     exner_in_w3            Exner pressure field in density space
  !> @param[in]     exner_in_wth           Exner pressure field in wth space
  !> @param[in]     u_in_w3                'Zonal' wind in density space
  !> @param[in]     v_in_w3                'Meridional' wind in density space
  !> @param[in]     w_in_wth               'Vertical' wind in theta space
  !> @param[in]     velocity_w2v           Velocity normal to cell top
  !> @param[in]     m_v_n                  Vapour mixing ratio at time level n
  !> @param[in]     m_cl_n                 Cloud liquid mixing ratio at time level n
  !> @param[in]     m_ci_n                 Cloud ice mixing ratio at time level n
  !> @param[in]     height_w3              Height of density space above surface
  !> @param[in]     height_wth             Height of theta space above surface
  !> @param[in]     shear                  3D wind shear on wtheta points
  !> @param[in]     delta                  Edge length on wtheta points
  !> @param[in]     max_diff_smag          Maximum diffusion coefficient allowed
  !> @param[in,out] zh_2d                  Boundary layer depth
  !> @param[in,out] z0msea_2d              Roughness length (sea)
  !> @param[in,out] z0m_2d                 Cell roughness length
  !> @param[in,out] ntml_2d                Number of turbulently mixed levels
  !> @param[in,out] cumulus_2d             Cumulus flag (true/false)
  !> @param[in]     tile_fraction          Surface tile fractions
  !> @param[in]     leaf_area_index        Leaf Area Index
  !> @param[in]     canopy_height          Canopy height
  !> @param[in]     sd_orog_2d             Standard deviation of orography
  !> @param[in]     peak_to_trough_orog    Half of peak-to-trough height over root(2) of orography
  !> @param[in]     silhouette_area_orog   Silhouette area of orography
  !> @param[in]     soil_albedo            Snow-free soil albedo
  !> @param[in]     soil_roughness         Bare soil surface roughness length
  !> @param[in]     soil_moist_wilt        Volumetric soil moisture at wilting point
  !> @param[in]     soil_moist_crit        Volumetric soil moisture at critical point
  !> @param[in]     soil_moist_sat         Volumetric soil moisture at saturation
  !> @param[in]     soil_thermal_cond      Soil thermal conductivity
  !> @param[in]     soil_suction_sat       Saturated soil water suction
  !> @param[in]     clapp_horn_b           Clapp and Hornberger b coefficient
  !> @param[in,out] soil_respiration       Soil respiration  (kg m-2 s-1)
  !> @param[in,out] thermal_cond_wet_soil  Thermal conductivity of wet soil (W m-1 K-1)
  !> @param[in]     sea_ice_temperature    Bulk temperature of sea-ice (K)
  !> @param[in]     sea_ice_conductivity   Sea ice thermal conductivity (W m-2 K-1)
  !> @param[in,out] tile_temperature       Surface tile temperatures
  !> @param[in]     tile_lw_grey_albedo    Surface tile longwave grey albedo
  !> @param[in]     tile_snow_mass         Snow mass on tiles (kg/m2)
  !> @param[in]     n_snow_layers          Number of snow layers on tiles
  !> @param[in]     snow_depth             Snow depth on tiles
  !> @param[in]     snow_layer_thickness   Thickness of snow layers (m)
  !> @param[in]     snow_layer_ice_mass    Mass of ice in snow layers (kg m-2)
  !> @param[in]     snow_layer_liq_mass    Mass of liquid in snow layers (kg m-2)
  !> @param[in]     snow_layer_temp        Temperature of snow layer (K)
  !> @param[in,out] surface_conductance    Surface conductance
  !> @param[in]     canopy_water           Canopy water on each tile
  !> @param[in]     soil_temperature       Soil temperature
  !> @param[in]     soil_moisture          Soil moisture content (kg m-2)
  !> @param[in]     unfrozen_soil_moisture Unfrozen soil moisture proportion
  !> @param[in]     frozen_soil_moisture   Frozen soil moisture proportion
  !> @param[in,out] tile_heat_flux         Surface heat flux
  !> @param[in,out] tile_moisture_flux     Surface moisture flux
  !> @param[in,out] gross_prim_prod        Gross Primary Productivity
  !> @param[in,out] net_prim_prod          Net Primary Productivity
  !> @param[in]     cos_zen_angle          Cosine of solar zenith angle
  !> @param[in]     skyview                Skyview / area enhancement factor
  !> @param[in]     sw_up_tile             Upwelling SW radiation on surface tiles
  !> @param[in]     sw_down_surf           Downwelling SW radiation at surface
  !> @param[in]     lw_down_surf           Downwelling LW radiation at surface
  !> @param[in]     sw_down_surf_blue      Photosynthetically active SW down
  !> @param[in]     dd_mf_cb               Downdraft massflux at cloud base (Pa/s)
  !> @param[in]     dtl_mphys              Microphysics liq temperature increment
  !> @param[in]     dmt_mphys              Microphysics total water increment
  !> @param[in]     sw_heating_rate        Shortwave radiation heating rate
  !> @param[in]     lw_heating_rate        Longwave radiation heating rate
  !> @param[in]     ozone                  Ozone field
  !> @param[in]     cf_bulk                Bulk cloud fraction
  !> @param[in]     cf_liquid              Liquid cloud fraction
  !> @param[in,out] rh_crit                Critical rel humidity
  !> @param[in,out] dsldzm                 Liquid potential temperature gradient in wth
  !> @param[in,out] wvar                   Vertical velocity variance in wth
  !> @param[in,out] visc_m_blend           Blended BL-Smag diffusion coefficient for momentum
  !> @param[in,out] visc_h_blend           Blended BL-Smag diffusion coefficient for scalars
  !> @param[in,out] dw_bl                  Vertical wind increment from BL scheme
  !> @param[in,out] rhokm_bl               Momentum eddy diffusivity on BL levels
  !> @param[in,out] surf_interp            Surface variables for regridding
  !> @param[in,out] rhokh_bl               Heat eddy diffusivity on BL levels
  !> @param[in,out] tke_bl                 Turbulent kinetic energy (m2 s-2)
  !> @param[in,out] ngstress_bl            Non-gradient stress function on BL levels
  !> @param[in,out] bq_bl                  Buoyancy parameter for moisture
  !> @param[in,out] bt_bl                  Buoyancy parameter for heat
  !> @param[in,out] moist_flux_bl          Vertical moisture flux on BL levels
  !> @param[in,out] heat_flux_bl           Vertical heat flux on BL levels
  !> @param[in,out] dtrdz_tq_bl            dt/(rho*r*r*dz) in wth
  !> @param[in,out] rdz_tq_bl              1/dz in w3
  !> @param[in,out] fd_taux                'Zonal' momentum stress from form drag
  !> @param[in,out] fd_tauy                'Meridional' momentum stress from form drag
  !> @param[in,out] lmix_bl                Turbulence mixing length in wth
  !> @param[in,out] gradrinr               Gradient Richardson number in wth
  !> @param[in,out] alpha1_tile            dqsat/dT in surface layer on tiles
  !> @param[in,out] ashtf_prime_tile       Heat flux coefficient on tiles
  !> @param[in,out] dtstar_tile            Change in surface temperature on tiles
  !> @param[in,out] fraca_tile             Fraction of moisture flux with only aerodynamic resistance
  !> @param[in,out] z0h_tile               Heat roughness length on tiles
  !> @param[in,out] z0m_tile               Momentum roughness length on tiles
  !> @param[in,out] rhokh_tile             Surface heat diffusivity on tiles
  !> @param[in,out] chr1p5m_tile           1.5m transfer coefficients on tiles
  !> @param[in,out] resfs_tile             Combined aerodynamic resistance
  !> @param[in,out] gc_tile                Stomatal conductance on tiles (m s-1)
  !> @param[in,out] canhc_tile             Canopy heat capacity on tiles
  !> @param[in,out] tile_water_extract     Extraction of water from each tile
  !> @param[in,out] blend_height_tq        Blending height for wth levels
  !> @param[in,out] z0m_eff                Grid mean effective roughness length
  !> @param[in,out] ustar                  Friction velocity
  !> @param[in,out] soil_moist_avail       Available soil moisture for evaporation
  !> @param[in,out] zh_nonloc              Depth of non-local BL scheme
  !> @param[in,out] zhsc_2d                Height of decoupled layer top
  !> @param[in,out] z_lcl                  Height of the LCL (wtheta levels)
  !> @param[in,out] inv_depth              Depth of BL top inversion layer
  !> @param[in,out] qcl_at_inv_top         Cloud water at top of inversion
  !> @param[in,out] shallow_flag           Indicator of shallow convection
  !> @param[in,out] uw0_flux               'Zonal' surface momentum flux
  !> @param[in,out] vw0 flux               'Meridional' surface momentum flux
  !> @param[in,out] lcl_height             Height of lifting condensation level (w3 levels)
  !> @param[in,out] parcel_top             Height of surface based parcel ascent
  !> @param[in,out] level_parcel_top       Model level of parcel_top
  !> @param[in,out] wstar_2d               BL velocity scale
  !> @param[in,out] thv_flux               Surface flux of theta_v
  !> @param[in,out] parcel_buoyancy        Integral of parcel buoyancy
  !> @param[in,out] qsat_at_lcl            Saturation specific hum at LCL
  !> @param[in,out] bl_weight_1dbl         Blending weight to 1D BL scheme in the BL
  !> @param[in,out] bl_type_ind            Diagnosed BL types
  !> @param[in,out] snow_unload_rate       Unloading of snow from PFTs by wind
  !> @param[in]     albedo_obs_scaling     Scaling factor to adjust albedos by
  !> @param[in,out] level_ent              Level of surface mixed layer inversion
  !> @param[in,out] level_ent_dsc          Level of decoupled stratocumulus inversion
  !> @param[in,out] ent_we_lim             Rho * entrainment rate at surface ML inversion (kg m-2 s-1)
  !> @param[in,out] ent_t_frac             Fraction of time surface ML inversion is above level
  !> @param[in,out] ent_zrzi               Level height as fraction of DSC inversion height above DSC ML base
  !> @param[in,out] ent_we_lim_dsc         Rho * entrainment rate at DSC inversion (kg m-2 s-1)
  !> @param[in,out] ent_t_frac_dsc         Fraction of time DSC inversion is above level
  !> @param[in,out] ent_zrzi_dsc           Level height as fraction of DSC inversion height above DSC ML base
  !> @param[in] soil_clay_2d               Soil clay fraction
  !> @param[in] soil_sand_2d               Soil sand fraction
  !> @param[in] dust_div_mrel              Relative soil mass in CLASSIC size divisions
  !> @param[in,out] dust_div_flux          Dust emission fluxes in CLASSIC size divisions (kg m-2 s-1)
  !> @param[in]     day_of_year            The day of the year
  !> @param[in]     urbwrr                 Urban repeating width ratio
  !> @param[in]     urbhwr                 Urban height-to-width ratio
  !> @param[in]     urbhgt                 Urban building height
  !> @param[in]     urbztm                 Urban effective roughness length
  !> @param[in]     urbdisp                Urban displacement height
  !> @param[in]     urbemisw               Urban wall emissivity
  !> @param[in]     urbemisr               Urban road emissivity
  !> @param[in,out] zht                    Diagnostic: turb mixing height
  !> @param[in,out] z0h_eff                Diagnostic: Gridbox mean effective roughness length for scalars
  !> @param[in,out] oblen                  Diagnostic: Obukhov length
  !> @param[in]     ndf_wth                Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth               Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth                Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_w3                 Number of DOFs per cell for density space
  !> @param[in]     undf_w3                Number of unique DOFs for density space
  !> @param[in]     map_w3                 Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_2d                 Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d                Number of unique DOFs for 2D fields
  !> @param[in]     map_2d                 Dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_tile               Number of DOFs per cell for tiles
  !> @param[in]     undf_tile              Number of total DOFs for tiles
  !> @param[in]     map_tile               Dofmap for cell for surface tiles
  !> @param[in]     ndf_pft                Number of DOFs per cell for PFTs
  !> @param[in]     undf_pft               Number of total DOFs for PFTs
  !> @param[in]     map_pft                Dofmap for cell for PFTs
  !> @param[in]     ndf_sice               Number of DOFs per cell for sice
  !> @param[in]     undf_sice              Number of total DOFs for sice
  !> @param[in]     map_sice               Dofmap for cell for sice
  !> @param[in]     ndf_snow               Number of DOFs per cell for snow
  !> @param[in]     undf_snow              Number of total DOFs for snow
  !> @param[in]     map_snow               Dofmap for cell for snow
  !> @param[in]     ndf_soil               Number of DOFs per cell for soil levels
  !> @param[in]     undf_soil              Number of total DOFs for soil levels
  !> @param[in]     map_soil               Dofmap for cell for soil levels
  !> @param[in]     ndf_surf               Number of DOFs per cell for surface variables
  !> @param[in]     undf_surf              Number of unique DOFs for surface variables
  !> @param[in]     map_surf               Dofmap for the cell at the base of the column for surface variables
  !> @param[in]     ndf_smtile             Number of DOFs per cell for soil levels and tiles
  !> @param[in]     undf_smtile            Number of total DOFs for soil levels and tiles
  !> @param[in]     map_smtile             Dofmap for cell for soil levels and tiles
  !> @param[in]     ndf_bl                 Number of DOFs per cell for BL types
  !> @param[in]     undf_bl                Number of total DOFs for BL types
  !> @param[in]     map_bl                 Dofmap for cell for BL types
  !> @param[in]     ndf_scal               Number of DOFs per cell for albedo scaling
  !> @param[in]     undf_scal              Number of total DOFs for albedo scaling
  !> @param[in]     map_scal               Dofmap for cell at the base of the column
  !> @param[in]     ndf_ent                Number of DOFs per cell for entrainment levels
  !> @param[in]     undf_ent               Number of total DOFs for entrainment levels
  !> @param[in]     map_ent                Dofmap for cell for entrainment levels
  !> @param[in]     ndf_dust               Number of DOFs per cell for dust divisions
  !> @param[in]     undf_dust              Number of total DOFs for dust divisions
  !> @param[in]     map_dust               Dofmap for cell for dust divisions

  subroutine bl_exp_code(nlayers,                               &
                         theta_in_wth,                          &
                         rho_in_w3,                             &
                         wetrho_in_wth,                         &
                         exner_in_w3,                           &
                         exner_in_wth,                          &
                         u_in_w3,                               &
                         u_w3_stencil_size, u_w3_stencil,       &
                         v_in_w3,                               &
                         v_w3_stencil_size, v_w3_stencil,       &
                         w_in_wth,                              &
                         velocity_w2v,                          &
                         m_v_n,                                 &
                         m_cl_n,                                &
                         m_ci_n,                                &
                         height_w3,                             &
                         height_wth,                            &
                         shear,                                 &
                         delta,                                 &
                         max_diff_smag,                         &
                         zh_2d,                                 &
                         z0msea_2d,                             &
                         z0m_2d,                                &
                         ntml_2d,                               &
                         cumulus_2d,                            &
                         tile_fraction,                         &
                         tile_stencil_size, tile_stencil,       &
                         leaf_area_index,                       &
                         canopy_height,                         &
                         sd_orog_2d,                            &
                         peak_to_trough_orog,                   &
                         silhouette_area_orog,                  &
                         soil_albedo,                           &
                         soil_roughness,                        &
                         soil_moist_wilt,                       &
                         soil_moist_crit,                       &
                         soil_moist_sat,                        &
                         soil_thermal_cond,                     &
                         soil_suction_sat,                      &
                         clapp_horn_b,                          &
                         soil_respiration,                      &
                         thermal_cond_wet_soil,                 &
                         sea_ice_temperature,                   &
                         sea_ice_conductivity,                  &
                         tile_temperature,                      &
                         tile_lw_grey_albedo,                   &
                         tile_snow_mass,                        &
                         n_snow_layers,                         &
                         snow_depth,                            &
                         snow_layer_thickness,                  &
                         snow_layer_ice_mass,                   &
                         snow_layer_liq_mass,                   &
                         snow_layer_temp,                       &
                         surface_conductance,                   &
                         canopy_water,                          &
                         soil_temperature,                      &
                         soil_moisture,                         &
                         unfrozen_soil_moisture,                &
                         frozen_soil_moisture,                  &
                         tile_heat_flux,                        &
                         tile_moisture_flux,                    &
                         gross_prim_prod,                       &
                         net_prim_prod,                         &
                         cos_zen_angle,                         &
                         skyview,                               &
                         sw_up_tile,                            &
                         sw_down_surf,                          &
                         lw_down_surf,                          &
                         sw_down_surf_blue,                     &
                         dd_mf_cb,                              &
                         dtl_mphys,                             &
                         dmt_mphys,                             &
                         sw_heating_rate,                       &
                         lw_heating_rate,                       &
                         ozone,                                 &
                         cf_bulk,                               &
                         cf_liquid,                             &
                         rh_crit,                               &
                         dsldzm,                                &
                         wvar,                                  &
                         visc_m_blend,                          &
                         visc_h_blend,                          &
                         dw_bl,                                 &
                         rhokm_bl,                              &
                         surf_interp,                           &
                         rhokh_bl,                              &
                         tke_bl,                                &
                         ngstress_bl,                           &
                         bq_bl,                                 &
                         bt_bl,                                 &
                         moist_flux_bl,                         &
                         heat_flux_bl,                          &
                         dtrdz_tq_bl,                           &
                         rdz_tq_bl,                             &
                         fd_taux,                               &
                         fd_tauy,                               &
                         lmix_bl,                               &
                         gradrinr,                              &
                         alpha1_tile,                           &
                         ashtf_prime_tile,                      &
                         dtstar_tile,                           &
                         fraca_tile,                            &
                         z0h_tile,                              &
                         z0m_tile,                              &
                         rhokh_tile,                            &
                         chr1p5m_tile,                          &
                         resfs_tile,                            &
                         gc_tile,                               &
                         canhc_tile,                            &
                         tile_water_extract,                    &
                         blend_height_tq,                       &
                         z0m_eff,                               &
                         ustar,                                 &
                         soil_moist_avail,                      &
                         zh_nonloc,                             &
                         zhsc_2d,                               &
                         z_lcl,                                 &
                         inv_depth,                             &
                         qcl_at_inv_top,                        &
                         shallow_flag,                          &
                         uw0_flux,                              &
                         vw0_flux,                              &
                         lcl_height,                            &
                         parcel_top,                            &
                         level_parcel_top,                      &
                         wstar_2d,                              &
                         thv_flux,                              &
                         parcel_buoyancy,                       &
                         qsat_at_lcl,                           &
                         bl_weight_1dbl,                        &
                         bl_type_ind,                           &
                         snow_unload_rate,                      &
                         albedo_obs_scaling,                    &
                         level_ent,                             &
                         level_ent_dsc,                         &
                         ent_we_lim,                            &
                         ent_t_frac,                            &
                         ent_zrzi,                              &
                         ent_we_lim_dsc,                        &
                         ent_t_frac_dsc,                        &
                         ent_zrzi_dsc,                          &
                         soil_clay_2d,                          &
                         soil_sand_2d,                          &
                         dust_div_mrel,                         &
                         dust_div_flux,                         &
                         day_of_year,                           &
                         urbwrr,                                &
                         urbhwr,                                &
                         urbhgt,                                &
                         urbztm,                                &
                         urbdisp,                               &
                         urbemisw,                              &
                         urbemisr,                              &
                         zht,                                   &
                         z0h_eff,                               &
                         oblen,                                 &
                         ndf_wth,                               &
                         undf_wth,                              &
                         map_wth,                               &
                         ndf_w3,                                &
                         undf_w3,                               &
                         map_w3,                                &
                         ndf_2d,                                &
                         undf_2d,                               &
                         map_2d,                                &
                         ndf_tile, undf_tile, map_tile,         &
                         ndf_pft, undf_pft, map_pft,            &
                         ndf_sice, undf_sice, map_sice,         &
                         ndf_snow, undf_snow, map_snow,         &
                         ndf_soil, undf_soil, map_soil,         &
                         ndf_surf, undf_surf, map_surf,         &
                         ndf_smtile, undf_smtile, map_smtile,   &
                         ndf_bl, undf_bl, map_bl,               &
                         ndf_scal, undf_scal, map_scal,         &
                         ndf_ent, undf_ent, map_ent,            &
                         ndf_dust, undf_dust, map_dust)

    !---------------------------------------
    ! LFRic modules
    !---------------------------------------
    use jules_control_init_mod, only: n_land_tile, n_sea_ice_tile, &
         first_sea_tile, first_sea_ice_tile

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use ancil_info, only: ssi_pts, sea_pts, sice_pts, sice_pts_ncat,          &
                          rad_nband, nsoilt, nsurft, dim_cslayer, nmasst
    use atm_fields_bounds_mod, only: tdims, pdims_s, pdims_l, pdims
    use atm_step_local, only: dim_cs1, co2_dim_len, co2_dim_row
    use bl_option_mod, only: alpha_cd
    use c_kappai, only: kappai, de
    use cv_run_mod, only: i_convection_vn, i_convection_vn_6a,                &
                          cldbase_opt_dp, cldbase_opt_md
    use dust_parameters_mod, only: ndiv, ndivh
    use jules_deposition_mod, only: l_deposition
    use jules_irrig_mod,          only: irr_crop, irr_crop_doell
    use jules_radiation_mod,      only: l_albedo_obs
    use jules_sea_seaice_mod, only: nice, nice_use
    use jules_snow_mod, only: nsmax, cansnowtile
    use jules_soil_mod, only: ns_deep, l_bedrock
    use jules_soil_biogeochem_mod, only: dim_ch4layer, soil_bgc_model,        &
                                        soil_model_ecosse, l_layeredc
    use jules_surface_types_mod, only: npft, ntype, ncpft, nnpft
    use jules_surface_mod,        only: l_urban2t, l_flake_model
    use jules_vegetation_mod,     only: l_crop, l_triffid, l_phenol,          &
                                        l_use_pft_psi, can_rad_mod, l_acclim
    use nlsizes_namelist_mod, only: row_length, rows, land_field,             &
         sm_levels, ntiles, bl_levels
    use planet_constants_mod, only: p_zero, kappa, planet_radius
    use rad_input_mod, only: co2_mmr
    use jules_urban_mod,           only: l_moruses
    use timestep_mod, only: timestep
    use theta_field_sizes,        only: t_i_length, t_j_length,               &
                                        u_i_length,u_j_length,                &
                                        v_i_length,v_j_length
    use veg3_parm_mod,            only: l_veg3

    ! spatially varying fields used from modules
    use level_heights_mod, only: r_theta_levels, r_rho_levels
    use turb_diff_ctl_mod, only: visc_m, visc_h, max_diff, delta_smag
    use solinc_data, only: sky

    ! subroutines used
    use atmos_physics2_save_restore_mod, only: ap2_init_conv_diag
    use bl_diags_mod, only: BL_diag, dealloc_bl_imp, dealloc_bl_expl
    use conv_diag_6a_mod, only: conv_diag_6a
    use ni_bl_ctl_mod, only: ni_bl_ctl
    use sf_diags_mod, only: sf_diag, dealloc_sf_expl, dealloc_sf_imp
    use sparm_mod, only: sparm
    use tilepts_mod, only: tilepts
    use tr_mix_mod, only: tr_mix

    !---------------------------------------
    ! JULES modules
    !---------------------------------------
    use crop_vars_mod,            only: crop_vars_type, crop_vars_data_type,   &
                                        crop_vars_alloc, crop_vars_assoc,      &
                                        crop_vars_dealloc, crop_vars_nullify
    use prognostics,              only: progs_data_type, progs_type,           &
                                        prognostics_alloc, prognostics_assoc,  &
                                        prognostics_dealloc, prognostics_nullify
    use jules_vars_mod,           only: jules_vars_type, jules_vars_data_type, &
                                        jules_vars_alloc, jules_vars_assoc,    &
                                        jules_vars_dealloc, jules_vars_nullify
    use p_s_parms,                only: psparms_type, psparms_data_type,       &
                                        psparms_alloc, psparms_assoc,          &
                                        psparms_dealloc, psparms_nullify
    use trif_vars_mod,            only: trif_vars_type, trif_vars_data_type,   &
                                        trif_vars_assoc, trif_vars_alloc,      &
                                        trif_vars_dealloc, trif_vars_nullify
    use aero,                     only: aero_type, aero_data_type,             &
                                        aero_assoc, aero_alloc,                &
                                        aero_dealloc, aero_nullify
    use urban_param_mod,          only: urban_param_type,                      &
                                        urban_param_data_type,                 &
                                        urban_param_assoc, urban_param_alloc,  &
                                        urban_param_dealloc, urban_param_nullify
    use trifctl,                  only: trifctl_type, trifctl_data_type,       &
                                        trifctl_assoc, trifctl_alloc,          &
                                        trifctl_dealloc, trifctl_nullify
    use coastal,                  only: coastal_type, coastal_data_type,       &
                                        coastal_assoc, coastal_alloc,          &
                                        coastal_dealloc, coastal_nullify
    use lake_mod,                 only: lake_type, lake_data_type,             &
                                        lake_assoc, lake_alloc,                &
                                        lake_dealloc, lake_nullify
    use ancil_info,               only: ainfo_type, ainfo_data_type,           &
                                        ancil_info_assoc, ancil_info_alloc,    &
                                        ancil_info_dealloc, ancil_info_nullify
    use jules_forcing_mod,        only: forcing_type, forcing_data_type,       &
                                        forcing_assoc, forcing_alloc,          &
                                        forcing_dealloc, forcing_nullify
    use fluxes_mod,               only: fluxes_type, fluxes_data_type,         &
                                        fluxes_alloc, fluxes_assoc,            &
                                        fluxes_nullify, fluxes_dealloc
    use jules_chemvars_mod,       only: chemvars_type, chemvars_data_type,     &
                                        chemvars_alloc, chemvars_assoc,        &
                                        chemvars_dealloc, chemvars_nullify

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, day_of_year
    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3
    integer(kind=i_def), intent(in) :: map_wth(ndf_wth)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)

    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile)
    integer(kind=i_def), intent(in) :: ndf_pft, undf_pft
    integer(kind=i_def), intent(in) :: map_pft(ndf_pft)
    integer(kind=i_def), intent(in) :: ndf_soil, undf_soil
    integer(kind=i_def), intent(in) :: map_soil(ndf_soil)
    integer(kind=i_def), intent(in) :: ndf_sice, undf_sice
    integer(kind=i_def), intent(in) :: map_sice(ndf_sice)
    integer(kind=i_def), intent(in) :: ndf_snow, undf_snow
    integer(kind=i_def), intent(in) :: map_snow(ndf_snow)
    integer(kind=i_def), intent(in) :: ndf_smtile, undf_smtile
    integer(kind=i_def), intent(in) :: map_smtile(ndf_smtile)

    integer(kind=i_def), intent(in) :: ndf_surf, undf_surf, ndf_bl, undf_bl
    integer(kind=i_def), intent(in) :: map_surf(ndf_surf)
    integer(kind=i_def), intent(in) :: map_bl(ndf_bl)
    integer(kind=i_def), intent(in) :: ndf_scal, undf_scal
    integer(kind=i_def), intent(in) :: map_scal(ndf_scal)
    integer(kind=i_def), intent(in) :: ndf_ent, undf_ent
    integer(kind=i_def), intent(in) :: map_ent(ndf_ent)
    integer(kind=i_def), intent(in) :: ndf_dust, undf_dust
    integer(kind=i_def), intent(in) :: map_dust(ndf_dust)

    integer(kind=i_def), intent(in) :: u_w3_stencil_size, v_w3_stencil_size
    integer(kind=i_def), dimension(ndf_w3,u_w3_stencil_size), intent(in) :: u_w3_stencil
    integer(kind=i_def), dimension(ndf_w3,v_w3_stencil_size), intent(in) :: v_w3_stencil

    integer(kind=i_def), intent(in) :: tile_stencil_size
    integer(kind=i_def), dimension(ndf_tile,tile_stencil_size), intent(in) :: tile_stencil

    real(kind=r_def), dimension(undf_wth), intent(inout):: rh_crit,            &
                                                           dsldzm,             &
                                                           wvar, dw_bl

    real(kind=r_def), dimension(undf_wth), intent(inout):: visc_h_blend,       &
                                                           visc_m_blend,       &
                                                           rhokm_bl,           &
                                                           tke_bl,             &
                                                           ngstress_bl,        &
                                                           bq_bl, bt_bl,       &
                                                           dtrdz_tq_bl,        &
                                                           lmix_bl,            &
                                                           gradrinr
    real(kind=r_def), dimension(undf_w3),  intent(inout):: rhokh_bl,           &
                                                           moist_flux_bl,      &
                                                           heat_flux_bl,       &
                                                           rdz_tq_bl, fd_taux, &
                                                           fd_tauy
    real(kind=r_def), dimension(undf_w3),  intent(in)   :: rho_in_w3,          &
                                                           exner_in_w3,        &
                                                           u_in_w3, v_in_w3, &
                                                           height_w3
    real(kind=r_def), dimension(undf_wth), intent(in)   :: theta_in_wth,       &
                                                           wetrho_in_wth,      &
                                                           exner_in_wth,       &
                                                           w_in_wth,           &
                                                           velocity_w2v,       &
                                                           m_v_n, m_cl_n,      &
                                                           m_ci_n,             &
                                                           height_wth,         &
                                                           shear,              &
                                                           delta,              &
                                                           max_diff_smag,      &
                                                           dtl_mphys,dmt_mphys,&
                                                           sw_heating_rate,    &
                                                           lw_heating_rate,    &
                                                           cf_bulk, cf_liquid, &
                                                           ozone

    real(kind=r_def), dimension(undf_2d), intent(inout) :: zh_2d,              &
                                                           zhsc_2d,            &
                                                           z0msea_2d,          &
                                                           z0m_2d

    integer(kind=i_def), dimension(undf_2d), intent(inout) :: ntml_2d,         &
                                                              cumulus_2d,      &
                                                              blend_height_tq, &
                                                              shallow_flag,    &
                                                              level_parcel_top


    real(kind=r_def), dimension(undf_2d), intent(inout) :: z0m_eff,            &
                                                           ustar,              &
                                                           soil_moist_avail,   &
                                                           zh_nonloc,          &
                                                           z_lcl,              &
                                                           inv_depth,          &
                                                           qcl_at_inv_top,     &
                                                           uw0_flux,           &
                                                           vw0_flux,           &
                                                           lcl_height,         &
                                                           parcel_top,         &
                                                           wstar_2d,           &
                                                           thv_flux,           &
                                                           parcel_buoyancy,    &
                                                           qsat_at_lcl,        &
                                                           bl_weight_1dbl

    real(kind=r_def), intent(in) :: tile_fraction(undf_tile)
    real(kind=r_def), intent(inout) :: tile_temperature(undf_tile)
    real(kind=r_def), intent(in) :: tile_lw_grey_albedo(undf_tile)
    real(kind=r_def), intent(in) :: tile_snow_mass(undf_tile)
    integer(kind=i_def), intent(in) :: n_snow_layers(undf_tile)
    real(kind=r_def), intent(in) :: snow_depth(undf_tile)
    real(kind=r_def), intent(in) :: canopy_water(undf_tile)
    real(kind=r_def), intent(inout) :: tile_heat_flux(undf_tile)
    real(kind=r_def), intent(inout) :: tile_moisture_flux(undf_tile)
    real(kind=r_def), intent(in) :: sw_up_tile(undf_tile)
    real(kind=r_def), intent(in) :: albedo_obs_scaling(undf_scal)

    real(kind=r_def), intent(in) :: leaf_area_index(undf_pft)
    real(kind=r_def), intent(in) :: canopy_height(undf_pft)

    real(kind=r_def), intent(in) :: sea_ice_temperature(undf_sice)
    real(kind=r_def), intent(in) :: sea_ice_conductivity(undf_sice)

    real(kind=r_def), intent(in) :: sd_orog_2d(undf_2d)
    real(kind=r_def), intent(in) :: peak_to_trough_orog(undf_2d)
    real(kind=r_def), intent(in) :: silhouette_area_orog(undf_2d)
    real(kind=r_def), intent(in) :: soil_albedo(undf_2d)
    real(kind=r_def), intent(in) :: soil_roughness(undf_2d)
    real(kind=r_def), intent(in) :: soil_thermal_cond(undf_2d)
    real(kind=r_def), intent(inout) :: surface_conductance(undf_2d)
    real(kind=r_def), intent(in) :: cos_zen_angle(undf_2d)
    real(kind=r_def), intent(in) :: skyview(undf_2d)
    real(kind=r_def), intent(in) :: sw_down_surf(undf_2d)
    real(kind=r_def), intent(in) :: lw_down_surf(undf_2d)
    real(kind=r_def), intent(in) :: sw_down_surf_blue(undf_2d)
    real(kind=r_def), intent(in) :: dd_mf_cb(undf_2d)
    real(kind=r_def), intent(inout) :: gross_prim_prod(undf_2d)
    real(kind=r_def), intent(inout) :: net_prim_prod(undf_2d)
    real(kind=r_def), intent(inout) :: soil_respiration(undf_2d)
    real(kind=r_def), intent(inout) :: thermal_cond_wet_soil(undf_2d)

    ! Urban morphology fields (surface_fields)
    real(kind=r_def), intent(in) :: urbwrr(undf_2d)
    real(kind=r_def), intent(in) :: urbhwr(undf_2d)
    real(kind=r_def), intent(in) :: urbhgt(undf_2d)
    real(kind=r_def), intent(in) :: urbztm(undf_2d)
    real(kind=r_def), intent(in) :: urbdisp(undf_2d)
    real(kind=r_def), intent(in) :: urbemisw(undf_2d)
    real(kind=r_def), intent(in) :: urbemisr(undf_2d)

    real(kind=r_def), intent(in) :: soil_moist_wilt(undf_2d)
    real(kind=r_def), intent(in) :: soil_moist_crit(undf_2d)
    real(kind=r_def), intent(in) :: soil_moist_sat(undf_2d)
    real(kind=r_def), intent(in) :: soil_suction_sat(undf_2d)
    real(kind=r_def), intent(in) :: clapp_horn_b(undf_2d)
    real(kind=r_def), intent(in) :: soil_temperature(undf_soil)
    real(kind=r_def), intent(in) :: soil_moisture(undf_soil)
    real(kind=r_def), intent(in) :: unfrozen_soil_moisture(undf_soil)
    real(kind=r_def), intent(in) :: frozen_soil_moisture(undf_soil)

    real(kind=r_def), intent(in) :: snow_layer_thickness(undf_snow)
    real(kind=r_def), intent(in) :: snow_layer_ice_mass(undf_snow)
    real(kind=r_def), intent(in) :: snow_layer_liq_mass(undf_snow)
    real(kind=r_def), intent(in) :: snow_layer_temp(undf_snow)

    real(kind=r_def), intent(inout) :: tile_water_extract(undf_smtile)

    integer(kind=i_def), dimension(undf_bl), intent(inout) :: bl_type_ind
    real(kind=r_def), dimension(undf_pft),  intent(inout)  :: snow_unload_rate
    real(kind=r_def), dimension(undf_surf), intent(inout)  :: surf_interp
    real(kind=r_def), dimension(undf_tile), intent(inout):: alpha1_tile,      &
                                                            ashtf_prime_tile, &
                                                            dtstar_tile,      &
                                                            fraca_tile,       &
                                                            z0h_tile,         &
                                                            z0m_tile,         &
                                                            rhokh_tile,       &
                                                            chr1p5m_tile,     &
                                                            resfs_tile,       &
                                                            gc_tile,          &
                                                            canhc_tile

    integer(kind=i_def), dimension(undf_2d),   intent(inout)  :: level_ent
    integer(kind=i_def), dimension(undf_2d),   intent(inout)  :: level_ent_dsc
    real(kind=r_def), dimension(undf_ent),  intent(inout)  :: ent_we_lim
    real(kind=r_def), dimension(undf_ent),  intent(inout)  :: ent_t_frac
    real(kind=r_def), dimension(undf_ent),  intent(inout)  :: ent_zrzi
    real(kind=r_def), dimension(undf_ent),  intent(inout)  :: ent_we_lim_dsc
    real(kind=r_def), dimension(undf_ent),  intent(inout)  :: ent_t_frac_dsc
    real(kind=r_def), dimension(undf_ent),  intent(inout)  :: ent_zrzi_dsc

    real(kind=r_def), dimension(undf_2d),   intent(in)     :: soil_clay_2d
    real(kind=r_def), dimension(undf_2d),   intent(in)     :: soil_sand_2d
    real(kind=r_def), dimension(undf_dust), intent(in)     :: dust_div_mrel
    real(kind=r_def), dimension(undf_dust), intent(inout)  :: dust_div_flux

    real(kind=r_def), pointer, intent(inout) :: zht(:)
    real(kind=r_def), pointer, intent(inout) :: z0h_eff(:)
    real(kind=r_def), pointer, intent(inout) :: oblen(:)
    !-----------------------------------------------------------------------
    ! Local variables for the kernel
    !-----------------------------------------------------------------------
    ! loop counters etc
    integer(i_def) :: k, i, i_tile, i_sice, n, i_snow, j

    ! local switches and scalars
    integer(i_um) :: error_code
    real(r_um) :: weight1, weight2, weight3
    logical :: l_aero_classic, l_spec_z0, l_extra_call, l_jules_call,        &
         l_cape_opt

    ! profile fields from level 1 upwards
    real(r_um), dimension(row_length,rows,nlayers) ::                        &
         rho_wet, rho_dry, z_rho, z_theta, bulk_cloud_fraction,              &
         liquid_cloud_fraction, rho_wet_tq, u_p, v_p, rhcpt

    real(r_um), dimension(0:row_length+1,0:rows+1,nlayers) :: rho_wet_rsq,   &
         p_rho_levels, u_px, v_px, exner_rho_levels

    ! profile field on boundary layer levels
    real(r_um), dimension(row_length,rows,bl_levels) ::                      &
         fqw, ftl, rhokh, bq_gb, bt_gb, dtrdz_charney_grid,                  &
         dtrdz_u, rdz_charney_grid, rhokm_mix, w_mixed, w_flux

    real(r_um), dimension(0:row_length+1,0:rows+1,bl_levels) :: rhokm,       &
         tau_fd_x, tau_fd_y

    ! profile fields from level 2 upwards
    real(r_um), dimension(row_length,rows,2:bl_levels) :: rdz_u

    real(r_um), dimension(0:row_length+1,0:rows+1,2:bl_levels) :: f_ngstress,&
         rhogamu, rhogamv

    ! profile fields from level 0 upwards
    real(r_um), dimension(row_length,rows,0:nlayers) ::                      &
         p_theta_levels, etadot, w, q, qcl, qcf, theta, exner_theta_levels

    real(r_um), dimension(co2_dim_len,co2_dim_row) :: co2

    ! profile fields with a hard-wired 2
    real(r_um), dimension(row_length,rows,2,bl_levels) :: rad_hr, micro_tends

    ! single level real fields
    real(r_um), dimension(row_length,rows) ::                                &
         p_star, lw_down, cos_zenith_angle, tstar, zh_prev, ddmfx,           &
         zlcl, zhpar, flux_e, flux_h, z0msea, photosynth_act_rad,            &
         z0m_sice_fmd, z0m_sice_skin, tstar_sea,                             &
         zh, dzh, tstar_land, tstar_ssi, dtstar_sea, wstar, wthvs, ice_fract,&
         tstar_sice, u_0_p, v_0_p, zlcl_uv, qsat_lcl, delthvu,               &
         alpha1_sea, ashtf_prime_sea, bl_type_1, bl_type_2, bl_type_3,       &
         bl_type_4, bl_type_5, bl_type_6, bl_type_7, chr1p5m_sice, rhokh_sea,&
         u_s, uw0, vw0, z0hssi, z0mssi, zhnl, ti_sice, zeroes, surf_dep_flux

    real(r_um), dimension(0:row_length+1,0:rows+1) :: u_0_px, v_0_px, flandg,&
         flandfac, fseafac, cdr10m, rhokm_land, rhokm_ssi

    real(r_um), dimension(row_length,rows,3) :: t_frac, t_frac_dsc, we_lim,  &
         we_lim_dsc, zrzi, zrzi_dsc

    ! single level integer fields
    integer(i_um), dimension(row_length,rows) :: ntml, ntpar, k_blend_tq,    &
         kent, kent_dsc

    integer(i_um), dimension(0:row_length+1,0:rows+1) :: k_blend_uv

    ! single level logical fields
    logical, dimension(row_length,rows) :: land_sea_mask, cumulus, l_shallow

    ! fields on sea-ice categories
    real(r_um), dimension(row_length,rows,nice_use) ::                       &
         tstar_sice_ncat, ice_fract_ncat, alpha1_sice,                       &
         ashtf_prime, rhokh_sice, k_sice_ncat, ti_sice_ncat, dtstar_sice

    ! field on land points and soil levels
    real(r_um), dimension(land_field,sm_levels) :: soil_layer_moisture,      &
         smvccl_soilt, smvcwt_soilt, smvcst_soilt, sthf_soilt,               &
         sthu_soilt, t_soil_soilt

    ! fields on land points
    real(r_um), dimension(land_field) :: hcon_soilt, sil_orog_land_gb,       &
         ho2r2_orog_gb, sd_orog, z0m_soil_gb, albsoil_soilt, gs_gb,          &
         fland, npp_gb, smc_soilt, hcons_soilt

    ! integer fields on land points
    integer, dimension(land_field) :: land_index

    ! integer fields on land points and tile types
    integer, dimension(land_field, ntype) :: surft_index

    ! integer fields on number of tile types
    integer, dimension(ntype) :: surft_pts

    ! fields on land points and tiles
    real(r_um), dimension(land_field,ntiles) :: canopy_surft, catch_surft,   &
         catch_snow_surft, snow_surft, z0_surft, z0h_bare_surft,             &
         tstar_surft, frac_surft, dtstar_surft,                              &
         alpha1, ashtf_prime_surft, chr1p5m, fraca, resfs, rhokh_surft,      &
         canhc_surft

    ! fields on land points and pfts
    real(r_um), dimension(land_field,npft) :: canht_pft, lai_pft

    ! field on surface tiles and soil levels
    real(r_um), dimension(land_field,sm_levels,ntiles) :: wt_ext_surft

    ! Fields which are not used and only required for subroutine argument list,
    ! hence are unset in the kernel
    ! if they become set, please move up to be with other variables
    integer(i_um) :: asteps_since_triffid
    integer(i_um), parameter :: nscmdpkgs=15
    logical,       parameter :: l_scmdiags(nscmdpkgs)=.false.

    real(r_um), dimension(row_length,rows,nlayers) :: tgrad_bm
    real(r_um), dimension(row_length,rows,2:nlayers+1) :: bl_w_var

    real(r_um), dimension(row_length,rows,nlayers) :: tnuc_new

    real(r_um), dimension(row_length,rows,bl_levels) ::                      &
         e_trb, tsq_trb, qsq_trb, cov_trb, dtrdz_v

    real(r_um), dimension(row_length,rows,2:bl_levels) :: rdz_v

    real(r_um), dimension(row_length,rows,0:bl_levels-1) :: taux_p, tauy_p

    real(r_um), dimension(row_length,rows,0:nlayers) :: conv_prog_precip

    real(r_um), dimension(row_length,rows) ::                                &
         z0h_scm, z0m_scm, soil_clay, soil_sand, dust_mrel1,                 &
         dust_mrel2, dust_mrel3, dust_mrel4, dust_mrel5, dust_mrel6,         &
         zhpar_shcu, t1_sd, q1_sd, qcl_inv_top, w_max,                       &
         ql_ad, cin_undilute, cape_undilute,                                 &
         entrain_coef, ustar_in, g_ccp, h_ccp,  ccp_strength, fb_surf,       &
         charnock_w, aresist, cu_over_orog, resist_b, rho_aresist, rib_gb,   &
         shallowc, vshr, z0m_eff_gb, zhsc, tnuc_nlcl

    integer(i_um), dimension(row_length,rows) :: nlcl, conv_type, nbdsc, ntdsc

    logical, dimension(row_length,rows) :: no_cumulus, l_congestus,          &
                                           l_congestus2, l_mid

    real(r_um), dimension(row_length,rows,nice_use) :: radnet_sice

    real(r_um), dimension(row_length,rows,ndiv) :: dust_flux, r_b_dust

    real(r_um), dimension(land_field) :: emis_soil

    real(r_um), dimension(land_field) :: resp_s_tot_soilt

    real(r_um), dimension(land_field,dim_cs1) :: resp_s_acc_gb_um

    real(r_um), dimension(land_field,ntiles) :: tsurf_elev_surft, epot_surft,&
         aresist_surft, dust_emiss_frac, flake, gc_surft, resft,             &
         resist_b_surft, rho_aresist_surft, tile_frac, u_s_std_surft

    real(r_um), dimension(land_field,ntiles,ndivh) :: u_s_t_dry_tile,        &
         u_s_t_tile

    real(r_um), dimension(land_field,npft) :: resp_w_pft, g_leaf_acc_pft,    &
         npp_acc_pft, resp_w_acc_pft


    !-----------------------------------------------------------------------
    ! JULES Types
    !-----------------------------------------------------------------------
    type(crop_vars_type) :: crop_vars
    type(crop_vars_data_type) :: crop_vars_data
    type(progs_type) :: progs
    type(progs_data_type) :: progs_data
    type(psparms_type) :: psparms
    type(psparms_data_type) :: psparms_data
    type(trif_vars_type) :: trif_vars
    type(trif_vars_data_type) :: trif_vars_data
    type(aero_type) :: aerotype
    type(aero_data_type) :: aero_data
    type(urban_param_type) :: urban_param
    type(urban_param_data_type) :: urban_param_data
    type(trifctl_type) :: trifctltype
    type(trifctl_data_type) :: trifctl_data
    type(coastal_type) :: coast
    type(coastal_data_type) :: coastal_data
    type(lake_type) :: lake_vars
    type(lake_data_type) :: lake_data
    type(ainfo_type) :: ainfo
    type(ainfo_data_type) :: ainfo_data
    type(forcing_type) :: forcing
    type(forcing_data_type) :: forcing_data
    type(fluxes_type) :: fluxes
    type(fluxes_data_type) :: fluxes_data
    type(chemvars_type) :: chemvars
    type(chemvars_data_type) :: chemvars_data
    type(jules_vars_type) :: jules_vars
    type(jules_vars_data_type), TARGET :: jules_vars_data

    integer(i_um) :: ndry_dep_species  ! For Dry deposition, dummy for now

    !-----------------------------------------------------------------------
    ! Initialisation of JULES data and pointer types
    !-----------------------------------------------------------------------
    call crop_vars_alloc(land_field, t_i_length, t_j_length,                   &
                     nsurft, ncpft,nsoilt, sm_levels, l_crop, irr_crop,        &
                     irr_crop_doell, crop_vars_data)
    call crop_vars_assoc(crop_vars, crop_vars_data)

    call prognostics_alloc(land_field, t_i_length, t_j_length,                 &
                      nsurft, npft, nsoilt, sm_levels, ns_deep, nsmax,         &
                      dim_cslayer, dim_cs1, dim_ch4layer,                      &
                      nice, nice_use, soil_bgc_model, soil_model_ecosse,       &
                      l_layeredc, l_triffid, l_phenol, l_bedrock, l_veg3,      &
                      nmasst, nnpft, l_acclim, progs_data)
    call prognostics_assoc(progs,progs_data)

    call psparms_alloc(land_field,t_i_length,t_j_length,                       &
                   nsoilt,sm_levels,dim_cslayer,nsurft,npft,                   &
                   soil_bgc_model,soil_model_ecosse,l_use_pft_psi,             &
                   psparms_data)
    call psparms_assoc(psparms, psparms_data)

    call trif_vars_alloc(land_field,                                           &
                     npft,dim_cslayer,nsoilt,dim_cs1,                          &
                     l_triffid, l_phenol, trif_vars_data)
    call trif_vars_assoc(trif_vars, trif_vars_data)

    call aero_alloc(land_field,t_i_length,t_j_length,                          &
                nsurft,ndiv, aero_data)
    call aero_assoc(aerotype, aero_data)

    call urban_param_alloc(land_field, l_urban2t, l_moruses, urban_param_data)
    call urban_param_assoc(urban_param, urban_param_data)

    call trifctl_alloc(land_field,                                             &
                   npft,dim_cslayer,dim_cs1,nsoilt,trifctl_data)
    call trifctl_assoc(trifctltype, trifctl_data)

    call coastal_alloc(land_field,t_i_length,t_j_length,                       &
                   u_i_length,u_j_length,                                      &
                   v_i_length,v_j_length,                                      &
                   nice_use,nice,coastal_data)
    call coastal_assoc(coast, coastal_data)

    call lake_alloc(land_field, l_flake_model, lake_data)
    call lake_assoc(lake_vars, lake_data)

    call ancil_info_alloc(land_field,t_i_length,t_j_length,                    &
                      nice,nsoilt,ntype,                                       &
                      ainfo_data)
    call ancil_info_assoc(ainfo, ainfo_data)

    call forcing_alloc(t_i_length,t_j_length,u_i_length, u_j_length,           &
                        v_i_length, v_j_length, forcing_data)
    call forcing_assoc(forcing, forcing_data)

    call fluxes_alloc(land_field, t_i_length, t_j_length,                      &
                      nsurft, npft, nsoilt, sm_levels,                         &
                      nice, nice_use,                                          &
                      fluxes_data)
    call fluxes_assoc(fluxes, fluxes_data)

    ! Set num dry dep species to fixed = 1 for now
    ndry_dep_species = 1
    call chemvars_alloc(land_field, t_i_length, t_j_length, npft, ntype,       &
                        l_deposition, ndry_dep_species, chemvars_data)
    call chemvars_assoc(chemvars, chemvars_data)

    ! Note, jules_vars needs setting up after the change to pdims_s below so is
    ! not with the rest of these allocations.

    !-----------------------------------------------------------------------
    ! Initialisation of variables and arrays
    !-----------------------------------------------------------------------
    ! diagnostic flags
    error_code=0
    ! other logicals
    l_aero_classic=.false.
    l_extra_call=.false.
    ! surface forcing
    if ( flux_bc_opt == flux_bc_opt_specified_scalars ) then
      flux_e(:,:)=fixed_flux_e
      flux_h(:,:)=fixed_flux_h
    end if
    if ( sea_surf_alg == sea_surf_alg_fixed_roughness ) then
      l_spec_z0 = .true.
      z0m_scm = 0.01_r_um
      z0h_scm = 0.001_r_um
    else
      l_spec_z0 = .false.
    end if

    ! Size this with stencil for use in UM routines called
    pdims_s%i_start=0
    pdims_s%i_end=row_length+1
    pdims_s%j_start=0
    pdims_s%j_end=rows+1
    pdims_l%halo_i=1
    pdims_l%halo_j=1
    deallocate(r_theta_levels)
    deallocate(r_rho_levels)
    allocate(r_theta_levels(0:row_length+1,0:rows+1,0:nlayers), source=rmdi)
    allocate(r_rho_levels(0:row_length+1,0:rows+1,nlayers), source=rmdi)
    rho_wet_rsq = 1.0_r_um

    ! Initialise those fields whose contents will not be fully set
    fraca            = 0.0_r_um
    fluxes%fqw_surft = 0.0_r_um
    epot_surft       = 0.0_r_um

    call jules_vars_alloc(land_field,ntype,nsurft,rad_nband,nsoilt,sm_levels,  &
            t_i_length, t_j_length, npft, bl_levels, pdims_s, pdims,      &
            l_albedo_obs, cansnowtile, l_deposition,                      &
            jules_vars_data)
    call jules_vars_assoc(jules_vars,jules_vars_data)

    ! surf_hgt_surft needs to be initialised in this kernel so it has the
    ! correct value for this grid-point when used within Jules
    jules_vars%surf_hgt_surft = 0.0_r_um

    if (can_rad_mod == 6) then
      jules_vars%diff_frac = 0.4_r_um
    else
      jules_vars%diff_frac = 0.0_r_um
    end if
    !-----------------------------------------------------------------------
    ! Mapping of LFRic fields into UM variables
    !-----------------------------------------------------------------------

    ! Land tile fractions
    flandg = 0.0_r_um
    do i = 1, n_land_tile
      flandg(1,1) = flandg(1,1) + real(tile_fraction(tile_stencil(1,1)+i-1), r_um)
      frac_surft(1, i) = real(tile_fraction(tile_stencil(1,1)+i-1), r_um)
    end do
    if (buddy_sea == buddy_sea_on) then
      do i = 1, n_land_tile
        flandg(0,1) = flandg(0,1) + real(tile_fraction(tile_stencil(1,2)+i-1), r_um)
        flandg(0,0) = flandg(0,0) + real(tile_fraction(tile_stencil(1,3)+i-1), r_um)
        flandg(1,0) = flandg(1,0) + real(tile_fraction(tile_stencil(1,4)+i-1), r_um)
        flandg(2,0) = flandg(2,0) + real(tile_fraction(tile_stencil(1,5)+i-1), r_um)
        flandg(2,1) = flandg(2,1) + real(tile_fraction(tile_stencil(1,6)+i-1), r_um)
        flandg(2,2) = flandg(2,2) + real(tile_fraction(tile_stencil(1,7)+i-1), r_um)
        flandg(1,2) = flandg(1,2) + real(tile_fraction(tile_stencil(1,8)+i-1), r_um)
      end do
      ! Set winds in stencil region
      u_px(0,1,1) = u_in_w3(u_w3_stencil(1,2))
      u_px(0,0,1) = u_in_w3(u_w3_stencil(1,3))
      u_px(1,0,1) = u_in_w3(u_w3_stencil(1,4))
      u_px(2,0,1) = u_in_w3(u_w3_stencil(1,5))
      u_px(2,1,1) = u_in_w3(u_w3_stencil(1,6))
      u_px(2,2,1) = u_in_w3(u_w3_stencil(1,7))
      u_px(1,2,1) = u_in_w3(u_w3_stencil(1,8))
      v_px(0,1,1) = v_in_w3(v_w3_stencil(1,2))
      v_px(0,0,1) = v_in_w3(v_w3_stencil(1,3))
      v_px(1,0,1) = v_in_w3(v_w3_stencil(1,4))
      v_px(2,0,1) = v_in_w3(v_w3_stencil(1,5))
      v_px(2,1,1) = v_in_w3(v_w3_stencil(1,6))
      v_px(2,2,1) = v_in_w3(v_w3_stencil(1,7))
      v_px(1,2,1) = v_in_w3(v_w3_stencil(1,8))
      ! Deal with corners of cubed-sphere
      if (tile_stencil_size == 9) then
        do i = 1, n_land_tile
          flandg(0,2) = flandg(0,2) + real(tile_fraction(tile_stencil(1,9)+i-1), r_um)
        end do
        u_px(0,2,1) = u_in_w3(u_w3_stencil(1,9))
        v_px(0,2,1) = v_in_w3(v_w3_stencil(1,9))
      else
        flandg(0,2) = 1.0_r_um
      end if
    else
      flandfac = 1.0_r_def
      fseafac = 1.0_r_def
    end if

    ! Sea-ice fraction
    i_sice = 0
    ice_fract = 0.0_r_um
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      ice_fract = ice_fract + real(tile_fraction(tile_stencil(1,1)+i-1), r_um)
      ice_fract_ncat(1, 1, i_sice) = real(tile_fraction(tile_stencil(1,1)+i-1), r_um)
    end do

    ! Because Jules tests on flandg < 1, we need to ensure this is exactly
    ! 1 when no sea or sea-ice is present
    if ( tile_fraction(tile_stencil(1,1)+first_sea_tile-1) == 0.0_r_def .and. &
         ice_fract(1,1) == 0.0_r_um) then
      flandg(1,1) = 1.0_r_um
    end if

    fland(1) = flandg(1,1)

    ! Jules requires fractions with respect to the land area
    if (flandg(1, 1) > 0.0_r_um) then
      land_field = 1
      land_index = 1
      frac_surft(1, 1:n_land_tile) = frac_surft(1, 1:n_land_tile) / flandg(1, 1)
      land_sea_mask = .true.
    else
      land_field = 0
      land_index = 0
      land_sea_mask = .false.
    end if

    ! Set type_pts and type_index
    call tilepts(land_field, frac_surft, surft_pts, surft_index,ainfo%l_lice_point)

    ! Jules requires sea-ice fractions with respect to the sea area
    if (ice_fract(1, 1) > 0.0_r_um) then
      ice_fract(1, 1) = ice_fract(1, 1) / (1.0_r_um - flandg(1, 1))
      ice_fract_ncat(1, 1, 1:n_sea_ice_tile) &
           = ice_fract_ncat(1, 1, 1:n_sea_ice_tile) / (1.0_r_um - flandg(1, 1))
    end if

    ! combined sea and sea-ice index
    ssi_pts = 1
    if (flandg(1, 1) < 1.0_r_um) then
      ainfo%ssi_index = 1
    else
      ainfo%ssi_index = 0
    end if
    ainfo%fssi_ij = 1.0_r_um - flandg(1, 1)

    ! individual sea and sea-ice indices
    ! first set defaults
    sice_pts = 0
    ainfo%sice_index = 0
    ainfo%sice_frac = 0.0_r_um
    sea_pts = 0
    ainfo%sea_index = 0
    ainfo%sea_frac = 0.0_r_um
    ! then calculate based on state
    if (ainfo%ssi_index(1) > 0) then
      if (ice_fract(1, 1) > 0.0_r_um) then
        sice_pts = 1
        ainfo%sice_index = 1
        ainfo%sice_frac = ice_fract(1, 1)
      end if
      if (ice_fract(1, 1) < 1.0_r_um) then
        sea_pts = 1
        ainfo%sea_index = 1
        ainfo%sea_frac = 1.0_r_um - ainfo%sice_frac
      end if
    end if

    ! multi-category sea-ice index
    do n = 1, nice_use
      if (ainfo%ssi_index(1) > 0 .and. ice_fract_ncat(1, 1, n) > 0.0_r_um) then
        sice_pts_ncat(n) = 1
        ainfo%sice_index_ncat(1, n) = 1
        ainfo%sice_frac_ncat(1, n) = ice_fract_ncat(1, 1, n)
      else
        sice_pts_ncat(n) = 0
        ainfo%sice_index_ncat(1, n) = 0
        ainfo%sice_frac_ncat(1, n) = 0.0_r_um
      end if
    end do

    ! Land tile temperatures
    tstar_land = 0.0_r_um
    do i = 1, n_land_tile
      tstar_surft(1, i) = real(tile_temperature(map_tile(1)+i-1), r_um)
      tstar_land = tstar_land + frac_surft(1, i) * tstar_surft(1, i)
    end do

    ! LW grey emissivities
    do i = 1, n_land_tile
      fluxes%emis_surft(1, i) = 1.0_r_um - real(tile_lw_grey_albedo(map_tile(1)+i-1), r_um)
    end do

    if (emis_method_soil /= emis_method_soil_fixed) fluxes%l_emis_surft_set(soil)=.TRUE.

    ! Sea temperature
    ! Default to temperature over frozen sea as the initialisation
    ! that follows does not initialise sea points if they are fully
    ! frozen
    tstar_sea = tfs
    if (tile_fraction(tile_stencil(1,1)+first_sea_tile-1) > 0.0_r_def) then
      tstar_sea = real(tile_temperature(map_tile(1)+first_sea_tile-1), r_um)
    end if

    ! Sea-ice temperatures
    i_sice = 0
    tstar_sice = 0.0_r_um
    tstar_sice_ncat = 0.0_r_um
    if (ice_fract(1, 1) > 0.0_r_um) then
      do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
        i_sice = i_sice + 1
        tstar_sice_ncat(1, 1, i_sice) = real(tile_temperature(map_tile(1)+i-1), r_um)
        tstar_sice = tstar_sice &
                   + ice_fract_ncat(1,1,i_sice) * tstar_sice_ncat(1,1,i_sice) &
                   / ice_fract
      end do
    end if

    ! Sea & Sea-ice temperature
    tstar_ssi = (1.0_r_um - ice_fract) * tstar_sea + ice_fract * tstar_sice

    ! Grid-box mean surface temperature
    tstar = flandg(1,1) * tstar_land + (1.0_r_um - flandg(1,1)) * tstar_ssi

    ! Sea-ice conductivity and bulk temperature
    ti_sice = 0.0_r_um
    if (ice_fract(1, 1) > 0.0_r_um) then
      do i = 1, n_sea_ice_tile
        k_sice_ncat(1, 1, i) = real(sea_ice_conductivity(map_sice(1)+i-1), r_um)
        ti_sice_ncat(1, 1, i) = real(sea_ice_temperature(map_sice(1)+i-1), r_um)
        ti_sice = ti_sice &
                + ice_fract_ncat(1,1,i) * ti_sice_ncat(1,1,i) / ice_fract
      end do
    end if

    do n = 1, npft
      ! Leaf area index
      lai_pft(1, n) = real(leaf_area_index(map_pft(1)+n-1), r_um)
      ! Canopy height
      canht_pft(1, n) = real(canopy_height(map_pft(1)+n-1), r_um)
    end do

    ! Urban ancillaries
    urban_param%wrr_gb   = real(urbwrr(map_2d(1)), r_um)
    urban_param%hwr_gb   = real(urbhwr(map_2d(1)), r_um)
    urban_param%hgt_gb   = real(urbhgt(map_2d(1)), r_um)
    urban_param%ztm_gb   = real(urbztm(map_2d(1)), r_um)
    urban_param%disp_gb  = real(urbdisp(map_2d(1)), r_um)
    urban_param%emisw_gb = real(urbemisw(map_2d(1)), r_um)
    urban_param%emisr_gb = real(urbemisr(map_2d(1)), r_um)

    ! Roughness length (z0_tile)
    z0m_soil_gb = real(soil_roughness(map_2d(1)), r_um)
    call sparm(land_field, n_land_tile, surft_pts, surft_index,         &
               frac_surft, canht_pft, lai_pft, z0m_soil_gb,             &
               catch_snow_surft, catch_surft, z0_surft, z0h_bare_surft, &
               urban_param%ztm_gb)

    ! Snow-free soil albedo
    albsoil_soilt = real(soil_albedo(map_2d(1)), r_um)

    do i = 1, sm_levels
      ! Volumetric soil moisture at wilting point (smvcwt_soilt)
      smvcwt_soilt(1, i) = real(soil_moist_wilt(map_2d(1)), r_um)
      ! Volumetric soil moisture at critical point (smvccl_soilt)
      smvccl_soilt(1, i) = real(soil_moist_crit(map_2d(1)), r_um)
      ! Volumetric soil moisture at saturation (smvcst_soilt)
      smvcst_soilt(1, i) = real(soil_moist_sat(map_2d(1)), r_um)
      ! Saturated soil water suction (sathh_soilt)
      psparms%sathh_soilt(1, 1, i) = real(soil_suction_sat(map_2d(1)), r_um)
      ! Clapp and Hornberger b coefficient (bexp_soilt)
      psparms%bexp_soilt(1, 1, i) = real(clapp_horn_b(map_2d(1)), r_um)
      ! Soil temperature (t_soil_soilt)
      t_soil_soilt(1, i) = real(soil_temperature(map_soil(1)+i-1), r_um)
      ! Soil moisture content (kg m-2, soil_layer_moisture)
      soil_layer_moisture(1, i) = real(soil_moisture(map_soil(1)+i-1), r_um)
      ! Unfrozen soil moisture proportion (sthu_soilt)
      sthu_soilt(1, i) = real(unfrozen_soil_moisture(map_soil(1)+i-1), r_um)
      ! Frozen soil moisture proportion (sthf_soilt)
      sthf_soilt(1, i) = real(frozen_soil_moisture(map_soil(1)+i-1), r_um)
    end do

    ! Soil thermal conductivity (hcon_soilt)
    hcon_soilt = real(soil_thermal_cond(map_2d(1)), r_um)

    ! Soil ancils dependant on smvcst_soilt (soil moisture saturation limit)
    if ( smvcst_soilt(1,1) > 0.0_r_um ) then
      ainfo%l_soil_point = .true.
    else
      ainfo%l_soil_point = .false.
    end if

    ! Scaling factors needed for use in surface exchange code
    if (albedo_obs .and. flandg(1, 1) > 0.0_r_um) then
      i_tile = 0
      do n = 1, rad_nband
        do i = 1, n_land_tile
          jules_vars%albobs_scaling_surft(1,i,n) = &
               albedo_obs_scaling(map_scal(1)+i_tile)
          ! Counting from 0 so increment index here
          i_tile = i_tile + 1
        end do
      end do
    end if

    ! Cosine of the solar zenith angle
    cos_zenith_angle = real(cos_zen_angle(map_2d(1)), r_um)

    if (topography == topography_horizon) then
      ! Set skyview factor used internally by JULES
      sky = real(skyview(map_2d(1)), r_um)
    end if

    ! Downwelling LW radiation at surface
    lw_down = real(lw_down_surf(map_2d(1)), r_um)

    ! Net SW radiation on tiles
    do i = 1, n_land_tile
      fluxes%sw_surft(1, i) = real(sw_down_surf(map_2d(1)) - &
                            sw_up_tile(map_tile(1)+i-1), r_um)
    end do

    ! Net SW on open sea
    fluxes%sw_sea(1) = real(sw_down_surf(map_2d(1)) - &
                     sw_up_tile(map_tile(1)+first_sea_tile-1), r_um)

    ! Net SW on sea-ice
    i_sice = 0
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      fluxes%sw_sicat(1, i_sice) = real(sw_down_surf(map_2d(1)) - &
                                 sw_up_tile(map_tile(1)+i-1), r_um)
    end do

    ! photosynthetically active downwelling SW radiation
    photosynth_act_rad = real(sw_down_surf_blue(map_2d(1)), r_um)

    ! Ozone
    chemvars%o3_gb = real(ozone(map_wth(1)), r_um)

    ! Carbon dioxide
    co2 = co2_mmr

    ! Standard deviation of orography
    sd_orog = real(sd_orog_2d(map_2d(1)), r_um)

    ! Half of peak-to-trough height over root(2) of orography (ho2r2_orog_gb)
    ho2r2_orog_gb = real(peak_to_trough_orog(map_2d(1)), r_um)

    sil_orog_land_gb = real(silhouette_area_orog(map_2d(1)), r_um)

    ! Surface conductance (gs_gb)
    gs_gb = real(surface_conductance(map_2d(1)), r_um)

    ! Canopy water on each tile (canopy_surft)
    do i = 1, n_land_tile
      canopy_surft(1, i) = real(canopy_water(map_tile(1)+i-1), r_um)
    end do

    i_snow = 0
    do i = 1, n_land_tile
      snow_surft(1, i) = real(tile_snow_mass(map_tile(1)+i-1), r_um)
      progs%nsnow_surft(1, i) = n_snow_layers(map_tile(1)+i-1)
      progs%snowdepth_surft(1, i) = real(snow_depth(map_tile(1)+i-1), r_um)
      do j = 1, nsmax
        progs%ds_surft(1, i, j) = real(snow_layer_thickness(map_snow(1)+i_snow), r_um)
        progs%sice_surft(1, i, j) = real(snow_layer_ice_mass(map_snow(1)+i_snow), r_um)
        progs%sliq_surft(1, i, j) = real(snow_layer_liq_mass(map_snow(1)+i_snow), r_um)
        progs%tsnow_surft(1, i, j) = real(snow_layer_temp(map_snow(1)+i_snow), r_um)
        ! Counting from 0 so increment here
        i_snow = i_snow + 1
      end do
    end do

    ! Dust source information
    soil_clay(1, 1) = real(soil_clay_2d(map_2d(1)), r_um)
    soil_sand(1, 1) = real(soil_sand_2d(map_2d(1)), r_um)
    if (ndiv > 0) dust_mrel1(1, 1) = real(dust_div_mrel(map_dust(1) + 0), r_um)
    if (ndiv > 1) dust_mrel2(1, 1) = real(dust_div_mrel(map_dust(1) + 1), r_um)
    if (ndiv > 2) dust_mrel3(1, 1) = real(dust_div_mrel(map_dust(1) + 2), r_um)
    if (ndiv > 3) dust_mrel4(1, 1) = real(dust_div_mrel(map_dust(1) + 3), r_um)
    if (ndiv > 4) dust_mrel5(1, 1) = real(dust_div_mrel(map_dust(1) + 4), r_um)
    if (ndiv > 5) dust_mrel6(1, 1) = real(dust_div_mrel(map_dust(1) + 5), r_um)
    dust_flux = 0.0_r_um

    !-----------------------------------------------------------------------
    ! For the initial implementation we pass each individual column
    ! of data to an array sized (1,1,k) to match the UMs (i,j,k) data
    ! layout.
    ! assuming map_wth(1) points to level 0
    ! and map_w3(1) points to level 1
    !-----------------------------------------------------------------------
    do k = 1, nlayers
      ! potential temperature on theta levels
      theta(1,1,k) = theta_in_wth(map_wth(1) + k)
      ! wet density on theta and rho levels
      rho_wet_tq(1,1,k) = wetrho_in_wth(map_wth(1) + k)
      ! dry density on rho levels
      rho_dry(1,1,k) = rho_in_w3(map_w3(1) + k-1)
      ! pressure on rho and theta levels
      p_rho_levels(1,1,k) = p_zero*(exner_in_w3(map_w3(1) + k-1))**(1.0_r_def/kappa)
      p_theta_levels(1,1,k) = p_zero*(exner_in_wth(map_wth(1) + k))**(1.0_r_def/kappa)
      ! exner pressure on rho and theta levels
      exner_rho_levels(1,1,k) = exner_in_w3(map_w3(1) + k-1)
      exner_theta_levels(1,1,k) = exner_in_wth(map_wth(1) + k)
      ! u wind on rho levels
      u_px(1,1,k) = u_in_w3(u_w3_stencil(1,1) + k-1)
      ! v wind on rho levels
      v_px(1,1,k) = v_in_w3(v_w3_stencil(1,1) + k-1)
      ! w wind on theta levels
      w(1,1,k) = w_in_wth(map_wth(1) + k)
      ! height of rho levels from centre of planet
      r_rho_levels(:,:,k) = height_w3(map_w3(1) + k-1) + planet_radius
      ! height of theta levels from centre of planet
      r_theta_levels(:,:,k) = height_wth(map_wth(1) + k) + planet_radius
      ! water vapour mixing ratio
      q(1,1,k) = m_v_n(map_wth(1) + k)
      ! cloud liquid mixing ratio
      qcl(1,1,k) = m_cl_n(map_wth(1) + k)
      ! cloud ice mixing ratio
      qcf(1,1,k) = m_ci_n(map_wth(1) + k)
      ! cloud fraction variables
      bulk_cloud_fraction(1,1,k) = cf_bulk(map_wth(1) + k)
      liquid_cloud_fraction(1,1,k) = cf_liquid(map_wth(1) + k)
    end do

    if ( smagorinsky ) then
      delta_smag(1,1) = delta(map_wth(1))
      max_diff(1,1) = max_diff_smag(map_wth(1))
      do k = 1, nlayers
        visc_m(1,1,k) = shear(map_wth(1) + k)
        visc_h(1,1,k) = shear(map_wth(1) + k)
      end do
    end if

    ! surface pressure
    p_theta_levels(1,1,0) = p_zero*(exner_in_wth(map_wth(1) + 0))**(1.0_r_def/kappa)
    p_star(1,1) = p_theta_levels(1,1,0)
    exner_theta_levels(1,1,0) = exner_in_wth(map_wth(1) + 0)
    ! near surface potential temperature
    theta(1,1,0) = theta_in_wth(map_wth(1) + 0)
    ! non-halo u and v winds
    u_p(1,1,:) = u_px(1,1,:)
    v_p(1,1,:) = v_px(1,1,:)
    ! surface currents
    u_0_px = 0.0
    v_0_px = 0.0
    u_0_p = 0.0
    v_0_p = 0.0
    ! near surface moisture fields
    q(1,1,0) = m_v_n(map_wth(1) + 0)
    qcl(1,1,0) = m_cl_n(map_wth(1) + 0)
    qcf(1,1,0) = m_ci_n(map_wth(1) + 0)
    ! surface height
    r_theta_levels(:,:,0) = height_wth(map_wth(1) + 0) + planet_radius
    ! height of levels above surface
    z_rho(1,1,:) = r_rho_levels(1,1,:)-r_theta_levels(1,1,0)
    z_theta(1,1,:) = r_theta_levels(1,1,1:nlayers)-r_theta_levels(1,1,0)
    ! vertical velocity
    w(1,1,0) = w_in_wth(map_wth(1) + 0)
    etadot(1,1,:) = velocity_w2v(map_wth(1):map_wth(1)+nlayers) / z_theta(1,1,nlayers)

    !-----------------------------------------------------------------------
    ! Things saved from one timestep to the next
    !-----------------------------------------------------------------------
    ! previous BL height
    zh(1,1) = zh_2d(map_2d(1))
    zh_prev = zh
    ! surface roughness
    z0msea(1,1) = z0msea_2d(map_2d(1))
    ! downdraft at cloud base
    ddmfx = dd_mf_cb(map_2d(1))
    ! needed to ensure neutral diagnostics can be calculated
    sf_diag%suv10m_n  = .true.
    sf_diag%l_u10m_n  = sf_diag%suv10m_n
    sf_diag%l_v10m_n  = sf_diag%suv10m_n
    sf_diag%l_mu10m_n = sf_diag%suv10m_n
    sf_diag%l_mv10m_n = sf_diag%suv10m_n
    ! needed to ensure z0h_eff is saved if wanted
    sf_diag%l_z0h_eff_gb = .not. associated(z0h_eff, empty_real_data)
    ! needed to ensure zht is saved if wanted
    bl_diag%l_zht      = .not. associated(zht, empty_real_data)
    bl_diag%l_oblen    = .not. associated(oblen, empty_real_data)
    bl_diag%l_weight1d = .true.
    !-----------------------------------------------------------------------
    ! Things saved from other parametrization schemes on this timestep
    !-----------------------------------------------------------------------
    do k = 1, bl_levels
      ! microphysics tendancy terms
      micro_tends(1,1,1,k) = dtl_mphys(map_wth(1)+k)/timestep
      micro_tends(1,1,2,k) = dmt_mphys(map_wth(1)+k)/timestep
      ! radiation tendancy terms
      rad_hr(1,1,1,k) = lw_heating_rate(map_wth(1)+k)
      rad_hr(1,1,2,k) = sw_heating_rate(map_wth(1)+k)
    end do

    !-----------------------------------------------------------------------
    ! code below here should mimic the call from the UMs atmos_physics2
    !-----------------------------------------------------------------------
    if (use_jules_flux) then
      l_jules_call=.true.
      call NI_bl_ctl (                                                         &
    !     IN parameters for SISL scheme
         outer_iterations, l_jules_call,                                       &
    !     IN time stepping information
         day_of_year,                                                          &
    !     IN switches
         L_aero_classic,                                                       &
    !     IN data fields.
         p_rho_levels, p_theta_levels, rho_wet_rsq,rho_wet,rho_dry, u_p, v_p,  &
         u_px, v_px, u_0_px, v_0_px,                                           &
         land_sea_mask, q, qcl, qcf, p_star, theta, exner_theta_levels, rad_hr,&
         micro_tends, soil_layer_moisture, rho_wet_tq, z_rho, z_theta,         &
    !     IN ancillary fields and fields needed to be kept from tstep to tstep
         hcon_soilt, smvccl_soilt, smvcwt_soilt, smvcst_soilt,                 &
         sthf_soilt, sthu_soilt, sil_orog_land_gb,                             &
         ho2r2_orog_gb, sd_orog, ice_fract_ncat, k_sice_ncat,                  &
         land_index, photosynth_act_rad, z0m_sice_fmd, z0m_sice_skin,          &
         soil_clay,soil_sand,dust_mrel1,dust_mrel2,                            &
         dust_mrel3,dust_mrel4,dust_mrel5,dust_mrel6,                          &
    !     IN additional variables for JULES
         canopy_surft, catch_surft, catch_snow_surft, snow_surft,              &
         z0_surft, z0h_bare_surft,                                             &
         z0m_soil_gb, lw_down, tstar_surft, tsurf_elev_surft,                  &
         co2,                                                                  &
         asteps_since_triffid,                                                 &
         frac_surft,canht_pft,lai_pft,fland,flandg,                            &
         albsoil_soilt, cos_zenith_angle,                                      &
    !     IN: input from the wave model
         charnock_w,                                                           &
    !     IN everything not covered so far
         t_soil_soilt, ti_sice, ti_sice_ncat, tstar, zh_prev,                  &
         ddmfx, bulk_cloud_fraction, liquid_cloud_fraction, zhpar, zlcl,       &
    !     IN SCM namelist data
         L_spec_z0, z0m_scm, z0h_scm, flux_e, flux_h, ustar_in,                &
    !     SCM diagnostics and STASH
         nSCMDpkgs, L_SCMDiags, BL_diag, sf_diag,                              &
    !     INOUT data
         gs_gb,z0msea,w,etadot,tstar_sea,tstar_sice_ncat,zh,dzh,               &
         cumulus, ntml,ntpar,l_shallow,                                        &
    !     INOUT additional variables for JULES
         g_leaf_acc_pft,npp_acc_pft,resp_w_acc_pft,resp_s_acc_gb_um,           &
    !     JULES TYPES (IN OUT)
         crop_vars, psparms, ainfo, trif_vars, aerotype, urban_param,          &
         progs, trifctltype, coast, jules_vars,                                &
         fluxes,                                                               &
         lake_vars,                                                            &
         forcing,                                                              &
         chemvars,                                                             &
    !     INOUT variables for TKE based turbulence schemes
         e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                         &
    !     INOUT variables from bdy_expl1 needed elsewhere
         bq_gb, bt_gb, dtrdz_charney_grid,rdz_charney_grid,                    &
         dtrdz_u, dtrdz_v, rdz_u, rdz_v, k_blend_tq, k_blend_uv,               &
    !     INOUT variables from Jules needed elsewhere
         flandfac,fseafac,rhokm_land,rhokm_ssi,cdr10m,                         &
         fqw, ftl, rib_gb, vshr, z0m_eff_gb, r_b_dust,                         &
         rho_aresist,aresist,resist_b, rhokm,rhokh,                            &
    !     INOUT variables required in IMP_SOLVER
         alpha1_sea, alpha1_sice, ashtf_prime_sea, ashtf_prime, u_s,           &
    !     INOUT additional variables for JULES
         radnet_sice,rho_aresist_surft,                                        &
         aresist_surft, resist_b_surft, alpha1, ashtf_prime_surft,             &
         epot_surft,                                                           &
         fraca,resfs,resft,rhokh_surft,rhokh_sice,rhokh_sea,                   &
         z0hssi,z0mssi,chr1p5m,chr1p5m_sice,smc_soilt,                         &
         npp_gb, resp_s_tot_soilt,                                             &
         resp_w_pft, gc_surft, canhc_surft, wt_ext_surft, flake,               &
         surft_index, surft_pts,                                               &
         tile_frac, tstar_land, tstar_ssi, dtstar_surft,                       &
         dtstar_sea, dtstar_sice, hcons_soilt, emis_soil,                      &
         t1_sd, q1_sd, fb_surf,                                                &
    !     OUT variables for message passing
         tau_fd_x, tau_fd_y, rhogamu, rhogamv, f_ngstress,                     &
    !     OUT diagnostics (done after implicit solver)
         zhnl, shallowc,cu_over_orog,bl_type_1,bl_type_2,bl_type_3,            &
         bl_type_4,bl_type_5,bl_type_6, bl_type_7, bl_w_var,                   &
    !     OUT variables required for mineral dust scheme
         dust_flux,dust_emiss_frac, u_s_t_tile,u_s_t_dry_tile,                 &
         u_s_std_surft, kent, we_lim, t_frac, zrzi,                            &
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,                     &
    !     OUT fields
         nbdsc,ntdsc,wstar,wthvs,uw0,vw0,taux_p,tauy_p,rhcpt,tgrad_bm          &
       )
    end if

    ! Use  convection switches to decide the value of  L_cape_opt
    if (i_convection_vn == i_convection_vn_6a ) then
      L_cape_opt = ( (cldbase_opt_dp == 3) .or. (cldbase_opt_md == 3) .or. &
                     (cldbase_opt_dp == 4) .or. (cldbase_opt_md == 4) .or. &
                     (cldbase_opt_dp == 5) .or. (cldbase_opt_md == 5) .or. &
                     (cldbase_opt_dp == 6) .or. (cldbase_opt_md == 6) )
    else
      L_cape_opt = .false.
    end if

    call ap2_init_conv_diag( rows, row_length, ntml, ntpar, nlcl, cumulus,  &
        l_shallow, l_mid, delthvu, ql_ad, zhpar, dzh, qcl_inv_top,          &
        zlcl, zlcl_uv, conv_type, no_cumulus, w_max, w, L_cape_opt)

    call conv_diag_6a(                                                  &
    !     IN Parallel variables
            row_length, rows                                            &
    !     IN model dimensions.
          , bl_levels                                                   &
          , p_rho_levels, p_theta_levels(1,1,1),exner_rho_levels        &
          , rho_wet, rho_wet_tq, z_theta, z_rho                         &
    !     IN Model switches
          , l_extra_call                                                &
          , no_cumulus                                                  &
    !     IN cloud data
          , qcf(1:row_length,1:rows,1:tdims%k_end)                      &
          , qcl(1:row_length,1:rows,1:tdims%k_end), bulk_cloud_fraction &
    !     IN everything not covered so far :
          , p_star, q(1:row_length,1:rows,1:tdims%k_end)                &
          , theta(tdims%i_start:tdims%i_end,                            &
              tdims%j_start:tdims%j_end,1:tdims%k_end)                  &
          , exner_theta_levels(tdims%i_start:tdims%i_end,               &
              tdims%j_start:tdims%j_end, 1:tdims%k_end)                 &
          , u_p, v_p, u_0_p, v_0_p                                      &
          , tstar_land, tstar_sea, tstar_sice, z0msea                   &
          , flux_e, flux_h, ustar_in, L_spec_z0, z0m_scm, z0h_scm       &
          , tstar, land_sea_mask, flandg, ice_fract                     &
          , w, w_max                                                    &
          , conv_prog_precip                                            &
          , g_ccp, h_ccp, ccp_strength                                  &
    !     IN surface fluxes
          , fb_surf, u_s                                                &
    !     SCM Diagnostics (dummy values in full UM)
          , nSCMDpkgs,L_SCMDiags                                        &
    !     OUT data required elsewhere in UM system :
          , zh,zhpar,dzh,qcl_inv_top,zlcl,zlcl_uv,delthvu,ql_ad         &
          , ntml,ntpar,nlcl                                             &
          , cumulus,l_shallow,l_congestus,l_congestus2                  &
          , conv_type                                                   &
          , CIN_undilute,CAPE_undilute, wstar, wthvs                    &
          , entrain_coef, qsat_lcl                                      &
          , Error_code                                                  &
          , tnuc_new, tnuc_nlcl  )

    l_jules_call=.false.
    call NI_bl_ctl (                                                    &
    !     IN parameters for SISL scheme
         outer_iterations, l_jules_call,                                &
    !     IN time stepping information
         day_of_year,                                                   &
    !     IN switches
         L_aero_classic,                                                &
    !     IN data fields.
         p_rho_levels, p_theta_levels, rho_wet_rsq,rho_wet,rho_dry, u_p,&
         v_p,                                                           &
         u_px, v_px, u_0_px, v_0_px,                                    &
         land_sea_mask, q, qcl, qcf, p_star, theta, exner_theta_levels, &
         rad_hr,                                                        &
         micro_tends, soil_layer_moisture, rho_wet_tq, z_rho, z_theta,  &
    !     IN ancillary fields and fields needed to be kept from tstep to tstep
         hcon_soilt, smvccl_soilt, smvcwt_soilt, smvcst_soilt,          &
         sthf_soilt, sthu_soilt, sil_orog_land_gb,                      &
    !-------------------------------------------------------------------------
         ho2r2_orog_gb, sd_orog, ice_fract_ncat, k_sice_ncat,           &
         land_index, photosynth_act_rad, z0m_sice_fmd, z0m_sice_skin,   &
         soil_clay,soil_sand,dust_mrel1,dust_mrel2,                     &
         dust_mrel3,dust_mrel4,dust_mrel5,dust_mrel6,                   &
    !     IN additional variables for JULES
         canopy_surft, catch_surft, catch_snow_surft, snow_surft,       &
         z0_surft, z0h_bare_surft,                                      &
         z0m_soil_gb, lw_down, tstar_surft, tsurf_elev_surft, &
         co2,                                                           &
         asteps_since_triffid,                                          &
         frac_surft,canht_pft,lai_pft,fland,flandg,                     &
         albsoil_soilt, cos_zenith_angle,                               &
    !     IN: input from the wave model
         charnock_w,                                                    &
    !     IN everything not covered so far
         t_soil_soilt, ti_sice, ti_sice_ncat, tstar, zh_prev,           &
         ddmfx, bulk_cloud_fraction, liquid_cloud_fraction, zhpar, zlcl,&
    !     IN SCM namelist data
         L_spec_z0, z0m_scm, z0h_scm, flux_e, flux_h, ustar_in,         &
    !     SCM diagnostics and STASH
         nSCMDpkgs, L_SCMDiags, BL_diag, sf_diag,                       &
    !     INOUT data
         gs_gb,z0msea,w,etadot,tstar_sea,tstar_sice_ncat,zh,dzh,        &
         cumulus, ntml,ntpar,l_shallow,                                 &
    !     INOUT additional variables for JULES
         g_leaf_acc_pft,npp_acc_pft,resp_w_acc_pft,resp_s_acc_gb_um,    &
    !     JULES TYPES (IN OUT)
         crop_vars, psparms, ainfo, trif_vars, aerotype, urban_param,   &
         progs, trifctltype, coast, jules_vars,                         &
         fluxes,                                                        &
         lake_vars,                                                     &
         forcing,                                                       &
         chemvars,                                                      &
    !     INOUT variables for TKE based turbulence schemes
         e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                  &
      ! INOUT variables from bdy_expl1 needed elsewhere
        bq_gb, bt_gb, dtrdz_charney_grid,rdz_charney_grid,              &
        dtrdz_u, dtrdz_v, rdz_u, rdz_v, k_blend_tq, k_blend_uv,         &
      ! INOUT variables from Jules needed elsewhere
        flandfac,fseafac,rhokm_land,rhokm_ssi,cdr10m,                   &
        fqw, ftl, rib_gb, vshr, z0m_eff_gb, r_b_dust,                   &
        rho_aresist,aresist,resist_b, rhokm,rhokh,                      &
      ! INOUT variables required in IMP_SOLVER
        alpha1_sea, alpha1_sice, ashtf_prime_sea, ashtf_prime, u_s,     &
      ! INOUT additional variables for JULES
        radnet_sice,rho_aresist_surft,                                  &
        aresist_surft, resist_b_surft, alpha1, ashtf_prime_surft,       &
        epot_surft,                                                     &
        fraca,resfs,resft,rhokh_surft,rhokh_sice,rhokh_sea,             &
        z0hssi,z0mssi,chr1p5m,chr1p5m_sice,smc_soilt,                   &
        npp_gb, resp_s_tot_soilt,                                       &
        resp_w_pft, gc_surft, canhc_surft, wt_ext_surft, flake,         &
        surft_index, surft_pts,                                         &
        tile_frac, tstar_land, tstar_ssi, dtstar_surft,                 &
        dtstar_sea, dtstar_sice, hcons_soilt, emis_soil,                &
        t1_sd, q1_sd, fb_surf,                                          &
      ! OUT variables for message passing
        tau_fd_x, tau_fd_y, rhogamu, rhogamv, f_ngstress,               &
      ! OUT diagnostics (done after implicit solver)
        zhnl, shallowc,cu_over_orog,bl_type_1,bl_type_2,bl_type_3,      &
        bl_type_4,bl_type_5,bl_type_6, bl_type_7, bl_w_var,             &
      ! OUT variables required for mineral dust scheme
        dust_flux,dust_emiss_frac, u_s_t_tile,u_s_t_dry_tile,           &
        u_s_std_surft, kent, we_lim, t_frac, zrzi,                      &
        kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,               &
      ! OUT fields
        nbdsc,ntdsc,wstar,wthvs,uw0,vw0,taux_p,tauy_p,rhcpt,tgrad_bm    &
     )

    if (bl_mix_w) then

      ! Interpolate rhokm from theta-levels to rho-levels, as is done for
      ! rhokh in bdy_expl2.
      ! NOTE: rhokm is defined on theta-levels with the k-indexing offset by
      ! 1 compared to the rest of the UM (k=1 is the surface).

      ! Bottom model-level is surface in both arrays, so no interp needed
      ! (for rhokh_mix, this is done in the JULES routine sf_impl2_jls).
      rhokm_mix(1,1,1) = rhokm(1,1,1)
      do k = 2, bl_levels-1
        weight1 = r_theta_levels(1,1,k) - r_theta_levels(1,1,k-1)
        weight2 = r_theta_levels(1,1,k) - r_rho_levels(1,1,k)
        weight3 = r_rho_levels(1,1,k)   - r_theta_levels(1,1,k-1)
        rhokm_mix(1,1,k) = (weight3/weight1) * rhokm(1,1,k+1) &
                         + (weight2/weight1) * rhokm(1,1,k)
        ! Scale exchange coefficients by 1/dz factor, as is done for
        ! rhokh_mix in bdy_impl4
        ! (note this doesn't need to be done for the surface exchange coef)
        rhokm_mix(1,1,k) = rhokm_mix(1,1,k) * rdz_charney_grid(1,1,k)
      end do
      k = bl_levels
      weight1 = r_theta_levels(1,1,k) - r_theta_levels(1,1,k-1)
      weight2 = r_theta_levels(1,1,k) - r_rho_levels(1,1,k)
      ! Assume rhokm(BL_LEVELS+1) is zero
      rhokm_mix(1,1,k) = (weight2/weight1) * rhokm(1,1,k)
      ! Scale exchange coefficients by 1/dz factor, as is done for
      ! rhokh_mix in bdy_impl4
      rhokm_mix(1,1,k) = rhokm_mix(1,1,k) * rdz_charney_grid(1,1,k)

      zeroes = 0.0_r_um
      w_mixed(:,:,:) = w(:,:,1:bl_levels)

      call  tr_mix (                                                 &
           ! IN fields
           bl_levels, alpha_cd, rhokm_mix(:,:,2:), rhokm_mix(:,:,1), &
           dtrdz_charney_grid, zeroes, zeroes,                       &
           kent, we_lim, t_frac, zrzi,                               &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,               &
           zhnl, zhsc, z_rho,                                        &
           ! INOUT / OUT fields
           w_mixed, w_flux, surf_dep_flux                            &
           )

      do k = 1, bl_levels
        dw_bl(map_wth(1)+k) = w_mixed(1,1,k) - w(1,1,k)
      end do

    end if

    ! 2D variables that need interpolating to cell faces
    surf_interp(map_surf(1)+0) = flandg(1,1)
    surf_interp(map_surf(1)+1) = rhokm_land(1,1)
    surf_interp(map_surf(1)+2) = rhokm_ssi(1,1)
    surf_interp(map_surf(1)+3) = flandfac(1,1)
    surf_interp(map_surf(1)+4) = fseafac(1,1)
    surf_interp(map_surf(1)+5) = zhnl(1,1)
    surf_interp(map_surf(1)+6) = fb_surf(1,1)
    surf_interp(map_surf(1)+7) = real(k_blend_uv(1,1)-1, r_def)
    surf_interp(map_surf(1)+8) = cdr10m(1,1)
    surf_interp(map_surf(1)+9) = real(sf_diag%cd10m_n(1,1), r_def)
    surf_interp(map_surf(1)+10)= real(sf_diag%cdr10m_n(1,1), r_def)
    do k=1,bl_levels
      rhokm_bl(map_wth(1) + k-1) = rhokm(1,1,k)
      rhokh_bl(map_w3(1) + k-1) = rhokh(1,1,k)
      bq_bl(map_wth(1) + k-1) = bq_gb(1,1,k)
      bt_bl(map_wth(1) + k-1) = bt_gb(1,1,k)
      moist_flux_bl(map_w3(1) + k-1) = fqw(1,1,k)
      heat_flux_bl(map_w3(1) + k-1) = ftl(1,1,k)
      dtrdz_tq_bl(map_wth(1) + k) = dtrdz_charney_grid(1,1,k)
      rdz_tq_bl(map_w3(1) + k-1) = rdz_charney_grid(1,1,k)
    end do

    level_ent(map_2d(1)) = int( kent(1,1), i_def )
    level_ent_dsc(map_2d(1)) = int( kent_dsc(1,1), i_def )
    do k = 1, 3
       ent_we_lim(map_ent(1) + k - 1) = real( we_lim(1,1,k), r_def )
       ent_t_frac(map_ent(1) + k - 1) = real( t_frac(1,1,k), r_def )
       ent_zrzi(map_ent(1) + k - 1) = real( zrzi(1,1,k), r_def )
       ent_we_lim_dsc(map_ent(1) + k - 1) = real( we_lim_dsc(1,1,k), r_def )
       ent_t_frac_dsc(map_ent(1) + k - 1) = real( t_frac_dsc(1,1,k), r_def )
       ent_zrzi_dsc(map_ent(1) + k - 1) = real( zrzi_dsc(1,1,k), r_def )
    end do

    if (formdrag == formdrag_dist_drag) then
      do k=1,bl_levels
        ! These fields will be passed to set wind, which maps w3 (cell centre)
        ! to w2 (cell face) vectors. However, they are actually defined in
        ! wtheta (cell top centre) and need mapping to fd1 (cell top edge).
        ! Set wind will therefore work correctly, but the indexing is shifted
        ! by half a level in the vertical for the input & output
        fd_taux(map_w3(1) + k-1) = tau_fd_x(1,1,k)
        fd_tauy(map_w3(1) + k-1) = tau_fd_y(1,1,k)
      end do
      do k=bl_levels+1,nlayers
        fd_taux(map_w3(1) + k-1) = 0.0_r_def
        fd_tauy(map_w3(1) + k-1) = 0.0_r_def
      end do
    end if
    gradrinr(map_wth(1)) = rib_gb(1,1)
    do k=2,bl_levels
      gradrinr(map_wth(1) + k-1) = BL_diag%gradrich(1,1,k)
      lmix_bl(map_wth(1) + k-1)  = BL_diag%elm3d(1,1,k)
      ngstress_bl(map_wth(1) + k-1) = f_ngstress(1,1,k)
      if (BL_diag%l_tke)                                                       &
         tke_bl(map_wth(1) + k-1) = BL_diag%tke(1,1,k)
    end do
    ! take BL weight from level 2 (first level above the surface)
    bl_weight_1dbl(map_2d(1)) = BL_diag%weight1d(1,1,2)
    do i = 1, n_land_tile
      alpha1_tile(map_tile(1)+i-1) = alpha1(1, i)
      ashtf_prime_tile(map_tile(1)+i-1) = ashtf_prime_surft(1, i)
      dtstar_tile(map_tile(1)+i-1) = dtstar_surft(1, i)
      fraca_tile(map_tile(1)+i-1) = fraca(1, i)
      z0h_tile(map_tile(1)+i-1) = fluxes%z0h_surft(1, i)
      z0m_tile(map_tile(1)+i-1) = fluxes%z0m_surft(1, i)
      rhokh_tile(map_tile(1)+i-1) = rhokh_surft(1, i)
      chr1p5m_tile(map_tile(1)+i-1) = chr1p5m(1, i)
      resfs_tile(map_tile(1)+i-1) = resfs(1, i)
      gc_tile(map_tile(1)+i-1) = gc_surft(1, i)
      canhc_tile(map_tile(1)+i-1) = canhc_surft(1, i)
    end do

    i_tile = 0
    do i = 1, n_land_tile
      do n = 1, sm_levels
        tile_water_extract(map_smtile(1)+i_tile) = wt_ext_surft(1,n,i)
        i_tile = i_tile + 1
      end do
    end do

    alpha1_tile(map_tile(1)+first_sea_tile-1) = alpha1_sea(1,1)
    ashtf_prime_tile(map_tile(1)+first_sea_tile-1) = ashtf_prime_sea(1,1)
    dtstar_tile(map_tile(1)+first_sea_tile-1) = dtstar_sea(1,1)
    rhokh_tile(map_tile(1)+first_sea_tile-1) = rhokh_sea(1,1)
    z0m_tile(map_tile(1)+first_sea_tile-1) = z0msea(1,1)

    i_sice = 0
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      alpha1_tile(map_tile(1)+i-1) = alpha1_sice(1,1,i_sice)
      ashtf_prime_tile(map_tile(1)+i-1) = ashtf_prime(1,1,i_sice)
      rhokh_tile(map_tile(1)+i-1) = rhokh_sice(1,1,i_sice)
      dtstar_tile(map_tile(1)+i-1) = dtstar_sice(1,1,i_sice)
    end do
    z0h_tile(map_tile(1)+first_sea_ice_tile-1) = z0hssi(1,1)
    z0m_tile(map_tile(1)+first_sea_ice_tile-1) = z0mssi(1,1)
    chr1p5m_tile(map_tile(1)+first_sea_ice_tile-1) = chr1p5m_sice(1,1)

    blend_height_tq(map_2d(1)) = k_blend_tq(1,1)
    z0m_eff(map_2d(1)) = z0m_eff_gb(1,1)
    ustar(map_2d(1)) = u_s(1,1)
    soil_moist_avail(map_2d(1)) = smc_soilt(1)
    zh_nonloc(map_2d(1)) = zhnl(1,1)
    z_lcl(map_2d(1)) = real(zlcl(1,1), r_def)
    inv_depth(map_2d(1)) = real(dzh(1,1), r_def)
    qcl_at_inv_top(map_2d(1)) = real(qcl_inv_top(1,1), r_def)
    if ( l_shallow(1,1) ) then
      shallow_flag(map_2d(1)) = 1_i_def
    else
      shallow_flag(map_2d(1)) = 0_i_def
    end if
    uw0_flux(map_2d(1)) = uw0(1,1)
    vw0_flux(map_2d(1)) = vw0(1,1)
    lcl_height(map_2d(1)) = zlcl_uv(1,1)
    parcel_top(map_2d(1)) = zhpar(1,1)
    level_parcel_top(map_2d(1)) = ntpar(1,1)
    wstar_2d(map_2d(1)) = wstar(1,1)
    thv_flux(map_2d(1)) = wthvs(1,1)
    parcel_buoyancy(map_2d(1)) = delthvu(1,1)
    qsat_at_lcl(map_2d(1)) = qsat_lcl(1,1)

    bl_type_ind(map_bl(1)+0) = bl_type_1(1,1)
    bl_type_ind(map_bl(1)+1) = bl_type_2(1,1)
    bl_type_ind(map_bl(1)+2) = bl_type_3(1,1)
    bl_type_ind(map_bl(1)+3) = bl_type_4(1,1)
    bl_type_ind(map_bl(1)+4) = bl_type_5(1,1)
    bl_type_ind(map_bl(1)+5) = bl_type_6(1,1)
    bl_type_ind(map_bl(1)+6) = bl_type_7(1,1)

    do n = 1, npft
      ! Unloading rate of snow from plant functional types
      snow_unload_rate(map_pft(1)+n-1) = real(jules_vars%unload_backgrnd_pft(1, n), r_def)
    end do

    do i = 1, n_land_tile
      ! Land tile temperatures
      tile_temperature(map_tile(1)+i-1) = real(tstar_surft(1, i), r_def)
      ! sensible heat flux
      tile_heat_flux(map_tile(1)+i-1) = real(fluxes%ftl_surft(1, i), r_def)
      ! moisture flux
      tile_moisture_flux(map_tile(1)+i-1) = real(fluxes%fqw_surft(1, i), r_def)
    end do

    ! Sea temperature
    tile_temperature(map_tile(1)+first_sea_tile-1) = real(tstar_sea(1,1), r_def)

    i_sice = 0
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      ! sea-ice temperature
      tile_temperature(map_tile(1)+i-1) = real(tstar_sice_ncat(1,1,i_sice), r_def)
      ! sea-ice heat flux
      tile_heat_flux(map_tile(1)+i-1) = real(fluxes%ftl_sicat(1,1,i_sice), r_def)
      ! sea-ice moisture flux
      tile_moisture_flux(map_tile(1)+i-1) = real(fluxes%fqw_sicat(1,1,i_sice), r_def)
    end do


    ! Cell roughness length for momentum
    z0m_2d(map_2d(1)) = real(sf_diag%z0m_gb(1,1), r_def)

    ! Dust fluxes
    do i = 1, ndiv
      dust_div_flux(map_dust(1) + i - 1) = real(dust_flux(1,1,i), r_def)
    end do

    ! update blended Smagorinsky diffusion coefficients only if using Smagorinsky scheme
    if ( smagorinsky ) then
      visc_m_blend(map_wth(1)) = visc_m(1,1,1)
      visc_h_blend(map_wth(1)) = visc_h(1,1,1)
      do k = 1, bl_levels-1
        visc_m_blend(map_wth(1) + k) = visc_m(1,1,k)
        visc_h_blend(map_wth(1) + k) = visc_h(1,1,k)
      end do
      visc_m_blend(map_wth(1) + bl_levels) = visc_m(1,1,bl_levels-1)
      visc_h_blend(map_wth(1) + bl_levels) = visc_h(1,1,bl_levels-1)
    endif

    ! write BL diagnostics
    if (flandg(1, 1) > 0.0_r_um) then
      gross_prim_prod(map_2d(1)) = real(sf_diag%gpp(1), r_def)
      net_prim_prod(map_2d(1)) = real(npp_gb(1), r_def)
      surface_conductance(map_2d(1)) = real(gs_gb(1), r_def)
      thermal_cond_wet_soil(map_2d(1)) = hcons_soilt(1)
      if (dim_cs1 == 4) then     ! l_triffid=TRUE (4 pool model)
        soil_respiration(map_2d(1)) = resp_s_tot_soilt(1)
      else      ! l_triffid=FALSE (1 pool model, no soil layers: dim_cslayer=1)
        soil_respiration(map_2d(1)) = trifctltype%resp_s_soilt(1,1,1,1)
      endif
    else
      gross_prim_prod(map_2d(1)) = 0.0_r_def
      net_prim_prod(map_2d(1)) = 0.0_r_def
      surface_conductance(map_2d(1)) = 0.0_r_def
      thermal_cond_wet_soil(map_2d(1)) = 0.0_r_def
      soil_respiration(map_2d(1)) = 0.0_r_def
    end if

    ! update BL prognostics
    zh_2d(map_2d(1))     = zh(1,1)
    zhsc_2d(map_2d(1))   = zhsc(1,1)
    z0msea_2d(map_2d(1)) = z0msea(1,1)
    ntml_2d(map_2d(1))   = ntml(1,1)
    if (cumulus(1,1)) then
      cumulus_2d(map_2d(1)) = 1_i_def
    else
      cumulus_2d(map_2d(1)) = 0_i_def
    endif

    if (rh_crit_opt == rh_crit_opt_tke) then
      do k = 1, nlayers
        rh_crit(map_wth(1)+k) = real(rhcpt(1,1,k), r_def)
      end do
      rh_crit(map_wth(1)) = rh_crit(map_wth(1)+1)
    end if

    ! Liquid temperature gradient for bimodal cloud scheme
    if (scheme == scheme_bimodal .or. &
         (scheme == scheme_pc2 .and. pc2ini == pc2ini_bimodal ) ) then
      do k = 1, nlayers
        dsldzm(map_wth(1)+k) = tgrad_bm(1,1,k)
      end do
    end if
    if (scheme == scheme_bimodal .or. turb_gen_mixph .or. &
         (scheme == scheme_pc2 .and. pc2ini == pc2ini_bimodal ) ) then
      do k = 2, nlayers
        wvar(map_wth(1)+k-1) = bl_w_var(1,1,k)
      end do
    end if

    if (.not. associated(zht, empty_real_data) ) then
      zht(map_2d(1)) = BL_diag%zht(1,1)
    end if
    if (.not. associated(oblen, empty_real_data) ) then
      oblen(map_2d(1)) = BL_diag%oblen(1,1)
    end if
    if (.not. associated(z0h_eff, empty_real_data) ) then
      z0h_eff(map_2d(1)) = sf_diag%z0h_eff_gb(1,1)
    end if

    ! deallocate diagnostics deallocated in atmos_physics2
    call dealloc_bl_expl(bl_diag)
    call dealloc_sf_expl(sf_diag)
    deallocate(BL_diag%gradrich)
    deallocate(BL_diag%elm3d)
    deallocate(BL_diag%tke)

    ! set this back to 1 before exit
    land_field = 1

    ! Set back to SCM for use elsewhere
    pdims_s%i_start=1
    pdims_s%i_end=row_length
    pdims_s%j_start=1
    pdims_s%j_end=rows
    pdims_l%halo_i=0
    pdims_l%halo_j=0
    deallocate(r_theta_levels)
    deallocate(r_rho_levels)
    allocate(r_theta_levels(row_length,rows,0:nlayers), source=rmdi)
    allocate(r_rho_levels(row_length,rows,nlayers), source=rmdi)

    call ancil_info_nullify(ainfo)
    call ancil_info_dealloc(ainfo_data)

    call forcing_nullify(forcing)
    call forcing_dealloc(forcing_data)

    call crop_vars_nullify(crop_vars)
    call crop_vars_dealloc(crop_vars_data)

    call lake_nullify(lake_vars)
    call lake_dealloc(lake_data)

    call coastal_nullify(coast)
    call coastal_dealloc(coastal_data)

    call trifctl_nullify(trifctltype)
    call trifctl_dealloc(trifctl_data)

    call urban_param_nullify(urban_param)
    call urban_param_dealloc(urban_param_data)

    call aero_nullify(aerotype)
    call aero_dealloc(aero_data)

    call trif_vars_nullify(trif_vars)
    call trif_vars_dealloc(trif_vars_data)

    call psparms_nullify(psparms)
    call psparms_dealloc(psparms_data)

    call jules_vars_dealloc(jules_vars_data)
    call jules_vars_nullify(jules_vars)

    call prognostics_nullify(progs)
    call prognostics_dealloc(progs_data)

    call fluxes_nullify(fluxes)
    call fluxes_dealloc(fluxes_data)

    call chemvars_nullify(chemvars)
    call chemvars_dealloc(chemvars_data)

  end subroutine bl_exp_code

end module bl_exp_kernel_mod
