!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the explicit UM boundary layer scheme.
!> @todo Doxygen comments are limited due to fparser bug - fix will be
!>       available with PSyclone 2.0.0
module bl_exp_kernel_mod

  use argument_mod,           only : arg_type,                   &
                                     GH_FIELD, GH_INTEGER,       &
                                     GH_READ, GH_WRITE, GH_INC,  &
                                     GH_READWRITE, CELLS,        &
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
                                     STENCIL, CROSS
  use constants_mod,          only : i_def, i_um, r_def, r_um, rmdi
  use fs_continuity_mod,      only : W3, Wtheta, W2
  use kernel_mod,             only : kernel_type
  use blayer_config_mod,      only : fixed_flux_e, fixed_flux_h, flux_bc_opt, &
                                     flux_bc_opt_specified_scalars, bl_mix_w
  use cloud_config_mod,       only : rh_crit_opt, rh_crit_opt_tke
  use convection_config_mod,  only : use_jules_flux
  use mixing_config_mod,      only : smagorinsky
  use surface_config_mod,     only : albedo_obs, sea_surf_alg, &
                                     sea_surf_alg_fixed_roughness, &
                                     formdrag, formdrag_dist_drag
  use timestepping_config_mod, only: outer_iterations

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: bl_exp_kernel_type
    private
    type(arg_type) :: meta_args(119) = (/                           &
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! theta_in_wth
        arg_type(GH_FIELD, GH_READ,      W3),                       &! rho_in_w3
        arg_type(GH_FIELD, GH_READ,      W3),                       &! wetrho_in_w3
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! wetrho_in_wth
        arg_type(GH_FIELD, GH_READ,      W3),                       &! exner_in_w3
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! exner_in_wth
        arg_type(GH_FIELD, GH_READ,      W3, STENCIL(CROSS)),       &! u1_in_w3
        arg_type(GH_FIELD, GH_READ,      W3, STENCIL(CROSS)),       &! u2_in_w3
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! u3_in_wth
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! m_v_n
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! m_cl_n
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! m_ci_n
        arg_type(GH_FIELD, GH_READ,      W3),                       &! height_w3
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! height_wth
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! shear
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! delta
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! max_diff_smag
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! zh_2d
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! z0msea_2d
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! ntml_2d
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cumulus_2d
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_2, STENCIL(CROSS)),&! tile_fraction
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_3),&! leaf_area_index
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_3),&! canopy_height
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sd_orog_2d
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! peak_to_trough_orog
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! silhouette_area_orog
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_albedo
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_roughness
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_wilt
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_crit
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_sat
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_thermal_cond
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_suction_sat
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! clapp_horn_b
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_carbon_content
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! soil_respiration
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! thermal_cond_wet_soil
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_temperature
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! tile_temperature
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_snow_mass
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! n_snow_layers
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! snow_depth
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! snow_layer_thickness
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! snow_layer_ice_mass
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! snow_layer_liq_mass
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! snow_layer_temp
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! surface_conductance
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! canopy_water
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! soil_temperature
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! soil_moisture
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! unfrozen_soil_moisture
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! frozen_soil_moisture
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! tile_heat_flux
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! tile_moisture_flux
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! gross_prim_prod
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! net_prim_prod
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! cos_zen_angle
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! sw_up_tile
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sw_down_surf
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! lw_down_surf
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sw_down_surf_blue
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! dd_mf_cb
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! dtl_mphys
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! dmt_mphys
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! sw_heating_rate
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! lw_heating_rate
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! ozone
        arg_type(GH_FIELD, GH_READ,      WTHETA),                   &! cf_bulk
        arg_type(GH_FIELD, GH_READWRITE, WTHETA),                   &! rh_crit
        arg_type(GH_FIELD, GH_WRITE,     WTHETA),                   &! visc_m_blend
        arg_type(GH_FIELD, GH_WRITE,     WTHETA),                   &! visc_h_blend
        arg_type(GH_FIELD, GH_INC,       W2),                       &! du_bl
        arg_type(GH_FIELD, GH_WRITE,     WTHETA),                   &! rhokm_bl
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_7),&! rhokm_surf
        arg_type(GH_FIELD, GH_WRITE,     W3),                       &! rhokh_bl
        arg_type(GH_FIELD, GH_WRITE,     WTHETA),                   &! ngstress_bl
        arg_type(GH_FIELD, GH_WRITE,     WTHETA),                   &! bq_bl
        arg_type(GH_FIELD, GH_WRITE,     WTHETA),                   &! bt_bl
        arg_type(GH_FIELD, GH_WRITE,     W3),                       &! moist_flux_bl
        arg_type(GH_FIELD, GH_WRITE,     W3),                       &! heat_flux_bl
        arg_type(GH_FIELD, GH_WRITE,     W3),                       &! dtrdz_uv_bl
        arg_type(GH_FIELD, GH_WRITE,     WTHETA),                   &! dtrdz_tq_bl
        arg_type(GH_FIELD, GH_WRITE,     W3),                       &! rdz_tq_bl
        arg_type(GH_FIELD, GH_WRITE,     WTHETA),                   &! rdz_uv_bl
        arg_type(GH_FIELD, GH_WRITE,     WTHETA),                   &! fd_taux
        arg_type(GH_FIELD, GH_WRITE,     WTHETA),                   &! fd_tauy
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! alpha1_tile
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! ashtf_prime_tile
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! dtstar_tile
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! fraca_tile
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! z0h_tile
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! z0m_tile
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! rhokh_tile
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! chr1p5m_tile
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! resfs_tile
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! canhc_tile
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_8),&! tile_water_extract
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! blend_height_tq
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! blend_height_uv
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! ustar
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_avail
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! zh_nonloc
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! z_lcl
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! inv_depth
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! qcl_at_inv_top
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! shallow_flag
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! uw0_flux
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! vw0_flux
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! lcl_height
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! parcel_top
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! level_parcel_top
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! wstar_2d
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! thv_flux
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! parcel_buoyancy
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! qsat_at_lcl
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_9),&! bl_type_ind
        arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_3),&! snow_unload_rate
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_10)&! albedo_obs_scaling
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass :: bl_exp_code
  end type

  public bl_exp_code

contains

  !> @brief Interface to the UM BL scheme
  !> @details The UM Boundary Layer scheme does:
  !>             vertical mixing of heat, momentum and moisture,
  !>             as documented in UMDP24
  !>          NB This version uses winds in w3 space (i.e. A-grid)
  !> @param[in]     nlayers              Number of layers
  !> @param[in]     theta_in_wth         Potential temperature field
  !> @param[in]     rho_in_w3            Density field in density space
  !> @param[in]     wetrho_in_w3         Wet density field in density space
  !> @param[in]     wetrho_in_wth        Wet density field in wth space
  !> @param[in]     exner_in_w3          Exner pressure field in density space
  !> @param[in]     exner_in_wth         Exner pressure field in wth space
  !> @param[in]     u1_in_w3             'Zonal' wind in density space
  !> @param[in]     u2_in_w3             'Meridional' wind in density space
  !> @param[in]     u3_in_wth            'Vertical' wind in theta space
  !> @param[in]     m_v_n                Vapour mixing ratio at time level n
  !> @param[in]     m_cl_n               Cloud liq mixing ratio at time level n
  !> @param[in]     m_ci_n               Cloud ice mixing ratio at time level n
  !> @param[in]     height_w3            Height of density space above surface
  !> @param[in]     height_wth           Height of theta space above surface
  !> @param[in]     shear                3D wind shear on wtheta points
  !> @param[in]     delta                Edge length on wtheta points
  !> @param[in]     max_diff_smag        Maximum diffusion coefficient allowed
  !> @param[in,out] zh_2d                Boundary layer depth
  !> @param[in,out] z0msea_2d            Roughness length
  !> @param[out]    ntml_2d              Number of turbulently mixed levels
  !> @param[out]    cumulus_2d           Cumulus flag (true/false)
  !> @param[in]     tile_fraction        Surface tile fractions
  !> @param[in]     leaf_area_index      Leaf Area Index
  !> @param[in]     canopy_height        Canopy height
  !> @param[in]     sd_orog_2d           Standard deviation of orography
  !> @param[in]     peak_to_trough_orog  Half of peak-to-trough height over root(2) of orography
  !> @param[in]     silhouette_area_orog Silhouette area of orography
  !> @param[in]     soil_albedo          Snow-free soil albedo
  !> @param[in]     soil_roughness       Bare soil surface roughness length
  !> @param[in]     soil_moist_wilt      Volumetric soil moisture at wilting point
  !> @param[in]     soil_moist_crit      Volumetric soil moisture at critical point
  !> @param[in]     soil_moist_sat       Volumetric soil moisture at saturation
  !> @param[in]     soil_thermal_cond    Soil thermal conductivity
  !> @param[in]     soil_suction_sat     Saturated soil water suction
  !> @param[in]     clapp_horn_b         Clapp and Hornberger b coefficient
  !> @param[in]     soil_carbon_content  Soil carbon content
  !> @param[out]    soil_respiration     Soil respiration  (kg m-2 s-1)
  !> @param[out]    thermal_cond_wet_soil Thermal conductivity of wet soil (W m-1 K-1)
  !> @param[in]     sea_ice_temperature  Bulk temperature of sea-ice (K)
  !> @param[in,out] tile_temperature     Surface tile temperatures
  !> @param[in]     tile_snow_mass       Snow mass on tiles (kg/m2)
  !> @param[in]     n_snow_layers        Number of snow layers on tiles
  !> @param[in]     snow_depth           Snow depth on tiles
  !> @param[in]     snow_layer_thickness Thickness of snow layers (m)
  !> @param[in]     snow_layer_ice_mass  Mass of ice in snow layers (kg m-2)
  !> @param[in]     snow_layer_liq_mass  Mass of liquid in snow layers (kg m-2)
  !> @param[in]     snow_layer_temp      Temperature of snow layer (K)
  !> @param[in,out] surface_conductance  Surface conductance
  !> @param[in]     canopy_water         Canopy water on each tile
  !> @param[in]     soil_temperature     Soil temperature
  !> @param[in]     soil_moisture        Soil moisture content (kg m-2)
  !> @param[in]     unfrozen_soil_moisture Unfrozen soil moisture proportion
  !> @param[in]     frozen_soil_moisture Frozen soil moisture proportion
  !> @param[out]    tile_heat_flux       Surface heat flux
  !> @param[out]    tile_moisture_flux   Surface moisture flux
  !> @param[out]    gross_prim_prod      Gross Primary Productivity
  !> @param[out]    net_prim_prod        Net Primary Productivity
  !> @param[in]     cos_zen_angle        Cosing of solar zenith angle
  !> @param[in]     sw_up_tile           Upwelling SW radiation on surface tiles
  !> @param[in]     sw_down_surf         Downwelling SW radiation at surface
  !> @param[in]     lw_down_surf         Downwelling LW radiation at surface
  !> @param[in]     sw_down_surf_blue    Photosynthetically active SW down
  !> @param[in]     dd_mf_cb             Downdraft massflux at cloud base (Pa/s)
  !> @param[in]     dtl_mphys            Microphysics liq temperature increment
  !> @param[in]     dmt_mphys            Microphysics total water increment
  !> @param[in]     sw_heating_rate      Shortwave radiation heating rate
  !> @param[in]     lw_heating_rate      Longwave radiation heating rate
  !> @param[in]     ozone                Ozone field
  !> @param[in]     cf_bulk              Bulk cloud fraction
  !> @param[in,out] rh_crit              Critical rel humidity
  !> @param[out]    visc_m_blend         Blended BL-Smag diffusion coefficient for momentum
  !> @param[out]    visc_h_blend         Blended BL-Smag diffusion coefficient for scalars
  !> @param[in,out] du_bl                Wind increment from BL scheme
  !> @param[out]    rhokm_bl             Momentum eddy diffusivity on BL levels
  !> @param[out]    rhokm_surf           Momentum eddy diffusivity for coastal tiling
  !> @param[out]    rhokh_bl             Heat eddy diffusivity on BL levels
  !> @param[out]    ngstress_bl          Non-gradient stress function on BL levels
  !> @param[out]    bq_bl                Buoyancy parameter for moisture
  !> @param[out]    bt_bl                Buoyancy parameter for heat
  !> @param[out]    moist_flux_bl        Vertical moisture flux on BL levels
  !> @param[out]    heat_flux_bl         Vertical heat flux on BL levels
  !> @param[out]    dtrdz_uv_bl          dt/(rho*r*r*dz) in w3 space
  !> @param[out]    dtrdz_tq_bl          dt/(rho*r*r*dz) in wth
  !> @param[out]    rdz_tq_bl            1/dz in w3
  !> @param[out]    rdz_uv_bl            1/dz in wth space
  !> @param[out]    alpha1_tile          dqsat/dT in surface layer on tiles
  !> @param[out]    ashtf_prime_tile     Heat flux coefficient on tiles
  !> @param[out]    dtstar_tile          Change in surface temperature on tiles
  !> @param[out]    fraca_tile           Fraction of moisture flux with only aerodynamic resistance
  !> @param[out]    z0h_tile             Heat roughness length on tiles
  !> @param[out]    z0m_tile             Momentum roughness length on tiles
  !> @param[out]    rhokh_tile           Surface heat diffusivity on tiles
  !> @param[out]    chr1p5m_tile         1.5m transfer coefficients on tiles
  !> @param[out]    resfs_tile           Combined aerodynamic resistance
  !> @param[out]    canhc_tile           Canopy heat capacity on tiles
  !> @param[out]    tile_water_extract   Extraction of water from each tile
  !> @param[out]    blend_height_tq      Blending height for wth levels
  !> @param[out]    blend_height_uv      Blending height for w3 levels
  !> @param[out]    ustar                Friction velocity
  !> @param[out]    soil_moist_avail     Available soil moisture for evaporation
  !> @param[out]    zh_nonloc            Depth of non-local BL scheme
  !> @param[out]    z_lcl                Height of the LCL (wtheta levels)
  !> @param[out]    inv_depth            Depth of BL top inversion layer
  !> @param[out]    qcl_at_inv_top       Cloud water at top of inversion
  !> @param[out]    shallow_flag         Indicator of shallow convection
  !> @param[out]    uw0_flux             'Zonal' surface momentum flux
  !> @param[out]    vw0 flux             'Meridional' surface momentum flux
  !> @param[out]    lcl_height           Height of lifting condensation level (w3 levels)
  !> @param[out]    parcel_top           Height of surface based parcel ascent
  !> @param[out]    level_parcel_top     Model level of parcel_top
  !> @param[out]    wstar_2d             BL velocity scale
  !> @param[out]    thv_flux             Surface flux of theta_v
  !> @param[out]    parcel_buoyancy      Integral of parcel buoyancy
  !> @param[out]    qsat_at_lcl          Saturation specific hum at LCL
  !> @param[out]    bl_type_ind          Diagnosed BL types
  !> @param[out]    snow_unload_rate     Unloading of snow from PFTs by wind
  !> @param[in]     albedo_obs_scaling   Scaling factor to adjust albedos by
  !> @param[in]     ndf_wth              Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth             Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth              Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_w3               Number of DOFs per cell for density space
  !> @param[in]     undf_w3              Number of unique DOFs for density space
  !> @param[in]     map_w3               Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_2d               Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d              Number of unique DOFs for 2D fields
  !> @param[in]     map_2d               Dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_tile             Number of DOFs per cell for tiles
  !> @param[in]     undf_tile            Number of total DOFs for tiles
  !> @param[in]     map_tile             Dofmap for cell for surface tiles
  !> @param[in]     ndf_pft              Number of DOFs per cell for PFTs
  !> @param[in]     undf_pft             Number of total DOFs for PFTs
  !> @param[in]     map_pft              Dofmap for cell for PFTs
  !> @param[in]     ndf_sice             Number of DOFs per cell for sice
  !> @param[in]     undf_sice            Number of total DOFs for sice
  !> @param[in]     map_sice             Dofmap for cell for sice
  !> @param[in]     ndf_snow             Number of DOFs per cell for snow
  !> @param[in]     undf_snow            Number of total DOFs for snow
  !> @param[in]     map_snow             Dofmap for cell for snow
  !> @param[in]     ndf_soil             Number of DOFs per cell for soil levels
  !> @param[in]     undf_soil            Number of total DOFs for soil levels
  !> @param[in]     map_soil             Dofmap for cell for soil levels
  !> @param[in]     ndf_surf             Number of DOFs per cell for surface variables
  !> @param[in]     undf_surf            Number of unique DOFs for surface variables
  !> @param[in]     map_surf             Dofmap for the cell at the base of the column for surface variables
  !> @param[in]     ndf_smtile           Number of DOFs per cell for soil levels and tiles
  !> @param[in]     undf_smtile          Number of total DOFs for soil levels and tiles
  !> @param[in]     map_smtile           Dofmap for cell for soil levels and tiles
  !> @param[in]     ndf_bl               Number of DOFs per cell for BL types
  !> @param[in]     undf_bl              Number of total DOFs for BL types
  !> @param[in]     map_bl               Dofmap for cell for BL types
  !> @param[in]     ndf_scal             Number of DOFs per cell for albedo scaling
  !> @param[in]     undf_scal            Number of total DOFs for albedo scaling
  !> @param[in]     map_scal             Dofmap for cell at the base of the column
  subroutine bl_exp_code(nlayers,                               &
                         theta_in_wth,                          &
                         rho_in_w3,                             &
                         wetrho_in_w3,                          &
                         wetrho_in_wth,                         &
                         exner_in_w3,                           &
                         exner_in_wth,                          &
                         u1_in_w3,                              &
                         u1_w3_stencil_size, u1_w3_stencil,     &
                         u2_in_w3,                              &
                         u2_w3_stencil_size, u2_w3_stencil,     &
                         u3_in_wth,                             &
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
                         soil_carbon_content,                   &
                         soil_respiration,                      &
                         thermal_cond_wet_soil,                 &
                         sea_ice_temperature,                   &
                         tile_temperature,                      &
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
                         rh_crit,                               &
                         visc_m_blend,                          &
                         visc_h_blend,                          &
                         du_bl,                                 &
                         rhokm_bl,                              &
                         rhokm_surf,                            &
                         rhokh_bl,                              &
                         ngstress_bl,                           &
                         bq_bl,                                 &
                         bt_bl,                                 &
                         moist_flux_bl,                         &
                         heat_flux_bl,                          &
                         dtrdz_uv_bl,                           &
                         dtrdz_tq_bl,                           &
                         rdz_tq_bl,                             &
                         rdz_uv_bl,                             &
                         fd_taux,                               &
                         fd_tauy,                               &
                         alpha1_tile,                           &
                         ashtf_prime_tile,                      &
                         dtstar_tile,                           &
                         fraca_tile,                            &
                         z0h_tile,                              &
                         z0m_tile,                              &
                         rhokh_tile,                            &
                         chr1p5m_tile,                          &
                         resfs_tile,                            &
                         canhc_tile,                            &
                         tile_water_extract,                    &
                         blend_height_tq,                       &
                         blend_height_uv,                       &
                         ustar,                                 &
                         soil_moist_avail,                      &
                         zh_nonloc,                             &
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
                         bl_type_ind,                           &
                         snow_unload_rate,                      &
                         albedo_obs_scaling,                    &
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
                         ndf_w2, undf_w2, map_w2,               &
                         ndf_surf, undf_surf, map_surf,         &
                         ndf_smtile, undf_smtile, map_smtile,   &
                         ndf_bl, undf_bl, map_bl,               &
                         ndf_scal, undf_scal, map_scal)

    !---------------------------------------
    ! LFRic modules
    !---------------------------------------
    use jules_control_init_mod, only: n_land_tile, n_sea_ice_tile, &
         first_sea_tile, first_sea_ice_tile

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use ancil_info, only: ssi_pts, sea_pts, sice_pts, sice_pts_ncat, rad_nband
    use atm_fields_bounds_mod, only: tdims, pdims_s, pdims_l
    use atm_step_local, only: dim_cs1, dim_cs2, co2_dim_len, co2_dim_row
    use c_kappai, only: kappai, de
    use cv_run_mod, only: i_convection_vn, i_convection_vn_6a,               &
                          cldbase_opt_dp, cldbase_opt_md
    use dust_parameters_mod, only: ndiv, ndivh
    use jules_sea_seaice_mod, only: nice_use
    use jules_snow_mod, only: nsmax
    use jules_surface_types_mod, only: npft, ntype
    use nlsizes_namelist_mod, only: row_length, rows, land_field,            &
         sm_levels, ntiles, bl_levels
    use planet_constants_mod, only: p_zero, kappa, planet_radius
    use rad_input_mod, only: co2_mmr
    use timestep_mod, only: timestep

    ! spatially varying fields used from modules
    use fluxes, only: sw_sea, sw_sicat
    use level_heights_mod, only: r_theta_levels, r_rho_levels
    use ozone_vars, only: o3_gb
    use turb_diff_ctl_mod, only: visc_m, visc_h, max_diff, delta_smag

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
    use jules_fields_mod, only : crop_vars, psparms, ainfo, trif_vars,         &
                                 aerotype, urban_param, progs, trifctltype,    &
                                 coast, jules_vars

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
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
    integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
    integer(kind=i_def), intent(in) :: map_w2(ndf_w2)
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

    integer(kind=i_def), intent(in) :: u1_w3_stencil_size, u2_w3_stencil_size
    integer(kind=i_def), dimension(ndf_w3,u1_w3_stencil_size), intent(in) :: u1_w3_stencil
    integer(kind=i_def), dimension(ndf_w3,u2_w3_stencil_size), intent(in) :: u2_w3_stencil

    integer(kind=i_def), intent(in) :: tile_stencil_size
    integer(kind=i_def), dimension(ndf_tile,tile_stencil_size), intent(in) :: tile_stencil

    real(kind=r_def), dimension(undf_wth), intent(inout):: rh_crit
    real(kind=r_def), dimension(undf_w2),  intent(inout):: du_bl

    real(kind=r_def), dimension(undf_wth), intent(inout):: visc_h_blend,       &
                                                           visc_m_blend,       &
                                                           rhokm_bl,           &
                                                           ngstress_bl,        &
                                                           bq_bl, bt_bl,       &
                                                           dtrdz_tq_bl,        &
                                                           rdz_uv_bl,          &
                                                           fd_taux, fd_tauy
    real(kind=r_def), dimension(undf_w3),  intent(inout):: rhokh_bl,           &
                                                           moist_flux_bl,      &
                                                           heat_flux_bl,       &
                                                           dtrdz_uv_bl,        &
                                                           rdz_tq_bl
    real(kind=r_def), dimension(undf_w3),  intent(in)   :: rho_in_w3,          &
                                                           wetrho_in_w3,       &
                                                           exner_in_w3,        &
                                                           u1_in_w3, u2_in_w3, &
                                                           height_w3
    real(kind=r_def), dimension(undf_wth), intent(in)   :: theta_in_wth,       &
                                                           wetrho_in_wth,      &
                                                           exner_in_wth,       &
                                                           u3_in_wth,          &
                                                           m_v_n, m_cl_n,      &
                                                           m_ci_n,             &
                                                           height_wth,         &
                                                           shear,              &
                                                           delta,              &
                                                           max_diff_smag,      &
                                                           dtl_mphys,dmt_mphys,&
                                                           sw_heating_rate,    &
                                                           lw_heating_rate,    &
                                                           cf_bulk, ozone

    real(kind=r_def), dimension(undf_2d), intent(inout) :: zh_2d,              &
                                                           z0msea_2d

    real(kind=r_def), dimension(undf_2d), intent(inout) :: ntml_2d,            &
                                                           cumulus_2d,         &
                                                           blend_height_tq,    &
                                                           blend_height_uv,    &
                                                           ustar,              &
                                                           soil_moist_avail,   &
                                                           zh_nonloc,          &
                                                           z_lcl,              &
                                                           inv_depth,          &
                                                           qcl_at_inv_top,     &
                                                           shallow_flag,       &
                                                           uw0_flux,           &
                                                           vw0_flux,           &
                                                           lcl_height,         &
                                                           parcel_top,         &
                                                           level_parcel_top,   &
                                                           wstar_2d,           &
                                                           thv_flux,           &
                                                           parcel_buoyancy,    &
                                                           qsat_at_lcl

    real(kind=r_def), intent(in) :: tile_fraction(undf_tile)
    real(kind=r_def), intent(inout) :: tile_temperature(undf_tile)
    real(kind=r_def), intent(in) :: tile_snow_mass(undf_tile)
    real(kind=r_def), intent(in) :: n_snow_layers(undf_tile)
    real(kind=r_def), intent(in) :: snow_depth(undf_tile)
    real(kind=r_def), intent(in) :: canopy_water(undf_tile)
    real(kind=r_def), intent(inout) :: tile_heat_flux(undf_tile)
    real(kind=r_def), intent(inout) :: tile_moisture_flux(undf_tile)
    real(kind=r_def), intent(in) :: sw_up_tile(undf_tile)
    real(kind=r_def), intent(in) :: albedo_obs_scaling(undf_scal)

    real(kind=r_def), intent(in) :: leaf_area_index(undf_pft)
    real(kind=r_def), intent(in) :: canopy_height(undf_pft)

    real(kind=r_def), intent(in) :: sea_ice_temperature(undf_sice)

    real(kind=r_def), intent(in) :: sd_orog_2d(undf_2d)
    real(kind=r_def), intent(in) :: peak_to_trough_orog(undf_2d)
    real(kind=r_def), intent(in) :: silhouette_area_orog(undf_2d)
    real(kind=r_def), intent(in) :: soil_albedo(undf_2d)
    real(kind=r_def), intent(in) :: soil_roughness(undf_2d)
    real(kind=r_def), intent(in) :: soil_thermal_cond(undf_2d)
    real(kind=r_def), intent(in) :: soil_carbon_content(undf_2d)
    real(kind=r_def), intent(inout) :: surface_conductance(undf_2d)
    real(kind=r_def), intent(in) :: cos_zen_angle(undf_2d)
    real(kind=r_def), intent(in) :: sw_down_surf(undf_2d)
    real(kind=r_def), intent(in) :: lw_down_surf(undf_2d)
    real(kind=r_def), intent(in) :: sw_down_surf_blue(undf_2d)
    real(kind=r_def), intent(in) :: dd_mf_cb(undf_2d)
    real(kind=r_def), intent(inout) :: gross_prim_prod(undf_2d)
    real(kind=r_def), intent(inout) :: net_prim_prod(undf_2d)
    real(kind=r_def), intent(inout) :: soil_respiration(undf_2d)
    real(kind=r_def), intent(inout) :: thermal_cond_wet_soil(undf_2d)

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

    real(kind=r_def), dimension(undf_bl),   intent(inout)  :: bl_type_ind
    real(kind=r_def), dimension(undf_pft),  intent(inout)  :: snow_unload_rate
    real(kind=r_def), dimension(undf_surf), intent(inout)  :: rhokm_surf
    real(kind=r_def), dimension(undf_tile), intent(inout):: alpha1_tile,      &
                                                            ashtf_prime_tile, &
                                                            dtstar_tile,      &
                                                            fraca_tile,       &
                                                            z0h_tile,         &
                                                            z0m_tile,         &
                                                            rhokh_tile,       &
                                                            chr1p5m_tile,     &
                                                            resfs_tile,       &
                                                            canhc_tile

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
         rho_wet, rho_dry, z_rho, z_theta, bulk_cloud_fraction, rho_wet_tq,  &
         u_p, v_p, zeros, rhcpt

    real(r_um), dimension(0:row_length+1,0:rows+1,nlayers) :: rho_wet_rsq,   &
         p_rho_levels, u_px, v_px, exner_rho_levels

    ! profile field on boundary layer levels
    real(r_um), dimension(row_length,rows,bl_levels) ::                      &
         fqw, ftl, rhokh, bq_gb, bt_gb, dtrdz_charney_grid,                  &
         dtrdz_u, rdz_charney_grid, rhokm_mix, w_mixed, w_flux

    real(r_um), dimension(bl_levels) :: gamma

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
         tstar_sice_ncat, fqw_ice, ftl_ice, ice_fract_ncat, alpha1_sice,     &
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

    ! fields on land points and carbon pools
    real(r_um), dimension(land_field,dim_cs1) :: cs_pool_gb_um, resp_s_gb_um

    ! fields on land points and tiles
    real(r_um), dimension(land_field,ntiles) :: canopy_surft, catch_surft,   &
         catch_snow_surft, snow_surft, z0_surft, z0h_bare_surft, sw_surft,   &
         tstar_surft, frac_surft, ftl_surft, fqw_surft, dtstar_surft,        &
         alpha1, ashtf_prime_surft, chr1p5m, fraca, resfs, rhokh_surft,      &
         z0h_surft, z0m_surft, canhc_surft

    ! fields on land points and pfts
    real(r_um), dimension(land_field,npft) :: canht_pft, lai_pft

    ! field on surface tiles and soil levels
    real(r_um), dimension(land_field,sm_levels,ntiles) :: wt_ext_surft

    ! Fields which are not used and only required for subroutine argument list,
    ! hence are unset in the kernel
    ! if they become set, please move up to be with other variables
    integer(i_um) :: asteps_since_triffid
    integer(i_um) :: curr_year, curr_day_number, curr_hour, curr_minute,     &
                     curr_second
    integer(i_um), parameter :: nscmdpkgs=15
    logical,       parameter :: l_scmdiags(nscmdpkgs)=.false.

    real(r_um), dimension(row_length,rows,nlayers) :: bl_w_var, tnuc_new

    real(r_um), dimension(row_length,rows,bl_levels) ::                      &
         e_trb, tsq_trb, qsq_trb, cov_trb, dtrdz_v

    real(r_um), dimension(row_length,rows,2:bl_levels) :: rdz_v

    real(r_um), dimension(row_length,rows,0:bl_levels-1) :: taux_p, tauy_p

    real(r_um), dimension(row_length,rows,0:nlayers) :: conv_prog_precip

    real(r_um), dimension(row_length,rows) ::                                &
         z0h_scm, z0m_scm, soil_clay, soil_sand, dust_mrel1,                 &
         dust_mrel2, dust_mrel3, dust_mrel4, dust_mrel5, dust_mrel6,         &
         zhpar_shcu, t1_sd, q1_sd, qcl_inv_top, w_max, deep_flag,            &
         past_precip, past_conv_ht, ql_ad, cin_undilute, cape_undilute,      &
         entrain_coef, ustar_in, g_ccp, h_ccp,  ccp_strength, fb_surf,       &
         charnock_w, aresist, cu_over_orog, resist_b, rho_aresist, rib_gb,   &
         shallowc, vshr, z0m_eff_gb, zhsc, tnuc_nlcl

    integer(i_um), dimension(row_length,rows) :: nlcl, conv_type, nbdsc, ntdsc

    logical, dimension(row_length,rows) :: no_cumulus, l_congestus,          &
                                           l_congestus2, l_mid

    real(r_um), dimension(row_length,rows,nice_use) :: radnet_sice

    real(r_um), dimension(row_length,rows,ndiv) :: dust_flux, r_b_dust

    real(r_um), dimension(land_field) :: emis_soil

    real(r_um), dimension(dim_cs2) :: resp_s_tot_soilt

    real(r_um), dimension(land_field,dim_cs2) :: resp_s_acc_gb_um

    real(r_um), dimension(land_field,ntiles) :: tsurf_elev_surft, epot_surft,&
         aresist_surft, dust_emiss_frac, emis_surft, flake, gc_surft, resft, &
         resist_b_surft, rho_aresist_surft, tile_frac, u_s_std_surft

    real(r_um), dimension(land_field,ntiles,ndivh) :: u_s_t_dry_tile,        &
         u_s_t_tile

    real(r_um), dimension(land_field,npft) :: resp_w_pft, g_leaf_acc_pft,    &
         npp_acc_pft, resp_w_acc_pft

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

    ! surf_hgt_surft needs to be initialised in this kernel so it has the
    ! correct value for this grid-point when used within Jules
    jules_vars%surf_hgt_surft = 0.0_r_um

    zeros=0.0_r_um

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

    !-----------------------------------------------------------------------
    ! Mapping of LFRic fields into UM variables
    !-----------------------------------------------------------------------

    ! Land tile fractions
    flandg = 0.0_r_um
    do i = 1, n_land_tile
      flandg(1,1) = flandg(1,1) + real(tile_fraction(tile_stencil(1,1)+i-1), r_um)
      flandg(0,1) = flandg(0,1) + real(tile_fraction(tile_stencil(1,2)+i-1), r_um)
      flandg(1,0) = flandg(1,0) + real(tile_fraction(tile_stencil(1,3)+i-1), r_um)
      flandg(2,1) = flandg(2,1) + real(tile_fraction(tile_stencil(1,4)+i-1), r_um)
      flandg(1,2) = flandg(1,2) + real(tile_fraction(tile_stencil(1,5)+i-1), r_um)
      frac_surft(1, i) = real(tile_fraction(tile_stencil(1,1)+i-1), r_um)
    end do
    ! Set these to land to exclude them until full stencil is available
    flandg(0,0) = 1.0_r_um
    flandg(0,2) = 1.0_r_um
    flandg(2,0) = 1.0_r_um
    flandg(2,2) = 1.0_r_um
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
    if (flandg(1,1) < 1.0_r_um .and. &
         tile_fraction(tile_stencil(1,1)+first_sea_tile-1) == 0.0_r_def .and. &
         ice_fract(1,1) == 0.0_r_um) then
      flandg(1,1) = 1.0_r_um
    end if

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

    ! Sea temperature
    tstar_sea = 0.0_r_um
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
        k_sice_ncat(1, 1, i) = 2.0_r_um * kappai / de
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

    ! Soil carbon content (cs_pool_gb_um)
    cs_pool_gb_um = real(soil_carbon_content(map_2d(1)), r_um)

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

    ! Downwelling LW radiation at surface
    lw_down = real(lw_down_surf(map_2d(1)), r_um)

    ! Net SW radiation on tiles
    do i = 1, n_land_tile
      sw_surft(1, i) = real(sw_down_surf(map_2d(1)) - &
                            sw_up_tile(map_tile(1)+i-1), r_um)
    end do

    ! Net SW on open sea
    sw_sea(1) = real(sw_down_surf(map_2d(1)) - &
                     sw_up_tile(map_tile(1)+first_sea_tile-1), r_um)

    ! Net SW on sea-ice
    i_sice = 0
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      sw_sicat(1, i_sice) = real(sw_down_surf(map_2d(1)) - &
                                 sw_up_tile(map_tile(1)+i-1), r_um)
    end do

    ! photosynthetically active downwelling SW radiation
    photosynth_act_rad = real(sw_down_surf_blue(map_2d(1)), r_um)

    ! Ozone
    o3_gb = real(ozone(map_wth(1)), r_um)

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
      progs%nsnow_surft(1, i) = int(n_snow_layers(map_tile(1)+i-1), i_um)
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
      rho_wet(1,1,k) = wetrho_in_w3(map_w3(1) + k-1)
      ! dry density on rho levels
      rho_dry(1,1,k) = rho_in_w3(map_w3(1) + k-1)
      ! pressure on rho and theta levels
      p_rho_levels(1,1,k) = p_zero*(exner_in_w3(map_w3(1) + k-1))**(1.0_r_def/kappa)
      p_theta_levels(1,1,k) = p_zero*(exner_in_wth(map_wth(1) + k))**(1.0_r_def/kappa)
      ! exner pressure on rho and theta levels
      exner_rho_levels(1,1,k) = exner_in_w3(map_w3(1) + k-1)
      exner_theta_levels(1,1,k) = exner_in_wth(map_wth(1) + k)
      ! u wind on rho levels
      u_px(1,1,k) = u1_in_w3(u1_w3_stencil(1,1) + k-1)
      u_px(0,1,k) = u1_in_w3(u1_w3_stencil(1,2) + k-1)
      u_px(1,0,k) = u1_in_w3(u1_w3_stencil(1,3) + k-1)
      u_px(2,1,k) = u1_in_w3(u1_w3_stencil(1,4) + k-1)
      u_px(1,2,k) = u1_in_w3(u1_w3_stencil(1,5) + k-1)
      ! v wind on rho levels
      v_px(1,1,k) = u2_in_w3(u2_w3_stencil(1,1) + k-1)
      v_px(0,1,k) = u2_in_w3(u2_w3_stencil(1,2) + k-1)
      v_px(1,0,k) = u2_in_w3(u2_w3_stencil(1,3) + k-1)
      v_px(2,1,k) = u2_in_w3(u2_w3_stencil(1,4) + k-1)
      v_px(1,2,k) = u2_in_w3(u2_w3_stencil(1,5) + k-1)
      ! w wind on theta levels
      w(1,1,k) = u3_in_wth(map_wth(1) + k)
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
    ! wet density multiplied by planet radius squared on rho levs
    rho_wet_rsq(1,1,:) = rho_wet(1,1,:) * r_rho_levels(1,1,:)**2
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
    w(1,1,0) = u3_in_wth(map_wth(1) + 0)
    etadot = w / z_theta(1,1,nlayers)

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
         curr_year, curr_day_number, curr_hour, curr_minute, curr_second,      &
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
         z0m_soil_gb, lw_down, sw_surft, tstar_surft, tsurf_elev_surft,        &
         co2,                                                                  &
         asteps_since_triffid,                                                 &
         cs_pool_gb_um,frac_surft,canht_pft,lai_pft,fland,flandg,              &
         albsoil_soilt, cos_zenith_angle,                                      &
    !     IN: input from the wave model
         charnock_w,                                                           &
    !     IN everything not covered so far
         t_soil_soilt, ti_sice,                                                &
         ti_sice_ncat,tstar,zh_prev,ddmfx,bulk_cloud_fraction,zhpar,zlcl,      &
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
         ftl_surft,radnet_sice,rho_aresist_surft,                              &
         aresist_surft, resist_b_surft, alpha1, ashtf_prime_surft,             &
         fqw_surft, epot_surft,                                                &
         fqw_ice,ftl_ice,fraca,resfs,resft,rhokh_surft,rhokh_sice,rhokh_sea,   &
         z0hssi,z0h_surft,z0mssi,z0m_surft,chr1p5m,chr1p5m_sice,smc_soilt,     &
         npp_gb, resp_s_gb_um, resp_s_tot_soilt,                               &
         resp_w_pft, gc_surft, canhc_surft, wt_ext_surft, flake,               &
         surft_index, surft_pts,                                               &
         tile_frac, tstar_land, tstar_ssi, dtstar_surft,                       &
         dtstar_sea, dtstar_sice, hcons_soilt, emis_surft, emis_soil,          &
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
         nbdsc,ntdsc,wstar,wthvs,uw0,vw0,taux_p,tauy_p,rhcpt,zeros             &
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
          , w, w_max, deep_flag, past_precip, past_conv_ht              &
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
         curr_year, curr_day_number, curr_hour, curr_minute, curr_second,&
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
         z0m_soil_gb, lw_down, sw_surft, tstar_surft, tsurf_elev_surft, &
         co2,                                                           &
         asteps_since_triffid,                                          &
         cs_pool_gb_um,frac_surft,canht_pft,lai_pft,fland,flandg,       &
         albsoil_soilt, cos_zenith_angle,                               &
    !     IN: input from the wave model
         charnock_w,                                                    &
    !     IN everything not covered so far
         t_soil_soilt, ti_sice,                                           &
         ti_sice_ncat,tstar,zh_prev,ddmfx,bulk_cloud_fraction,zhpar,zlcl, &
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
        ftl_surft,radnet_sice,rho_aresist_surft,                        &
        aresist_surft, resist_b_surft, alpha1, ashtf_prime_surft,       &
        fqw_surft, epot_surft,                                          &
        fqw_ice,ftl_ice,fraca,resfs,resft,rhokh_surft,rhokh_sice,rhokh_sea, &
        z0hssi,z0h_surft,z0mssi,z0m_surft,chr1p5m,chr1p5m_sice,smc_soilt, &
        npp_gb, resp_s_gb_um, resp_s_tot_soilt,                         &
        resp_w_pft, gc_surft, canhc_surft, wt_ext_surft, flake,         &
        surft_index, surft_pts,                                         &
        tile_frac, tstar_land, tstar_ssi, dtstar_surft,                 &
        dtstar_sea, dtstar_sice, hcons_soilt, emis_surft, emis_soil,    &
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
        nbdsc,ntdsc,wstar,wthvs,uw0,vw0,taux_p,tauy_p,rhcpt,zeros       &
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
      gamma = 1.0_r_um

      call  tr_mix (                                              &
           ! IN fields
           bl_levels, gamma, rhokm_mix(:,:,2:), rhokm_mix(:,:,1), &
           dtrdz_charney_grid, zeroes, zeroes,                    &
           kent, we_lim, t_frac, zrzi,                            &
           kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,            &
           zhnl, zhsc, z_rho,                                     &
           ! INOUT / OUT fields
           w_mixed, w_flux, surf_dep_flux                         &
           )

      do k = 1, bl_levels
        du_bl(map_w2(5)+k) = w_mixed(1,1,k) - w(1,1,k)
      end do

    end if

    rhokm_surf(map_surf(1)+0) = rhokm_land(1,1)
    rhokm_surf(map_surf(1)+1) = rhokm_ssi(1,1)
    rhokm_surf(map_surf(1)+2) = flandg(1,1)
    rhokm_surf(map_surf(1)+3) = flandfac(1,1)
    rhokm_surf(map_surf(1)+4) = fseafac(1,1)
    do k=1,bl_levels
      rhokm_bl(map_wth(1) + k) = rhokm(1,1,k)
      rhokh_bl(map_w3(1) + k) = rhokh(1,1,k)
      bq_bl(map_wth(1) + k) = bq_gb(1,1,k)
      bt_bl(map_wth(1) + k) = bt_gb(1,1,k)
      moist_flux_bl(map_w3(1) + k) = fqw(1,1,k)
      heat_flux_bl(map_w3(1) + k) = ftl(1,1,k)
      dtrdz_uv_bl(map_w3(1) + k) = dtrdz_u(1,1,k)
      dtrdz_tq_bl(map_wth(1) + k) = dtrdz_charney_grid(1,1,k)
      rdz_tq_bl(map_w3(1) + k) = rdz_charney_grid(1,1,k)
    end do
    if (formdrag == formdrag_dist_drag) then
      do k=1,bl_levels
        fd_taux(map_wth(1) + k) = tau_fd_x(1,1,k)
        fd_tauy(map_wth(1) + k) = tau_fd_y(1,1,k)
      end do
    end if
    do k=2,bl_levels
      rdz_uv_bl(map_wth(1) + k) = rdz_u(1,1,k)
      ngstress_bl(map_wth(1) + k) = f_ngstress(1,1,k)
    end do

    do i = 1, n_land_tile
      alpha1_tile(map_tile(1)+i-1) = alpha1(1, i)
      ashtf_prime_tile(map_tile(1)+i-1) = ashtf_prime_surft(1, i)
      dtstar_tile(map_tile(1)+i-1) = dtstar_surft(1, i)
      fraca_tile(map_tile(1)+i-1) = fraca(1, i)
      z0h_tile(map_tile(1)+i-1) = z0h_surft(1, i)
      z0m_tile(map_tile(1)+i-1) = z0m_surft(1, i)
      rhokh_tile(map_tile(1)+i-1) = rhokh_surft(1, i)
      chr1p5m_tile(map_tile(1)+i-1) = chr1p5m(1, i)
      resfs_tile(map_tile(1)+i-1) = resfs(1, i)
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

    blend_height_tq(map_2d(1)) = real(k_blend_tq(1,1), r_def)
    blend_height_uv(map_2d(1)) = real(k_blend_uv(1,1), r_def)
    ustar(map_2d(1)) = u_s(1,1)
    soil_moist_avail(map_2d(1)) = smc_soilt(1)
    zh_nonloc(map_2d(1)) = zhnl(1,1)
    z_lcl(map_2d(1)) = real(zlcl(1,1), r_def)
    inv_depth(map_2d(1)) = real(dzh(1,1), r_def)
    qcl_at_inv_top(map_2d(1)) = real(qcl_inv_top(1,1), r_def)
    if ( l_shallow(1,1) ) then
      shallow_flag(map_2d(1)) = 1.0_r_def
    else
      shallow_flag(map_2d(1)) = 0.0_r_def
    end if
    uw0_flux(map_2d(1)) = uw0(1,1)
    vw0_flux(map_2d(1)) = vw0(1,1)
    lcl_height(map_2d(1)) = zlcl_uv(1,1)
    parcel_top(map_2d(1)) = zhpar(1,1)
    level_parcel_top(map_2d(1)) = real(ntpar(1,1), r_def)
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
      tile_heat_flux(map_tile(1)+i-1) = real(ftl_surft(1, i), r_def)
      ! moisture flux
      tile_moisture_flux(map_tile(1)+i-1) = real(fqw_surft(1, i), r_def)
    end do

    ! Sea temperature
    tile_temperature(map_tile(1)+first_sea_tile-1) = real(tstar_sea(1,1), r_def)
    ! heat flux - NB actually grid box mean but okay for aquaplanet
    tile_heat_flux(map_tile(1)+first_sea_tile-1) = real(ftl(1,1,1), r_def)
    ! moisture flux - NB actually grid box mean but okay for aquaplanet
    tile_moisture_flux(map_tile(1)+first_sea_tile-1) = real(fqw(1,1,1), r_def)

    i_sice = 0
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      ! sea-ice temperature
      tile_temperature(map_tile(1)+i-1) = real(tstar_sice_ncat(1,1,i_sice), r_def)
      ! sea-ice heat flux
      tile_heat_flux(map_tile(1)+i-1) = real(ftl_ice(1,1,i_sice), r_def)
      ! sea-ice moisture flux
      tile_moisture_flux(map_tile(1)+i-1) = real(fqw_ice(1,1,i_sice), r_def)
    end do

    ! update blended Smagorinsky diffusion coefficients only if using Smagorinsky scheme
    if ( smagorinsky ) then
      do k = 1, nlayers
        visc_m_blend(map_wth(1) + k) = visc_m(1,1,k)
        visc_h_blend(map_wth(1) + k) = visc_h(1,1,k)
      end do
    endif

    ! write BL diagnostics
    if (flandg(1, 1) > 0.0_r_um) then
      gross_prim_prod(map_2d(1)) = real(sf_diag%gpp(1), r_def)
      net_prim_prod(map_2d(1)) = real(npp_gb(1), r_def)
      surface_conductance(map_2d(1)) = real(gs_gb(1), r_def)
      thermal_cond_wet_soil(map_2d(1)) = hcons_soilt(1)
      soil_respiration(map_2d(1)) = resp_s_gb_um(1, 1)
    else
      gross_prim_prod(map_2d(1)) = 0.0_r_def
      net_prim_prod(map_2d(1)) = 0.0_r_def
      surface_conductance(map_2d(1)) = 0.0_r_def
      thermal_cond_wet_soil(map_2d(1)) = 0.0_r_def
      soil_respiration(map_2d(1)) = 0.0_r_def
    end if

    ! update BL prognostics
    zh_2d(map_2d(1))     = zh(1,1)
    z0msea_2d(map_2d(1)) = z0msea(1,1)
    ntml_2d(map_2d(1))   = real(ntml(1,1))
    if (cumulus(1,1)) then
      cumulus_2d(map_2d(1)) = 1.0_r_def
    else
      cumulus_2d(map_2d(1)) = 0.0_r_def
    endif

    if (rh_crit_opt == rh_crit_opt_tke) then
      do k = 1, nlayers
        rh_crit(map_wth(1)+k) = real(rhcpt(1,1,k), r_def)
      end do
      rh_crit(map_wth(1)) = rh_crit(map_wth(1)+1)
    end if

    ! deallocate diagnostics deallocated in atmos_physics2
    call dealloc_bl_expl(bl_diag)
    call dealloc_sf_expl(sf_diag)
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

  end subroutine bl_exp_code

end module bl_exp_kernel_mod
