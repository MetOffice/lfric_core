!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the UM boundary layer scheme.
!>
module bl_kernel_mod

  use argument_mod,           only : arg_type,                     &
                                     GH_FIELD, GH_READ, GH_WRITE,  &
                                     GH_READWRITE, CELLS, GH_INC,  &
                                     GH_INTEGER, ANY_SPACE_1,      &
                                     ANY_SPACE_2, ANY_SPACE_3,     &
                                     ANY_SPACE_4
  use section_choice_config_mod, only : cloud, &
                                        cloud_um
  use constants_mod,          only : i_def, i_um, r_def, r_um
  use fs_continuity_mod,      only : W3, Wtheta
  use kernel_mod,             only : kernel_type
  use blayer_config_mod,      only : fixed_flux_e, fixed_flux_h
  use mixing_config_mod,      only : smagorinsky
  use timestepping_config_mod, only: outer_iterations

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: bl_kernel_type
    private
    type(arg_type) :: meta_args(122) = (/            &
        arg_type(GH_INTEGER,  GH_READ),              &! outer
        arg_type(GH_FIELD,   GH_READ,   WTHETA),     &! theta_in_wth
        arg_type(GH_FIELD,   GH_READ,   W3),         &! rho_in_w3
        arg_type(GH_FIELD,   GH_READ,   W3),         &! wetrho_in_w3
        arg_type(GH_FIELD,   GH_READ,   WTHETA),     &! wetrho_in_wth
        arg_type(GH_FIELD,   GH_READ,   W3),         &! exner_in_w3
        arg_type(GH_FIELD,   GH_READ,   WTHETA),     &! exner_in_wth
        arg_type(GH_FIELD,   GH_READ,   W3),         &! u1_in_w3
        arg_type(GH_FIELD,   GH_READ,   W3),         &! u2_in_w3
        arg_type(GH_FIELD,   GH_READ,   WTHETA),     &! u3_in_wth
        arg_type(GH_FIELD,   GH_READ,   WTHETA),     &! m_v_n
        arg_type(GH_FIELD,   GH_READ,   WTHETA),     &! m_cl_n
        arg_type(GH_FIELD,   GH_READ,   WTHETA),     &! m_ci_n
        arg_type(GH_FIELD,   GH_READ,   WTHETA),     &! theta_star
        arg_type(GH_FIELD,   GH_READ,   W3),         &! u1_star
        arg_type(GH_FIELD,   GH_READ,   W3),         &! u2_star
        arg_type(GH_FIELD,   GH_READ,   WTHETA),     &! u3_star
        arg_type(GH_FIELD,   GH_READ,   W3),         &! height_w3
        arg_type(GH_FIELD,   GH_READ,   WTHETA),     &! height_wth
        arg_type(GH_FIELD,   GH_READ,   WTHETA),     &! shear
        arg_type(GH_FIELD,   GH_READ,   WTHETA),     &! delta
        arg_type(GH_FIELD,   GH_READ,   WTHETA),     &! max_diff_smag
        arg_type(GH_FIELD,   GH_INC,  ANY_SPACE_1),  &! zh_2d
        arg_type(GH_FIELD,   GH_INC,  ANY_SPACE_1),  &! z0msea_2d
        arg_type(GH_FIELD,   GH_INC,  ANY_SPACE_1),  &! ntml_2d
        arg_type(GH_FIELD,   GH_INC,  ANY_SPACE_1),  &! cumulus_2d
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_2), & ! tile_fraction
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_3), & ! leaf_area_index
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_3), & ! canopy_height
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_1), & ! sd_orog_2d
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_1), & ! peak_to_trough_orog
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_1), & ! silhouette_area_orog
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_1), & ! soil_albedo
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_1), & ! soil_roughness
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_4), & ! soil_moist_wilt
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_4), & ! soil_moist_crit
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_4), & ! soil_moist_sat
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_4), & ! soil_cond_sat
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_4), & ! soil_thermal_cap
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_1), & ! soil_thermal_cond
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_4), & ! soil_suction_sat
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_4), & ! clapp_horn_b
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_1), & ! soil_carbon_content
        arg_type(GH_FIELD,   GH_INC,   ANY_SPACE_2), & ! tile_temperature
        arg_type(GH_FIELD,   GH_INC,   ANY_SPACE_2), & ! tile_snow_mass
        arg_type(GH_FIELD,   GH_INC,   ANY_SPACE_2), & ! n_snow_layers
        arg_type(GH_FIELD,   GH_INC,   ANY_SPACE_2), & ! snow_depth
        arg_type(GH_FIELD,   GH_INC,   ANY_SPACE_1), & ! surface_conductance
        arg_type(GH_FIELD,   GH_INC,   ANY_SPACE_2), & ! canopy_water
        arg_type(GH_FIELD,   GH_INC,   ANY_SPACE_4), & ! soil_temperature
        arg_type(GH_FIELD,   GH_INC,   ANY_SPACE_4), & ! soil_moisture
        arg_type(GH_FIELD,   GH_INC,   ANY_SPACE_4), & ! unfrozen_soil_moisture
        arg_type(GH_FIELD,   GH_INC,   ANY_SPACE_4), & ! frozen_soil_moisture
        arg_type(GH_FIELD,   GH_WRITE, ANY_SPACE_2), & ! tile_heat_flux
        arg_type(GH_FIELD,   GH_WRITE, ANY_SPACE_2), & ! tile_moisture_flux
        arg_type(GH_FIELD,   GH_WRITE, ANY_SPACE_1), & ! gross_prim_prod
        arg_type(GH_FIELD,   GH_WRITE, ANY_SPACE_1), & ! net_prim_prod
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_1), & ! cos_zen_angle
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_2), & ! sw_up_tile
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_1), & ! sw_down_surf
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_1), & ! lw_down_surf
        arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_1), & ! sw_down_surf_blue
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! dtheta_bl
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! dt_bl
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! dmv_bl
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! dt_conv
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! dmv_conv
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! dmcl_conv
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! dmcf_conv
        arg_type(GH_FIELD,   GH_READWRITE,  W3),     &! du_conv
        arg_type(GH_FIELD,   GH_READWRITE,  W3),     &! dv_conv
        arg_type(GH_FIELD,   GH_READ,  WTHETA),      &! dtl_mphys
        arg_type(GH_FIELD,   GH_READ,  WTHETA),      &! dmt_mphys
        arg_type(GH_FIELD,   GH_READ,  WTHETA),      &! sw_heating_rate
        arg_type(GH_FIELD,   GH_READ,  WTHETA),      &! lw_heating_rate
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! m_v
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! m_cl
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! m_ci
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! m_r
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! m_g
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! cf_area
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! cf_ice
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! cf_liq
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! cf_bulk
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! rh_crit_wth
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! visc_m_blend
        arg_type(GH_FIELD,   GH_READWRITE,  WTHETA), &! visc_h_blend
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! cca
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! ccw
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! deep_in_col
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! shallow_in_col
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! mid_in_col
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! freeze_level
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! deep_prec
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! shallow_prec
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! mid_prec
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! deep_term
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! cape_timescale
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! conv_rain
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! conv_snow
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! lowest_cv_base
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! lowest_cv_top
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! cv_base
        arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1),&! cv_top
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! massflux_up
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! massflux_down
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! entrain_up
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! entrain_down
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! detrain_up
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! detrain_down
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! dd_dt
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! dd_dq
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! deep_dt
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! deep_dq
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! deep_massflux
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! deep_tops
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! shallow_dt
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! shallow_dq
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! shallow_massflux
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! mid_dt
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA),     &! mid_dq
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA)      &! mid_massflux
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass ::bl_code
  end type

  !-----------------------------------------------------------------------------
  ! Constructors
  !-----------------------------------------------------------------------------

  ! Overload the default structure constructor for function space
  interface bl_kernel_type
    module procedure bl_kernel_constructor
  end interface

  public bl_code

contains

  type(bl_kernel_type) function bl_kernel_constructor() result(self)
    implicit none
    return
  end function bl_kernel_constructor

  !> @brief Interface to the UM BL scheme
  !> @details The UM Boundary Layer scheme does:
  !>             vertical mixing of heat, momentum and moisture,
  !>             as documented in UMDP24
  !>          NB This version uses winds in w3 space (i.e. A-grid)
  !>            and doesn't currently feed-back wind increments
  !> @param[in]     nlayers       Number of layers
  !> @param[in]     outer         Outer loop counter
  !> @param[in]     theta_in_wth  Potential temperature field
  !> @param[in]     rho_in_w3     Density field in density space
  !> @param[in]     wetrho_in_w3  Wet density field in density space
  !> @param[in]     wetrho_in_wth Wet density field in wth space
  !> @param[in]     exner_in_w3   Exner pressure field in density space
  !> @param[in]     exner_in_wth  Exner pressure field in wth space
  !> @param[in]     u1_in_w3      'Zonal' wind in density space
  !> @param[in]     u2_in_w3      'Meridional' wind in density space
  !> @param[in]     u3_in_wth     'Vertical' wind in theta space
  !> @param[in]     m_v_n         Vapour mixing ratio at time level n
  !> @param[in]     m_cl_n        Cloud liquid mixing ratio at time level n
  !> @param[in]     m_ci_n        Cloud ice mixing ratio at time level n
  !> @param[in]     theta_star    Potential temperature predictor after advection
  !> @param[in]     u1_star       'Zonal' wind predictor after advection
  !> @param[in]     u2_star       'Meridional' wind predictor after advection
  !> @param[in]     u3_star       'Vertical' wind predictor after advection
  !> @param[in]     height_w3     Height of density space levels above surface
  !> @param[in]     height_wth    Height of theta space levels above surface
  !> @param[in]     shear         3D wind shear on wtheta points
  !> @param[in]     delta         Edge length on wtheta points
  !> @param[in]     max_diff_smag maximum diffusion coefficient allowed in this run
  !> @param[in,out] zh_2d         Boundary layer depth
  !> @param[in,out] z0msea_2d     Roughness length
  !> @param[in,out] ntml_2d       Number of turbulently mixed levels
  !> @param[in,out] cumulus_2d    Cumulus flag (true/false)
  !> @param[in] tile_fraction      Surface tile fractions
  !> @param[in] leaf_area_index    Leaf Area Index
  !> @param[in] canopy_height      Canopy height
  !> @param[in] sd_orog_2d         Standard deviation of orography
  !> @param[in] peak_to_trough_orog  Half of peak-to-trough height over root(2) of orography
  !> @param[in] silhouette_area_orog Silhouette area of orography
  !> @param[in] soil_albedo        Snow-free soil albedo
  !> @param[in] soil_roughness     Bare soil surface roughness length
  !> @param[in] soil_moist_wilt    Volumetric soil moisture at wilting point
  !> @param[in] soil_moist_crit    Volumetric soil moisture at critical point
  !> @param[in] soil_moist_sat     Volumetric soil moisture at saturation
  !> @param[in] soil_cond_sat      Saturated soil conductivity
  !> @param[in] soil_thermal_cap   Soil thermal capacity
  !> @param[in] soil_thermal_cond  Soil thermal conductivity
  !> @param[in] soil_suction_sat   Saturated soil water suction
  !> @param[in] clapp_horn_b       Clapp and Hornberger b coefficient
  !> @param[in] soil_carbon_content Soil carbon content
  !> @param[in,out] tile_temperature   Surface tile temperatures
  !> @param[in,out] tile_snow_mass     Snow mass on tiles (kg/m2)
  !> @param[in,out] n_snow_layers      Number of snow layers on tiles
  !> @param[in,out] snow_depth         Snow depth on tiles
  !> @param[in,out] surface_conductance Surface conductance
  !> @param[in,out] canopy_water       Canopy water on each tile
  !> @param[in,out] soil_temperature   Soil temperature
  !> @param[in,out] soil_moisture      Soil moisture content (kg m-2)
  !> @param[in,out] unfrozen_soil_moisture Unfrozen soil moisture proportion
  !> @param[in,out] frozen_soil_moisture Frozen soil moisture proportion
  !> @param[out]    tile_heat_flux      Surface heat flux
  !> @param[out]    tile_moisture_flux  Surface moisture flux
  !> @param[out]    gross_prim_prod     Gross Primary Productivity
  !> @param[out]    net_prim_prod       Net Primary Productivity
  !> @param[in]     cos_zen_angle       Cosing of solar zenith angle
  !> @param[in]     sw_up_tile          Upwelling SW radiation on surface tiles
  !> @param[in]     sw_down_surf        Downwelling SW radiation at surface
  !> @param[in]     lw_down_surf        Downwelling LW radiation at surface
  !> @param[in]     sw_down_surf_blue   Photosynthetically active SW down
  !> @param[out]    dtheta_bl     BL theta increment
  !> @param[out]    dt_bl         BL temperature increment
  !> @param[out]    dmv_bl        BL vapour increment
  !> @param[in,out] dt_conv       Convection temperature increment
  !> @param[in,out] dmv_conv      Convection vapour increment
  !> @param[in,out] dmcl_conv     Convection cloud liquid increment
  !> @param[in,out] dmcf_conv     Convection cloud ice increment
  !> @param[in,out] du_conv       Convection du increment
  !> @param[in,out] dv_conv       Convection dv increment
  !> @param[in]     dtl_mphys       Microphysics liquid temperature increment
  !> @param[in]     dmt_mphys       Microphysics total water increment
  !> @param[in]     sw_heating_rate Shortwave radiation heating rate
  !> @param[in]     lw_heating_rate Longwave radiation heating rate
  !> @param[in,out] m_v             Vapour mixing ration after advection
  !> @param[in,out] m_cl            Cloud liquid mixing ratio after advection
  !> @param[in,out] m_ci            Cloud ice mixing ratio after advection
  !> @param[in,out] m_r             Rain mixing ratio after advection
  !> @param[in,out] m_g             Graupel mixing ratio after advection
  !> @param[in,out] cf_area         Area cloud fraction
  !> @param[in,out] cf_ice          Ice cloud fraction
  !> @param[in,out] cf_liq          Liquid cloud fraction
  !> @param[in,out] cf_bulk         Bulk cloud fraction
  !> @param[in,out] rh_crit_wth   Critical rel humidity in pot temperature space
  !> @param[in,out] visc_m_blend  Blended BL-Smag diffusion coefficient for momentum
  !> @param[in,out] visc_h_blend  Blended BL-Smag diffusion coefficient for scalars
  !> @param[out]    cca           convective cloud amount (fraction)
  !> @param[out]    ccw           convective cloud water (kg/kg) (can be ice or liquid)
  !> @param[out]    deep_in_col   indicator of deep in column
  !> @param[out]    shallow_in_col indicator of shallow in column
  !> @param[out]    mid_in_col    indicator of mid in column
  !> @param[out]    freeze_level  level number of freezing level
  !> @param[out]    deep_prec     precipitation rate from deep convection(kg/m2/s)
  !> @param[out]    shallow_prec  precipitation rate from shallow convection(kg/m2/s)
  !> @param[out]    mid_prec      precipitation rate from mid convection(kg/m2/s)
  !> @param[out]    deep_term     termination level number of deep convection
  !> @param[out]    cape_timescale cape timescale (s)
  !> @param[out]    conv_rain     surface rainfall rate from convection (kg/m2/s)
  !> @param[out]    conv_snow     surface snowfall rate from convection (kg/m2/s)
  !> @param[out]    lowest_cv_base level number for start of convection in column
  !> @param[out]    lowest_cv_top level number for end of lowest convection in column
  !> @param[out]    cv_base       level number of base of highest convection in column
  !> @param[out]    cv_top        level number for end of highest convection in column
  !> @param[out]    massflux_up   convective upwards mass flux (Pa/s)
  !> @param[out]    massflux_down convective downwards mass flux (Pa/s)
  !> @param[out]    entrain_up    convective upwards entrainment
  !> @param[out]    entrain_down  convective downwards entrainment
  !> @param[out]    detrain_up    convective upwards detrainment
  !> @param[out]    detrain_down  convective downwards detrainment
  !> @param[out]    dd_dt         temperature increment from downdraughts per timestep
  !> @param[out]    dd_dq         vapour increment from downdraughts per timestep
  !> @param[out]    deep_dt       temperature increment from deep convection per timestep
  !> @param[out]    deep_dq       vapour increment from deep convection per timestep
  !> @param[out]    deep_massflux upward mass flux from deep convection
  !> @param[out]    deep_tops     set to 1.0 if deep stops at this model level
  !> @param[out]    shallow_dt    temperature increment from shallow convection per timestep
  !> @param[out]    shallow_dq    vapour increment from shallow convection per timestep
  !> @param[out]    shallow_massflux  upward mass flux from shallow convection
  !> @param[out]    mid_dt        temperature increment from mid convection per timestep
  !> @param[out]    mid_dq        vapour increment from mid convection per timestep
  !> @param[out]    mid_massflux  upward mass flux from mid convection
  !> @param[in]     ndf_wth       Number of degrees of freedom per cell for potential temperature space
  !> @param[in]     undf_wth      Number unique of degrees of freedom  for potential temperature space
  !> @param[in]     map_wth       dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_w3        Number of degrees of freedom per cell for density space
  !> @param[in]     undf_w3       Number unique of degrees of freedom  for density space
  !> @param[in]     map_w3        dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_2d        Number of degrees of freedom per cell for 2D fields
  !> @param[in]     undf_2d       Number unique of degrees of freedom  for 2D fields
  !> @param[in]     map_2d        dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_tile      Number of DOFs per cell for tiles
  !> @param[in]     undf_tile     Number of total DOFs for tiles
  !> @param[in]     map_tile      Dofmap for cell for surface tiles
  !> @param[in]     ndf_pft       Number of DOFs per cell for PFTs
  !> @param[in]     undf_pft      Number of total DOFs for PFTs
  !> @param[in]     map_pft       Dofmap for cell for PFTs
  !> @param[in]     ndf_soil      Number of DOFs per cell for soil levels
  !> @param[in]     undf_soil     Number of total DOFs for soil levels
  !> @param[in]     map_soil      Dofmap for cell for soil levels
  subroutine bl_code(nlayers,                               &
                     outer,                                 &
                     theta_in_wth,                          &
                     rho_in_w3,                             &
                     wetrho_in_w3,                          &
                     wetrho_in_wth,                         &
                     exner_in_w3,                           &
                     exner_in_wth,                          &
                     u1_in_w3,                              &
                     u2_in_w3,                              &
                     u3_in_wth,                             &
                     m_v_n,                                 &
                     m_cl_n,                                &
                     m_ci_n,                                &
                     theta_star,                            &
                     u1_star,                               &
                     u2_star,                               &
                     u3_star,                               &
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
                     soil_cond_sat,                         &
                     soil_thermal_cap,                      &
                     soil_thermal_cond,                     &
                     soil_suction_sat,                      &
                     clapp_horn_b,                          &
                     soil_carbon_content,                   &
                     tile_temperature,                      &
                     tile_snow_mass,                        &
                     n_snow_layers,                         &
                     snow_depth,                            &
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
                     dtheta_bl,                             &
                     dt_bl,                                 &
                     dmv_bl,                                &
                     dt_conv,                               &
                     dmv_conv,                              &
                     dmcl_conv,                             &
                     dmcf_conv,                             &
                     du_conv,                               &
                     dv_conv,                               &
                     dtl_mphys,                             &
                     dmt_mphys,                             &
                     sw_heating_rate,                       &
                     lw_heating_rate,                       &
                     m_v,                                   &
                     m_cl,                                  &
                     m_ci,                                  &
                     m_r,                                   &
                     m_g,                                   &
                     cf_area,                               &
                     cf_ice,                                &
                     cf_liq,                                &
                     cf_bulk,                               &
                     rh_crit_wth,                           &
                     visc_m_blend,                          &
                     visc_h_blend,                          &
                     cca,                                   &
                     ccw,                                   &
                     deep_in_col,                           &
                     shallow_in_col,                        &
                     mid_in_col,                            &
                     freeze_level,                          &
                     deep_prec,                             &
                     shallow_prec,                          &
                     mid_prec,                              &
                     deep_term,                             &
                     cape_timescale,                        &
                     conv_rain,                             &
                     conv_snow,                             &
                     lowest_cv_base,                        &
                     lowest_cv_top,                         &
                     cv_base,                               &
                     cv_top,                                &
                     massflux_up,                           &
                     massflux_down,                         &
                     entrain_up,                            &
                     entrain_down,                          &
                     detrain_up,                            &
                     detrain_down,                          &
                     dd_dt,                                 &
                     dd_dq,                                 &
                     deep_dt,                               &
                     deep_dq,                               &
                     deep_massflux,                         &
                     deep_tops,                             &
                     shallow_dt,                            &
                     shallow_dq,                            &
                     shallow_massflux,                      &
                     mid_dt,                                &
                     mid_dq,                                &
                     mid_massflux,                          &
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
                     ndf_soil, undf_soil, map_soil)

    !---------------------------------------
    ! LFRic modules
    !---------------------------------------
    use init_jules_alg_mod, only: &
         n_surf_tile, n_land_tile, n_sea_tile, n_sea_ice_tile, &
         first_land_tile, first_sea_tile, first_sea_ice_tile,  &
         n_soil_levs

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    ! structures holding diagnostic arrays - not used
    use bl_diags_mod, only: BL_diag, dealloc_bl_imp
    use sf_diags_mod, only: sf_diag, dealloc_sf_expl, dealloc_sf_imp
    ! other modules containing stuff passed to BL
    use ancil_info, only: l_soil_point, ssi_pts,                            &
         sea_pts, sice_pts, ssi_index, sea_index, sice_index, fssi_ij,      &
         sea_frac, sice_frac, sice_pts_ncat, sice_index_ncat, sice_frac_ncat
    use atm_fields_bounds_mod, only: tdims, udims, vdims, udims_s, vdims_s, &
         pdims
    use atm_step_local, only: dim_cs1, dim_cs2
    use atmos_physics2_alloc_mod !everything
    use bdy_expl3_mod, only: bdy_expl3
    use bl_option_mod, only: flux_bc_opt, specified_fluxes_only
    use conv_diag_6a_mod, only: conv_diag_6a
    use gen_phys_inputs_mod, only: l_mr_physics
    use jules_sea_seaice_mod, only: nice_use
    use jules_surface_types_mod, only: npft, ntype
    use level_heights_mod, only: r_theta_levels, r_rho_levels, eta_theta_levels
    use ni_bl_ctl_mod, only: ni_bl_ctl
    use ni_imp_ctl_mod, only: ni_imp_ctl
    use nlsizes_namelist_mod, only: row_length, rows, land_field,&
         sm_levels, ntiles, model_levels, bl_levels, tr_vars, n_cca_lev
    use ozone_vars, only: o3_gb
    use planet_constants_mod, only: p_zero, kappa, planet_radius
    use prognostics, only: snowdepth_surft, nsnow_surft
    use p_s_parms, only: bexp_soilt, sathh_soilt, hcap_soilt, satcon_soilt
    use rad_input_mod, only: co2_mmr
    use timestep_mod, only: timestep
    use turb_diff_ctl_mod, only: visc_m, visc_h, max_diff, delta_smag

    ! Jules related subroutines
    use sparm_mod, only: sparm
    use tilepts_mod, only: tilepts

    ! Flags controlling the saving of fields
    use jules_surface_mod, ONLY: ISrfExCnvGust, IP_SrfExWithCnv

    ! Convection scheme modules 
    use mphys_inputs_mod, only: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
    use cloud_inputs_mod, only: i_cld_vn
    use pc2_constants_mod, only: i_cld_pc2

    use cv_run_mod, only: l_param_conv, i_convection_vn, i_convection_vn_6a, &
                          n_conv_calls, iconv_deep, iconv_shallow, l_mom,    &
                          qmin_conv, l_safe_conv, cldbase_opt_dp,            &
                          cldbase_opt_md

    use glue_conv_6a_mod,                only: glue_conv_6a
    use scm_convss_dg_mod,               only: scm_convss_dg_type
    use atmos_physics2_save_restore_mod, only: ap2_init_conv_diag

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: outer
    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3

    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile)
    integer(kind=i_def), intent(in) :: ndf_pft, undf_pft
    integer(kind=i_def), intent(in) :: map_pft(ndf_pft)
    integer(kind=i_def), intent(in) :: ndf_soil, undf_soil
    integer(kind=i_def), intent(in) :: map_soil(ndf_soil)

    integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth
    integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
    integer(kind=i_def), dimension(ndf_2d),  intent(in) :: map_2d

    real(kind=r_def), dimension(undf_wth), intent(out)  :: dtheta_bl, dt_bl,   &
                                                           dmv_bl
    real(kind=r_def), dimension(undf_wth), intent(inout):: m_v, m_cl, m_ci,    &
                                                           m_r, m_g,           &
                                                           cf_area, cf_ice,    &
                                                           cf_liq, cf_bulk,    &
                                                           rh_crit_wth,        &
                                                           visc_h_blend,       &
                                                           visc_m_blend
    real(kind=r_def), dimension(undf_w3),  intent(in)   :: rho_in_w3,          &
                                                           wetrho_in_w3,       &
                                                           exner_in_w3,        &
                                                           u1_in_w3, u2_in_w3, &
                                                           u1_star, u2_star,   &
                                                           height_w3
    real(kind=r_def), dimension(undf_wth), intent(in)   :: theta_in_wth,       &
                                                           wetrho_in_wth,      &
                                                           exner_in_wth,       &
                                                           u3_in_wth,          &
                                                           m_v_n, m_cl_n,      &
                                                           m_ci_n,             &
                                                           theta_star,         &
                                                           u3_star,            &
                                                           height_wth,         &
                                                           shear,              &
                                                           delta,              &
                                                           max_diff_smag,      &
                                                           dtl_mphys,dmt_mphys,&
                                                           sw_heating_rate,    &
                                                           lw_heating_rate
    real(kind=r_def), dimension(undf_wth), intent(inout) :: dt_conv, dmv_conv, &
                                                          dmcl_conv, dmcf_conv
    real(kind=r_def), dimension(undf_wth), intent(out) :: cca, ccw,            &
                          massflux_up, massflux_down, entrain_up,entrain_down, &
                          detrain_up, detrain_down, dd_dt, dd_dq,              &
                          deep_dt, deep_dq, deep_massflux, deep_tops,          &
                          shallow_dt, shallow_dq, shallow_massflux,            &
                          mid_dt, mid_dq, mid_massflux 

    real(kind=r_def), dimension(undf_w3),  intent(inout) :: du_conv, dv_conv

    real(kind=r_def), dimension(undf_2d), intent(inout) :: zh_2d,              &
                                                           z0msea_2d, ntml_2d, &
                                                           cumulus_2d

    real(kind=r_def), dimension(undf_2d), intent(out) ::  deep_in_col,         &
                           shallow_in_col, mid_in_col,                         &
                           freeze_level, deep_prec, shallow_prec, mid_prec,    &
                           deep_term, cape_timescale, conv_rain, conv_snow,    &
                           lowest_cv_base, lowest_cv_top, cv_base, cv_top

    real(kind=r_def), intent(in) :: tile_fraction(undf_tile)
    real(kind=r_def), intent(inout) :: tile_temperature(undf_tile)
    real(kind=r_def), intent(inout) :: tile_snow_mass(undf_tile)
    real(kind=r_def), intent(inout) :: n_snow_layers(undf_tile)
    real(kind=r_def), intent(inout) :: snow_depth(undf_tile)
    real(kind=r_def), intent(inout) :: canopy_water(undf_tile)
    real(kind=r_def), intent(out) :: tile_heat_flux(undf_tile)
    real(kind=r_def), intent(out) :: tile_moisture_flux(undf_tile)
    real(kind=r_def), intent(in) :: sw_up_tile(undf_tile)

    real(kind=r_def), intent(in) :: leaf_area_index(undf_pft)
    real(kind=r_def), intent(in) :: canopy_height(undf_pft)

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
    real(kind=r_def), intent(out) :: gross_prim_prod(undf_2d)
    real(kind=r_def), intent(out) :: net_prim_prod(undf_2d)

    real(kind=r_def), intent(in) :: soil_moist_wilt(undf_soil)
    real(kind=r_def), intent(in) :: soil_moist_crit(undf_soil)
    real(kind=r_def), intent(in) :: soil_moist_sat(undf_soil)
    real(kind=r_def), intent(in) :: soil_cond_sat(undf_soil)
    real(kind=r_def), intent(in) :: soil_thermal_cap(undf_soil)
    real(kind=r_def), intent(in) :: soil_suction_sat(undf_soil)
    real(kind=r_def), intent(in) :: clapp_horn_b(undf_soil)
    real(kind=r_def), intent(inout) :: soil_temperature(undf_soil)
    real(kind=r_def), intent(inout) :: soil_moisture(undf_soil)
    real(kind=r_def), intent(inout) :: unfrozen_soil_moisture(undf_soil)
    real(kind=r_def), intent(inout) :: frozen_soil_moisture(undf_soil)

    ! Local variables for the kernel
    integer(i_def) :: k, i, i_tile, i_sice, i_pft, n

    ! switches and model parameters/dimensions/time etc
    integer(i_um) :: cycleno, error_code
    integer(i_um) :: curr_year, curr_day_number, curr_hour, curr_minute, &
                     curr_second
    integer(i_um) :: co2_dim_len, co2_dim_row
    integer(i_um) :: asteps_since_triffid
    integer(i_um), parameter :: nscmdpkgs=15
    logical,       parameter :: l_scmdiags(nscmdpkgs)=.false.
    logical :: l_aero_classic, l_spec_z0, &
               l_extra_call, l_calc_at_p, l_jules_call

    ! profile fields from level 1 upwards
    real(r_um), dimension(row_length,rows,nlayers) ::                        &
         p_rho_levels, rho_wet_rsq, rho_wet, rho_dry, z_rho, z_theta,        &
         bulk_cloud_fraction, bl_w_var, rhcpt, t_latest, q_latest,           &
         qcl_latest, qcf_latest, cf_latest, cfl_latest, cff_latest, cca_3d,  &
         cca0, ccw0, area_cloud_fraction, cloud_fraction_liquid,             &
         cloud_fraction_frozen
    real(r_um), dimension(row_length,rows,bl_levels) ::                 &
         e_trb, tsq_trb, qsq_trb, cov_trb, tau_fd_x, tau_fd_y, rhogamu, &
         rhogamv, f_ngstress, fqw, ftl, rhokh
    real(r_um), dimension(row_length,rows,nlayers) :: rho_wet_tq
                             !IB removed -1 from nlayers here
    real(r_um), dimension(row_length,rows,nlayers) :: exner_rho_levels
                             !IB removed +1 from nlayers here
    ! profile fields on u/v points
    real(r_um), dimension(row_length,rows,nlayers) :: u, v, r_u, r_v
    real(r_um), dimension(row_length,rows,bl_levels) ::            &
         taux_fd_u, tauy_fd_v, rhogamu_u, rhogamv_v, f_ngstress_u, &
         f_ngstress_v, taux, tauy
    ! profile fields from level 0 upwards
    real(r_um), dimension(row_length,rows,0:nlayers) ::                 &
         p_theta_levels, w_copy, etadot_copy, R_w, p_rho_minus_one, w,  &
         exner_rho_minus_one

    ! profile fields with 0 level which isn't used
    real(r_um), dimension(row_length,rows,0:nlayers) ::                      &
         q, qcl, qcf, theta, exner_theta_levels, aerosol, dust_div1,         &
         dust_div2, dust_div3, dust_div4, dust_div5, dust_div6, so2, dms,    &
         so4_aitken, so4_accu, so4_diss, nh3, soot_new, soot_aged, soot_cld, &
         bmass_new, bmass_aged, bmass_cld, ocff_new, ocff_aged, ocff_cld,    &
         nitr_acc, nitr_diss, ozone_tracer, conv_prog_precip
    ! multi tracer variables
    real(r_um), dimension(row_length,rows,0:nlayers,tr_vars) :: &
         free_tracers
    ! profile fields with a hard-wired 2
    real(r_um), dimension(row_length,rows,2,bl_levels) :: rad_hr, micro_tends
    ! single level fields
    real(r_um), dimension(row_length,rows) ::                                &
         p_star, lw_down, cos_zenith_angle, ti_sicat, tstar, zh_prev, ddmfx, &
         zlcl, zhpar, z0h_scm, z0m_scm, flux_e, flux_h, z0msea,              &
         photosynth_act_rad, soil_clay, soil_sand, dust_mrel1, dust_mrel2,   &
         dust_mrel3, dust_mrel4, dust_mrel5, dust_mrel6, tstar_sea, zh, dzh, &
         zhpar_shcu, flandfac, fseafac, rhokm_land, rhokm_ssi, cdr10m,       &
         tstar_land, tstar_ssi, dtstar_sea, t1_sd, q1_sd, wstar, wthvs,      &
         xx_cos_theta_latitude, ice_fract, ls_rain, ls_snow, conv_rain_copy, &
         conv_snow_copy, qcl_inv_top, co2_emits, co2flux, tscrndcl_ssi,      &
         tstbtrans, sum_eng_fluxes, sum_moist_flux, drydep2, olr,            &
         surf_ht_flux_land, zlcl_mixed, theta_star_surf, qv_star_surf,       &
         snowmelt, tstar_sice, u_0_p, v_0_p, w_max, deep_flag, past_precip,  &
         past_conv_ht, zlcl_uv, ql_ad, cin_undilute, cape_undilute,          &
         entrain_coef, qsat_lcl, delthvu, dtstar_sice, ustar_in, g_ccp,      &
         h_ccp, fb_surf, charnock_w, uwind_wav, vwind_wav, sstfrz
    ! single level fields on u/v points
    real(r_um), dimension(row_length,rows) :: u_0, v_0, rhokm_u_land,    &
         rhokm_u_ssi, rhokm_v_land, rhokm_v_ssi, flandfac_u, flandfac_v, &
         fseafac_u, fseafac_v, taux_land, tauy_land, taux_ssi, tauy_ssi
    ! single level fields
    integer(i_um), dimension(row_length,rows) :: nlcl, ntml, ntpar, lcbase, &
                                                 ccb0, cct0, conv_type,     &
                                                 nbdsc, ntdsc
    ! single level fields
    logical, dimension(row_length,rows) :: land_sea_mask, cumulus,       &
                                           l_shallow,                    &
                                           no_cumulus, l_congestus,      &
                                           l_congestus2
    ! fields on ice categories
    real(r_um), dimension(row_length,rows,nice_use) ::                    &
         ice_fract_ncat_sicat, k_sice_sicat, co2, ti_cat_sicat,           &
         tstar_sice_sicat, radnet_sice, fqw_ice, ftl_ice, di_ncat_sicat,  &
         ice_fract_ncat
    ! field on land points and soil levels
    real(r_um), dimension(land_field,sm_levels) :: soil_layer_moisture, &
         smvccl_soilt, smvcwt_soilt, smvcst_soilt, sthf_soilt,          &
         sthu_soilt, ext, t_soil_soilt
    ! fields on land points
    real(r_um), dimension(land_field) :: hcon_soilt, sil_orog_land_gb,   &
         ho2r2_orog_gb, sd_orog, z0m_soil_gb, albsoil_soilt, gs_gb
    integer, dimension(land_field) :: land_index
    real(r_um), dimension(land_field,dim_cs1) :: cs_pool_gb_um
    real(r_um), dimension(land_field,dim_cs2) :: resp_s_acc_gb_um
    ! fields on land points and tiles
    real(r_um), dimension(land_field,ntiles) :: canopy_surft, catch_surft, &
         catch_snow_surft, snow_surft, z0_surft, z0h_bare_surft, sw_surft, &
         tstar_surft, frac_surft, ftl_surft, fqw_surft, epot_surft,        &
         dtstar_surft, tsurf_elev_surft, tscrndcl_surft, ei_surft,         &
         ecan_surft, melt_surft, surf_htf_surft
    ! fields on land points and pfts
    real(r_um), dimension(land_field,npft) :: canht_pft, lai_pft,          &
         g_leaf_acc_pft, npp_acc_pft, resp_w_acc_pft
    ! stashwork arrays
    real(r_um) :: stashwork3(1), stashwork9(1)

    !-----------------------------------------------------------------------
    ! Extra Convection local arrays
    !-----------------------------------------------------------------------

    integer(i_um) :: n_cumulus, n_deep, n_shallow, n_congestus, n_mid,   &
                    ntra_flds, ntra_lev, segments, n_conv_levels,        &
                    call_number, ntra_fld, seg_num 

    ! single level integer fields
    integer(i_um), dimension(row_length,rows) ::             &
        it_lcbase,it_lctop, it_ccb, it_cct,                  &
        it_ccb0, it_cct0, it_kterm_deep, it_kterm_shall,     &
        it_cg_term, it_lcbase0, freeze_lev, ccb, cct, lctop

    logical :: l_tracer, l_calc_dxek, l_q_interact, l_scm_convss_dg,   &
               l_cape_opt

    logical, dimension(row_length,rows) :: it_mid_level,l_mid

    real(r_um)  ::                                                        &
        timestep_conv, one_over_conv_calls, orig_value

    ! multi level BUT not level 0 - arrays
    real(r_um), dimension(row_length,rows,nlayers) ::                      &
        theta_conv, q_conv, qcl_conv, qcf_conv, qrain_conv, qcf2_conv,     &
        qgraup_conv, cf_liquid_conv, cf_frozen_conv, bulk_cf_conv,         &
        u_conv,  v_conv, dq_add, ccw_3d,                                   &
        dthbydt, dqbydt, dqclbydt, dqcfbydt, dcflbydt, dcffbydt, dbcfbydt, &
        dubydt_p, dvbydt_p, rho_dry_theta

    real(r_um), dimension(row_length,rows,nlayers) ::                      &
        it_ccw, it_ccw0, it_conv_rain_3d, it_conv_snow_3d, it_cca,         &
        it_cca0, it_w2p, it_cca0_dp, it_cca0_md, it_cca0_sh,               &
        it_up_flux_half, it_up_flux, it_dwn_flux,                          &
        it_entrain_up, it_detrain_up, it_entrain_dwn, it_detrain_dwn,      &
        it_uw_dp, it_vw_dp, it_uw_shall, it_vw_shall, it_uw_mid, it_vw_mid,&
        it_wqt_flux, it_wthetal_flux, it_wthetav_flux, it_wql_flux,        &
        it_mf_deep, it_mf_congest, it_mf_shall, it_mf_midlev,              &
        it_dt_deep, it_dt_congest, it_dt_shall, it_dt_midlev,              &
        it_dq_deep, it_dq_congest, it_dq_shall, it_dq_midlev,              &
        it_du_deep, it_du_congest, it_du_shall, it_du_midlev,              &
        it_dv_deep, it_dv_congest, it_dv_shall, it_dv_midlev,              &
        it_dt_dd, it_dq_dd, it_du_dd, it_dv_dd, it_area_ud, it_area_dd

    ! Increments for PC2 from convection
    real(r_um), dimension(row_length,rows,nlayers) ::                  &
        cf_liquid_inc, cf_frozen_inc, bulk_cf_inc

    ! single level fields
    real(r_um), dimension(row_length,rows) ::                          &
        it_lcca,it_cca_2d, it_cclwp,                                   &
        it_cclwp0, it_conv_rain, it_conv_snow,                         &
        it_precip_dp, it_precip_sh, it_precip_md, it_cape_out,         &
        it_dp_cfl_limited, it_md_cfl_limited,                          &
        ind_cape_reduced, cape_ts_used, it_ind_deep, it_ind_shall,     &
        it_precip_cg, it_wstar_dn, it_wstar_up, it_mb1, it_mb2

    ! total tracer - has to be allocatable
    real(r_um), dimension(:,:,:,:), allocatable :: &
        tot_tracer

    ! Only required for argument list to glue_conv6a
    type(scm_convss_dg_type), allocatable :: scm_convss_dg(:,:)

    !-----------------------------------------------------------------------
    ! Initialisation of variables and arrays
    ! These will need to be configured/passed in by LFRic, but for
    ! the initial implementation these will be set here...
    !-----------------------------------------------------------------------
    ! dynamics loop counter
    cycleno=outer
    ! model dimensions
    co2_dim_len=1
    co2_dim_row=1
    ! diagnostic flags
    error_code=0
    ! other logicals
    l_aero_classic=.false.
    l_extra_call=.false.
    l_jules_call=.false.
    ! surface forcing
    if ( flux_bc_opt == specified_fluxes_only ) then
      flux_e(:,:)=fixed_flux_e
      flux_h(:,:)=fixed_flux_h
    end if
    l_spec_z0=.false.


    if ( .not. allocated(uhalo) ) then
      call atmos_physics2_alloc( land_field, &
                                 ntiles,     &
                                 1_i_um,     &
                                 1_i_um,     &
                                 npft,ntype, &
                                 1_i_um,     &
                                 1_i_um,     &
                                 sm_levels,  &
                                 bl_levels,  &
                                 nice_use )
    end if

    !----------------------------------------------------------------------
    ! Surface fields as needed by Jules
    !----------------------------------------------------------------------
    ! Land tile fractions
    i_tile = 0
    flandg = 0.0_r_um
    do i = first_land_tile, first_land_tile + n_land_tile - 1
      i_tile = i_tile + 1
      flandg = flandg + real(tile_fraction(map_tile(i)), r_um)
      frac_surft(1, i_tile) = real(tile_fraction(map_tile(i)), r_um)
    end do
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
    call tilepts(land_field, frac_surft, surft_pts, surft_index)

    ! Sea-ice fraction
    i_sice = 0
    ice_fract = 0.0_r_um
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      ice_fract = ice_fract + real(tile_fraction(map_tile(i)), r_um)
      ice_fract_ncat(1, 1, i_sice) = real(tile_fraction(map_tile(i)), r_um)
    end do

    ! Jules requires sea-ice fractions with respect to the sea area
    if (ice_fract(1, 1) > 0.0_r_um) then
      ice_fract(1, 1) = ice_fract(1, 1) / (1.0_r_um - flandg(1, 1))
      ice_fract_ncat(1, 1, 1:n_sea_ice_tile) &
           = ice_fract_ncat(1, 1, 1:n_sea_ice_tile) / (1.0_r_um - flandg(1, 1))
    end if
    ice_fract_ncat_sicat = ice_fract_ncat

    ! combined sea and sea-ice index
    ssi_pts = 1
    if (flandg(1, 1) < 1.0_r_um) then
      ssi_index = 1
    else
      ssi_index = 0
    end if
    fssi_ij = 1.0_r_um - flandg(1, 1)

    ! individual sea and sea-ice indices
    if (ssi_index(1) > 0) then
      if (ice_fract(1, 1) > 0.0_r_um) then
        sice_pts = 1
        sice_index = 1
        sice_frac = ice_fract(1, 1)
      else
        sice_pts = 0
        sice_index = 0
        sice_frac = 0.0_r_um
      end if
      if (ice_fract(1, 1) < 1.0_r_um) then
        sea_pts = 1
        sea_index = 1
        sea_frac = 1.0_r_um - sice_frac
      else
        sea_pts = 0
        sea_index = 0
        sea_frac = 0.0_r_um
      end if
    end if

    ! multi-category sea-ice index
    do n = 1, nice_use
      if (ssi_index(1) > 0 .and. ice_fract_ncat(1, 1, n) > 0.0_r_um) then
        sice_pts_ncat(n) = 1
        sice_index_ncat(1, n) = 1
        sice_frac_ncat(1, n) = ice_fract_ncat(1, 1, n)
      else
        sice_pts_ncat(n) = 0
        sice_index_ncat(1, n) = 0
        sice_frac_ncat(1, n) = 0.0_r_um
      end if
    end do

    ! Land tile temperatures
    i_tile = 0
    tstar_land = 0.0_r_um
    do i = first_land_tile, first_land_tile + n_land_tile - 1
      i_tile = i_tile + 1
      tstar_surft(1, i_tile) = real(tile_temperature(map_tile(i)), r_um)
      tstar_land = tstar_land + frac_surft(1, i_tile) * tstar_surft(1, i_tile)
    end do

    ! Sea temperature
    tstar_sea = real(tile_temperature(map_tile(first_sea_tile)), r_um)

    ! Sea-ice temperatures
    i_sice = 0
    tstar_sice = 0.0_r_um
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      tstar_sice_sicat(1, 1, i_sice) = real(tile_temperature(map_tile(i)), r_um)
      tstar_sice = tstar_sice &
                 + ice_fract_ncat(1,1,i_sice) * tstar_sice_sicat(1,1,i_sice)
    end do

    ! Sea & Sea-ice temperature
    tstar_ssi = (1.0_r_um - ice_fract) * tstar_sea + ice_fract * tstar_sice

    ! Grid-box mean surface temperature
    tstar = flandg * tstar_land + (1.0_r_um - flandg) * tstar_ssi

    do i_pft = 1, npft
      ! Leaf area index
      lai_pft(1, i_pft) = real(leaf_area_index(map_pft(i_pft)), r_um)
      ! Canopy height
      canht_pft(1, i_pft) = real(canopy_height(map_pft(i_pft)), r_um)
    end do

    ! Roughness length (z0_tile)
    z0m_soil_gb = real(soil_roughness(map_2d(1)), r_um)
    call sparm(land_field, n_land_tile, surft_pts, surft_index,         &
               frac_surft, canht_pft, lai_pft, z0m_soil_gb,             &
               catch_snow_surft, catch_surft, z0_surft, z0h_bare_surft)

    ! Snow-free soil albedo
    albsoil_soilt = real(soil_albedo(map_2d(1)), r_um)

    do i = 1, n_soil_levs
      ! Volumetric soil moisture at wilting point (smvcwt_soilt)
      smvcwt_soilt(1, i) = real(soil_moist_wilt(map_soil(i)), r_um)
      ! Volumetric soil moisture at critical point (smvccl_soilt)
      smvccl_soilt(1, i) = real(soil_moist_crit(map_soil(i)), r_um)
      ! Volumetric soil moisture at saturation (smvcst_soilt)
      smvcst_soilt(1, i) = real(soil_moist_sat(map_soil(i)), r_um)
      ! Saturated soil conductivity (satcon_soilt)
      satcon_soilt(1, 1, i) = real(soil_cond_sat(map_soil(i)), r_um)
      ! Soil thermal capacity (hcap_soilt)
      hcap_soilt(1, 1, i) = real(soil_thermal_cap(map_soil(i)), r_um)
      ! Saturated soil water suction (sathh_soilt)
      sathh_soilt(1, 1, i) = real(soil_suction_sat(map_soil(i)), r_um)
      ! Clapp and Hornberger b coefficient (bexp_soilt)
      bexp_soilt(1, 1, i) = real(clapp_horn_b(map_soil(i)), r_um)
      ! Soil temperature (t_soil_soilt)
      t_soil_soilt(1, i) = real(soil_temperature(map_soil(i)), r_um)
      ! Soil moisture content (kg m-2, soil_layer_moisture)
      soil_layer_moisture(1, i) = real(soil_moisture(map_soil(i)), r_um)
      ! Unfrozen soil moisture proportion (sthu_soilt)
      sthu_soilt(1, i) = real(unfrozen_soil_moisture(map_soil(i)), r_um)
      ! Frozen soil moisture proportion (sthf_soilt)
      sthf_soilt(1, i) = real(frozen_soil_moisture(map_soil(i)), r_um)
    end do
    ! this needs level 0, but set equal to level 1 in all applications anyway
    satcon_soilt(1, 1, 0) = satcon_soilt(1, 1, 1)

    ! Soil thermal conductivity (hcon_soilt)
    hcon_soilt = real(soil_thermal_cond(map_2d(1)), r_um)

    ! Soil carbon content (cs_pool_gb_um)
    cs_pool_gb_um = real(soil_carbon_content(map_2d(1)), r_um)

    ! Soil ancils dependant on smvcst_soilt (soil moisture saturation limit)
    if( smvcst_soilt(1,1) > 0.0_r_um )then
      l_soil_point = .true.
    else
      l_soil_point = .false.
    end if

    ! Cosine of the solar zenith angle
    cos_zenith_angle = real(cos_zen_angle(map_2d(1)), r_um)

    ! Downwelling LW radiation at surface
    lw_down = real(lw_down_surf(map_2d(1)), r_um)

    ! Net SW radiation on tiles
    i_tile = 0
    do i = first_land_tile, first_land_tile + n_land_tile - 1
      i_tile = i_tile + 1
      sw_surft(1, i_tile) = real(sw_down_surf(map_2d(1)) - &
                                 sw_up_tile(map_tile(i)), r_um)
    end do

    ! photosynthetically active downwelling SW radiation
    photosynth_act_rad = real(sw_down_surf_blue(map_2d(1)), r_um)

    ! Ozone
    o3_gb = 0.0_r_um

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
    i_tile = 0
    do i = first_land_tile, first_land_tile + n_land_tile - 1
      i_tile = i_tile + 1
      canopy_surft(1, i_tile) = real(canopy_water(map_tile(i)), r_um)
    end do

    i_tile = 0
    do i = first_land_tile, first_land_tile + n_land_tile - 1
      i_tile = i_tile + 1
      ! Lying snow mass on land tiles
      snow_surft(1, i_tile) = real(tile_snow_mass(map_tile(i)), r_um)
      ! Number of snow layers on tiles (nsnow_surft)
      nsnow_surft(1, i_tile) = int(n_snow_layers(map_tile(i)), i_um)
      ! Snow depth on tiles (snowdepth_surft)
      snowdepth_surft(1, i_tile) = real(snow_depth(map_tile(i)), r_um)
    end do

    !-----------------------------------------------------------------------
    ! Main model prognostics - passed in from Gungho
    ! For the initial implementation we pass each individual column
    ! of data to an array sized (1,1,k) to match the UMs (i,j,k) data 
    ! layout.
    !-----------------------------------------------------------------------

    ! assuming map_wth(1) points to level 0
    ! and map_w3(1) points to level 1

    do k = 1, nlayers
      ! potential temperature on theta levels
      theta(1,1,k) = theta_in_wth(map_wth(1) + k)
      ! wet density on theta levels
      rho_wet_tq(1,1,k) = wetrho_in_wth(map_wth(1) + k)
      ! wet density on rho levels
      rho_wet(1,1,k) = wetrho_in_w3(map_w3(1) + k-1)
      ! dry density on rho levels
      rho_dry(1,1,k) = rho_in_w3(map_w3(1) + k-1)
      ! pressure on rho levels
      p_rho_levels(1,1,k) = p_zero*(exner_in_w3(map_w3(1) + k-1))**(1.0_r_def/kappa)
      ! pressure on theta levels
      p_theta_levels(1,1,k) = p_zero*(exner_in_wth(map_wth(1) + k))**(1.0_r_def/kappa)
      ! exner pressure on rho and theta levels
      exner_rho_levels(1,1,k) = exner_in_w3(map_w3(1) + k-1)
      exner_theta_levels(1,1,k) = exner_in_wth(map_wth(1) + k)
      ! u wind on rho levels
      u_p(1,1,k) = u1_in_w3(map_w3(1) + k-1)
      ! v wind on rho levels
      v_p(1,1,k) = u2_in_w3(map_w3(1) + k-1)
      ! w wind on theta levels
      w(1,1,k) = u3_in_wth(map_wth(1) + k)
      ! height of rho levels from centre of planet
      r_rho_levels(1,1,k) = height_w3(map_w3(1) + k-1)                  &
                          + planet_radius
      ! height of theta levels from centre of planet
      r_theta_levels(1,1,k) = height_wth(map_wth(1) + k)                &
                            + planet_radius
      ! water vapour mixing ratio
      q(1,1,k) = m_v_n(map_wth(1) + k)
      ! cloud liquid mixing ratio
      qcl(1,1,k) = m_cl_n(map_wth(1) + k)
      ! cloud ice mixing ratio
      qcf(1,1,k) = m_ci_n(map_wth(1) + k)
      ! cloud fraction variables
      bulk_cloud_fraction(1,1,k) = cf_bulk(map_wth(1) + k)
      area_cloud_fraction(1,1,k) = cf_area(map_wth(1) + k)
      cloud_fraction_liquid(1,1,k) = cf_liq(map_wth(1) + k)
      cloud_fraction_frozen(1,1,k) = cf_ice(map_wth(1) + k)
      ! 3D RH_crit field
      rhcpt(1,1,k) = rh_crit_wth(map_wth(1) + k) 
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
    ! setup odd array which is on rho levels but without level 1
    p_rho_minus_one(1,1,0) = p_theta_levels(1,1,0)
    p_rho_minus_one(1,1,1:nlayers-1) = p_rho_levels(1,1,2:nlayers)
    p_rho_minus_one(1,1,nlayers) = 0.0_r_um
    ! near surface potential temperature
    theta(1,1,0) = theta_in_wth(map_wth(1) + 0)
    ! wet density multiplied by planet radius squared on rho levs
    rho_wet_rsq = rho_wet * r_rho_levels**2
    ! extended halo u and v winds
    u_px = u_p
    v_px = v_p
    ! u and v winds on their native grid
    u = u_p
    v = v_p
    ! surface currents
    u_0_px = 0.0
    v_0_px = 0.0
    u_0 = 0.0
    v_0 = 0.0
    u_0_p = 0.0
    v_0_p = 0.0
    ! near surface moisture fields
    q(1,1,0) = m_v_n(map_wth(1) + 0)
    qcl(1,1,0) = m_cl_n(map_wth(1) + 0)
    qcf(1,1,0) = m_ci_n(map_wth(1) + 0)
    ! surface height
    r_theta_levels(1,1,0) = height_wth(map_wth(1) + 0) + planet_radius
    ! eta space, 0-1 scaled height
    eta_theta_levels(:) = (r_theta_levels(1,1,:)-r_theta_levels(1,1,0))  &
                        /(r_theta_levels(1,1,nlayers)-planet_radius)
    ! height of levels above surface
    z_rho = r_rho_levels-r_theta_levels(1,1,0)
    z_theta(1,1,:) = r_theta_levels(1,1,1:nlayers)-r_theta_levels(1,1,0)
    ! vertical velocity
    w(1,1,0) = u3_in_wth(map_wth(1) + 0)
    w_copy = w
    etadot_copy = w_copy / r_theta_levels(1,1,nlayers)

    if (i_convection_vn == i_convection_vn_6a ) then
      do k = 1, nlayers-1
        exner_rho_minus_one(1,1,k) = exner_rho_levels(1,1,k+1)
      end do 

      exner_rho_minus_one(1,1,0)   = exner_theta_levels(1,1,0)
      exner_rho_minus_one(1,1,nlayers) = 0.0_r_um
    end if

    !-----------------------------------------------------------------------
    ! Things saved from one timestep to the next
    !-----------------------------------------------------------------------
    ! previous BL height
    zh(1,1) = zh_2d(map_2d(1))
    zh_prev = zh
    ! surface roughness
    z0msea(1,1) = z0msea_2d(map_2d(1))
    ! downdraft at cloud base
    ddmfx = 0.0

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

    ! Use  convection switches to decide the value of  L_cape_opt
    IF (i_convection_vn == i_convection_vn_6a ) THEN
      L_cape_opt = ( (cldbase_opt_dp == 3) .OR. (cldbase_opt_md == 3) .OR. &
                     (cldbase_opt_dp == 4) .OR. (cldbase_opt_md == 4) .OR. &
                     (cldbase_opt_dp == 5) .OR. (cldbase_opt_md == 5) .OR. &
                     (cldbase_opt_dp == 6) .OR. (cldbase_opt_md == 6) )
    ELSE
      L_cape_opt = .false.
    END IF
    CALL ap2_init_conv_diag( rows, row_length, ntml, ntpar, nlcl, cumulus,  &
        l_shallow, l_mid, delthvu, ql_ad, zhpar, dzh, qcl_inv_top,          &
        zlcl, zlcl_uv, conv_type, no_cumulus, w_max, w_copy, L_cape_opt)

    CALL conv_diag_6a(                                                  &

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
          , w_copy, w_max, deep_flag, past_precip, past_conv_ht         &
          , conv_prog_precip                                            &
          , g_ccp, h_ccp                                                &
    !
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
            )

    CALL NI_bl_ctl (                                                    &
    !     IN parameters for SISL scheme
         CycleNo, l_jules_call,                                         &
    !     IN time stepping information
         curr_year, curr_day_number, curr_hour, curr_minute, curr_second,&
    !     IN switches
         L_aero_classic,                                                &
    !     IN data fields.
         p_rho_levels, p_theta_levels, rho_wet_rsq,rho_wet,rho_dry, u_p, v_p,&
         u_px, v_px, u_0_px, v_0_px,                                    &
         land_sea_mask, q, qcl, qcf, p_star, theta, exner_theta_levels, rad_hr,&
         micro_tends, soil_layer_moisture, rho_wet_tq, z_rho, z_theta,  &
    !     IN ancillary fields and fields needed to be kept from tstep to tstep
         hcon_soilt, smvccl_soilt, smvcwt_soilt, smvcst_soilt,          &
         sthf_soilt, sthu_soilt, sil_orog_land_gb,                      &
    !-------------------------------------------------------------------------
         ho2r2_orog_gb, sd_orog, ice_fract_ncat_sicat, k_sice_sicat,    &
         land_index, photosynth_act_rad,                                &
         soil_clay,soil_sand,dust_mrel1,dust_mrel2,                     &
         dust_mrel3,dust_mrel4,dust_mrel5,dust_mrel6,                   &
    !     IN additional variables for JULES
         canopy_surft, catch_surft, catch_snow_surft, snow_surft,       &
         z0_surft, z0h_bare_surft,                                      &
         z0m_soil_gb, lw_down, sw_surft, tstar_surft, tsurf_elev_surft, &
         co2(1:co2_dim_len,1:co2_dim_row,1),                            &
         asteps_since_triffid,                                          &
         cs_pool_gb_um,frac_surft,canht_pft,lai_pft,fland,flandg,       &
         albsoil_soilt, cos_zenith_angle,                               &
    !     IN: input from the wave model
         charnock_w,                                                    &
    !     IN everything not covered so far
         t_soil_soilt, ti_sicat,                                        &
         ti_cat_sicat,tstar,zh_prev,ddmfx,bulk_cloud_fraction,zhpar,zlcl, &
    !     IN SCM namelist data
         L_spec_z0, z0m_scm, z0h_scm, flux_e, flux_h, ustar_in,         &
    !     SCM diagnostics and STASH
         nSCMDpkgs, L_SCMDiags, BL_diag, sf_diag,                       &
    !     INOUT data
         gs_gb,z0msea,w_copy,etadot_copy,tstar_sea,tstar_sice_sicat,zh,dzh,&
         cumulus, ntml,ntpar,l_shallow,                                 &
    !     INOUT additional variables for JULES
         g_leaf_acc_pft,npp_acc_pft,resp_w_acc_pft,resp_s_acc_gb_um,    &
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
        nbdsc,ntdsc,wstar,wthvs,uw0,vw0,taux_p,tauy_p,rhcpt             &
     )

    !-----------------------------------------------------------------------
    ! re-grid stuff onto u/v points
    ! The UM requires this regridding.  This will be dealt with 
    ! differently in LFRic (see #1059)
    !-----------------------------------------------------------------------
    flandg_u=flandg
    flandg_v=flandg
    rhokm_u=rhokm
    rhokm_v=rhokm
    rhokm_u_ssi=rhokm_ssi
    rhokm_v_ssi=rhokm_ssi
    rhokm_u_land=rhokm_land
    rhokm_v_land=rhokm_land
    f_ngstress_u=f_ngstress
    f_ngstress_v=f_ngstress
    taux_fd_u=tau_fd_x
    tauy_fd_v=tau_fd_y

    l_calc_at_p = .FALSE.  ! Flag for calling this routine on p-grid
      CALL bdy_expl3 (                                                 &
      ! IN grid related variables
           bl_levels, l_calc_at_p,                                     &
           udims, vdims, udims_s, vdims_s, pdims,                      &
           udims, vdims,                                               &
      ! IN SCM diags
           nSCMDpkgs,L_SCMDiags,                                       &
      ! IN variables used in flux calculations
           u, v, u_0, v_0, rhokm_u_land, rhokm_v_land, flandfac_u,     &
           flandfac_v, rhokm_u_ssi, rhokm_v_ssi, fseafac_u, fseafac_v, &
           flandg_u, flandg_v, zhnl, rdz_u, rdz_v, rhokm_u, rhokm_v,   &
           taux_fd_u, tauy_fd_v,                                       &
           rhogamu_u, rhogamv_v, f_ngstress_u, f_ngstress_v,           &
      ! OUT explicit momentum fluxes
           taux_land, tauy_land, taux_ssi, tauy_ssi, taux, tauy        &
                   )

    !-----------------------------------------------------------------------
    ! Needed by ni_imp_ctl whether convection called or not
    !-----------------------------------------------------------------------
    cca_3d = 0.0_r_um
    ccw_3d = 0.0_r_um
    lcbase(1,1) = 0     ! default for no convection and convective cloud

    !========================================================================
    ! Call to 6A Gregory-Rowntree convection scheme 
    !========================================================================

    if (i_convection_vn == i_convection_vn_6a ) then

      ! Current assumptions about setup based on GA7

      l_tracer = .false.   ! not allowing tracers as none so far
      l_calc_dxek  = ( i_cld_vn == i_cld_pc2 )
      l_q_interact = l_calc_dxek 

      segments = 1          ! i.e. one column
      seg_num = map_wth(1)  ! Only used by debugging error message from UM
                            ! This is probably most useful to indicate
                            ! which column has a problem

      n_conv_levels = tdims%k_end
      if (l_mom) then
        ! Limit convection calling levels to maximum of model_levels - 1
        ! This is because CMT increments to u and v exist on n_levels + 1
        if (n_conv_levels  >   nlayers - 1 ) then
          n_conv_levels = nlayers - 1
        end if
      end if

      timestep_conv=timestep/real(n_conv_calls)

      ! Sub-timestep scheme

      one_over_conv_calls = 1.0_r_um/real(n_conv_calls)

      l_mid(1,1) = .true.

      do k=1,tdims%k_end
        theta_conv(1,1,k) = theta_star(map_wth(1) + k)
        ! Pointing to _star values
        q_conv(1,1,k)   = m_v(map_wth(1) + k)
        qcl_conv(1,1,k) = m_cl(map_wth(1) + k)
        qcf_conv(1,1,k) = m_ci(map_wth(1) + k)

        cf_liquid_conv(1,1,k) = cf_liq(map_wth(1) + k)
        cf_frozen_conv(1,1,k) = cf_ice(map_wth(1) + k)
        bulk_cf_conv(1,1,k)   = cf_bulk(map_wth(1) + k)
        if (l_mcr_qrain) then 
          qrain_conv(1,1,k) = m_r(map_wth(1) + k)
        end if
        if (l_mcr_qgraup) then 
          qgraup_conv(1,1,k) = m_g(map_wth(1) + k)
        end if

        ! Note need to pass down qcf2 
        qcf2_conv(1,1,k)= 0.0_r_um
      end do

      ccb(1,1) = 0_i_um
      cct(1,1) = 0_i_um
      lctop(1,1) = 0_i_um

      ! Check for negative (less than a minimum) q being passed to convection
      if (l_safe_conv) then
        do k=1,tdims%k_end
          if (q_conv(1,1,k) < qmin_conv) then
            ! Should really be warning print statements if this is happening
            ! but log_event calls not allowed.
            ! store q added to ensure sensible profile
            dq_add(1,1,k) = qmin_conv - q_conv(1,1,k)
            ! Reset to qmin in non-conservative way
            q_conv(1,1,k) = qmin_conv
          else
            dq_add(1,1,k) = 0.0
          end if
        end do ! k
      end if
      if (l_mom) then
        do k=1,tdims%k_end
          u_conv(1,1,k) = u1_star(map_w3(1) + k-1)
          v_conv(1,1,k) = u2_star(map_w3(1) + k-1)
        end do ! k
      end if

      if ( cycleno == outer_iterations .AND. l_tracer ) then
        ! Need tracers no code for this at present
        ntra_fld = 1             ! can't have 0 sized arrays
        ntra_lev = 1
        ! allocate tot_tracer
        allocate( tot_tracer(1, 1, ntra_lev, ntra_fld) )
      else
        ntra_fld = 1             ! can't have 0 sized arrays
        ntra_lev = 1
        ! allocate dummy space in tot_tracer
        allocate( tot_tracer(1, 1, ntra_lev, ntra_fld) )
      end if

      l_scm_convss_dg = .false.  ! not allowing sub-timestep SCM diags

      n_congestus = 0  
      n_deep = 0
      n_shallow = 0

      if (cumulus(1,1)) then
        n_cumulus = 1
        if (iconv_deep >  0  .AND. .NOT. l_shallow(1,1) ) then
          n_deep = 1
        endif
        if (iconv_shallow >  0  .AND. l_shallow(1,1) ) then
          n_shallow = 1
        endif
        n_congestus = 1   ! as UM though not actually using scheme
      else
        n_cumulus = 0
      end if

      ! Loop over convection calls per model time step
      do call_number = 1, n_conv_calls

        if (l_mid(1,1)) then
          n_mid = 1
        else
          n_mid = 0
        end if

        ! Initialise convection work arrays holding information from each
        ! iteration (sub-step)
        it_cca(1,1,:n_cca_lev)  = 0.0_r_um
        it_cca0(1,1,:n_cca_lev) = 0.0_r_um
        it_cca0_dp(1,1,:n_cca_lev) = 0.0_r_um
        it_cca0_sh(1,1,:n_cca_lev) = 0.0_r_um
        it_cca0_md(1,1,:n_cca_lev) = 0.0_r_um
      
        it_ccw(1,1,:)  = 0.0_r_um
        it_ccw0(1,1,:) = 0.0_r_um
        it_conv_rain_3d(1,1,:) = 0.0_r_um
        it_conv_snow_3d(1,1,:) = 0.0_r_um
        it_w2p(1,1,:) = 0.0_r_um
        it_up_flux(1,1,:) = 0.0_r_um

        it_lcca(1,1)   = 0.0_r_um
        it_lcbase(1,1) = 0
        it_lctop(1,1)  = 0

        it_ccb(1,1)    = 0
        it_cct(1,1)    = 0
        it_cca_2d(1,1) = 0.0_r_um
        it_cclwp(1,1)  = 0.0_r_um

        it_ccb0(1,1)   = 0
        it_cct0(1,1)   = 0
        it_cclwp0(1,1) = 0.0_r_um

        it_conv_rain(1,1) = 0.0_r_um
        it_conv_snow(1,1) = 0.0_r_um
        it_precip_dp(1,1) = 0.0_r_um
        it_precip_sh(1,1) = 0.0_r_um
        it_precip_md(1,1) = 0.0_r_um
        it_cape_out(1,1)  = 0.0_r_um
        it_kterm_deep(1,1)  = 0
        it_kterm_shall(1,1) = 0
        it_mid_level(1,1) = .FALSE.
        it_dp_cfl_limited(1,1) = 0.0_r_um
        it_md_cfl_limited(1,1) = 0.0_r_um

        it_precip_cg(1,1) = 0.0_r_um
        it_wstar_up(1,1)  = 0.0_r_um
        it_mb1(1,1) = 0.0_r_um
        it_mb2(1,1) = 0.0_r_um
        it_cg_term(1,1) = 0


        call glue_conv_6a                                                     &
          ( rows*row_length, segments, n_conv_levels, bl_levels, call_number, &
           seg_num, theta_conv, q_conv, qcl_conv, qcf_conv                    &
          , qrain_conv, qgraup_conv, qcf2_conv                                &
          , cf_liquid_conv, cf_frozen_conv, bulk_cf_conv                      &
          , p_star, land_sea_mask                                             &
          , u_conv, v_conv, w(1,1,1)                                          &
          , tot_tracer, dthbydt, dqbydt,   dqclbydt, dqcfbydt                 &
          , dcflbydt, dcffbydt, dbcfbydt, dubydt_p, dvbydt_p                  &
          , it_conv_rain, it_conv_snow, it_conv_rain_3d, it_conv_snow_3d      &
          , it_cca0_dp, it_cca0_md, it_cca0_sh                                &
          , it_cca0,  it_ccb0, it_cct0, it_cclwp0, it_ccw0, it_lcbase0        &
          , it_lctop,  it_lcca                                                &
          , it_cca,   it_ccb,  it_cct,  it_cclwp,  it_ccw,  it_lcbase         &
          , it_cca_2d, freeze_lev, it_dp_cfl_limited, it_md_cfl_limited       &
          , it_mid_level, it_kterm_deep, it_kterm_shall                       &
          , it_precip_dp, it_precip_sh, it_precip_md, it_precip_cg            &
          , it_wstar_dn,  it_wstar_up                                         &
          , it_mb1, it_mb2, it_cg_term, n_cumulus                             &
          , uw0, vw0, w_max                                                   &
          , zlcl, zlcl_uv, zhpar, entrain_coef                                &
          , conv_prog_precip, deep_flag, past_precip, past_conv_ht            &
          , it_cape_out, n_deep, n_congestus, n_shallow                       &
          , n_mid, r_rho_levels(1,1,1), r_theta_levels(1,1,1)                 &
          , rho_wet, rho_wet_tq, rho_dry, rho_dry_theta, delta_smag           &
          , exner_rho_levels, exner_rho_minus_one, exner_theta_levels         &
          , p_rho_minus_one, p_theta_levels                                   &
          , z_theta, z_rho, timestep_conv                                     &
          , t1_sd, q1_sd, ntml, ntpar                                         &
          , conv_type, l_shallow                                              &
          , l_congestus, l_mid, cumulus                                       &
          , wstar, wthvs, delthvu, ql_ad, qsat_lcl, ftl, fqw                  &
          , l_tracer, ntra_fld, ntra_lev, n_cca_lev, l_mcr_qrain              &
          , l_mcr_qgraup, l_mcr_qcf2, l_calc_dxek , l_q_interact              &
          , it_up_flux_half, it_up_flux,      it_dwn_flux                     &
          , it_entrain_up,   it_detrain_up, it_entrain_dwn,  it_detrain_dwn   &
          , it_uw_dp,        it_vw_dp                                         &
          , it_uw_shall,     it_vw_shall, it_uw_mid,  it_vw_mid               &
          , it_wqt_flux,  it_wthetal_flux, it_wthetav_flux, it_wql_flux       &
          , it_mf_deep,      it_mf_congest, it_mf_shall,  it_mf_midlev        &
          , it_dt_deep,      it_dt_congest, it_dt_shall,  it_dt_midlev        &
          , it_dq_deep,      it_dq_congest, it_dq_shall,  it_dq_midlev        &
          , it_du_deep,      it_du_congest, it_du_shall,  it_du_midlev        &
          , it_dv_deep,      it_dv_congest, it_dv_shall,  it_dv_midlev        &
          , ind_cape_reduced,  cape_ts_used, it_ind_deep, it_ind_shall        &
          , it_w2p, it_dt_dd, it_dq_dd, it_du_dd, it_dv_dd, it_area_ud        &
          , it_area_dd, scm_convss_dg, l_scm_convss_dg                        &
           )

        ! Mid-level convection only possible on subsequent sub-steps if 
        ! occurs on first step or column is diagnosed as cumulus
        l_mid(1,1) = it_mid_level(1,1) .or. cumulus(1,1)

        ! Update cloud info from substep
        ! Highest convective layer properties - diagnostic
        !------------------------------------
        ! max cct across total number of calls to convection
        ! Note that diagnostic is a real not an integer
        cct(1,1) = max( cct(1,1),it_cct(1,1))
        ! min ccb across total number of calls to convection
        ! excluding ccb=0
        if (cv_base(map_2d(1)) > 0 .AND. it_ccb(1,1) > 0) then
          ccb(1,1) = min(ccb(1,1),it_ccb(1,1))
        else
          ccb(1,1) = max(ccb(1,1),it_ccb(1,1))
        end if

        ! Lowest convective layer properties
        !------------------------------------
        ! max lctop across total number of calls to convection
        lctop(1,1) = max(lctop(1,1),it_lctop(1,1))


        ! min lcbase across total number of calls to convection
        ! excluding lcbase=0
        if (lcbase(1,1) > 0 .AND. it_lcbase(1,1) > 0) then
          lcbase(1,1) = min(lcbase(1,1), it_lcbase(1,1))
        else
          lcbase(1,1) = max(lcbase(1,1), it_lcbase(1,1))
        end if

        do k=1, model_levels
          ccw_3d(1,1,k) = ccw_3d(1,1,k) + one_over_conv_calls*it_ccw0(1,1,k)
          cca_3d(1,1,k) = cca_3d(1,1,k) + one_over_conv_calls*it_cca0(1,1,k)
          ! Assuming lccrad = .true.
          cca_3d(1,1,k) = min(cca_3d(1,1,k), 1.0)
        end do

        ! single level convection diagnostics
        freeze_level(map_2d(1)) = freeze_level(map_2d(1)) +                   &
                                    real(freeze_lev(1,1)) *one_over_conv_calls

        if (it_mid_level(1,1)) then
          mid_in_col(map_2d(1)) = mid_in_col(map_2d(1)) + one_over_conv_calls
        end if 
        deep_in_col(map_2d(1)) = deep_in_col(map_2d(1)) +                     &
                                       it_ind_deep(1,1)*one_over_conv_calls
        shallow_in_col(map_2d(1)) = shallow_in_col(map_2d(1)) +               &
                                       it_ind_shall(1,1)*one_over_conv_calls

        deep_term(map_2d(1)) = deep_term(map_2d(1))+                          &
                                   real(it_kterm_deep(1,1)) *one_over_conv_calls

        cape_timescale(map_2d(1)) =cape_timescale(map_2d(1))     +            &
                                      cape_ts_used(1,1) *one_over_conv_calls

        conv_rain(map_2d(1)) = conv_rain(map_2d(1))+                          &
                                     it_conv_rain(1,1) *one_over_conv_calls
        conv_snow(map_2d(1)) = conv_snow(map_2d(1))+                          &
                                     it_conv_snow(1,1) *one_over_conv_calls
        deep_prec(map_2d(1)) = deep_prec(map_2d(1))+                          &
                                     it_precip_dp(1,1) *one_over_conv_calls
        shallow_prec(map_2d(1)) = shallow_prec(map_2d(1))+                    &
                                     it_precip_sh(1,1) *one_over_conv_calls
        mid_prec(map_2d(1)) = mid_prec(map_2d(1))+                            &
                                     it_precip_md(1,1) *one_over_conv_calls

        ! Frequency of deep convection terminating on level k
        if (it_ind_deep(1,1) == 1.0_r_um) then
          k = it_kterm_deep(1,1)
          if (k > 0) then  ! in case still get a zero value
            deep_tops(map_wth(1)+k) = deep_tops(map_wth(1)+k) + one_over_conv_calls
          end if
        end if

        ! update input fields *_conv for next substep 
        do k = 1, n_conv_levels
          theta_conv(1,1,k) = theta_conv(1,1,k)                         &
                                + dthbydt(1,1,k) * timestep_conv
          q_conv(1,1,k)     = q_conv(1,1,k)                             &
                                + dqbydt(1,1,k) * timestep_conv
          qcl_conv(1,1,k)   = qcl_conv(1,1,k)                           &
                                +(dqclbydt(1,1,k) * timestep_conv)
          qcf_conv(1,1,k)   = qcf_conv(1,1,k)                           &
                                +(dqcfbydt(1,1,k) * timestep_conv)
          cf_liquid_conv(1,1,k) = cf_liquid_conv(1,1,k)                 &
                                    +(dcflbydt(1,1,k) * timestep_conv)
          cf_frozen_conv(1,1,k) = cf_frozen_conv(1,1,k)                 &
                                    + (dcffbydt(1,1,k) * timestep_conv)
          bulk_cf_conv(1,1,k)   = bulk_cf_conv(1,1,k)                   &
                                    +(dbcfbydt(1,1,k) * timestep_conv)
          dt_conv(map_wth(1) + k)   = dt_conv(map_wth(1) + k)               &
                                       + dthbydt(1,1,k) * timestep_conv
          dmv_conv(map_wth(1) + k)  =  dmv_conv(map_wth(1) + k)             &
                                       + dqbydt(1,1,k) * timestep_conv
          dmcl_conv(map_wth(1) + k) =  dmcl_conv(map_wth(1) + k)            &
                                       + dqclbydt(1,1,k) * timestep_conv
          dmcf_conv(map_wth(1) + k) =  dmcf_conv(map_wth(1) + k)            &
                                       + dqcfbydt(1,1,k) * timestep_conv

          ! Update diagnostics   
          massflux_up(map_wth(1) + k) = massflux_up(map_wth(1) + k) +          &
                                        it_up_flux(1,1,k)*one_over_conv_calls
          massflux_down(map_wth(1) + k) = massflux_down(map_wth(1) + k)+       &
                                        it_dwn_flux(1,1,k)*one_over_conv_calls
          entrain_up(map_wth(1) + k) = entrain_up(map_wth(1) + k) +            &
                                        it_entrain_up(1,1,k)*one_over_conv_calls
          entrain_down(map_wth(1) + k) = entrain_down(map_wth(1) + k) +        &
                                        it_entrain_dwn(1,1,k)*one_over_conv_calls
          detrain_up(map_wth(1) + k) = detrain_up(map_wth(1) + k) +            &
                                        it_detrain_up(1,1,k)*one_over_conv_calls
          detrain_down(map_wth(1) + k) = detrain_down(map_wth(1) + k) +        &
                                        it_detrain_dwn(1,1,k)*one_over_conv_calls
          dd_dt(map_wth(1) + k) = dd_dt(map_wth(1) + k) +                      &
                                        it_dt_dd(1,1,k)*one_over_conv_calls
          dd_dq(map_wth(1) + k) = dd_dq(map_wth(1) + k) +                      &
                                        it_dq_dd(1,1,k)*one_over_conv_calls

          deep_massflux(map_wth(1) + k) = deep_massflux(map_wth(1) + k) +      &
                                        it_mf_deep(1,1,k)*one_over_conv_calls
          deep_dt(map_wth(1) + k) = deep_dt(map_wth(1) + k) +                  &
                                        it_dt_deep(1,1,k)*one_over_conv_calls
          deep_dq(map_wth(1) + k) = deep_dq(map_wth(1) + k) +                  &
                                        it_dq_deep(1,1,k)*one_over_conv_calls
          shallow_massflux(map_wth(1) + k) = shallow_massflux(map_wth(1) + k) +&
                                        it_mf_shall(1,1,k)*one_over_conv_calls
          shallow_dt(map_wth(1) + k) = shallow_dt(map_wth(1) + k) +            &
                                        it_dt_shall(1,1,k)*one_over_conv_calls
          shallow_dq(map_wth(1) + k) = shallow_dq(map_wth(1) + k) +            &
                                        it_dq_shall(1,1,k)*one_over_conv_calls
          mid_massflux(map_wth(1) + k) = mid_massflux(map_wth(1) + k) +        &
                                        it_mf_midlev(1,1,k)*one_over_conv_calls
          mid_dt(map_wth(1) + k) = mid_dt(map_wth(1) + k) +                    &
                                        it_dt_midlev(1,1,k)*one_over_conv_calls
          mid_dq(map_wth(1) + k) = mid_dq(map_wth(1) + k) +                    &
                                        it_dq_midlev(1,1,k)*one_over_conv_calls
        end do

        if (l_mom) then
          do k = 1, n_conv_levels
            u_conv(1,1,k) = u_conv(1,1,k)   + dubydt_p(1,1,k) * timestep_conv
            v_conv(1,1,k) = v_conv(1,1,k)   + dvbydt_p(1,1,k) * timestep_conv
            ! total increments
            du_conv(map_w3(1) + k -1) = du_conv(map_w3(1) + k -1) + dubydt_p(1,1,k) * timestep_conv
            dv_conv(map_w3(1) + k -1) = dv_conv(map_w3(1) + k -1) + dvbydt_p(1,1,k) * timestep_conv
          end do
        end if    !l_mom

        ! Would have PC2 checks here

      end do    ! loop over calls to convection

      if (l_safe_conv) then
        do k = 1, n_conv_levels
          dmv_conv(map_wth(1) + k) = dmv_conv(map_wth(1) + k) + dq_add(1,1,k)
        end do 
      end if 

      ! Convection/PC2 checks
      ! Protect against generation of inconsistently low cloud
      ! fraction implying very high in-cloud condensate amounts.
      ! In-cloud condensate amounts above 2.0e-3 lead to
      ! cloud fraction being increased (up to a value of 1.0)

      do k = 1, n_conv_levels
        ! Liquid cloud fraction
        if (cf_liquid_conv(1,1,k) > 0.0_r_um) then
          if ( (qcl_conv(1,1,k)/cf_liquid_conv(1,1,k) ) > 2.0e-3_r_um ) then
            orig_value = cf_liquid_conv(1,1,k)
            cf_liquid_conv(1,1,k) = min(1.0,qcl_conv(1,1,k)/2.0e-3_r_um)
            cf_liquid_inc(1,1,k)  = cf_liquid_inc(1,1,k)                   &
                                    + cf_liquid_conv(1,1,k) - orig_value
            bulk_cf_conv(1,1,k) = bulk_cf_conv(1,1,k)                      &
                                    + cf_liquid_conv(1,1,k) - orig_value
            bulk_cf_inc(1,1,k) = bulk_cf_inc(1,1,k)                        &
                                    + cf_liquid_conv(1,1,k) - orig_value
          end if
        end if

        ! Ice cloud fraction
        if (cf_frozen_conv(1,1,k) > 0.0_r_um) then
          if ( (qcf_conv(1,1,k)/cf_frozen_conv(1,1,k)) > 2.0e-3_r_um ) then
            orig_value = cf_frozen_conv(1,1,k)
            cf_frozen_conv(1,1,k) = min(1.0,qcf_conv(1,1,k)/2.0e-3_r_um)
            cf_frozen_inc(1,1,k)  = cf_frozen_inc(1,1,k)                  &
                                      + cf_frozen_conv(1,1,k) - orig_value
            bulk_cf_conv(1,1,k) = bulk_cf_conv(1,1,k)                     &
                                      + cf_frozen_conv(1,1,k) - orig_value
            bulk_cf_inc(1,1,k) = bulk_cf_inc(1,1,k)                       &
                                      + cf_frozen_conv(1,1,k) - orig_value
          end if
        end if
      end do

      ! Store convective downdraught mass fluxes at cloud base
      ! if required for surface exchange.
      if (isrfexcnvgust == ip_srfexwithcnv) then
        if (ccb(1,1) > 0) then
          ddmfx(1,1)=massflux_down( map_wth(1) + ccb(1,1))
        else
          ddmfx(1,1)=0.0
        end if
      end if

      ! copy convective cloud fraction into prognostic array
      do k = 1, n_conv_levels
        cca(map_wth(1) + k) =  min(cca_3d(1,1,k), 1.0)
        ccw(map_wth(1) + k) =  ccw_3d(1,1,k)
      end do

      ! Copy integers into real diagnostic arrays
      cv_top(map_2d(1))        = real(cct(1,1))
      cv_base(map_2d(1))       = real(ccb(1,1))
      lowest_cv_top(map_2d(1)) = real(lctop(1,1))

      ! Would need to copy tracers back!
      deallocate(tot_tracer)

    end if
    !========================================================================
    ! End of convection
    !========================================================================

    !-----------------------------------------------------------------------
    ! increments / fields with increments added
    !-----------------------------------------------------------------------
    do k = 1, nlayers
      t_latest(1,1,k) = theta_star(map_wth(1) + k) * exner_theta_levels(1,1,k) &
                      + dt_conv(map_wth(1) + k)
      r_u(1,1,k) = u1_star(map_w3(1) + k-1) - u1_in_w3(map_w3(1) + k-1) &
                        + du_conv(map_w3(1) + k-1)
      r_v(1,1,k) = u2_star(map_w3(1) + k-1) - u2_in_w3(map_w3(1) + k-1) &
                        + dv_conv(map_w3(1) + k-1)
      r_w(1,1,k) = u3_star(map_wth(1) + k) - u3_in_wth(map_wth(1) + k)
      q_latest(1,1,k)   = m_v(map_wth(1) + k)  + dmv_conv(map_wth(1) + k)
      qcl_latest(1,1,k) = m_cl(map_wth(1) + k)
      qcf_latest(1,1,k) = m_ci(map_wth(1) + k)
! When PC2 working
!      qcl_latest(1,1,k) = m_cl(map_wth(1) + k)  + dmcl_conv(map_wth(1) + k)
!      qcf_latest(1,1,k) = m_ci(map_wth(1) + k)  + dmcf_conv(map_wth(1) + k)
      ! diagnostic increments
      dt_bl(map_wth(1)+k) = t_latest(1,1,k)
      dmv_bl(map_wth(1)+k) = q_latest(1,1,k)
    end do
    r_w(1,1,0) = u3_star(map_wth(1) + 0) - u3_in_wth(map_wth(1) + 0)

    conv_rain_copy(1,1) = conv_rain(map_2d(1))
    conv_snow_copy(1,1) = conv_snow(map_2d(1))

    CALL NI_imp_ctl (                                                   &
    ! IN Model switches
            CycleNo                                                     &
    ! IN trig arrays
          , xx_cos_theta_latitude                                       &
    ! IN data fields.
          , p_theta_levels, p_rho_minus_one, rho_wet_rsq, rho_wet_tq    &
          , u, v, w                                                     &
          , land_sea_mask, q, qcl, qcf, p_star, theta, exner_theta_levels&
    ! IN ancillary fields and fields needed to be kept from tstep to tstep
          , sil_orog_land_gb, ho2r2_orog_gb                             &
          , ice_fract, di_ncat_sicat, ice_fract_ncat, k_sice_sicat      &
          , u_0, v_0, land_index, cca, lcbase, ccb0, cct0               &
          , ls_rain, ls_snow, conv_rain_copy, conv_snow_copy            &
    ! IN variables required from BDY_LAYR
          , alpha1_sea, alpha1_sice, ashtf_prime_sea, ashtf_prime, bq_gb, bt_gb&
          , dtrdz_charney_grid, rdz_charney_grid, dtrdz_u, dtrdz_v      &
          , rdz_u, rdz_v, cdr10m_u, cdr10m_v, z_theta                   &
          , k_blend_tq, k_blend_uv, u_s, rhokm, rhokm_u, rhokm_v        &
    ! IN diagnostics (started or from) BDY_LAYR
          , rib_gb,zlcl, zhnl, dzh, qcl_inv_top, zh                     &
          , bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6 &
          , bl_type_7, z0m_eff_gb, ntml, cumulus                        &
    ! IN data required for tracer mixing :
          , rho_aresist,aresist,r_b_dust                                &
          , kent, we_lim, t_frac, zrzi                                  &
          , kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                  &
          , zhsc,z_rho,dust_flux,dust_emiss_frac                        &
          , u_s_t_tile,u_s_t_dry_tile,u_s_std_surft                     &
     ! IN additional variables for JULES. Now includes lai_pft, canht_pft.
          , surft_pts,surft_index,tile_frac,canopy_surft                &
          , alpha1,fraca,rhokh_surft,smc_soilt,chr1p5m,resfs,z0hssi,z0mssi &
          , canhc_surft,flake,wt_ext_surft,lw_down,lai_pft,canht_pft    &
          , sw_surft,ashtf_prime_surft,gc_surft,aresist_surft           &
          , resft,rhokh_sice,rhokh_sea,z0h_surft,z0m_surft              &
          , chr1p5m_sice                                                &
          , fland, flandg, flandg_u,flandg_v                            &
          , emis_surft, t_soil_soilt, snow_surft, sstfrz                &
    ! IN JULES variables for STASH
          , gs_gb, npp_gb, resp_s_gb_um                                 &
          , resp_s_tot_soilt,cs_pool_gb_um                              &
          , catch_surft                                                 &
          , co2_emits, co2flux                                          &
    ! INOUT diagnostic info
          , STASHwork3, STASHwork9                                      &
    ! SCM Diagnostics (dummy in full UM) & bl diags
          , nSCMDpkgs, L_SCMDiags, bl_diag, sf_diag                     &
    ! INOUT (Note ti_cat_sicat and ti_sicat are IN only if l_sice_multilayers=T)
          , TScrnDcl_SSI, TScrnDcl_surft, tStbTrans                     &
          , cca0, fqw, ftl, taux, tauy, rhokh                           &
          , fqw_ice,ftl_ice,dtstar_surft,dtstar_sea,dtstar_sice,ti_cat_sicat &
          , area_cloud_fraction, bulk_cloud_fraction                    &
          , t_latest, q_latest, qcl_latest, qcf_latest                  &
          , cf_latest, cfl_latest, cff_latest                           &
          , R_u, R_v, R_w, cloud_fraction_liquid, cloud_fraction_frozen &
          , sum_eng_fluxes,sum_moist_flux, rhcpt                        &
    ! INOUT tracer fields
          , aerosol, free_tracers,  resist_b,  resist_b_surft           &
          , dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6 &
          , drydep2, so2, dms, so4_aitken, so4_accu, so4_diss, nh3      &
          , soot_new, soot_aged, soot_cld, bmass_new, bmass_aged        &
          , bmass_cld, ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss &
          , co2, ozone_tracer                                           &
    ! INOUT additional variables for JULES
          , tstar_surft, fqw_surft, epot_surft, ftl_surft               &
          , radnet_sice,olr,tstar_sice_sicat,tstar_ssi                  &
          , tstar_sea,taux_land,taux_ssi,tauy_land,tauy_ssi,Error_code  &
    ! OUT fields
          , surf_ht_flux_land, zlcl_mixed                               &
          , theta_star_surf, qv_star_surf                               &
    ! OUT additional variables for JULES
          , tstar, ti_sicat, ext, snowmelt,tstar_land,tstar_sice, ei_surft &
          , ecan_surft, melt_surft, surf_htf_surft                      &
    ! OUT fields for coupling to the wave model
          , uwind_wav, vwind_wav                                        &
            )

    !-----------------------------------------------------------------------
    ! update main model prognostics
    ! Currently we're only modifying temperature.
    !-----------------------------------------------------------------------

    do k = 1, nlayers
      ! potential temperature increment on theta levels
      dtheta_bl(map_wth(1) + k) = t_latest(1,1,k)                       &
                                   / exner_theta_levels(1,1,k)          &
                                - theta_star(map_wth(1)+k)
      ! water vapour on theta levels
      m_v(map_wth(1) + k)  = q_latest(1,1,k)
      ! cloud liquid and ice water on theta levels
      m_cl(map_wth(1) + k) = qcl_latest(1,1,k)
      m_ci(map_wth(1) + k) = qcf_latest(1,1,k)
      ! diagnostic increments
      dt_bl(map_wth(1)+k)  = t_latest(1,1,k) - dt_bl(map_wth(1)+k)
      dmv_bl(map_wth(1)+k) = q_latest(1,1,k) - dmv_bl(map_wth(1)+k)
    end do
    ! copy down lowest level to surface as done in UM
    dtheta_bl(map_wth(1) + 0) = dtheta_bl(map_wth(1) + 1)
    m_v(map_wth(1) + 0)  = m_v(map_wth(1) + 1)
    m_cl(map_wth(1) + 0) = m_cl(map_wth(1) + 1)
    m_ci(map_wth(1) + 0) = m_ci(map_wth(1) + 1)
    dt_bl(map_wth(1)+0)  = dt_bl(map_wth(1)+1)
    dmv_bl(map_wth(1)+0) = dmv_bl(map_wth(1)+1)

    ! update blended Smagorinsky diffusion coefficients only if using Smagorinsky scheme
    if ( smagorinsky ) then
      do k = 1, nlayers
        visc_m_blend(map_wth(1) + k) = visc_m(1,1,k)
        visc_h_blend(map_wth(1) + k) = visc_h(1,1,k)
      end do
    endif

    ! update cloud fractions only if using cloud scheme and only on last
    ! dynamics iteration
    if ( cloud == cloud_um .and. &
         outer == outer_iterations ) then
      do k = 1, nlayers
        cf_bulk(map_wth(1) + k) = bulk_cloud_fraction(1,1,k)
        cf_ice(map_wth(1) + k)  = cloud_fraction_frozen(1,1,k)
        cf_liq(map_wth(1) + k)  = cloud_fraction_liquid(1,1,k)
        cf_area(map_wth(1) + k) = area_cloud_fraction(1,1,k)
      end do
      cf_ice(map_wth(1) + 0)  = cf_ice(map_wth(1) + 1)
      cf_liq(map_wth(1) + 0)  = cf_liq(map_wth(1) + 1)
      cf_bulk(map_wth(1) + 0) = cf_bulk(map_wth(1) + 1)
      cf_area(map_wth(1) + 0) = cf_area(map_wth(1) + 1)
    endif

    ! update BL prognostics
    ! There is no need to loop over nlayers_2d, since here
    ! nlayers_2d=1
    if (outer == outer_iterations) then

      i_tile = 0
      do i = first_land_tile, first_land_tile + n_land_tile - 1
        i_tile = i_tile + 1
        ! Land tile temperatures
        tile_temperature(map_tile(i)) = real(tstar_surft(1, i_tile), r_def)
         ! sensible heat flux
        tile_heat_flux(map_tile(i)) = real(ftl_surft(1, i_tile), r_def)
        ! moisture flux
        tile_moisture_flux(map_tile(i)) = real(fqw_surft(1, i_tile), r_def)
      end do

      ! Surface conductance (gs_gb)
      surface_conductance(map_2d(1)) = real(gs_gb(1), r_def)

      ! Sea temperature
      tile_temperature(map_tile(first_sea_tile)) = real(tstar_sea(1,1), r_def)
      ! heat flux - NB actually grid box mean but okay for aquaplanet
      tile_heat_flux(map_tile(first_sea_tile)) = real(ftl(1,1,1), r_def)
      ! moisture flux - NB actually grid box mean but okay for aquaplanet
      tile_moisture_flux(map_tile(first_sea_tile)) = real(fqw(1,1,1), r_def)

      i_sice = 0
      do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
        i_sice = i_sice + 1
        ! sea-ice temperature
        tile_temperature(map_tile(i)) = real(tstar_sice_sicat(1,1,i_sice), r_def)
        ! sea-ice heat flux
        tile_heat_flux(map_tile(i)) = real(ftl_ice(1,1,i_sice), r_def)
        ! sea-ice moisture flux
        tile_moisture_flux(map_tile(i)) = real(fqw_ice(1,1,i_sice), r_def)
      end do

      if (flandg(1, 1) > 0.0_r_um) then
        ! Gross primary productivity diagnostic
        gross_prim_prod(map_2d(1)) = real(sf_diag%gpp(1), r_def)
        ! Net primary productivity diagnostic
        net_prim_prod(map_2d(1)) = real(npp_gb(1), r_def)
      end if

      zh_2d(map_2d(1))     = zh(1,1)
      z0msea_2d(map_2d(1)) = z0msea(1,1)
      ntml_2d(map_2d(1))   = real(ntml(1,1))
      if (cumulus(1,1)) then
        cumulus_2d(map_2d(1)) = 1.0_r_def
      else
        cumulus_2d(map_2d(1)) = 0.0_r_def
      endif
      ! from convection
      lowest_cv_base(map_2d(1)) = lcbase(1,1)
    endif

    ! deallocate diagnostics deallocated in atmos_physics2
    CALL dealloc_bl_imp(bl_diag)
    CALL dealloc_sf_imp(sf_diag)
    CALL dealloc_sf_expl(sf_diag)

    ! set this back to 1 before exit
    land_field = 1

  end subroutine bl_code

end module bl_kernel_mod
