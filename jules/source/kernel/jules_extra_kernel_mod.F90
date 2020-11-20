!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to JULES surf_couple_extra.
!>
module jules_extra_kernel_mod

  use argument_mod,            only : arg_type,                               &
                                      GH_FIELD, GH_READ, GH_READWRITE, CELLS, &
                                      ANY_DISCONTINUOUS_SPACE_1,              &
                                      ANY_DISCONTINUOUS_SPACE_2,              &
                                      ANY_DISCONTINUOUS_SPACE_3,              &
                                      ANY_DISCONTINUOUS_SPACE_4,              &
                                      ANY_DISCONTINUOUS_SPACE_5
  use constants_mod,           only : i_def, i_um, r_def, r_um
  use kernel_mod,              only : kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: jules_extra_kernel_type
    private
    type(arg_type) :: meta_args(50) = (/                             &
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! ls_rain
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! conv_rain
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! ls_snow
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! conv_snow
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! tile_fraction
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! leaf_area_index
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! canopy_height
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! snow_unload_rate
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_moist_wilt
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_moist_crit
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_moist_sat
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_cond_sat
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_thermal_cap
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_thermal_cond
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_suction_sat
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! clapp_horn_b
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_carbon_content
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! mean_topog_index
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! a_sat_frac
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! c_sat_frac
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! a_wet_frac
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! c_wet_frac
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! tile_temperature
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! net_prim_prod
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! snow_sublimation
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! surf_heat_flux
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! canopy_evap
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_4), & ! water_extraction
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! thermal_cond_wet_soil
        arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! soil_respiration
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4), & ! soil_temperature
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4), & ! soil_moisture
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4), & ! unfrozen_soil_moisture
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4), & ! frozen_soil_moisture
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! canopy_water
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! tile_snow_mass
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! tile_snow_rgrain
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! n_snow_layers
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! snow_depth
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! snow_under_canopy
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! snowpack_density
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_5), & ! snow_layer_thickness
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_5), & ! snow_layer_ice_mass
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_5), & ! snow_layer_liq_mass
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_5), & ! snow_layer_temp
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_5), & ! snow_layer_rgrain
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! total_snowmelt
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! soil_sat_frac
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! water_table
        arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1)  & ! wetness_under_soil
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass ::jules_extra_code
  end type

  public jules_extra_code

contains

  !> @brief Interface to JULES surf_couple_extra
  !> @details JULES surf_couple_extra calculates the surface and soil fluxes
  !> and stores of Water (via the hydrol{ogy}, snow and river_control
  !> modules) and Carbon (via triffid and inferno/fire modules). Note: only the
  !> hydrol and snow components are implemented so far.

  !> @param[in]     nlayers                Number of layers
  !> @param[in]     ls_rain                Large-scale rainfall rate (kg m-2 s-1)
  !> @param[in]     conv_rain              Convective rainfall rate (kg m-2 s-1)
  !> @param[in]     ls_snow                Large-scale snowfall rate (kg m-2 s-1)
  !> @param[in]     conv_snow              Convective snowfall rate (kg m-2 s-1)
  !> @param[in]     tile_fraction          Surface tile fractions
  !> @param[in]     leaf_area_index        Leaf Area Index
  !> @param[in]     canopy_height          Canopy height (m)
  !> @param[in]     snow_unload_rate       Unloading of snow from PFTs by wind
  !> @param[in]     soil_moist_wilt        Volumetric soil moist at wilting pt
  !> @param[in]     soil_moist_crit        Volumetric soil moist at critical pt
  !> @param[in]     soil_moist_sat         Volumetric soil moist at saturation
  !> @param[in]     soil_cond_sat          Saturated soil thermal conductivity (kg m-2 s-1)
  !> @param[in]     soil_thermal_cap       Soil thermal capacity (J m-3 K-1)
  !> @param[in]     soil_thermal_cond      Soil thermal conductivity (W m-1 K-1)
  !> @param[in]     soil_suction_sat       Saturated soil water suction (m)
  !> @param[in]     clapp_horn_b           Clapp and Hornberger b coefficient
  !> @param[in]     soil_carbon_content    Soil carbon content (kg m-2)
  !> @param[in]     mean_topog_index       Mean topographic index
  !> @param[in]     a_sat_frac             a gridbox saturated fraction
  !> @param[in]     c_sat_frac             c gridbox saturated fraction
  !> @param[in]     a_wet_frac             a gridbox wet fraction
  !> @param[in]     c_wet_frac             c gridbox wet fraction
  !> @param[in]     tile_temperature       Surface tile temperatures (K)
  !> @param[in]     net_prim_prod          Net Primary Productivity (kg m-2 s-1)
  !> @param[in]     snow_sublimation       Sublimation of snow (kg m-2 s-1)
  !> @param[in]     surf_heat_flux         Surface heat flux (W m-2)
  !> @param[in]     canopy_evap            Canopy evaporation from land tiles (kg m-2 s-1)
  !> @param[in]     water_extraction       Extraction of water from each soil layer (kg m-2 s-1)
  !> @param[in]     thermal_cond_wet_soil  Thermal conductivity of soil (W m-1 K-1)
  !> @param[in]     soil_respiration       Soil respiration (kg m-2 s-1)
  !> @param[in,out] soil_temperature       Soil temperature (K)
  !> @param[in,out] soil_moisture          Soil moisture content (kg m-2)
  !> @param[in,out] unfrozen_soil_moisture Unfrozen soil moisture proportion
  !> @param[in,out] frozen_soil_moisture   Frozen soil moisture proportion
  !> @param[in,out] canopy_water           Canopy water on each tile (kg m-2)
  !> @param[in,out] tile_snow_mass         Snow mass on tiles (kg m-2)
  !> @param[in,out] tile_snow_rgrain       Snow grain radius on tiles (microns)
  !> @param[in,out] n_snow_layers          Number of snow layers on tiles
  !> @param[in,out] snow_depth             Snow depth on tiles (m)
  !> @param[in,out] snow_under_canopy      Amount of snow under canopy (kg m-2)
  !> @param[in,out] snowpack_density       Density of snow on ground (kg m-3)
  !> @param[in,out] snow_layer_thickness   Thickness of snow layers (m)
  !> @param[in,out] snow_layer_ice_mass    Mass of ice in snow layers (kg m-2)
  !> @param[in,out] snow_layer_liq_mass    Mass of liquid in snow layers (kg m-2)
  !> @param[in,out] snow_layer_temp        Temperature of snow layer (K)
  !> @param[in,out] snow_layer_rgrain      Grain radius of snow layer (microns)
  !> @param[in,out] total_snowmelt         Surface plus canopy snowmelt rate (kg m-2 s-1)
  !> @param[in,out] soil_sat_frac          Soil saturated fraction
  !> @param[in,out] water_table            Water table depth (m)
  !> @param[in,out] wetness_under_soil     Soil wetness below soil column
  !> @param[in]     ndf_2d                 Total DOFs per cell for 2D fields
  !> @param[in]     undf_2d                Unique DOFs per cell for 2D fields
  !> @param[in]     map_2d                 DOFmap for cells for 2D fields
  !> @param[in]     ndf_tile               Total DOFs per cell for surface tiles
  !> @param[in]     undf_tile              Unique DOFs per cell for surface tiles
  !> @param[in]     map_tile               DOFmap for cells for surface tiles
  !> @param[in]     ndf_pft                Total DOFs per cell for plant types
  !> @param[in]     undf_pft               Unique DOFs per cell for plant types
  !> @param[in]     map_pft                DOFmap for cells for plant types
  !> @param[in]     ndf_soil               Total DOFs per cell for soil levels
  !> @param[in]     undf_soil              Unique DOFs per cell for soil levels
  !> @param[in]     map_soil               DOFmap for cells for soil levels
  !> @param[in]     ndf_snow               Total DOFs per cell for snow layers
  !> @param[in]     undf_snow              Unique DOFs per cell for snow layers
  !> @param[in]     map_snow               DOFmap for cells for snow layers
  subroutine jules_extra_code(             &
               nlayers,                    &
               ls_rain,                    &
               conv_rain,                  &
               ls_snow,                    &
               conv_snow,                  &
               tile_fraction,              &
               leaf_area_index,            &
               canopy_height,              &
               snow_unload_rate,           &
               soil_moist_wilt,            &
               soil_moist_crit,            &
               soil_moist_sat,             &
               soil_cond_sat,              &
               soil_thermal_cap,           &
               soil_thermal_cond,          &
               soil_suction_sat,           &
               clapp_horn_b,               &
               soil_carbon_content,        &
               mean_topog_index,           &
               a_sat_frac,                 &
               c_sat_frac,                 &
               a_wet_frac,                 &
               c_wet_frac,                 &
               tile_temperature,           &
               net_prim_prod,              &
               snow_sublimation,           &
               surf_heat_flux,             &
               canopy_evap,                &
               water_extraction,           &
               thermal_cond_wet_soil,      &
               soil_respiration,           &
               soil_temperature,           &
               soil_moisture,              &
               unfrozen_soil_moisture,     &
               frozen_soil_moisture,       &
               canopy_water,               &
               tile_snow_mass,             &
               tile_snow_rgrain,           &
               n_snow_layers,              &
               snow_depth,                 &
               snow_under_canopy,          &
               snowpack_density,           &
               snow_layer_thickness,       &
               snow_layer_ice_mass,        &
               snow_layer_liq_mass,        &
               snow_layer_temp,            &
               snow_layer_rgrain,          &
               total_snowmelt,             &
               soil_sat_frac,              &
               water_table,                &
               wetness_under_soil,         &
               ndf_2d,                     &
               undf_2d,                    &
               map_2d,                     &
               ndf_tile,                   &
               undf_tile,                  &
               map_tile,                   &
               ndf_pft,                    &
               undf_pft,                   &
               map_pft,                    &
               ndf_soil,                   &
               undf_soil,                  &
               map_soil,                   &
               ndf_snow,                   &
               undf_snow,                  &
               map_snow                    &
                            )

    !---------------------------------------
    ! LFRic modules
    !---------------------------------------
    use jules_control_init_mod, only: &
         nsurft => n_land_tile, n_land_tile
    use jules_physics_init_mod, only: decrease_sath_cond

    ! Module imports for surf_couple_extra JULESvn5.4
    use ancil_info, only: nsoilt, l_soil_point, soil_pts, lice_pts, l_lice_point
    use atm_step_local, only: dim_cs1, land_pts_trif, npft_trif
    use cderived_mod, only: delta_lambda, delta_phi
    use jules_surface_types_mod, only: npft, ntype
    use jules_snow_mod, only: nsmax
    use nlsizes_namelist_mod, only: row_length, rows, land_pts => land_field, &
                                    sm_levels, ntiles
    use UM_ParCore, only: nproc

    ! Spatially varying fields used from modules
    use jules_internal, only: unload_backgrnd_pft
    use p_s_parms, only: bexp_soilt, sathh_soilt, hcap_soilt, hcon_soilt,     &
                         satcon_soilt, smvccl_soilt, smvcwt_soilt,            &
                         smvcst_soilt
    use prognostics, only: snowdepth_surft, nsnow_surft, rgrainl_surft,       &
                           tsnow_surft, sliq_surft, sice_surft, ds_surft,     &
                           rho_snow_grnd_surft
    use trignometric_mod, only: cos_theta_latitude

    ! Jules related subroutines
    use sparm_mod, only             : sparm
    use tilepts_mod, only           : tilepts
    use surf_couple_extra_mod, only : surf_couple_extra
    use infiltration_rate_mod, only : infiltration_rate

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d, ndf_tile, undf_tile,  &
                                       ndf_pft, undf_pft, ndf_soil,           &
                                       undf_soil, ndf_snow, undf_snow

    integer(kind=i_def), dimension(ndf_2d),    intent(in) :: map_2d
    integer(kind=i_def), dimension(ndf_tile),  intent(in) :: map_tile
    integer(kind=i_def), dimension(ndf_pft),   intent(in) :: map_pft
    integer(kind=i_def), dimension(ndf_soil),  intent(in) :: map_soil
    integer(kind=i_def), dimension(ndf_snow),  intent(in) :: map_snow

    real(kind=r_def), dimension(undf_2d), intent(in) :: ls_rain, conv_rain,   &
                                                        ls_snow, conv_snow

    real(kind=r_def), intent(in)    :: tile_fraction(undf_tile)
    real(kind=r_def), intent(in)    :: snow_sublimation(undf_tile)
    real(kind=r_def), intent(in)    :: surf_heat_flux(undf_tile)
    real(kind=r_def), intent(in)    :: canopy_evap(undf_tile)
    real(kind=r_def), intent(in)    :: tile_temperature(undf_tile)

    real(kind=r_def), intent(in)    :: leaf_area_index(undf_pft)
    real(kind=r_def), intent(in)    :: canopy_height(undf_pft)
    real(kind=r_def), intent(in)    :: snow_unload_rate(undf_pft)

    real(kind=r_def), intent(in)    :: soil_moist_wilt(undf_2d)
    real(kind=r_def), intent(in)    :: soil_moist_crit(undf_2d)
    real(kind=r_def), intent(in)    :: soil_moist_sat(undf_2d)
    real(kind=r_def), intent(in)    :: soil_cond_sat(undf_2d)
    real(kind=r_def), intent(in)    :: soil_thermal_cap(undf_2d)
    real(kind=r_def), intent(in)    :: soil_suction_sat(undf_2d)
    real(kind=r_def), intent(in)    :: clapp_horn_b(undf_2d)
    real(kind=r_def), intent(in)    :: water_extraction(undf_soil)

    real(kind=r_def), intent(in)    :: soil_thermal_cond(undf_2d)
    real(kind=r_def), intent(in)    :: soil_carbon_content(undf_2d)
    real(kind=r_def), intent(in)    :: mean_topog_index(undf_2d)
    real(kind=r_def), intent(in)    :: a_sat_frac(undf_2d)
    real(kind=r_def), intent(in)    :: c_sat_frac(undf_2d)
    real(kind=r_def), intent(in)    :: a_wet_frac(undf_2d)
    real(kind=r_def), intent(in)    :: c_wet_frac(undf_2d)
    real(kind=r_def), intent(in)    :: net_prim_prod(undf_2d)
    real(kind=r_def), intent(in)    :: thermal_cond_wet_soil(undf_2d)
    real(kind=r_def), intent(in)    :: soil_respiration(undf_2d)

    real(kind=r_def), intent(inout) :: canopy_water(undf_tile)
    real(kind=r_def), intent(inout) :: tile_snow_mass(undf_tile)
    real(kind=r_def), intent(inout) :: tile_snow_rgrain(undf_tile)
    real(kind=r_def), intent(inout) :: n_snow_layers(undf_tile)
    real(kind=r_def), intent(inout) :: snow_depth(undf_tile)
    real(kind=r_def), intent(inout) :: snow_under_canopy(undf_tile)
    real(kind=r_def), intent(inout) :: snowpack_density(undf_tile)
    real(kind=r_def), intent(inout) :: total_snowmelt(undf_tile)

    real(kind=r_def), intent(inout) :: snow_layer_thickness(undf_snow)
    real(kind=r_def), intent(inout) :: snow_layer_ice_mass(undf_snow)
    real(kind=r_def), intent(inout) :: snow_layer_liq_mass(undf_snow)
    real(kind=r_def), intent(inout) :: snow_layer_temp(undf_snow)
    real(kind=r_def), intent(inout) :: snow_layer_rgrain(undf_snow)

    real(kind=r_def), intent(inout) :: soil_temperature(undf_soil)
    real(kind=r_def), intent(inout) :: soil_moisture(undf_soil)
    real(kind=r_def), intent(inout) :: unfrozen_soil_moisture(undf_soil)
    real(kind=r_def), intent(inout) :: frozen_soil_moisture(undf_soil)

    real(kind=r_def), intent(inout) :: soil_sat_frac(undf_2d)
    real(kind=r_def), intent(inout) :: water_table(undf_2d)
    real(kind=r_def), intent(inout) :: wetness_under_soil(undf_2d)

    ! Local variables for the kernel
    integer :: i,j, i_pft, i_snow

!------------------------------------------------------------------------------
    ! JULES surf_couple_extra subroutine arguments declared using JULESvn5.4
    !       order and naming

    ! Integer parameters
    integer(i_um), parameter :: river_row_length = 1, river_rows = 1,         &
                                aocpl_row_length = 1, aocpl_p_rows = 1

    ! Integer indices (module intent=in)
    integer(i_um) :: a_step, g_p_field, g_r_field, global_row_length,         &
         global_rows, global_river_row_length, global_river_rows

    integer(i_um), dimension(land_pts) :: land_index, lice_index, soil_index

    integer(i_um), dimension(nsurft) :: surft_pts

    integer(i_um), dimension(land_pts,nsurft) :: surft_index

    ! Logical (module intent=in)
    logical :: smlt, stf_sub_surf_roff

    logical, dimension(row_length,rows) :: land_sea_mask

    ! Real variables (module intent=in)
    ! Driving data
    real(r_um), dimension(row_length, rows) :: ls_rain_ij, con_rain_ij,       &
         ls_snow_ij, con_snow_ij, pstar_ij, ls_graup_ij, tl_1_ij,             &
         lw_down_ij, qw_1_ij, u_1_ij, v_1_ij, cca_2d, soil_clay_ij,           &
         flash_rate_ancil, pop_den_ancil, flandg

    ! State
    real(r_um), dimension(land_pts, ntiles) :: tstar_surft,                   &
    ! Fluxes
         ei_surft, surf_htf_surft, ecan_surft, tile_frac, sw_surft,           &
         fqw_surft

    real(r_um), dimension(land_pts, nsoilt, sm_levels) :: ext_soilt

    real(r_um), dimension(land_pts, nsoilt) :: fexp_soilt, gamtot_soilt,      &
         ti_mean_soilt, ti_sig_soilt, a_fsat_soilt, c_fsat_soilt,             &
         a_fwet_soilt, c_fwet_soilt

    real(r_um), dimension(land_pts) :: npp_gb, frac_agr_gb, z0m_soil_gb

    ! River routing
    real(r_um), dimension(aocpl_row_length) :: xpa
    real(r_um), dimension(0:aocpl_row_length) :: xua
    real(r_um), dimension(aocpl_row_length+1) :: xva
    real(r_um), dimension(aocpl_p_rows) :: ypa, yua
    real(r_um), dimension(0:aocpl_p_rows) :: yva
    real(r_um), dimension(river_row_length, river_rows) :: trivdir, trivseq
    real(r_um), dimension(row_length, rows) :: r_area, slope, flowobs1,       &
         r_inext, r_jnext, r_land

    ! Integers (module intent = in out)
    integer(i_um) :: a_steps_since_riv, asteps_since_triffid

    ! Real variables (module intent = in out)
    real(r_um), dimension(land_pts, nsurft) :: melt_surft, catch_snow_surft,  &
         tsurf_elev_surft, rgrain_surft, snow_grnd_surft, snow_surft,         &
         infil_surft, catch_surft, canopy_surft, frac_surft

    real(r_um), dimension(row_length, rows) :: substore, surfstore, flowin,   &
         bflowin, acc_lake_evap

    real(r_um), dimension(river_row_length, river_rows) :: twatstor

    real(r_um), dimension(land_pts, npft) :: g_leaf_acc_pft,                  &
         g_leaf_phen_acc_pft, lai_pft, canht_pft

    real(r_um), dimension(land_pts_trif, npft_trif) :: npp_acc_pft,           &
         resp_w_acc_pft

    real(r_um), dimension(land_pts_trif, dim_cs1) :: resp_s_acc_gb_um,        &
         cs_pool_gb_um

    real(r_um), dimension(land_pts) :: snomlt_sub_htf, dhf_surf_minus_soil,   &
         ls_rainfrac_gb, tot_surf_runoff, tot_sub_runoff, inlandout_atm_gb

    real(r_um), dimension(land_pts, nsoilt) :: hcons_soilt, fsat_soilt,       &
         fwetl_soilt, zw_soilt, sthzw_soilt, cs_ch4_soilt

    real(r_um), dimension(land_pts, nsoilt, 1, dim_cs1) :: resp_s_soilt

    real(r_um), dimension(land_pts, nsoilt, sm_levels) :: t_soil_soilt,       &
         smcl_soilt, sthf_soilt, sthu_soilt

    ! Real variables (module intent = out)
    real(r_um), dimension(land_pts) :: sub_surf_roff, surf_roff, tot_tfall,   &
         snowmelt_gb, rrun, rflow, canopy_gb

    real(r_um), dimension(row_length, rows) :: snomlt_surf_htf, snowmelt_ij,  &
         snow_mass_ij

    real(r_um), dimension(land_pts, nsurft) :: snow_soil_htf, smc_soilt,      &
         z0_surft, z0h_bare_surft

    !-------------------------------------------------------------------

    ! Data from 2D fields
    ls_rain_ij(1,1)  = ls_rain(map_2d(1))   ! Large-scale rainfall rate
    con_rain_ij(1,1) = conv_rain(map_2d(1)) ! Convective rainfall rate
    ls_snow_ij(1,1)  = ls_snow(map_2d(1))   ! Large-scale snowfall rate
    con_snow_ij(1,1) = conv_snow(map_2d(1)) ! Convective snowfallfall rate

    !----------------------------------------------------------------------
    ! Surface fields as needed by Jules
    !----------------------------------------------------------------------

    ! Ancillaries:
    ! Land tile fractions (frac_surft)
    flandg = 0.0_r_um
    do i = 1, n_land_tile
      flandg = flandg + real(tile_fraction(map_tile(i)), r_um)
      frac_surft(1, i)     = real(tile_fraction(map_tile(i)), r_um)
      ei_surft(1, i)       = real(snow_sublimation(map_tile(i)), r_um)
      surf_htf_surft(1, i) = real(surf_heat_flux(map_tile(i)), r_um)
      ecan_surft(1, i)     = real(canopy_evap(map_tile(i)), r_um)
    end do

    ! Jules requires fractions with respect to the land area
    if (flandg(1, 1) > 0.0_r_um) then
      land_pts = 1
      land_index = 1
      frac_surft(1, 1:n_land_tile) = frac_surft(1, 1:n_land_tile) / flandg(1, 1)
      land_sea_mask = .true.
    else
      land_pts = 0
      land_index = 0
      land_sea_mask = .false.
    end if

    ! Logical flags controlling diagnostic calculations
    smlt = .false.
    stf_sub_surf_roff = .false.

    ! GPW fudge: set tile_frac = frac_surft
    tile_frac = frac_surft

    ! Set type_pts and type_index
    call tilepts(land_pts, frac_surft, surft_pts, surft_index)

    ! Vegetation prognostics (for TRIFFID)
    do i_pft = 1, npft
      ! Leaf area index
      lai_pft(1, i_pft) = real(leaf_area_index(map_pft(i_pft)), r_um)
      ! Canopy height
      canht_pft(1, i_pft) = real(canopy_height(map_pft(i_pft)), r_um)
      ! Unloading rate of snow from plant functional types
      unload_backgrnd_pft(1, i_pft) = real(snow_unload_rate(map_pft(i_pft)), r_um)
    end do

    ! Get catch_snow_surft and catch_surft from call to sparm
    call sparm(land_pts, n_land_tile, surft_pts, surft_index,                 &
               frac_surft, canht_pft, lai_pft, z0m_soil_gb,                   &
               catch_snow_surft, catch_surft, z0_surft, z0h_bare_surft)

    ! Soil carbon content (cs_pool_gb_um, )
    cs_pool_gb_um = real(soil_carbon_content(map_2d(1)), r_um)

    ! GPW fudge while TRIFFID not implemented: set cs_ch4_soilt = cs_pool_gb_um
    cs_ch4_soilt = cs_pool_gb_um

    ! Decrease in saturated conductivity of soil with depth
    fexp_soilt = decrease_sath_cond

    ! Mean grid box topographic index
    ti_mean_soilt = real(mean_topog_index(map_2d(1)), r_um)

    ! a saturated soil fraction
    a_fsat_soilt = real(a_sat_frac(map_2d(1)), r_um)

    ! c saturated soil fraction
    c_fsat_soilt = real(c_sat_frac(map_2d(1)), r_um)

    ! a soil wetness fraction
    a_fwet_soilt = real(a_wet_frac(map_2d(1)), r_um)

    ! c soil wetness fraction
    c_fwet_soilt = real(c_wet_frac(map_2d(1)), r_um)

    ! Net primary productivity diagnostic
    npp_gb = real(net_prim_prod(map_2d(1)), r_um)

    ! Prognostics:
    ! Land tile temperatures
    do i = 1, n_land_tile
      tstar_surft(1, i) = real(tile_temperature(map_tile(i)), r_um)
    end do

    ! Soil ancillaries and prognostics
    satcon_soilt(1, 1, 0) = real(soil_cond_sat(map_2d(1)), r_um)
    do i = 1, sm_levels
    ! Volumetric soil moisture at wilting point (smvcwt_soilt)
      smvcwt_soilt(1, 1, i) = real(soil_moist_wilt(map_2d(1)), r_um)
    ! Volumetric soil moisture at critical point (smvccl_soilt)
      smvccl_soilt(1, 1, i) = real(soil_moist_crit(map_2d(1)), r_um)
    ! Volumetric soil moisture at saturation (smvcst_soilt)
      smvcst_soilt(1, 1, i) = real(soil_moist_sat(map_2d(1)), r_um)
    ! Saturated soil conductivity (satcon_soilt)
      satcon_soilt(1, 1, i) = real(soil_cond_sat(map_2d(1)), r_um)
    ! Soil thermal capacity (hcap_soilt)
      hcap_soilt(1, 1, i) = real(soil_thermal_cap(map_2d(1)), r_um)
    ! Saturated soil water suction (sathh_soilt)
      sathh_soilt(1, 1, i) = real(soil_suction_sat(map_2d(1)), r_um)
    ! Clapp and Hornberger b coefficient (bexp_soilt)
      bexp_soilt(1, 1, i) = real(clapp_horn_b(map_2d(1)), r_um)
    ! Soil temperature (t_soil_soilt)
      t_soil_soilt(1, 1, i) = real(soil_temperature(map_soil(i)), r_um)
    ! Soil moisture content (kg m-2, soil_layer_moisture/smcl_soilt)
      smcl_soilt(1, 1, i) = real(soil_moisture(map_soil(i)), r_um)
    ! Unfrozen soil moisture proportion (sthu_soilt)
      sthu_soilt(1, 1, i) = real(unfrozen_soil_moisture(map_soil(i)), r_um)
    ! Frozen soil moisture proportion (sthf_soilt)
      sthf_soilt(1, 1, i) = real(frozen_soil_moisture(map_soil(i)), r_um)
    ! Water extraction (ext_soilt [=ext in bl_kernel))
      ext_soilt(1, 1, i) = real(water_extraction(map_soil(i)), r_um)
    end do

    ! Soil thermal conductivity (hcon_soilt, hcons_soilt)
    hcon_soilt = real(soil_thermal_cond(map_2d(1)), r_um)
    hcons_soilt = real(thermal_cond_wet_soil(map_2d(1)), r_um)

    ! Soil and land ice ancils dependant on smvcst_soilt
    ! (soil moisture saturation limit)
    if ( smvcst_soilt(1, 1, 1) > 0.0_r_um ) then
      soil_pts = 1
      soil_index = 1
      l_soil_point = .true.
      lice_pts = 0
      lice_index = 0
      l_lice_point = .false.
    else
      soil_pts = 0
      soil_index = 0
      l_soil_point = .false.
      lice_pts = 1
      lice_index = 1
      l_lice_point = .true.
    end if

    ! Calculate the infiltration rate
    call infiltration_rate(land_pts, ntiles, surft_pts, surft_index,          &
                           satcon_soilt, frac_surft, infil_surft)

    ! Canopy water on each tile (canopy_surft)
    do i = 1, n_land_tile
      canopy_surft(1, i) = real(canopy_water(map_tile(i)), r_um)
    end do

    ! Snow prognostics
    i_snow = 0
    do i = 1, n_land_tile
      ! Lying snow mass on land tiles
      snow_surft(1, i) = real(tile_snow_mass(map_tile(i)), r_um)
      ! Snow grain size on tiles (microns)
      rgrain_surft(1, i) = real(tile_snow_rgrain(map_tile(i)), r_um)
      ! Number of snow layers on tiles (nsnow_surft)
      nsnow_surft(1, i) = int(n_snow_layers(map_tile(i)), i_um)
      ! Snow depth on tiles (snowdepth_surft)
      snowdepth_surft(1, i) = real(snow_depth(map_tile(i)), r_um)
      ! Snow mass under canopy
      snow_grnd_surft(1, i) = real(snow_under_canopy(map_tile(i)), r_um)
      ! Snowpack density (rho_snow_grnd_surft)
      rho_snow_grnd_surft(1, i) = real(snowpack_density(map_tile(i)), r_um)
      do j = 1, nsmax
        i_snow = i_snow + 1
        ! Thickness of snow layers
        ds_surft(1, i, j) = real(snow_layer_thickness(map_snow(i_snow)), r_um)
        ! Mass of ice in snow layers
        sice_surft(1, i, j) = real(snow_layer_ice_mass(map_snow(i_snow)), r_um)
        ! Mass of liquid in snow layers
        sliq_surft(1, i, j) = real(snow_layer_liq_mass(map_snow(i_snow)), r_um)
        ! Temperature of snow layers
        tsnow_surft(1, i, j) = real(snow_layer_temp(map_snow(i_snow)), r_um)
        ! Grain size of snow layers
        rgrainl_surft(1, i, j) = real(snow_layer_rgrain(map_snow(i_snow)), r_um)
      end do
    end do

    ! Snow melt
    do i = 1, n_land_tile
      melt_surft(1, i) = real(total_snowmelt(map_tile(i)), r_um)
    end do

    ! Soil saturated fraction
    fsat_soilt = real(soil_sat_frac(map_2d(1)), r_um)

    ! Water table depth
    zw_soilt = real(water_table(map_2d(1)), r_um)

    ! Soil wetness below soil column
    sthzw_soilt = real(wetness_under_soil(map_2d(1)), r_um)

    ! Soil respiration (resp_s_soilt [=resp_s_soilt in bl_kernel])
    resp_s_soilt = real(soil_respiration(map_2d(1)), r_um)

  !----------------------------------------------------------------------------
  ! Call to surf_couple_extra using JULESvn5.4 standalone variable names

    call surf_couple_extra(                                                   &
    !Driving data and associated INTENT(IN)
    ls_rain_ij, con_rain_ij, ls_snow_ij, con_snow_ij, tl_1_ij, lw_down_ij,    &
    qw_1_ij, u_1_ij, v_1_ij, pstar_ij,                                        &

    !Fluxes INTENT(IN)
    ei_surft, surf_htf_surft, ecan_surft, ext_soilt, sw_surft,                &

    !Misc INTENT(IN)
    a_step, smlt, tile_frac, hcons_soilt,                                     &

    !Fluxes INTENT(IN OUT)
    melt_surft,                                                               &

    !Fluxes INTENT(OUT)
    snomlt_surf_htf, snowmelt_ij, snomlt_sub_htf, sub_surf_roff, surf_roff,   &
    tot_tfall, snowmelt_gb, rrun, rflow, snow_soil_htf,                       &

    !IN
    land_pts, row_length, rows, river_row_length, river_rows, land_index,     &
    ls_graup_ij, cca_2d, nsurft, surft_pts, surft_index, tstar_surft,         &
    lice_pts, lice_index, soil_pts, soil_index, stf_sub_surf_roff,fexp_soilt, &
    gamtot_soilt, ti_mean_soilt, ti_sig_soilt, cs_ch4_soilt, flash_rate_ancil,&
    pop_den_ancil, a_fsat_soilt, c_fsat_soilt, a_fwet_soilt, c_fwet_soilt,    &
    ntype, fqw_surft, delta_lambda, delta_phi,                                &
    cos_theta_latitude, aocpl_row_length, aocpl_p_rows, xpa, xua, xva, ypa,   &
    yua, yva, g_p_field, g_r_field, nproc, global_row_length, global_rows,    &
    global_river_row_length, global_river_rows, flandg, trivdir, trivseq,     &
    r_area, slope, flowobs1, r_inext, r_jnext, r_land, frac_agr_gb,           &
    soil_clay_ij, resp_s_soilt, npp_gb, z0m_soil_gb,                          &

    !IN OUT
    a_steps_since_riv, t_soil_soilt, tsurf_elev_surft, rgrain_surft,          &
    snow_grnd_surft, snow_surft, smcl_soilt, sthf_soilt, sthu_soilt,          &
    canopy_surft, fsat_soilt, fwetl_soilt, zw_soilt, sthzw_soilt,             &
    snow_mass_ij, ls_rainfrac_gb, substore, surfstore, flowin, bflowin,       &
    tot_surf_runoff, tot_sub_runoff, acc_lake_evap, twatstor,                 &
    asteps_since_triffid, g_leaf_acc_pft, g_leaf_phen_acc_pft,                &
    npp_acc_pft, resp_s_acc_gb_um, resp_w_acc_pft, cs_pool_gb_um, frac_surft, &
    lai_pft, canht_pft, catch_snow_surft, catch_surft, infil_surft,           &
    inlandout_atm_gb,                                                         &

    !OUT- mostly for SCM diagnostics output below
    dhf_surf_minus_soil, canopy_gb, smc_soilt, z0_surft, z0h_bare_surft,      &

    ! IN
    land_sea_mask                                                             &
    )

  !---------------------------------------------------------------------------
  ! Return the updated prognostic values to jules_prognostics.

    i_snow = 0
    do i = 1, n_land_tile
      ! Canopy water on each tile (canopy_surft)
      canopy_water(map_tile(i)) = real(canopy_surft(1, i), r_def)
      ! Lying snow mass on land tiles
      tile_snow_mass(map_tile(i)) = real(snow_surft(1, i), r_def)
      ! Number of snow layers on tiles (nsnow_surft)
      n_snow_layers(map_tile(i)) = real(nsnow_surft(1, i), r_def)
      ! Snow depth on tiles (snowdepth_surft)
      snow_depth(map_tile(i)) = real(snowdepth_surft(1, i), r_def)
      ! Snow grain size on tiles (microns)
      tile_snow_rgrain(map_tile(i)) = real(rgrain_surft(1, i), r_def)
      ! Snow mass under canopy
      snow_under_canopy(map_tile(i)) = real(snow_grnd_surft(1, i), r_def)
      ! Snowpack density (rho_snow_grnd_surft)
      snowpack_density(map_tile(i)) = real(rho_snow_grnd_surft(1, i), r_def)
      ! Total snowmelt
      total_snowmelt(map_tile(i)) = real(melt_surft(1, i), r_def)
      do j = 1, nsmax
        i_snow = i_snow + 1
        ! Thickness of snow layers
        snow_layer_thickness(map_snow(i_snow)) = real(ds_surft(1, i, j), r_def)
        ! Mass of ice in snow layers
        snow_layer_ice_mass(map_snow(i_snow)) = real(sice_surft(1, i, j), r_def)
        ! Mass of liquid in snow layers
        snow_layer_liq_mass(map_snow(i_snow)) = real(sliq_surft(1, i, j), r_def)
        ! Temperature of snow layers
        snow_layer_temp(map_snow(i_snow)) = real(tsnow_surft(1, i, j), r_def)
        ! Grain size of snow layers
        snow_layer_rgrain(map_snow(i_snow)) = real(rgrainl_surft(1, i, j), r_def)
      end do
    end do

    do i = 1, sm_levels
      ! Soil temperature (t_soil_soilt)
      soil_temperature(map_soil(i)) = real(t_soil_soilt(1, 1, i), r_def)
      ! Soil moisture content (kg m-2, soil_layer_moisture)
      soil_moisture(map_soil(i)) = real(smcl_soilt(1, 1, i), r_def)
      ! Unfrozen soil moisture proportion (sthu_soilt)
      unfrozen_soil_moisture(map_soil(i)) = real(sthu_soilt(1, 1, i), r_def)
      ! Frozen soil moisture proportion (sthf_soilt)
      frozen_soil_moisture(map_soil(i)) = real(sthf_soilt(1, 1, i), r_def)
    end do

    ! Soil saturated fraction
    soil_sat_frac(map_2d(1)) = real(fsat_soilt(1,1), r_def)

    ! Water table depth
    water_table(map_2d(1)) = real(zw_soilt(1,1), r_def)

    ! Wetness below soil column
    wetness_under_soil(map_2d(1)) = real(sthzw_soilt(1,1), r_def)

    ! Reset land_pts to 1 before exit
    ! This is required because prior to calling the kernel, it is unknown
    ! whether this is a land point or not. Local variables are dimensioned
    ! using this as their size, hence it needs to be at least 1 to allow
    ! for the possibility this is a land point. However it needs to be
    ! set correctly in the kernel for use in Jules. Therefore it needs setting
    ! back to 1 so that variables are dimensioned correctly on the next
    ! kernel call.
    land_pts = 1

  end subroutine jules_extra_code

end module jules_extra_kernel_mod
