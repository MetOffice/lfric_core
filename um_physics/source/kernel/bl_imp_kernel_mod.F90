!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the implicit UM boundary layer scheme.
!>
module bl_imp_kernel_mod

  use argument_mod,           only : arg_type,                  &
                                     GH_FIELD, GH_INTEGER,      &
                                     GH_READ, GH_WRITE, GH_INC, &
                                     GH_READWRITE, CELLS,       &
                                     ANY_DISCONTINUOUS_SPACE_1, &
                                     ANY_DISCONTINUOUS_SPACE_2, &
                                     ANY_DISCONTINUOUS_SPACE_3, &
                                     ANY_DISCONTINUOUS_SPACE_4, &
                                     ANY_DISCONTINUOUS_SPACE_5, &
                                     ANY_DISCONTINUOUS_SPACE_6, &
                                     ANY_DISCONTINUOUS_SPACE_7
  use section_choice_config_mod, only : cloud, cloud_um
  use cloud_config_mod,       only : scheme, scheme_smith, scheme_pc2
  use constants_mod,          only : i_def, i_um, r_def, r_um
  use fs_continuity_mod,      only : W3, Wtheta, W2
  use kernel_mod,             only : kernel_type
  use timestepping_config_mod, only: outer_iterations
  use physics_config_mod,     only : lowest_level,          &
                                     lowest_level_constant, &
                                     lowest_level_gradient, &
                                     lowest_level_flux
  use planet_config_mod,      only : cp

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: bl_imp_kernel_type
    private
    type(arg_type) :: meta_args(84) = (/                              &
        arg_type(GH_INTEGER, GH_READ),                                &! outer
        arg_type(GH_FIELD,   GH_READ,      WTHETA),                   &! theta_in_wth
        arg_type(GH_FIELD,   GH_READ,      W3),                       &! wetrho_in_w3
        arg_type(GH_FIELD,   GH_READ,      WTHETA),                   &! wetrho_in_wth
        arg_type(GH_FIELD,   GH_READ,      W3),                       &! exner_in_w3
        arg_type(GH_FIELD,   GH_READ,      WTHETA),                   &! exner_in_wth
        arg_type(GH_FIELD,   GH_READ,      W2),                       &! u_physics
        arg_type(GH_FIELD,   GH_READ,      WTHETA),                   &! m_v_n
        arg_type(GH_FIELD,   GH_READ,      WTHETA),                   &! m_cl_n
        arg_type(GH_FIELD,   GH_READ,      WTHETA),                   &! m_ci_n
        arg_type(GH_FIELD,   GH_READ,      WTHETA),                   &! theta_star
        arg_type(GH_FIELD,   GH_READ,      W2),                       &! u_physics_star
        arg_type(GH_FIELD,   GH_READ,      W3),                       &! height_w3
        arg_type(GH_FIELD,   GH_READ,      WTHETA),                   &! height_wth
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! ntml_2d
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! cumulus_2d
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_fraction
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_3),&! leaf_area_index
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_3),&! canopy_height
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! peak_to_trough_orog
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! silhouette_area_orog
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_thickness
        arg_type(GH_FIELD,   GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_temperature
        arg_type(GH_FIELD,   GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! tile_temperature
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_snow_mass
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! n_snow_layers
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! snow_depth
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! canopy_water
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! soil_temperature
        arg_type(GH_FIELD,   GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! tile_heat_flux
        arg_type(GH_FIELD,   GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! tile_moisture_flux
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! sw_up_tile
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sw_down_surf
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! lw_down_surf
        arg_type(GH_FIELD,   GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! snow_sublimation  (kg m-2 s-1)
        arg_type(GH_FIELD,   GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! surf_heat_flux (W m-2)
        arg_type(GH_FIELD,   GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! canopy_evap (kg m-2 s-1)
        arg_type(GH_FIELD,   GH_WRITE,     ANY_DISCONTINUOUS_SPACE_5),&! water_extraction (kg m-2 s-1)
        arg_type(GH_FIELD,   GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! total_snowmelt (kg m-2 s-1)
        arg_type(GH_FIELD,   GH_WRITE,     WTHETA),                   &! dtheta_bl
        arg_type(GH_FIELD,   GH_INC,       W2),                       &! du_bl_w2
        arg_type(GH_FIELD,   GH_READ,      WTHETA),                   &! dt_conv
        arg_type(GH_FIELD,   GH_READ,      W2),                       &! du_conv_w2
        arg_type(GH_FIELD,   GH_READWRITE, WTHETA),                   &! m_v
        arg_type(GH_FIELD,   GH_READWRITE, WTHETA),                   &! m_cl
        arg_type(GH_FIELD,   GH_READWRITE, WTHETA),                   &! m_ci
        arg_type(GH_FIELD,   GH_READWRITE, WTHETA),                   &! cf_area
        arg_type(GH_FIELD,   GH_READWRITE, WTHETA),                   &! cf_ice
        arg_type(GH_FIELD,   GH_READWRITE, WTHETA),                   &! cf_liq
        arg_type(GH_FIELD,   GH_READWRITE, WTHETA),                   &! cf_bulk
        arg_type(GH_FIELD,   GH_READ,      WTHETA),                   &! rh_crit_wth
        arg_type(GH_FIELD,   GH_READ,      W2),                       &! rhokm_w2
        arg_type(GH_FIELD,   GH_READ,      W2),                       &! rhokm_surf_w2
        arg_type(GH_FIELD,   GH_READ,      W3),                       &! rhokh_bl
        arg_type(GH_FIELD,   GH_READ,      W2),                       &! ngstress_w2
        arg_type(GH_FIELD,   GH_READ,      WTHETA),                   &! bq_bl
        arg_type(GH_FIELD,   GH_READ,      WTHETA),                   &! bt_bl
        arg_type(GH_FIELD,   GH_READ,      W3),                       &! moist_flux_bl
        arg_type(GH_FIELD,   GH_READ,      W3),                       &! heat_flux_bl
        arg_type(GH_FIELD,   GH_READ,      W2),                       &! dtrdz_w2
        arg_type(GH_FIELD,   GH_READ,      WTHETA),                   &! dtrdz_tq_bl
        arg_type(GH_FIELD,   GH_READ,      W3),                       &! rdz_tq_bl
        arg_type(GH_FIELD,   GH_READ,      W2),                       &! rdz_w2
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! alpha1_tile
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! ashtf_prime_tile
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! dtstar_tile
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! fraca_tile
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! z0h_tile
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! z0m_tile
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! rhokh_tile
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! chr1p5m_tile
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! resfs_tile
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! canhc_tile
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! tile_water_extract
        arg_type(GH_FIELD,   GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! z_lcl
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! inv_depth
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! qcl_at_inv_top
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! blend_height_tq
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! blend_height_uv
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! ustar
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_avail
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! zh_nonloc
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! zh_2d
        arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_7) &! bl_types
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass :: bl_imp_code
  end type

  public bl_imp_code

contains

  !> @brief Interface to the implicit UM BL scheme
  !> @details The UM Boundary Layer scheme does:
  !>             vertical mixing of heat, momentum and moisture,
  !>             as documented in UMDP24
  !>          NB This version uses winds in w3 space (i.e. A-grid)
  !> @param[in]     nlayers              Number of layers
  !> @param[in]     outer                Outer loop counter
  !> @param[in]     theta_in_wth         Potential temperature field
  !> @param[in]     wetrho_in_w3         Wet density field in density space
  !> @param[in]     wetrho_in_wth        Wet density field in wth space
  !> @param[in]     exner_in_w3          Exner pressure field in density space
  !> @param[in]     exner_in_wth         Exner pressure field in wth space
  !> @param[in]     u_physics            Wind in native space at time n
  !> @param[in]     m_v_n                Vapour mixing ratio at time level n
  !> @param[in]     m_cl_n               Cloud liq mixing ratio at time level n
  !> @param[in]     m_ci_n               Cloud ice mixing ratio at time level n
  !> @param[in]     theta_star           Potential temperature after advection
  !> @param[in]     u_physics_star       Wind in native space after advection
  !> @param[in]     height_w3            Height of density space above surface
  !> @param[in]     height_wth           Height of theta space above surface
  !> @param[in]     ntml_2d              Number of turbulently mixed levels
  !> @param[in]     cumulus_2d           Cumulus flag (true/false)
  !> @param[in]     tile_fraction        Surface tile fractions
  !> @param[in]     leaf_area_index      Leaf Area Index
  !> @param[in]     canopy_height        Canopy height
  !> @param[in]     peak_to_trough_orog  Half of peak-to-trough height over root(2) of orography
  !> @param[in]     silhouette_area_orog Silhouette area of orography
  !> @param[in]     sea_ice_thickness    Depth of sea-ice (m)
  !> @param[in,out] sea_ice_temperature  Bulk temperature of sea-ice (K)
  !> @param[in,out] tile_temperature     Surface tile temperatures
  !> @param[in]     tile_snow_mass       Snow mass on tiles (kg/m2)
  !> @param[in]     n_snow_layers        Number of snow layers on tiles
  !> @param[in]     snow_depth           Snow depth on tiles
  !> @param[in]     canopy_water         Canopy water on each tile
  !> @param[in]     soil_temperature     Soil temperature
  !> @param[in,out] tile_heat_flux       Surface heat flux
  !> @param[in,out] tile_moisture_flux   Surface moisture flux
  !> @param[in]     sw_up_tile           Upwelling SW radiation on surface tiles
  !> @param[in]     sw_down_surf         Downwelling SW radiation at surface
  !> @param[in]     lw_down_surf         Downwelling LW radiation at surface
  !> @param[out]    snow_sublimation     Sublimation of snow
  !> @param[out]    surf_heat_flux       Surface heat flux
  !> @param[out]    canopy_evap          Canopy evaporation from land tiles
  !> @param[out]    water_extraction     Extraction of water from each soil layer
  !> @param[out]    total_snowmelt       Surface plus canopy snowmelt rate
  !> @param[out]    dtheta_bl            BL theta increment
  !> @param[in,out] du_bl_w2             BL wind increment
  !> @param[in]     dt_conv              Convection temperature increment
  !> @param[in]     du_conv_w2           Convection wind increment
  !> @param[in,out] m_v                  Vapour mixing ration after advection
  !> @param[in,out] m_cl                 Cloud liq mixing ratio after advection
  !> @param[in,out] m_ci                 Cloud ice mixing ratio after advection
  !> @param[in,out] cf_area              Area cloud fraction
  !> @param[in,out] cf_ice               Ice cloud fraction
  !> @param[in,out] cf_liq               Liquid cloud fraction
  !> @param[in,out] cf_bulk              Bulk cloud fraction
  !> @param[in]     rh_crit_wth          Critical relative humidity
  !> @param[in]     rhokm_w2             Momentum eddy diffusivity mapped to cell faces
  !> @param[in]     rhokm_surf_w2        Surface eddy diffusivity mapped to cell faces
  !> @param[in]     rhokh_bl             Heat eddy diffusivity on BL levels
  !> @param[in]     ngstress_w2          NG stress function mapped to cell faces
  !> @param[in]     bq_bl                Buoyancy parameter for moisture
  !> @param[in]     bt_bl                Buoyancy parameter for heat
  !> @param[in]     moist_flux_bl        Vertical moisture flux on BL levels
  !> @param[in]     heat_flux_bl         Vertical heat flux on BL levels
  !> @param[in]     dtrdz_w2             dt/(rho*r*r*dz) mapped to cell faces
  !> @param[in]     dtrdz_tq_bl          dt/(rho*r*r*dz) in wth
  !> @param[in]     rdz_tq_bl            1/dz in w3
  !> @param[in]     rdz_w2               1/dz mapped to cell faces
  !> @param[in]     alpha1_tile          dqsat/dT in surface layer on tiles
  !> @param[in]     ashtf_prime_tile     Heat flux coefficient on tiles
  !> @param[in]     dtstar_tile          Change in surface temperature on tiles
  !> @param[in]     fraca_tile           Fraction of moisture flux with only aerodynamic resistance
  !> @param[in]     z0h_tile             Heat roughness length on tiles
  !> @param[in]     z0m_tile             Momentum roughness length on tiles
  !> @param[in]     rhokh_tile           Surface heat diffusivity on tiles
  !> @param[in]     chr1p5m_tile         1.5m transfer coefficients on tiles
  !> @param[in]     resfs_tile           Combined aerodynamic resistance
  !> @param[in]     canhc_tile           Canopy heat capacity on tiles
  !> @param[in]     tile_water_extract   Extraction of water from each tile
  !> @param[in,out] z_lcl                Height of the LCL
  !> @param[in]     inv_depth            Depth of BL top inversion layer
  !> @param[in]     qcl_at_inv_top       Cloud water at top of inversion
  !> @param[in]     blend_height_tq      Blending height for wth levels
  !> @param[in]     blend_height_uv      Blending height for w3 levels
  !> @param[in]     ustar                Friction velocity
  !> @param[in]     soil_moist_avail     Available soil moisture for evaporation
  !> @param[in]     zh_nonloc            Depth of non-local BL scheme
  !> @param[in]     zh_2d                Total BL depth
  !> @param[in]     bl_types             Diagnosed BL types
  !> @param[in]     ndf_wth              Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth             Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth              Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_w3               Number of DOFs per cell for density space
  !> @param[in]     undf_w3              Number of unique DOFs for density space
  !> @param[in]     map_w3               Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_w2               Number of DOFs per cell for w2 space
  !> @param[in]     undf_w2              Number of unique DOFs for w2 space
  !> @param[in]     map_w2               Dofmap for the cell at the base of the column for w2 space
  !> @param[in]     ndf_2d               Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d              Number of unique DOFs for 2D fields
  !> @param[in]     map_2d               Dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_tile             Number of DOFs per cell for tiles
  !> @param[in]     undf_tile            Number of total DOFs for tiles
  !> @param[in]     map_tile             Dofmap for cell for surface tiles
  !> @param[in]     ndf_pft              Number of DOFs per cell for PFTs
  !> @param[in]     undf_pft             Number of total DOFs for PFTs
  !> @param[in]     map_pft              Dofmap for cell for PFTs
  !> @param[in]     ndf_sice             Number of DOFs per cell for sice levels
  !> @param[in]     undf_sice            Number of total DOFs for sice levels
  !> @param[in]     map_sice             Dofmap for cell for sice levels
  !> @param[in]     ndf_soil             Number of DOFs per cell for soil levels
  !> @param[in]     undf_soil            Number of total DOFs for soil levels
  !> @param[in]     map_soil             Dofmap for cell for soil levels
  !> @param[in]     ndf_smtile           Number of DOFs per cell for soil levels and tiles
  !> @param[in]     undf_smtile          Number of total DOFs for soil levels and tiles
  !> @param[in]     map_smtile           Dofmap for cell for soil levels and tiles
  !> @param[in]     ndf_bl               Number of DOFs per cell for BL types
  !> @param[in]     undf_bl              Number of total DOFs for BL types
  !> @param[in]     map_bl               Dofmap for cell for BL types
  subroutine bl_imp_code(nlayers,                            &
                         outer,                              &
                         theta_in_wth,                       &
                         wetrho_in_w3,                       &
                         wetrho_in_wth,                      &
                         exner_in_w3,                        &
                         exner_in_wth,                       &
                         u_physics,                          &
                         m_v_n,                              &
                         m_cl_n,                             &
                         m_ci_n,                             &
                         theta_star,                         &
                         u_physics_star,                     &
                         height_w3,                          &
                         height_wth,                         &
                         ntml_2d,                            &
                         cumulus_2d,                         &
                         tile_fraction,                      &
                         leaf_area_index,                    &
                         canopy_height,                      &
                         peak_to_trough_orog,                &
                         silhouette_area_orog,               &
                         sea_ice_thickness,                  &
                         sea_ice_temperature,                &
                         tile_temperature,                   &
                         tile_snow_mass,                     &
                         n_snow_layers,                      &
                         snow_depth,                         &
                         canopy_water,                       &
                         soil_temperature,                   &
                         tile_heat_flux,                     &
                         tile_moisture_flux,                 &
                         sw_up_tile,                         &
                         sw_down_surf,                       &
                         lw_down_surf,                       &
                         snow_sublimation,                   &
                         surf_heat_flux,                     &
                         canopy_evap,                        &
                         water_extraction,                   &
                         total_snowmelt,                     &
                         dtheta_bl,                          &
                         du_bl_w2,                           &
                         dt_conv,                            &
                         du_conv_w2,                         &
                         m_v,                                &
                         m_cl,                               &
                         m_ci,                               &
                         cf_area,                            &
                         cf_ice,                             &
                         cf_liq,                             &
                         cf_bulk,                            &
                         rh_crit_wth,                        &
                         rhokm_w2,                           &
                         rhokm_surf_w2,                      &
                         rhokh_bl,                           &
                         ngstress_w2,                        &
                         bq_bl,                              &
                         bt_bl,                              &
                         moist_flux_bl,                      &
                         heat_flux_bl,                       &
                         dtrdz_w2,                           &
                         dtrdz_tq_bl,                        &
                         rdz_tq_bl,                          &
                         rdz_w2,                             &
                         alpha1_tile,                        &
                         ashtf_prime_tile,                   &
                         dtstar_tile,                        &
                         fraca_tile,                         &
                         z0h_tile,                           &
                         z0m_tile,                           &
                         rhokh_tile,                         &
                         chr1p5m_tile,                       &
                         resfs_tile,                         &
                         canhc_tile,                         &
                         tile_water_extract,                 &
                         z_lcl,                              &
                         inv_depth,                          &
                         qcl_at_inv_top,                     &
                         blend_height_tq,                    &
                         blend_height_uv,                    &
                         ustar,                              &
                         soil_moist_avail,                   &
                         zh_nonloc,                          &
                         zh_2d,                              &
                         bl_types,                           &
                         ndf_wth,                            &
                         undf_wth,                           &
                         map_wth,                            &
                         ndf_w3,                             &
                         undf_w3,                            &
                         map_w3,                             &
                         ndf_w2,                             &
                         undf_w2,                            &
                         map_w2,                             &
                         ndf_2d,                             &
                         undf_2d,                            &
                         map_2d,                             &
                         ndf_tile, undf_tile, map_tile,      &
                         ndf_pft, undf_pft, map_pft,         &
                         ndf_sice, undf_sice, map_sice,      &
                         ndf_soil, undf_soil, map_soil,      &
                         ndf_smtile, undf_smtile, map_smtile,&
                         ndf_bl, undf_bl, map_bl)

    !---------------------------------------
    ! LFRic modules
    !---------------------------------------
    use jules_control_init_mod, only: n_land_tile, n_sea_ice_tile, &
         first_sea_tile, first_sea_ice_tile

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use ancil_info, only: ssi_pts,                                          &
         sea_pts, sice_pts, ssi_index, sea_index, sice_index, fssi_ij,      &
         sea_frac, sice_frac, sice_pts_ncat, sice_index_ncat, sice_frac_ncat
    use atm_fields_bounds_mod, only: udims, vdims, udims_s, vdims_s,        &
         pdims
    use atm_step_local, only: dim_cs1, dim_cs2
    use c_kappai, only: kappai, de
    use dust_parameters_mod, only: ndiv, ndivh
    use jules_sea_seaice_mod, only: nice_use
    use jules_snow_mod, only: cansnowtile, rho_snow_const, l_snowdep_surf
    use jules_surface_types_mod, only: npft, ntype, lake, nnvg
    use jules_vegetation_mod, only: can_model
    use nlsizes_namelist_mod, only: row_length, rows, land_field, &
         sm_levels, ntiles, bl_levels, tr_vars
    use pftparm, only: emis_pft
    use planet_constants_mod, only: p_zero, kappa, planet_radius
    use nvegparm, only: emis_nvg
    use rad_input_mod, only: co2_mmr

    ! spatially varying fields used from modules
    use fluxes, only: sw_sicat
    use level_heights_mod, only: r_theta_levels, r_rho_levels
    use prognostics, only: nsnow_surft
    use jules_mod, only: snowdep_surft

    ! subroutines used
    use bdy_expl3_mod, only: bdy_expl3
    use bl_diags_mod, only: bl_diag, dealloc_bl_imp, alloc_bl_expl
    use sf_diags_mod, only: sf_diag, dealloc_sf_expl, dealloc_sf_imp,       &
                            alloc_sf_expl
    use ni_imp_ctl_mod, only: ni_imp_ctl
    use tilepts_mod, only: tilepts

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: outer

    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3, ndf_w2
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d, undf_w2
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3
    integer(kind=i_def), intent(in) :: map_wth(ndf_wth)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)
    integer(kind=i_def), intent(in) :: map_w2(ndf_w2)

    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile)
    integer(kind=i_def), intent(in) :: ndf_pft, undf_pft
    integer(kind=i_def), intent(in) :: map_pft(ndf_pft)
    integer(kind=i_def), intent(in) :: ndf_sice, undf_sice
    integer(kind=i_def), intent(in) :: map_sice(ndf_sice)
    integer(kind=i_def), intent(in) :: ndf_soil, undf_soil
    integer(kind=i_def), intent(in) :: map_soil(ndf_soil)
    integer(kind=i_def), intent(in) :: ndf_smtile, undf_smtile
    integer(kind=i_def), intent(in) :: map_smtile(ndf_smtile)

    integer(kind=i_def), intent(in) :: ndf_bl, undf_bl
    integer(kind=i_def), intent(in) :: map_bl(ndf_bl)

    real(kind=r_def), dimension(undf_w2),  intent(inout):: du_bl_w2
    real(kind=r_def), dimension(undf_wth), intent(out)  :: dtheta_bl
    real(kind=r_def), dimension(undf_wth), intent(inout):: m_v, m_cl, m_ci,    &
                                                           cf_area, cf_ice,    &
                                                           cf_liq, cf_bulk
    real(kind=r_def), dimension(undf_w3),  intent(in)   :: wetrho_in_w3,       &
                                                           exner_in_w3,        &
                                                           height_w3,          &
                                                           rhokh_bl,           &
                                                           moist_flux_bl,      &
                                                           heat_flux_bl,       &
                                                           rdz_tq_bl
    real(kind=r_def), dimension(undf_wth), intent(in)   :: theta_in_wth,       &
                                                           wetrho_in_wth,      &
                                                           exner_in_wth,       &
                                                           m_v_n, m_cl_n,      &
                                                           m_ci_n,             &
                                                           theta_star,         &
                                                           height_wth,         &
                                                           dt_conv,            &
                                                           rh_crit_wth,        &
                                                           bq_bl, bt_bl,       &
                                                           dtrdz_tq_bl

    real(kind=r_def), dimension(undf_w2),  intent(in)   :: u_physics,         &
                                                           u_physics_star,    &
                                                           rhokm_w2,          &
                                                           ngstress_w2,       &
                                                           dtrdz_w2,          &
                                                           rdz_w2,            &
                                                           rhokm_surf_w2,     &
                                                           du_conv_w2

    real(kind=r_def), dimension(undf_2d), intent(in) :: ntml_2d,              &
                                                        cumulus_2d, zh_2d,    &
                                                        blend_height_tq,      &
                                                        blend_height_uv,      &
                                                        ustar,                &
                                                        soil_moist_avail,     &
                                                        zh_nonloc

    real(kind=r_def), intent(in)    :: tile_fraction(undf_tile)
    real(kind=r_def), intent(inout) :: tile_temperature(undf_tile)
    real(kind=r_def), intent(in)    :: tile_snow_mass(undf_tile)
    real(kind=r_def), intent(in)    :: n_snow_layers(undf_tile)
    real(kind=r_def), intent(in)    :: snow_depth(undf_tile)
    real(kind=r_def), intent(in)    :: canopy_water(undf_tile)
    real(kind=r_def), intent(inout) :: tile_heat_flux(undf_tile)
    real(kind=r_def), intent(inout) :: tile_moisture_flux(undf_tile)
    real(kind=r_def), intent(in)    :: sw_up_tile(undf_tile)
    real(kind=r_def), intent(out)   :: snow_sublimation(undf_tile)
    real(kind=r_def), intent(out)   :: surf_heat_flux(undf_tile)
    real(kind=r_def), intent(out)   :: canopy_evap(undf_tile)
    real(kind=r_def), intent(out)   :: total_snowmelt(undf_tile)

    real(kind=r_def), intent(in) :: leaf_area_index(undf_pft)
    real(kind=r_def), intent(in) :: canopy_height(undf_pft)

    real(kind=r_def), intent(in) :: sea_ice_thickness(undf_sice)
    real(kind=r_def), intent(inout) :: sea_ice_temperature(undf_sice)

    real(kind=r_def), intent(in) :: peak_to_trough_orog(undf_2d)
    real(kind=r_def), intent(in) :: silhouette_area_orog(undf_2d)
    real(kind=r_def), intent(in) :: sw_down_surf(undf_2d)
    real(kind=r_def), intent(in) :: lw_down_surf(undf_2d)
    real(kind=r_def), intent(inout) :: z_lcl(undf_2d)
    real(kind=r_def), intent(in) :: inv_depth(undf_2d)
    real(kind=r_def), intent(in) :: qcl_at_inv_top(undf_2d)

    real(kind=r_def), intent(in) :: soil_temperature(undf_soil)
    real(kind=r_def), intent(out):: water_extraction(undf_soil)

    real(kind=r_def), intent(in) :: tile_water_extract(undf_smtile)

    real(kind=r_def), dimension(undf_bl),   intent(in)  :: bl_types
    real(kind=r_def), dimension(undf_tile), intent(in)  :: alpha1_tile,      &
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
    integer(i_def) :: k, i, i_tile, i_sice, i_pft, n

    ! local switches and scalars
    integer(i_um) :: error_code
    logical :: l_calc_at_p

    ! profile fields from level 1 upwards
    real(r_um), dimension(row_length,rows,nlayers) ::                        &
         p_rho_levels, rho_wet_rsq, rho_wet, z_rho, z_theta,                 &
         bulk_cloud_fraction, rhcpt, t_latest, q_latest, qcl_latest,         &
         qcf_latest, cca_3d, area_cloud_fraction,                            &
         cloud_fraction_liquid, cloud_fraction_frozen, rho_wet_tq,           &
         cf_latest, cfl_latest, cff_latest, zeros

    ! profile field on boundary layer levels
    real(r_um), dimension(row_length,rows,bl_levels) :: fqw, ftl, rhokh,     &
         bq_gb, bt_gb, dtrdz_charney_grid, rdz_charney_grid, rhokm

    ! profile fields on u/v points and all levels
    real(r_um), dimension(row_length,rows,nlayers) :: u, v, r_u, r_v

    ! profile fields on u/v points and BL levels
    real(r_um), dimension(row_length,rows,bl_levels) :: taux, tauy,          &
         dtrdz_u, dtrdz_v, rhokm_u, rhokm_v

    ! profile fields from level 2 upwards
    real(r_um), dimension(row_length,rows,2:bl_levels) :: rdz_u, rdz_v,      &
         f_ngstress_u, f_ngstress_v

    ! profile fields from level 0 upwards
    real(r_um), dimension(row_length,rows,0:nlayers) ::                      &
         p_theta_levels, R_w, p_rho_minus_one, w,                            &
         q, qcl, qcf, theta, exner_theta_levels, co2

    ! single level real fields
    real(r_um), dimension(row_length,rows) ::                                &
         p_star, lw_down, tstar, tstar_sea, zh, tstar_land, tstar_ssi,       &
         dtstar_sea, ice_fract, tstar_sice, alpha1_sea,                      &
         ashtf_prime_sea, bl_type_1, bl_type_2, bl_type_3, bl_type_4,        &
         bl_type_5, bl_type_6, bl_type_7, chr1p5m_sice, flandg, rhokh_sea,   &
         u_s, z0hssi, z0mssi, zhnl, zlcl_mix, zlcl, dzh, qcl_inv_top

    ! single level real fields on u/v points
    real(r_um), dimension(row_length,rows) :: u_0, v_0, rhokm_u_land,        &
         rhokm_u_ssi, rhokm_v_land, rhokm_v_ssi, taux_land, tauy_land,       &
         taux_ssi, tauy_ssi, flandg_u, flandg_v

    ! single level integer fields
    integer(i_um), dimension(row_length,rows) :: ntml, lcbase, k_blend_tq,   &
         k_blend_uv

    ! single level logical fields
    logical, dimension(row_length,rows) :: land_sea_mask, cumulus

    ! fields on sea-ice categories
    real(r_um), dimension(row_length,rows,nice_use) ::                       &
         tstar_sice_ncat, fqw_ice, ftl_ice, ice_fract_ncat, alpha1_sice,     &
         ashtf_prime, rhokh_sice, k_sice_ncat, ti_sice_ncat, di_sice_ncat,   &
         dtstar_sice

    ! field on land points and soil levels
    real(r_um), dimension(land_field,sm_levels) :: t_soil_soilt, ext

    ! real fields on land points
    real(r_um), dimension(land_field) :: sil_orog_land_gb, ho2r2_orog_gb,    &
         gs_gb, fland, smc_soilt

    ! integer fields on land points
    integer, dimension(land_field) :: land_index

    ! integer fields on land points and tile types
    integer, dimension(land_field, ntype) :: surft_index

    ! integer fields on number of tile types
    integer, dimension(ntype) :: surft_pts

    ! fields on land points and surface tiles
    real(r_um), dimension(land_field,ntiles) :: canopy_surft, snow_surft,    &
         sw_surft, tstar_surft, frac_surft, ftl_surft, fqw_surft,            &
         dtstar_surft, alpha1, ashtf_prime_surft, chr1p5m, fraca, resfs,     &
         rhokh_surft, z0h_surft, z0m_surft, resft, flake, emis_surft,        &
         canhc_surft, ei_surft, ecan_surft, melt_surft, surf_htf_surft

    ! fields on land points and plant functional types
    real(r_um), dimension(land_field,npft) :: canht_pft, lai_pft

    ! field on surface tiles and soil levels
    real(r_um), dimension(land_field,sm_levels,ntiles) :: wt_ext_surft

    ! Fields which are not used and only required for subroutine argument list,
    ! hence are unset in the kernel
    ! if they become set, please move up to be with other variables
    integer(i_um), parameter :: nscmdpkgs=15
    logical,       parameter :: l_scmdiags(nscmdpkgs)=.false.

    real(r_um), dimension(row_length,rows,nlayers) :: cca0

    real(r_um), dimension(row_length,rows,bl_levels) :: taux_fd_u, tauy_fd_v,&
         rhogamu_u, rhogamv_v

    real(r_um), dimension(row_length,rows,0:nlayers) ::                      &
         aerosol, dust_div1, dust_div2, dust_div3, dust_div4, dust_div5,     &
         dust_div6, so2, dms, so4_aitken, so4_accu, so4_diss, nh3, soot_new, &
         soot_aged, soot_cld, bmass_new, bmass_aged, bmass_cld, ocff_new,    &
         ocff_aged, ocff_cld, nitr_acc, nitr_diss, ozone_tracer, qrain

    real(r_um), dimension(row_length,rows,0:nlayers,tr_vars) :: free_tracers

    real(r_um), dimension(row_length,rows) :: ti_sice,                       &
         xx_cos_theta_latitude, ls_rain, ls_snow, conv_rain, conv_snow,      &
         co2_emits, co2flux, tscrndcl_ssi, tstbtrans,                        &
         sum_eng_fluxes, sum_moist_flux, drydep2, olr, surf_ht_flux_land,    &
         theta_star_surf, qv_star_surf, snowmelt, uwind_wav,                 &
         vwind_wav, sstfrz, aresist, resist_b, rho_aresist, rib_gb,          &
         z0m_eff_gb, zhsc, cdr10m_u, cdr10m_v

    real(r_um), dimension(row_length,rows,3) :: t_frac, t_frac_dsc, we_lim,  &
         we_lim_dsc, zrzi, zrzi_dsc

    real(r_um), dimension(row_length,rows) :: flandfac_u, flandfac_v,        &
         fseafac_u, fseafac_v

    integer(i_um), dimension(row_length,rows) :: ccb0, cct0, kent, kent_dsc

    real(r_um), dimension(row_length,rows,nice_use) :: radnet_sice

    real(r_um), dimension(row_length,rows,ndiv) :: dust_flux, r_b_dust

    real(r_um), dimension(dim_cs2) :: resp_s_tot_soilt

    real(r_um), dimension(land_field,dim_cs1) :: resp_s_gb_um, cs_pool_gb_um

    real(r_um), dimension(land_field) :: npp_gb

    real(r_um), dimension(land_field,ntiles) :: catch_surft,                 &
         tscrndcl_surft, epot_surft, aresist_surft, dust_emiss_frac,         &
         gc_surft, resist_b_surft, u_s_std_surft

    real(r_um), dimension(land_field,ntiles,ndivh) :: u_s_t_dry_tile,        &
         u_s_t_tile

    real(r_um) :: stashwork3(1), stashwork9(1)

    !-----------------------------------------------------------------------
    ! Initialisation of variables and arrays
    !-----------------------------------------------------------------------
    error_code=0

    zeros=0.0_r_um

    call alloc_sf_expl(sf_diag, outer == outer_iterations)
    call alloc_bl_expl(bl_diag, outer == outer_iterations)

    if (bl_diag%l_tke) then
      allocate(bl_diag%tke(pdims%i_start:pdims%i_end,                   &
                           pdims%j_start:pdims%j_end,bl_levels))
      do k = 1, bl_levels
        bl_diag%tke(:,:,k) = 0.0
      end do
    else
      allocate(bl_diag%tke(1,1,1))
    end if

    if (bl_diag%l_elm3d) then
      allocate(bl_diag%elm3d(pdims%i_start:pdims%i_end,                 &
                           pdims%j_start:pdims%j_end,bl_levels))
      do k = 1, bl_levels
        bl_diag%elm3d(:,:,k) = 0.0
      end do
    else
      allocate(bl_diag%elm3d(1,1,1))
    end if


    ! Needed to stop crash in tr_mix - will ultimately need passing from
    ! bl_exp_kernel to here, but science is currently unused
    kent = 2
    kent_dsc = 2

    !-----------------------------------------------------------------------
    ! Mapping of LFRic fields into UM variables
    !-----------------------------------------------------------------------

    ! Land tile fractions
    flandg = 0.0_r_um
    do i = 1, n_land_tile
      flandg = flandg + real(tile_fraction(map_tile(i)), r_um)
      frac_surft(1, i) = real(tile_fraction(map_tile(i)), r_um)
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
    tstar_land = 0.0_r_um
    do i = 1, n_land_tile
      tstar_surft(1, i) = real(tile_temperature(map_tile(i)), r_um)
      tstar_land = tstar_land + frac_surft(1, i) * tstar_surft(1, i)
      ! sensible heat flux
      ftl_surft(1, i) = real(tile_heat_flux(map_tile(i)), r_um)
      ! moisture flux
      fqw_surft(1, i) = real(tile_moisture_flux(map_tile(i)), r_um)
    end do

    ! Sea temperature
    tstar_sea = real(tile_temperature(map_tile(first_sea_tile)), r_um)

    ! Sea-ice temperatures
    i_sice = 0
    tstar_sice = 0.0_r_um
    if (ice_fract(1, 1) > 0.0_r_um) then
      do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
        i_sice = i_sice + 1
        tstar_sice_ncat(1, 1, i_sice) = real(tile_temperature(map_tile(i)), r_um)
        tstar_sice = tstar_sice &
                   + ice_fract_ncat(1,1,i_sice) * tstar_sice_ncat(1,1,i_sice) &
                   / ice_fract
        ! sea-ice heat flux
        ftl_ice(1,1,i_sice) = real(tile_heat_flux(map_tile(i)), r_um)
        ! sea-ice moisture flux
        fqw_ice(1,1,i_sice) = real(tile_moisture_flux(map_tile(i)), r_um)
      end do
    end if

    ! Sea & Sea-ice temperature
    tstar_ssi = (1.0_r_um - ice_fract) * tstar_sea + ice_fract * tstar_sice

    ! Grid-box mean surface temperature
    tstar = flandg * tstar_land + (1.0_r_um - flandg) * tstar_ssi

    ! Sea-ice conductivity, bulk temperature and thickness
    do i = 1, n_sea_ice_tile
      k_sice_ncat(1, 1, i) = 2.0_r_um * kappai / de
      ti_sice_ncat(1, 1, i) = real(sea_ice_temperature(map_sice(i)), r_um)
      di_sice_ncat(1, 1, i) = real(sea_ice_thickness(map_sice(i)), r_um)
    end do

    do i_pft = 1, npft
      ! Leaf area index
      lai_pft(1, i_pft) = real(leaf_area_index(map_pft(i_pft)), r_um)
      ! Canopy height
      canht_pft(1, i_pft) = real(canopy_height(map_pft(i_pft)), r_um)
    end do

    do i = 1, sm_levels
      ! Soil temperature (t_soil_soilt)
      t_soil_soilt(1, i) = real(soil_temperature(map_soil(i)), r_um)
    end do

    ! Downwelling LW radiation at surface
    lw_down = real(lw_down_surf(map_2d(1)), r_um)

    ! Net SW radiation on tiles
    do i = 1, n_land_tile
      sw_surft(1, i) = real(sw_down_surf(map_2d(1)) - &
                            sw_up_tile(map_tile(i)), r_um)
    end do

    ! Net SW on sea-ice
    i_sice = 0
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      sw_sicat(1, i_sice) = real(sw_down_surf(map_2d(1)) - &
                                 sw_up_tile(map_tile(i)), r_um)
    end do

    ! Carbon dioxide
    co2 = co2_mmr

    ! Half of peak-to-trough height over root(2) of orography (ho2r2_orog_gb)
    ho2r2_orog_gb = real(peak_to_trough_orog(map_2d(1)), r_um)

    sil_orog_land_gb = real(silhouette_area_orog(map_2d(1)), r_um)

    ! Canopy water on each tile (canopy_surft)
    do i = 1, n_land_tile
      canopy_surft(1, i) = real(canopy_water(map_tile(i)), r_um)
    end do

    do i = 1, n_land_tile
      ! Lying snow mass on land tiles
      snow_surft(1, i) = real(tile_snow_mass(map_tile(i)), r_um)
      ! Number of snow layers on tiles (nsnow_surft)
      nsnow_surft(1, i) = int(n_snow_layers(map_tile(i)), i_um)
      ! Equivalent snowdepth for surface calculations.
      ! code copied from jules_land_sf_explicit
      ! 4 is a magic number inherited from Jules, meaning radiative canopy
      ! with heat capacity and snow beneath
      if ( (can_model == 4) .and. cansnowtile(i) .and. l_snowdep_surf) then
        snowdep_surft(1, i) = snow_surft(1, i) / rho_snow_const
      else
        snowdep_surft(1, i) = real(snow_depth(map_tile(i)), r_um)
      end if
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
      ! pressure on rho and theta levels
      p_rho_levels(1,1,k) = p_zero*(exner_in_w3(map_w3(1) + k-1))**(1.0_r_def/kappa)
      p_theta_levels(1,1,k) = p_zero*(exner_in_wth(map_wth(1) + k))**(1.0_r_def/kappa)
      ! exner pressure on theta levels
      exner_theta_levels(1,1,k) = exner_in_wth(map_wth(1) + k)
      ! u wind on rho levels
      u(1,1,k) = u_physics(map_w2(1) + k-1)
      ! v wind on rho levels
      v(1,1,k) = u_physics(map_w2(2) + k-1)
      ! w wind on theta levels
      w(1,1,k) = u_physics(map_w2(5) + k)
      ! height of rho levels from centre of planet
      r_rho_levels(1,1,k) = height_w3(map_w3(1) + k-1) + planet_radius
      ! height of theta levels from centre of planet
      r_theta_levels(1,1,k) = height_wth(map_wth(1) + k) + planet_radius
      ! water vapour mixing ratio
      q(1,1,k) = m_v_n(map_wth(1) + k)
      ! cloud liquid mixing ratio
      qcl(1,1,k) = m_cl_n(map_wth(1) + k)
      ! cloud ice mixing ratio
      qcf(1,1,k) = m_ci_n(map_wth(1) + k)
      ! 3D RH_crit field
      rhcpt(1,1,k) = rh_crit_wth(map_wth(1) + k)
    end do

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
    ! surface currents
    u_0 = 0.0
    v_0 = 0.0
    ! near surface moisture fields
    q(1,1,0) = m_v_n(map_wth(1) + 0)
    qcl(1,1,0) = m_cl_n(map_wth(1) + 0)
    qcf(1,1,0) = m_ci_n(map_wth(1) + 0)
    ! surface height
    r_theta_levels(1,1,0) = height_wth(map_wth(1) + 0) + planet_radius
    ! height of levels above surface
    z_rho = r_rho_levels-r_theta_levels(1,1,0)
    z_theta(1,1,:) = r_theta_levels(1,1,1:nlayers)-r_theta_levels(1,1,0)
    ! vertical velocity
    w(1,1,0) = u_physics(map_w2(5) + 0)

    !-----------------------------------------------------------------------
    ! Things passed from other parametrization schemes on this timestep
    !-----------------------------------------------------------------------
    cumulus(1,1) = (cumulus_2d(map_2d(1)) > 0.5_r_def)
    ntml(1,1) = INT(ntml_2d(map_2d(1)))

    rhokm_u_land(1,1) = rhokm_surf_w2(map_w2(1) + 1)
    rhokm_u_ssi(1,1) = rhokm_surf_w2(map_w2(1) + 2)
    flandg_u(1,1) = rhokm_surf_w2(map_w2(1) + 3)
    rhokm_v_land(1,1) = rhokm_surf_w2(map_w2(2) + 1)
    rhokm_v_ssi(1,1) = rhokm_surf_w2(map_w2(2) + 2)
    flandg_v(1,1) = rhokm_surf_w2(map_w2(2) + 3)
    do k = 1, bl_levels
      rhokm_u(1,1,k) = rhokm_w2(map_w2(1) + k)
      rhokm_v(1,1,k) = rhokm_w2(map_w2(2) + k)
      rhokh(1,1,k) = rhokh_bl(map_w3(1) + k)
      bq_gb(1,1,k) = bq_bl(map_wth(1) + k)
      bt_gb(1,1,k) = bt_bl(map_wth(1) + k)
      fqw(1,1,k) = moist_flux_bl(map_w3(1) + k)
      ftl(1,1,k) = heat_flux_bl(map_w3(1) + k)
      dtrdz_u(1,1,k) = dtrdz_w2(map_w2(1) + k)
      dtrdz_v(1,1,k) = dtrdz_w2(map_w2(2) + k)
      dtrdz_charney_grid(1,1,k) = dtrdz_tq_bl(map_wth(1) + k)
      rdz_charney_grid(1,1,k) = rdz_tq_bl(map_w3(1) + k)
    end do
    do k = 2, bl_levels
      rdz_u(1,1,k) = rdz_w2(map_w2(1) + k)
      rdz_v(1,1,k) = rdz_w2(map_w2(2) + k)
      f_ngstress_u(1,1,k) = ngstress_w2(map_w2(1) + k)
      f_ngstress_v(1,1,k) = ngstress_w2(map_w2(2) + k)
    end do

    do i = 1, n_land_tile
      alpha1(1, i) = alpha1_tile(map_tile(i))
      ashtf_prime_surft(1, i) = ashtf_prime_tile(map_tile(i))
      dtstar_surft(1, i) = dtstar_tile(map_tile(i))
      fraca(1, i) = fraca_tile(map_tile(i))
      z0h_surft(1, i) = z0h_tile(map_tile(i))
      z0m_surft(1, i) = z0m_tile(map_tile(i))
      rhokh_surft(1, i) = rhokh_tile(map_tile(i))
      chr1p5m(1, i) = chr1p5m_tile(map_tile(i))
      resfs(1, i) = resfs_tile(map_tile(i))
      canhc_surft(1, i) = canhc_tile(map_tile(i))
    end do
    ! recalculate the total resistance factor
    resft = fraca + (1.0 - fraca) * resfs
    resft(1,lake) = 1.0
    ! recalculate the total lake fraction
    flake = 0.0
    flake(1,lake) = 1.0
    ! recalculate the surface emissivity
    emis_surft(1,1:npft) = emis_pft(1:npft)
    emis_surft(1,npft+1:npft+nnvg) = emis_nvg(1:nnvg)

    i_tile = 0
    do i = 1, n_land_tile
      do n = 1, sm_levels
        i_tile = i_tile + 1
        wt_ext_surft(1,n,i) = tile_water_extract(map_smtile(i_tile))
      end do
    end do

    alpha1_sea(1,1) = alpha1_tile(map_tile(first_sea_tile))
    ashtf_prime_sea(1,1) = ashtf_prime_tile(map_tile(first_sea_tile))
    dtstar_sea(1,1) = dtstar_tile(map_tile(first_sea_tile))
    rhokh_sea(1,1) = rhokh_tile(map_tile(first_sea_tile))

    i_sice = 0
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      alpha1_sice(1,1,i_sice) = alpha1_tile(map_tile(i))
      ashtf_prime(1,1,i_sice) = ashtf_prime_tile(map_tile(i))
      rhokh_sice(1,1,i_sice) = rhokh_tile(map_tile(i))
      dtstar_sice(1,1,i_sice) = dtstar_tile(map_tile(i))
    end do
    z0hssi(1,1) = z0h_tile(map_tile(first_sea_ice_tile))
    z0mssi(1,1) = z0m_tile(map_tile(first_sea_ice_tile))
    chr1p5m_sice(1,1) = chr1p5m_tile(map_tile(first_sea_ice_tile))

    k_blend_tq(1,1) = int(blend_height_tq(map_2d(1)), i_um)
    k_blend_uv(1,1) = int(blend_height_uv(map_2d(1)), i_um)
    u_s(1,1) = ustar(map_2d(1))
    smc_soilt(1) = soil_moist_avail(map_2d(1))
    zhnl(1,1) = zh_nonloc(map_2d(1))
    zh(1,1) = zh_2d(map_2d(1))
    zlcl(1,1) = real(z_lcl(map_2d(1)), r_um)
    dzh(1,1) = real(inv_depth(map_2d(1)), r_um)
    qcl_inv_top(1,1) = real(qcl_at_inv_top(map_2d(1)), r_um)
    zlcl_mix = 0.0_r_um

    bl_type_1(1,1) = bl_types(map_bl(1))
    bl_type_2(1,1) = bl_types(map_bl(2))
    bl_type_3(1,1) = bl_types(map_bl(3))
    bl_type_4(1,1) = bl_types(map_bl(4))
    bl_type_5(1,1) = bl_types(map_bl(5))
    bl_type_6(1,1) = bl_types(map_bl(6))
    bl_type_7(1,1) = bl_types(map_bl(7))

    l_calc_at_p = .FALSE.  ! Flag for calling this routine on p-grid
    call bdy_expl3 (                                                 &
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
    lcbase(1,1) = 0     ! default for no convection and convective cloud

    !-----------------------------------------------------------------------
    ! increments / fields with increments added
    !-----------------------------------------------------------------------
    do k = 1, nlayers
      t_latest(1,1,k) = theta_star(map_wth(1) + k) * exner_theta_levels(1,1,k) &
                      + dt_conv(map_wth(1) + k)
      r_u(1,1,k) = u_physics_star(map_w2(1) + k-1)                             &
                   - u_physics(map_w2(1) + k-1) + du_conv_w2(map_w2(1) + k-1)
      r_v(1,1,k) = u_physics_star(map_w2(2) + k-1)                             &
                   - u_physics(map_w2(2) + k-1) + du_conv_w2(map_w2(2) + k-1)
      r_w(1,1,k) = u_physics_star(map_w2(5) + k) - u_physics(map_w2(5) + k)
      q_latest(1,1,k)   = m_v(map_wth(1) + k)
      qcl_latest(1,1,k) = m_cl(map_wth(1) + k)
      qcf_latest(1,1,k) = m_ci(map_wth(1) + k)
    end do
    r_w(1,1,0) = u_physics_star(map_w2(5) + 0) - u_physics(map_w2(5) + 0)

    if (scheme == scheme_pc2) then
      do k = 1, nlayers
        ! Assign _latest with current updated values
        cfl_latest(1,1,k) = cf_liq(map_wth(1) + k)
        cff_latest(1,1,k) = cf_ice(map_wth(1) + k)
        cf_latest(1,1,k)  = cf_bulk(map_wth(1) + k)
      end do
    end if

    call NI_imp_ctl (                                                   &
    ! IN Model switches
            outer                                                       &
    ! IN trig arrays
          , xx_cos_theta_latitude                                       &
    ! IN data fields.
          , p_theta_levels, p_rho_minus_one, rho_wet_rsq, rho_wet_tq    &
          , u, v, w                                                     &
          , land_sea_mask, q, qcl, qcf, p_star, theta, qrain            &
          , exner_theta_levels                                          &
    ! IN ancillary fields and fields needed to be kept from tstep to tstep
          , sil_orog_land_gb, ho2r2_orog_gb                             &
          , ice_fract, di_sice_ncat, ice_fract_ncat, k_sice_ncat        &
          , u_0, v_0, land_index, cca_3d, lcbase, ccb0, cct0            &
          , ls_rain, ls_snow, conv_rain, conv_snow                      &
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
          , surft_pts,surft_index,frac_surft,canopy_surft               &
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
    ! SCM Diagnostics (dummy in full UM) & BL diags
          , nSCMDpkgs, L_SCMDiags, bl_diag, sf_diag                     &
    ! INOUT (Note ti_sice_ncat and ti_sice are IN if l_sice_multilayers=T)
          , TScrnDcl_SSI, TScrnDcl_surft, tStbTrans                     &
          , cca0, fqw, ftl, taux, tauy, rhokh                           &
          , fqw_ice,ftl_ice,dtstar_surft,dtstar_sea,dtstar_sice,ti_sice_ncat &
          , area_cloud_fraction, bulk_cloud_fraction                    &
          , t_latest, q_latest, qcl_latest, qcf_latest                  &
          , cf_latest, cfl_latest, cff_latest                           &
          , R_u, R_v, R_w, cloud_fraction_liquid, cloud_fraction_frozen &
          , sum_eng_fluxes,sum_moist_flux, rhcpt                        &
    ! Dummy arguments for bimodal scheme compliance
          , zeros,zeros,zeros,zeros,zeros                               &
    ! INOUT tracer fields
          , aerosol, free_tracers,  resist_b,  resist_b_surft           &
          , dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6 &
          , drydep2, so2, dms, so4_aitken, so4_accu, so4_diss, nh3      &
          , soot_new, soot_aged, soot_cld, bmass_new, bmass_aged        &
          , bmass_cld, ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss &
          , co2, ozone_tracer                                           &
    ! INOUT additional variables for JULES
          , tstar_surft, fqw_surft, epot_surft, ftl_surft               &
          , radnet_sice,olr,tstar_sice_ncat,tstar_ssi                   &
          , tstar_sea,taux_land,taux_ssi,tauy_land,tauy_ssi,Error_code  &
    ! OUT fields
          , surf_ht_flux_land, zlcl_mix                                 &
          , theta_star_surf, qv_star_surf                               &
    ! OUT additional variables for JULES
          , tstar, ti_sice, ext, snowmelt,tstar_land,tstar_sice, ei_surft &
          , ecan_surft, melt_surft, surf_htf_surft                      &
    ! OUT fields for coupling to the wave model
          , uwind_wav, vwind_wav                                        &
            )

    !--------------------------------------------------------------
    ! update the first two dofs of the wind increment
    !-------------------------------------------------------------
    do k = 1, nlayers
      du_bl_w2(map_w2(1) + k - 1) = r_u(1,1,k)                          &
         - (u_physics_star(map_w2(1) + k-1) - u_physics(map_w2(1) + k-1))
      du_bl_w2(map_w2(2) + k - 1) = r_v(1,1,k)                          &
         - (u_physics_star(map_w2(2) + k-1) - u_physics(map_w2(2) + k-1))
    end do

    !---------------------------------------------------------------
    ! now need to re-run using the 2nd two dofs
    !---------------------------------------------------------------

    ! First we need to re-set everythign which has been over-written or
    ! allocated/deallocated
    call dealloc_bl_imp(bl_diag)
    call dealloc_sf_imp(sf_diag)
    call alloc_bl_expl(bl_diag, outer == outer_iterations)

    if (bl_diag%l_tke) then
      allocate(bl_diag%tke(pdims%i_start:pdims%i_end,                   &
                           pdims%j_start:pdims%j_end,bl_levels))
      do k = 1, bl_levels
        bl_diag%tke(:,:,k) = 0.0
      end do
    else
      allocate(bl_diag%tke(1,1,1))
    end if

    if (bl_diag%l_elm3d) then
      allocate(bl_diag%elm3d(pdims%i_start:pdims%i_end,                &
                           pdims%j_start:pdims%j_end,bl_levels))
      do k = 1, bl_levels
        bl_diag%elm3d(:,:,k) = 0.0
      end do
    else
      allocate(bl_diag%elm3d(1,1,1))
    end if

    ! Land tile temperatures
    tstar_land = 0.0_r_um
    do i = 1, n_land_tile
      tstar_surft(1, i) = real(tile_temperature(map_tile(i)), r_um)
      tstar_land = tstar_land + frac_surft(1, i) * tstar_surft(1, i)
      ! sensible heat flux
      ftl_surft(1, i) = real(tile_heat_flux(map_tile(i)), r_um)
      ! moisture flux
      fqw_surft(1, i) = real(tile_moisture_flux(map_tile(i)), r_um)
    end do

    ! Sea temperature
    tstar_sea = real(tile_temperature(map_tile(first_sea_tile)), r_um)

    ! Sea-ice temperatures
    i_sice = 0
    tstar_sice = 0.0_r_um
    if (ice_fract(1, 1) > 0.0_r_um) then
      do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
        i_sice = i_sice + 1
        tstar_sice_ncat(1, 1, i_sice) = real(tile_temperature(map_tile(i)), r_um)
        tstar_sice = tstar_sice &
                   + ice_fract_ncat(1,1,i_sice) * tstar_sice_ncat(1,1,i_sice) &
                   / ice_fract
        ! sea-ice heat flux
        ftl_ice(1,1,i_sice) = real(tile_heat_flux(map_tile(i)), r_um)
        ! sea-ice moisture flux
        fqw_ice(1,1,i_sice) = real(tile_moisture_flux(map_tile(i)), r_um)
      end do
    end if

    ! Sea & Sea-ice temperature
    tstar_ssi = (1.0_r_um - ice_fract) * tstar_sea + ice_fract * tstar_sice

    ! Grid-box mean surface temperature
    tstar = flandg * tstar_land + (1.0_r_um - flandg) * tstar_ssi

    ! Sea-ice conductivity, bulk temperature and thickness
    do i = 1, n_sea_ice_tile
      ti_sice_ncat(1, 1, i) = real(sea_ice_temperature(map_sice(i)), r_um)
    end do

    ! Things passed from other parametrization schemes on this timestep

    ! Update variables with 2nd two DOFs
    do k = 1, nlayers
      ! u wind on rho levels
      u(1,1,k) = u_physics(map_w2(3) + k-1)
      ! v wind on rho levels
      v(1,1,k) = u_physics(map_w2(4) + k-1)
    end do

    rhokm_u_land(1,1) = rhokm_surf_w2(map_w2(3) + 1)
    rhokm_u_ssi(1,1) = rhokm_surf_w2(map_w2(3) + 2)
    flandg_u(1,1) = rhokm_surf_w2(map_w2(3) + 3)
    rhokm_v_land(1,1) = rhokm_surf_w2(map_w2(4) + 1)
    rhokm_v_ssi(1,1) = rhokm_surf_w2(map_w2(4) + 2)
    flandg_v(1,1) = rhokm_surf_w2(map_w2(4) + 3)
    do k = 1, bl_levels
      rhokm_u(1,1,k) = rhokm_w2(map_w2(3) + k)
      rhokm_v(1,1,k) = rhokm_w2(map_w2(4) + k)
      rhokh(1,1,k) = rhokh_bl(map_w3(1) + k)
      fqw(1,1,k) = moist_flux_bl(map_w3(1) + k)
      ftl(1,1,k) = heat_flux_bl(map_w3(1) + k)
      dtrdz_u(1,1,k) = dtrdz_w2(map_w2(3) + k)
      dtrdz_v(1,1,k) = dtrdz_w2(map_w2(4) + k)
    end do
    do k = 2, bl_levels
      rdz_u(1,1,k) = rdz_w2(map_w2(3) + k)
      rdz_v(1,1,k) = rdz_w2(map_w2(4) + k)
      f_ngstress_u(1,1,k) = ngstress_w2(map_w2(3) + k)
      f_ngstress_v(1,1,k) = ngstress_w2(map_w2(4) + k)
    end do

    do i = 1, n_land_tile
      dtstar_surft(1, i) = dtstar_tile(map_tile(i))
    end do

    dtstar_sea(1,1) = dtstar_tile(map_tile(first_sea_tile))

    i_sice = 0
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      dtstar_sice(1,1,i_sice) = dtstar_tile(map_tile(i))
    end do

    ! Re-call the scheme
    call bdy_expl3 (                                                 &
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
    ! increments / fields with increments added
    !-----------------------------------------------------------------------
    do k = 1, nlayers
      t_latest(1,1,k) = theta_star(map_wth(1) + k) * exner_theta_levels(1,1,k) &
                      + dt_conv(map_wth(1) + k)
      r_u(1,1,k) = u_physics_star(map_w2(3) + k-1)                             &
                   - u_physics(map_w2(3) + k-1) + du_conv_w2(map_w2(3) + k-1)
      r_v(1,1,k) = u_physics_star(map_w2(4) + k-1)                             &
                   - u_physics(map_w2(4) + k-1) + du_conv_w2(map_w2(4) + k-1)
      r_w(1,1,k) = u_physics_star(map_w2(5) + k) - u_physics(map_w2(5) + k)
      q_latest(1,1,k)   = m_v(map_wth(1) + k)
      qcl_latest(1,1,k) = m_cl(map_wth(1) + k)
      qcf_latest(1,1,k) = m_ci(map_wth(1) + k)
    end do
    r_w(1,1,0) = u_physics_star(map_w2(5) + 0) - u_physics(map_w2(5) + 0)

    if (scheme == scheme_pc2) then
      do k = 1, nlayers
        ! Assign _latest with current updated values
        cfl_latest(1,1,k) = cf_liq(map_wth(1) + k)
        cff_latest(1,1,k) = cf_ice(map_wth(1) + k)
        cf_latest(1,1,k)  = cf_bulk(map_wth(1) + k)
      end do
    end if

    call NI_imp_ctl (                                                   &
    ! IN Model switches
            outer                                                       &
    ! IN trig arrays
          , xx_cos_theta_latitude                                       &
    ! IN data fields.
          , p_theta_levels, p_rho_minus_one, rho_wet_rsq, rho_wet_tq    &
          , u, v, w                                                     &
          , land_sea_mask, q, qcl, qcf, p_star, theta, qrain            &
          , exner_theta_levels                                          &
    ! IN ancillary fields and fields needed to be kept from tstep to tstep
          , sil_orog_land_gb, ho2r2_orog_gb                             &
          , ice_fract, di_sice_ncat, ice_fract_ncat, k_sice_ncat        &
          , u_0, v_0, land_index, cca_3d, lcbase, ccb0, cct0            &
          , ls_rain, ls_snow, conv_rain, conv_snow                      &
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
          , surft_pts,surft_index,frac_surft,canopy_surft               &
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
    ! SCM Diagnostics (dummy in full UM) & BL diags
          , nSCMDpkgs, L_SCMDiags, bl_diag, sf_diag                     &
    ! INOUT (Note ti_sice_ncat and ti_sice are IN if l_sice_multilayers=T)
          , TScrnDcl_SSI, TScrnDcl_surft, tStbTrans                     &
          , cca0, fqw, ftl, taux, tauy, rhokh                           &
          , fqw_ice,ftl_ice,dtstar_surft,dtstar_sea,dtstar_sice,ti_sice_ncat &
          , area_cloud_fraction, bulk_cloud_fraction                    &
          , t_latest, q_latest, qcl_latest, qcf_latest                  &
          , cf_latest, cfl_latest, cff_latest                           &
          , R_u, R_v, R_w, cloud_fraction_liquid, cloud_fraction_frozen &
          , sum_eng_fluxes,sum_moist_flux, rhcpt                        &
    ! Dummy arguments for bimodal scheme compliance
          , zeros,zeros,zeros,zeros,zeros                               &
    ! INOUT tracer fields
          , aerosol, free_tracers,  resist_b,  resist_b_surft           &
          , dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6 &
          , drydep2, so2, dms, so4_aitken, so4_accu, so4_diss, nh3      &
          , soot_new, soot_aged, soot_cld, bmass_new, bmass_aged        &
          , bmass_cld, ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss &
          , co2, ozone_tracer                                           &
    ! INOUT additional variables for JULES
          , tstar_surft, fqw_surft, epot_surft, ftl_surft               &
          , radnet_sice,olr,tstar_sice_ncat,tstar_ssi                   &
          , tstar_sea,taux_land,taux_ssi,tauy_land,tauy_ssi,Error_code  &
    ! OUT fields
          , surf_ht_flux_land, zlcl_mix                                 &
          , theta_star_surf, qv_star_surf                               &
    ! OUT additional variables for JULES
          , tstar, ti_sice, ext, snowmelt,tstar_land,tstar_sice, ei_surft &
          , ecan_surft, melt_surft, surf_htf_surft                      &
    ! OUT fields for coupling to the wave model
          , uwind_wav, vwind_wav                                        &
            )

    !-----------------------------------------------------------------------
    ! update main model prognostics
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
      ! wind increments
      du_bl_w2(map_w2(3) + k - 1) = r_u(1,1,k)                          &
         - (u_physics_star(map_w2(3) + k-1) - u_physics(map_w2(3) + k-1))
      du_bl_w2(map_w2(4) + k - 1) = r_v(1,1,k)                          &
         - (u_physics_star(map_w2(4) + k-1) - u_physics(map_w2(4) + k-1))
    end do

    ! Update lowest-level values
    select case(lowest_level)
      case(lowest_level_constant)
        dtheta_bl(map_wth(1)) = t_latest(1,1,1) / exner_theta_levels(1,1,1)   &
                              - theta_star(map_wth(1))
        m_v(map_wth(1))  = m_v(map_wth(1) + 1)

      case(lowest_level_gradient)
        dtheta_bl(map_wth(1)) = t_latest(1,1,1) / exner_theta_levels(1,1,1)   &
                              - z_theta(1,1,1) * (                            &
                                t_latest(1,1,2) / exner_theta_levels(1,1,2)   &
                              - t_latest(1,1,1) / exner_theta_levels(1,1,1) ) &
                              / (z_theta(1,1,2) - z_theta(1,1,1))             &
                              - theta_star(map_wth(1))
        m_v(map_wth(1))  = m_v(map_wth(1) + 1)                                &
                         - z_theta(1,1,1) * (                                 &
                           m_v(map_wth(1) + 2) - m_v(map_wth(1) + 1) )        &
                         / (z_theta(1,1,2) - z_theta(1,1,1))
      case(lowest_level_flux)
        dtheta_bl(map_wth(1)) = t_latest(1,1,1) / exner_theta_levels(1,1,1)   &
                              + ftl(1,1,1) / (cp * rhokh(1,1,1))              &
                              - theta_star(map_wth(1))
        m_v(map_wth(1))  = m_v(map_wth(1) + 1)                                &
                         + fqw(1,1,1) / rhokh(1,1,1)
    end select
    m_cl(map_wth(1)) = m_cl(map_wth(1) + 1)
    m_ci(map_wth(1)) = m_ci(map_wth(1) + 1)

    ! update cloud fractions only if using cloud scheme
    if ( cloud == cloud_um ) then
      if ( scheme == scheme_smith ) then
        do k = 1, nlayers
          cf_bulk(map_wth(1) + k) = bulk_cloud_fraction(1,1,k)
          cf_ice(map_wth(1) + k)  = cloud_fraction_frozen(1,1,k)
          cf_liq(map_wth(1) + k)  = cloud_fraction_liquid(1,1,k)
          cf_area(map_wth(1) + k) = area_cloud_fraction(1,1,k)
        end do
      else if ( scheme == scheme_pc2 ) then
        do k = 1, nlayers
          cf_ice (map_wth(1) + k) = cff_latest(1,1,k)
          cf_liq (map_wth(1) + k) = cfl_latest(1,1,k)
          cf_bulk(map_wth(1) + k) = cf_latest(1,1,k)
          cf_area(map_wth(1) + k) = area_cloud_fraction(1,1,k)
        end do
      end if
      cf_ice(map_wth(1) + 0)  = cf_ice(map_wth(1) + 1)
      cf_liq(map_wth(1) + 0)  = cf_liq(map_wth(1) + 1)
      cf_bulk(map_wth(1) + 0) = cf_bulk(map_wth(1) + 1)
      cf_area(map_wth(1) + 0) = cf_area(map_wth(1) + 1)
    endif

    ! update BL prognostics
    if (outer == outer_iterations) then

      z_lcl(map_2d(1)) = zlcl_mix(1,1)

      ! Update land tiles
      do i = 1, n_land_tile
        tile_temperature(map_tile(i)) = real(tstar_surft(1, i), r_def)
        tile_heat_flux(map_tile(i)) = real(ftl_surft(1, i), r_def)
        tile_moisture_flux(map_tile(i)) = real(fqw_surft(1, i), r_def)
        snow_sublimation(map_tile(i)) = real(ei_surft(1, i), r_def)
        ! NB - net surface heat flux
        surf_heat_flux(map_tile(i)) = real(surf_htf_surft(1, i), r_def)
        canopy_evap(map_tile(i)) = real(ecan_surft(1, i), r_def)
        total_snowmelt(map_tile(i)) = real(melt_surft(1, i), r_def)
      end do

      ! Update sea tile
      tile_temperature(map_tile(first_sea_tile)) = real(tstar_sea(1,1), r_def)
      ! NB actually grid box mean but okay for aquaplanet
      tile_heat_flux(map_tile(first_sea_tile)) = real(ftl(1,1,1), r_def)
      ! NB actually grid box mean but okay for aquaplanet
      tile_moisture_flux(map_tile(first_sea_tile)) = real(fqw(1,1,1), r_def)

      ! Update sea-ice tiles
      i_sice = 0
      do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
        i_sice = i_sice + 1
        tile_temperature(map_tile(i)) = real(tstar_sice_ncat(1,1,i_sice), r_def)
        tile_heat_flux(map_tile(i)) = real(ftl_ice(1,1,i_sice), r_def)
        tile_moisture_flux(map_tile(i)) = real(fqw_ice(1,1,i_sice), r_def)
      end do

      ! Sea-ice bulk temperature
      do i = 1, n_sea_ice_tile
        sea_ice_temperature(map_sice(i)) = real(ti_sice_ncat(1, 1, i), r_def)
      end do

      do i = 1, sm_levels
        water_extraction(map_soil(i)) = real(ext(1, i), r_def)
      end do

    endif

    ! deallocate diagnostics deallocated in atmos_physics2
    call dealloc_bl_imp(bl_diag)
    call dealloc_sf_imp(sf_diag)
    call dealloc_sf_expl(sf_diag)

    ! set this back to 1 before exit
    land_field = 1

  end subroutine bl_imp_code

end module bl_imp_kernel_mod
