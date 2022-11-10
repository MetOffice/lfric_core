!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the implicit UM boundary layer scheme.
!>
module bl_imp_kernel_mod

  use argument_mod,           only : arg_type,                  &
                                     GH_FIELD, GH_SCALAR,       &
                                     GH_INTEGER, GH_REAL,       &
                                     GH_READ, GH_WRITE, GH_INC, &
                                     GH_READWRITE, CELL_COLUMN, &
                                     ANY_DISCONTINUOUS_SPACE_1, &
                                     ANY_DISCONTINUOUS_SPACE_2, &
                                     ANY_DISCONTINUOUS_SPACE_3, &
                                     ANY_DISCONTINUOUS_SPACE_4, &
                                     ANY_DISCONTINUOUS_SPACE_5, &
                                     ANY_DISCONTINUOUS_SPACE_6, &
                                     ANY_DISCONTINUOUS_SPACE_7
  use section_choice_config_mod, only : cloud, cloud_um
  use blayer_config_mod,         only : fric_heating
  use cloud_config_mod,          only : scheme, scheme_smith, scheme_pc2, &
                                        scheme_bimodal
  use radiation_config_mod,      only : topography, topography_horizon
  use constants_mod,             only : i_def, i_um, r_def, r_um
  use empty_data_mod,            only : empty_real_data
  use fs_continuity_mod,         only : W3, Wtheta
  use kernel_mod,                only : kernel_type
  use timestepping_config_mod,   only : outer_iterations
  use mixing_config_mod,         only : leonard_term
  use physics_config_mod,        only : lowest_level,          &
                                        lowest_level_constant, &
                                        lowest_level_gradient, &
                                        lowest_level_flux
  use planet_config_mod,         only : cp
  use jules_control_init_mod,    only : n_sea_ice_tile
  use water_constants_mod,       only : tfs, lc, lf
  use derived_config_mod,        only : l_esm_couple
  use surface_config_mod,        only : emis_method_soil, emis_method_soil_fixed
  use jules_surface_types_mod,   only: soil

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: bl_imp_kernel_type
    private
    type(arg_type) :: meta_args(106) = (/                                         &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                &! outer
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! theta_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! wetrho_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! wetrho_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! exner_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! exner_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_v_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_cl_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_ci_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! theta_star
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! height_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! height_wth
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! ntml_2d
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! cumulus_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_fraction
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_3),&! leaf_area_index
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_3),&! canopy_height
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! peak_to_trough_orog
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! silhouette_area_orog
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_thickness
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_temperature
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_4),&! sea_ice_conductivity
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! tile_temperature
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_lw_grey_albedo
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! screen_temperature
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! time_since_transition
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! latitude
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_snow_mass
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! n_snow_layers
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! snow_depth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! canopy_water
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_5),&! soil_temperature
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! tile_heat_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2),&! tile_moisture_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! sw_up_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! sw_down_surf
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! lw_down_surf
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! skyview
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! snowice_sublimation (kg m-2 s-1)
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! surf_heat_flux (W m-2)
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! canopy_evap (kg m-2 s-1)
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_5),&! water_extraction (kg m-2 s-1)
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! snowice_melt (kg m-2 s-1)
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     WTHETA),                   &! dtheta_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! diss_u
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! diss_v
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! dt_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! m_v
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! m_cl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! m_ci
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cf_area
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cf_ice
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cf_liq
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cf_bulk
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! rh_crit_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! dsldzm
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! wvar
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! gradrinr
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! tau_dec_bm
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! tau_hom_bm
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! tau_mph_bm
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! rhokh_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! bq_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! bt_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! moist_flux_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! heat_flux_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! dtrdz_tq_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! rdz_tq_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! thetal_inc_leonard_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! mt_inc_leonard_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! alpha1_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! ashtf_prime_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! dtstar_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! fraca_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! z0h_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! z0m_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! rhokh_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! chr1p5m_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! resfs_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! canhc_tile
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_6),&! tile_water_extract
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! z_lcl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! inv_depth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! qcl_at_inv_top
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! blend_height_tq
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! ustar
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! soil_moist_avail
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! zh_nonloc
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! zh_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! zhsc_2d
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_7),&! bl_type_ind
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! t1p5m_surft
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! q1p5m_surft
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! t1p5m
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! q1p5m
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! qcl1p5m
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! rh1p5m
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! latent_heat
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! snomlt_surf_htf
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! soil_evap
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! surf_ht_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! soil_surf_ht_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! surf_sw_net
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! surf_radnet
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2),&! surf_lw_up
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2) &! surf_lw_down
         /)

    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: bl_imp_code
  end type

  public :: bl_imp_code

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
  !> @param[in]     m_v_n                Vapour mixing ratio at time level n
  !> @param[in]     m_cl_n               Cloud liq mixing ratio at time level n
  !> @param[in]     m_ci_n               Cloud ice mixing ratio at time level n
  !> @param[in]     theta_star           Potential temperature after advection
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
  !> @param[in]     sea_ice_conductivity Sea ice thermal conductivity (W m-2 K-1)
  !> @param[in,out] tile_temperature     Surface tile temperatures
  !> @param[in]     tile_lw_grey_albedo  Surface tile longwave grey albedo
  !> @param[in,out] screen_temperature   Tiled screen level liquid temperature
  !> @param[in,out] time_since_transition Time since decoupled screen transition
  !> @param[in]     latitude             Latitude of cell centre
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
  !> @param[in]     skyview              Skyview / area enhancement factor
  !> @param[in,out] snowice_sublimation  Sublimation of snow and ice
  !> @param[in,out] surf_heat_flux       Surface heat flux
  !> @param[in,out] canopy_evap          Canopy evaporation from land tiles
  !> @param[in,out] water_extraction     Extraction of water from each soil layer
  !> @param[in,out] snowice_melt         Surface, canopy and sea ice, snow and ice melt rate
  !> @param[in,out] dtheta_bl            BL theta increment
  !> @param[in]     diss_u               Zonal Molecular dissipation rate
  !> @param[in]     diss_v               Meridional Molecular dissipation rate
  !> @param[in]     dt_conv              Convection temperature increment
  !> @param[in,out] m_v                  Vapour mixing ration after advection
  !> @param[in,out] m_cl                 Cloud liq mixing ratio after advection
  !> @param[in,out] m_ci                 Cloud ice mixing ratio after advection
  !> @param[in,out] cf_area              Area cloud fraction
  !> @param[in,out] cf_ice               Ice cloud fraction
  !> @param[in,out] cf_liq               Liquid cloud fraction
  !> @param[in,out] cf_bulk              Bulk cloud fraction
  !> @param[in]     rh_crit_wth          Critical relative humidity
  !> @param[in]     dsldzm               Liquid potential temperature gradient in wth
  !> @param[in]     wvar                 Vertical velocity variance in wth
  !> @param[in]     gradrinr             Gradient Richardson number in wth
  !> @param[in]     tau_dec_bm           Decorrelation time scale in wth
  !> @param[in]     tau_hom_bm           Homogenisation time scale in wth
  !> @param[in]     tau_mph_bm           Phase-relaxation time scale in wth
  !> @param[in]     rhokh_bl             Heat eddy diffusivity on BL levels
  !> @param[in]     bq_bl                Buoyancy parameter for moisture
  !> @param[in]     bt_bl                Buoyancy parameter for heat
  !> @param[in,out] moist_flux_bl        Vertical moisture flux on BL levels
  !> @param[in,out] heat_flux_bl         Vertical heat flux on BL levels
  !> @param[in]     dtrdz_tq_bl          dt/(rho*r*r*dz) in wth
  !> @param[in]     rdz_tq_bl            1/dz in w3
  !> @param[in]     thetal_inc_leonard_wth  Leonard term increment for water potential temperature
  !> @param[in]     mt_inc_leonard_wth   Leonard term increment for total moisture
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
  !> @param[in]     ustar                Friction velocity
  !> @param[in]     soil_moist_avail     Available soil moisture for evaporation
  !> @param[in]     zh_nonloc            Depth of non-local BL scheme
  !> @param[in]     zh_2d                Total BL depth
  !> @param[in]     zhsc_2d              Height of decoupled layer top
  !> @param[in]     bl_type_ind          Diagnosed BL types
  !> @param[in,out] t1p5m_surft          Diagnostic: 1.5m temperature for land tiles
  !> @param[in,out] q1p5m_surft          Diagnostic: 1.5m specific humidity for land tiles
  !> @param[in,out] t1p5m                Diagnostic: 1.5m temperature
  !> @param[in,out] q1p5m                Diagnostic: 1.5m specific humidity
  !> @param[in,out] qcl1p5m              Diagnostic: 1.5m specific cloud water
  !> @param[in,out] rh1p5m               Diagnostic: 1.5m relative humidity
  !> @param[in,out] latent_heat          Diagnostic: Surface latent heat flux
  !> @param[in,out] snomlt_surf_htf      Diagnostic: Grid mean suface snowmelt heat flux
  !> @param[in,out] soil_evap            Diagnostic: Grid mean evapotranspiration from the soil
  !> @param[in,out] surf_ht_flux         Diagnostic: Surface to sub-surface heat flux
  !> @param[in,out] soil_surf_ht_flux    Diagnostic: Grid mean surface soil heat flux
  !> @param[in,out] surf_sw_net          Diagnostic: Net surface shortwave radiation
  !> @param[in,out] surf_radnet          Diagnostic: Net surface radiation
  !> @param[in,out] surf_lw_up           Diagnostic: Upward surface longtwave radiation
  !> @param[in,out] surf_lw_down         Diagnostic: Downward surface longwave radiation
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
                         m_v_n,                              &
                         m_cl_n,                             &
                         m_ci_n,                             &
                         theta_star,                         &
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
                         sea_ice_conductivity,               &
                         tile_temperature,                   &
                         tile_lw_grey_albedo,                &
                         screen_temperature,                 &
                         time_since_transition,              &
                         latitude,                           &
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
                         skyview,                            &
                         snowice_sublimation,                &
                         surf_heat_flux,                     &
                         canopy_evap,                        &
                         water_extraction,                   &
                         snowice_melt,                       &
                         dtheta_bl,                          &
                         diss_u,                             &
                         diss_v,                             &
                         dt_conv,                            &
                         m_v,                                &
                         m_cl,                               &
                         m_ci,                               &
                         cf_area,                            &
                         cf_ice,                             &
                         cf_liq,                             &
                         cf_bulk,                            &
                         rh_crit_wth,                        &
                         dsldzm,                             &
                         wvar,                               &
                         gradrinr,                           &
                         tau_dec_bm,                         &
                         tau_hom_bm,                         &
                         tau_mph_bm,                         &
                         rhokh_bl,                           &
                         bq_bl,                              &
                         bt_bl,                              &
                         moist_flux_bl,                      &
                         heat_flux_bl,                       &
                         dtrdz_tq_bl,                        &
                         rdz_tq_bl,                          &
                         thetal_inc_leonard_wth,             &
                         mt_inc_leonard_wth,                 &
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
                         ustar,                              &
                         soil_moist_avail,                   &
                         zh_nonloc,                          &
                         zh_2d,                              &
                         zhsc_2d,                            &
                         bl_type_ind,                        &
                         t1p5m_surft, q1p5m_surft,           &
                         t1p5m, q1p5m, qcl1p5m, rh1p5m,      &
                         latent_heat, snomlt_surf_htf,       &
                         soil_evap,                          &
                         surf_ht_flux, soil_surf_ht_flux,    &
                         surf_sw_net, surf_radnet,           &
                         surf_lw_up, surf_lw_down,           &
                         ndf_wth,                            &
                         undf_wth,                           &
                         map_wth,                            &
                         ndf_w3,                             &
                         undf_w3,                            &
                         map_w3,                             &
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
    use ancil_info, only: ssi_pts, sea_pts, sice_pts, sice_pts_ncat,           &
                          nsurft, nsoilt, dim_cslayer, rad_nband, nmasst
    use atm_fields_bounds_mod, only: pdims, pdims_s
    use atm_step_local, only: dim_cs1
    use csigma, only: sbcon
    use dust_parameters_mod, only: ndiv, ndivh
    use jules_deposition_mod, only: l_deposition
    use jules_irrig_mod, only: irr_crop, irr_crop_doell
    use jules_sea_seaice_mod, only: nice, nice_use, emis_sea
    use jules_snow_mod, only: cansnowtile, rho_snow_const, l_snowdep_surf,nsmax
    use jules_surface_types_mod, only: npft, ntype, lake, nnvg, ncpft, nnpft
    use jules_surface_mod, only: l_flake_model
    use jules_vegetation_mod, only: can_model, l_crop, l_triffid, l_phenol,    &
                                    can_rad_mod, l_acclim
    use jules_radiation_mod, only: l_albedo_obs
    use jules_soil_mod, only: ns_deep, l_bedrock
    use jules_soil_biogeochem_mod, only: dim_ch4layer, soil_bgc_model,         &
                                         soil_model_ecosse, l_layeredc
    use nlsizes_namelist_mod, only: row_length, rows, land_field,              &
                                    sm_levels, ntiles, bl_levels, tr_vars
    use pftparm, only: emis_pft
    use planet_constants_mod, only: p_zero, kappa, planet_radius, two_omega
    use nvegparm, only: emis_nvg
    use rad_input_mod, only: co2_mmr
    use timestep_mod, only: timestep
    use theta_field_sizes, only: t_i_length, t_j_length, &
                                 u_i_length,u_j_length,  &
                                 v_i_length,v_j_length
    use veg3_parm_mod, only: l_veg3

    ! spatially varying fields used from modules
    use level_heights_mod, only: r_theta_levels, r_rho_levels
    use dyn_coriolis_mod, only: f3_at_u
    use leonard_incs_mod, only: thetal_inc_leonard, qw_inc_leonard
    use solinc_data, only: sky

    ! subroutines used
    use bl_diags_mod, only: bl_diag, dealloc_bl_imp, alloc_bl_expl
    use sf_diags_mod, only: sf_diag, dealloc_sf_expl, dealloc_sf_imp,       &
                            alloc_sf_expl
    use ni_imp_ctl_mod, only: ni_imp_ctl
    use tilepts_mod, only: tilepts
    use ls_cld_mod, only: ls_cld
    use qsat_mod, only: qsat

    !---------------------------------------
    ! JULES modules
    !---------------------------------------
    use crop_vars_mod,            only: crop_vars_type, crop_vars_data_type,   &
                                        crop_vars_alloc, crop_vars_assoc, &
                                        crop_vars_nullify, crop_vars_dealloc
    use prognostics,              only: progs_data_type, progs_type,           &
                                        prognostics_alloc, prognostics_assoc,  &
                                        prognostics_nullify, prognostics_dealloc
    use jules_vars_mod,           only: jules_vars_type, jules_vars_data_type, &
                                        jules_vars_alloc, jules_vars_assoc,    &
                                        jules_vars_dealloc, jules_vars_nullify
    use aero,                     only: aero_type, aero_data_type,             &
                                        aero_alloc, aero_assoc, &
                                        aero_nullify, aero_dealloc
    use coastal,                  only: coastal_type, coastal_data_type,       &
                                        coastal_assoc, coastal_alloc, &
                                        coastal_dealloc, coastal_nullify
    use lake_mod,                 only: lake_type, lake_data_type,             &
                                        lake_assoc, lake_alloc, &
                                        lake_dealloc, lake_nullify
    use jules_forcing_mod,        only: forcing_type, forcing_data_type,       &
                                        forcing_assoc, forcing_alloc,          &
                                        forcing_nullify, forcing_dealloc
    use ancil_info,               only: ainfo_type, ainfo_data_type,           &
                                        ancil_info_assoc, ancil_info_alloc,    &
                                        ancil_info_dealloc, ancil_info_nullify
    use fluxes_mod,               only: fluxes_type, fluxes_data_type,         &
                                        fluxes_alloc, fluxes_assoc,            &
                                        fluxes_nullify, fluxes_dealloc
    use trifctl,                  only: trifctl_type

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: outer

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
    integer(kind=i_def), intent(in) :: ndf_sice, undf_sice
    integer(kind=i_def), intent(in) :: map_sice(ndf_sice)
    integer(kind=i_def), intent(in) :: ndf_soil, undf_soil
    integer(kind=i_def), intent(in) :: map_soil(ndf_soil)
    integer(kind=i_def), intent(in) :: ndf_smtile, undf_smtile
    integer(kind=i_def), intent(in) :: map_smtile(ndf_smtile)

    integer(kind=i_def), intent(in) :: ndf_bl, undf_bl
    integer(kind=i_def), intent(in) :: map_bl(ndf_bl)

    real(kind=r_def), dimension(undf_w3),  intent(inout) :: moist_flux_bl,     &
                                                            heat_flux_bl
    real(kind=r_def), dimension(undf_wth), intent(inout) :: dtheta_bl,         &
                                                            m_v, m_cl, m_ci,   &
                                                            cf_area, cf_ice,   &
                                                            cf_liq, cf_bulk
    real(kind=r_def), dimension(undf_w3),  intent(in)   :: wetrho_in_w3,       &
                                                           exner_in_w3,        &
                                                           height_w3,          &
                                                           rhokh_bl,           &
                                                           rdz_tq_bl,          &
                                                           diss_u, diss_v
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
                                                           dtrdz_tq_bl,        &
                                                           dsldzm,             &
                                                           wvar,               &
                                                           gradrinr,           &
                                                           tau_dec_bm,         &
                                                           tau_hom_bm,         &
                                                           tau_mph_bm,         &
                                                           mt_inc_leonard_wth, &
                                                           thetal_inc_leonard_wth

    integer(kind=i_def), dimension(undf_2d), intent(in) :: ntml_2d,            &
                                                           cumulus_2d,         &
                                                           blend_height_tq

    real(kind=r_def), dimension(undf_2d), intent(in) :: zh_2d,                &
                                                        ustar,                &
                                                        soil_moist_avail,     &
                                                        zh_nonloc, zhsc_2d

    real(kind=r_def), intent(in)    :: tile_fraction(undf_tile)
    real(kind=r_def), intent(inout) :: tile_temperature(undf_tile)
    real(kind=r_def), intent(in)    :: tile_lw_grey_albedo(undf_tile)
    real(kind=r_def), intent(inout) :: screen_temperature(undf_tile)
    real(kind=r_def), intent(inout) :: time_since_transition(undf_2d)
    real(kind=r_def), intent(in)    :: latitude(undf_2d)
    real(kind=r_def), intent(in)    :: tile_snow_mass(undf_tile)
    integer(kind=i_def), intent(in) :: n_snow_layers(undf_tile)
    real(kind=r_def), intent(in)    :: snow_depth(undf_tile)
    real(kind=r_def), intent(in)    :: canopy_water(undf_tile)
    real(kind=r_def), intent(inout) :: tile_heat_flux(undf_tile)
    real(kind=r_def), intent(inout) :: tile_moisture_flux(undf_tile)
    real(kind=r_def), intent(in)    :: sw_up_tile(undf_tile)
    real(kind=r_def), intent(inout)   :: snowice_sublimation(undf_tile)
    real(kind=r_def), intent(inout)   :: surf_heat_flux(undf_tile)
    real(kind=r_def), intent(inout)   :: canopy_evap(undf_tile)
    real(kind=r_def), intent(inout)   :: snowice_melt(undf_tile)

    real(kind=r_def), intent(in) :: leaf_area_index(undf_pft)
    real(kind=r_def), intent(in) :: canopy_height(undf_pft)

    real(kind=r_def), intent(in) :: sea_ice_thickness(undf_sice)
    real(kind=r_def), intent(inout) :: sea_ice_temperature(undf_sice)
    real(kind=r_def), intent(inout) :: sea_ice_conductivity(undf_sice)

    real(kind=r_def), intent(in) :: peak_to_trough_orog(undf_2d)
    real(kind=r_def), intent(in) :: silhouette_area_orog(undf_2d)
    real(kind=r_def), intent(in) :: sw_down_surf(undf_2d)
    real(kind=r_def), intent(in) :: lw_down_surf(undf_2d)
    real(kind=r_def), intent(in) :: skyview(undf_2d)
    real(kind=r_def), intent(inout) :: z_lcl(undf_2d)
    real(kind=r_def), intent(in) :: inv_depth(undf_2d)
    real(kind=r_def), intent(in) :: qcl_at_inv_top(undf_2d)
    real(kind=r_def), pointer, intent(inout) :: t1p5m_surft(:)
    real(kind=r_def), pointer, intent(inout) :: q1p5m_surft(:)
    real(kind=r_def), pointer, intent(inout) :: t1p5m(:), q1p5m(:)
    real(kind=r_def), pointer, intent(inout) :: qcl1p5m(:), rh1p5m(:)
    real(kind=r_def), pointer, intent(inout) :: latent_heat(:)
    real(kind=r_def), pointer, intent(inout) :: snomlt_surf_htf(:)
    real(kind=r_def), pointer, intent(inout) :: soil_evap(:)
    real(kind=r_def), pointer, intent(inout) :: surf_ht_flux(:)
    real(kind=r_def), pointer, intent(inout) :: soil_surf_ht_flux(:)
    real(kind=r_def), pointer, intent(inout) :: surf_sw_net(:)
    real(kind=r_def), pointer, intent(inout) :: surf_radnet(:)
    real(kind=r_def), pointer, intent(inout) :: surf_lw_up(:)
    real(kind=r_def), pointer, intent(inout) :: surf_lw_down(:)

    real(kind=r_def), intent(in) :: soil_temperature(undf_soil)
    real(kind=r_def), intent(inout):: water_extraction(undf_soil)

    real(kind=r_def), intent(in) :: tile_water_extract(undf_smtile)

    integer(kind=i_def), dimension(undf_bl), intent(in) :: bl_type_ind
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

    real(kind=r_um) :: target_mass(pdims%i_start:pdims%i_end,                &
                                   pdims%j_start:pdims%j_end,1)

    !-----------------------------------------------------------------------
    ! Local variables for the kernel
    !-----------------------------------------------------------------------
    ! loop counters etc
    integer(i_def) :: k, i, i_tile, i_sice, n, icode

    ! local switches and scalars
    integer(i_um) :: error_code
    logical, parameter :: l_calc_at_p=.false.

    ! profile fields from level 1 upwards
    real(r_um), dimension(row_length,rows,nlayers) ::                        &
         p_rho_levels, rho_wet_rsq, rho_wet, z_rho, z_theta,                 &
         bulk_cloud_fraction, rhcpt, t_latest, q_latest, qcl_latest,         &
         qcf_latest, qcf2_latest, cca_3d, area_cloud_fraction,               &
         cloud_fraction_liquid, cloud_fraction_frozen, rho_wet_tq,           &
         cf_latest, cfl_latest, cff_latest

    ! profile field on boundary layer levels
    real(r_um), dimension(row_length,rows,bl_levels) :: fqw, ftl, rhokh,     &
         bq_gb, bt_gb, dtrdz_charney_grid, rdz_charney_grid, rhokm

    ! profile fields on u/v points and all levels
    real(r_um), dimension(row_length,rows,nlayers) :: u, v, r_u, r_v

    ! profile fields on u/v points and BL levels
    real(r_um), dimension(row_length,rows,bl_levels) :: taux, tauy,          &
         dtrdz_u, dtrdz_v, rhokm_u, rhokm_v, dissip_u, dissip_v

    ! profile fields from level 2 upwards
    real(r_um), dimension(row_length,rows,2:bl_levels) :: rdz_u, rdz_v

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
         u_s, z0hssi, z0mssi, zhnl, zlcl_mix, zlcl, dzh, qcl_inv_top,        &
         work_2d_1, work_2d_2, work_2d_3, qcl1p5m_loc

    ! single level real fields on u/v points
    real(r_um), dimension(row_length,rows) :: u_0, v_0, taux_land, tauy_land,&
         taux_ssi, tauy_ssi, flandg_u, flandg_v

    ! single level integer fields
    integer(i_um), dimension(row_length,rows) :: ntml, lcbase, k_blend_tq,   &
         k_blend_uv

    ! single level logical fields
    logical, dimension(row_length,rows) :: land_sea_mask, cumulus

    ! fields on sea-ice categories
    real(r_um), dimension(row_length,rows,nice_use) ::                       &
         tstar_sice_ncat, ice_fract_ncat, alpha1_sice,                       &
         ashtf_prime, rhokh_sice, k_sice_ncat, ti_sice_ncat, di_sice_ncat,   &
         dtstar_sice

    ! field on land points and soil levels
    real(r_um), dimension(land_field,sm_levels) :: t_soil_soilt

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
         tstar_surft, frac_surft,                                            &
         dtstar_surft, alpha1, ashtf_prime_surft, chr1p5m, fraca, resfs,     &
         rhokh_surft, resft, flake,                                          &
         canhc_surft

    ! fields on land points and plant functional types
    real(r_um), dimension(land_field,npft) :: canht_pft, lai_pft

    ! field on surface tiles and soil levels
    real(r_um), dimension(land_field,sm_levels,ntiles) :: wt_ext_surft

    ! Fields which are not used and only required for subroutine argument list,
    ! hence are unset in the kernel
    ! if they become set, please move up to be with other variables
    integer(i_um), parameter :: nscmdpkgs=15
    logical,       parameter :: l_scmdiags(nscmdpkgs)=.false.

    real(r_um), dimension(row_length,rows,nlayers) :: cca0, ccw0

    real(r_um), dimension(row_length,rows,0:nlayers) ::                      &
         aerosol, dust_div1, dust_div2, dust_div3, dust_div4, dust_div5,     &
         dust_div6, so2, dms, so4_aitken, so4_accu, so4_diss, nh3, soot_new, &
         soot_aged, soot_cld, bmass_new, bmass_aged, bmass_cld, ocff_new,    &
         ocff_aged, ocff_cld, nitr_acc, nitr_diss, ozone_tracer, qrain

    real(r_um), dimension(row_length,rows,nlayers) ::                        &
         tgrad_in, tau_dec_in, tau_hom_in, tau_mph_in

    real(r_um), dimension(row_length,rows,nlayers) :: wvar_in

    real(r_um), dimension(row_length,rows,0:nlayers,tr_vars) :: free_tracers

    real(r_um), dimension(row_length,rows) :: ti_sice,                       &
         xx_cos_theta_latitude, ls_rain, ls_snow, conv_rain, conv_snow,      &
         co2_emits, co2flux, co2_flux_tot, tscrndcl_ssi, tstbtrans,          &
         sum_eng_fluxes, sum_moist_flux, drydep2, olr, surf_ht_flux_land,    &
         theta_star_surf, qv_star_surf, uwind_wav,                           &
         vwind_wav, sstfrz, aresist, resist_b, rho_aresist, rib_gb,          &
         z0m_eff_gb, zhsc, cdr10m_u, cdr10m_v

    real(r_um), dimension(row_length,rows,3) :: t_frac, t_frac_dsc, we_lim,  &
         we_lim_dsc, zrzi, zrzi_dsc

    integer(i_um), dimension(row_length,rows) :: lcbase0, ccb0, cct0,        &
         kent, kent_dsc

    real(r_um), dimension(row_length,rows,nice_use) :: radnet_sice

    real(r_um), dimension(row_length,rows,ndiv) :: dust_flux, r_b_dust

    real(r_um), dimension(land_field) :: resp_s_tot_soilt

    real(r_um), dimension(land_field) :: npp_gb

    real(r_um), dimension(land_field,ntiles) :: catch_surft,                 &
         tscrndcl_surft, epot_surft, aresist_surft, dust_emiss_frac,         &
         gc_surft, resist_b_surft, u_s_std_surft

    real(r_um), dimension(land_field,ntiles,ndivh) :: u_s_t_dry_tile,        &
         u_s_t_tile

    real(r_um) :: stashwork3(1), stashwork9(1)


    !-----------------------------------------------------------------------
    ! JULES Types
    !-----------------------------------------------------------------------
    type(crop_vars_type) :: crop_vars
    type(crop_vars_data_type) :: crop_vars_data
    type(progs_type) :: progs
    type(progs_data_type) :: progs_data
    type(jules_vars_type) :: jules_vars
    type(jules_vars_data_type) :: jules_vars_data
    type(aero_type) :: aerotype
    type(aero_data_type) :: aero_data
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
    type(trifctl_type) :: trifctltype

    !-----------------------------------------------------------------------
    ! Initialisation of JULES data and pointer types
    !-----------------------------------------------------------------------
    call crop_vars_alloc(land_field, t_i_length, t_j_length,                  &
                     nsurft, ncpft,nsoilt, sm_levels, l_crop, irr_crop,       &
                     irr_crop_doell, crop_vars_data)

    call crop_vars_assoc(crop_vars, crop_vars_data)

    call prognostics_alloc(land_field, t_i_length, t_j_length,                &
                      nsurft, npft, nsoilt, sm_levels, ns_deep, nsmax,        &
                      dim_cslayer, dim_cs1, dim_ch4layer,                     &
                      nice, nice_use, soil_bgc_model, soil_model_ecosse,      &
                      l_layeredc, l_triffid, l_phenol, l_bedrock, l_veg3,     &
                      nmasst, nnpft, l_acclim, progs_data)
    call prognostics_assoc(progs,progs_data)

    call jules_vars_alloc(land_field,ntype,nsurft,rad_nband,nsoilt,sm_levels, &
                t_i_length, t_j_length, npft, bl_levels, pdims_s, pdims,      &
                l_albedo_obs, cansnowtile, l_deposition,                      &
                jules_vars_data)
    call jules_vars_assoc(jules_vars,jules_vars_data)

    if (can_rad_mod == 6) then
      jules_vars%diff_frac = 0.4_r_um
    else
      jules_vars%diff_frac = 0.0_r_um
    end if

    call aero_alloc(land_field,t_i_length,t_j_length,                         &
                nsurft,ndiv, aero_data)
    call aero_assoc(aerotype, aero_data)

    call coastal_alloc(land_field,t_i_length,t_j_length,                      &
                   u_i_length,u_j_length,                                     &
                   v_i_length,v_j_length,                                     &
                   nice_use,nice,coastal_data)
    call coastal_assoc(coast, coastal_data)

    call lake_alloc(land_field, l_flake_model, lake_data)
    call lake_assoc(lake_vars, lake_data)

    call ancil_info_alloc(land_field,t_i_length,t_j_length,                   &
                      nice,nsoilt,ntype,                                      &
                      ainfo_data)
    call ancil_info_assoc(ainfo, ainfo_data)

    call forcing_alloc(t_i_length,t_j_length,u_i_length, u_j_length,          &
                       v_i_length, v_j_length, forcing_data)
    call forcing_assoc(forcing, forcing_data)

    call fluxes_alloc(land_field, t_i_length, t_j_length,                     &
                      nsurft, npft, nsoilt, sm_levels,                        &
                      nice, nice_use,                                         &
                      fluxes_data)
    call fluxes_assoc(fluxes, fluxes_data)
    !-----------------------------------------------------------------------
    ! Initialisation of variables and arrays
    !-----------------------------------------------------------------------
    error_code=0

    ! neutral wind diagnostics are calculated in bl_imp_du_kernel so set
    ! flags to false here (to avoid divide by zero at end of imp_solver)
    sf_diag%suv10m_n = .false.
    sf_diag%l_u10m_n  = sf_diag%suv10m_n
    sf_diag%l_v10m_n  = sf_diag%suv10m_n
    sf_diag%l_mu10m_n = sf_diag%suv10m_n
    sf_diag%l_mv10m_n = sf_diag%suv10m_n
    call alloc_sf_expl(sf_diag, outer == outer_iterations)
    call alloc_bl_expl(bl_diag, outer == outer_iterations)

! Set logical flags for sf_diags
    sf_diag%smlt = .not. associated(snomlt_surf_htf, empty_real_data)
    sf_diag%l_lw_surft = .not. associated(surf_lw_up, empty_real_data)  &
                    .or. .not. associated(surf_lw_down, empty_real_data)
    sf_diag%l_lw_up_sice_weighted_cat = .not. associated(surf_lw_up, empty_real_data)


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

    if (bl_diag%l_gradrich) then
      allocate(bl_diag%gradrich(pdims%i_start:pdims%i_end,                   &
                           pdims%j_start:pdims%j_end,bl_levels))
      do k = 1, bl_levels-1
        bl_diag%gradrich(1,1,k+1) = gradrinr(map_wth(1) + k)
      end do
    else
      allocate(bl_diag%gradrich(1,1,1))
    end if

    ! Following variables need to be initialised to stop crashed in unused
    ! UM code
    kent = 2
    kent_dsc = 2
    olr = 300.0_r_um
    fluxes%fqw_surft = 0.0_r_um
    cdr10m_u       = 0.0_r_um
    cdr10m_v       = 0.0_r_um
    conv_rain      = 0.0_r_um
    conv_snow      = 0.0_r_um
    cca0           = 0.0_r_um
    epot_surft     = 0.0_r_um

    !-----------------------------------------------------------------------
    ! Mapping of LFRic fields into UM variables
    !-----------------------------------------------------------------------

    ! Fields for decoupled screen temperature diagnostic
    tstbtrans      = time_since_transition(map_2d(1))
    f3_at_u        = two_omega * sin(latitude(map_2d(1)))

    ! Land tile fractions
    flandg = 0.0_r_um
    do i = 1, n_land_tile
      flandg = flandg + real(tile_fraction(map_tile(1)+i-1), r_um)
      frac_surft(1, i) = real(tile_fraction(map_tile(1)+i-1), r_um)
    end do

    ! Sea-ice fraction
    i_sice = 0
    ice_fract = 0.0_r_um
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      ice_fract = ice_fract + real(tile_fraction(map_tile(1)+i-1), r_um)
      ice_fract_ncat(1, 1, i_sice) = real(tile_fraction(map_tile(1)+i-1), r_um)
    end do

    ! Because Jules tests on flandg < 1, we need to ensure this is exactly
    ! 1 when no sea or sea-ice is present
    if ( tile_fraction(map_tile(1)+first_sea_tile-1) == 0.0_r_def .and. &
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
      tscrndcl_surft(1, i) = real(screen_temperature(map_tile(1)+i-1), r_um)
      tstar_land = tstar_land + frac_surft(1, i) * tstar_surft(1, i)
      ! sensible heat flux
      fluxes%ftl_surft(1, i) = real(tile_heat_flux(map_tile(1)+i-1), r_um)
      ! moisture flux
      fluxes%fqw_surft(1, i) = real(tile_moisture_flux(map_tile(1)+i-1), r_um)
    end do

    ! Sea temperature
    ! Default to temperature over frozen sea as the initialisation
    ! that follows does not initialise sea points if they are fully
    ! frozen
    tstar_sea = tfs
    if (tile_fraction(map_tile(1)+first_sea_tile-1) > 0.0_r_def) then
      tstar_sea = real(tile_temperature(map_tile(1)+first_sea_tile-1), r_um)
    end if
    tscrndcl_ssi=real(screen_temperature(map_tile(1)+first_sea_tile-1), r_um)

    ! Sea-ice temperatures
    fluxes%ftl_sicat = 0.0_r_um
    fluxes%fqw_sicat = 0.0_r_um
    tstar_sice       = 0.0_r_um
    tstar_sice_ncat  = 0.0_r_um

    i_sice = 0
    if (ice_fract(1, 1) > 0.0_r_um) then
      do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
        i_sice = i_sice + 1
        tstar_sice_ncat(1, 1, i_sice) = real(tile_temperature(map_tile(1)+i-1), r_um)
        tstar_sice = tstar_sice &
                   + ice_fract_ncat(1,1,i_sice) * tstar_sice_ncat(1,1,i_sice) &
                   / ice_fract
        ! sea-ice heat flux
        fluxes%ftl_sicat(1,1,i_sice) = real(tile_heat_flux(map_tile(1)+i-1), r_um)
        ! sea-ice moisture flux
        fluxes%fqw_sicat(1,1,i_sice) = real(tile_moisture_flux(map_tile(1)+i-1), r_um)
      end do
    end if

    ! Sea & Sea-ice temperature
    tstar_ssi = (1.0_r_um - ice_fract) * tstar_sea + ice_fract * tstar_sice

    ! Grid-box mean surface temperature
    tstar = flandg * tstar_land + (1.0_r_um - flandg) * tstar_ssi

    ! Sea-ice conductivity, bulk temperature and thickness
    do i = 1, n_sea_ice_tile
      k_sice_ncat(1, 1, i) = real(sea_ice_conductivity(map_sice(1)+i-1), r_um)
      ti_sice_ncat(1, 1, i) = real(sea_ice_temperature(map_sice(1)+i-1), r_um)
      di_sice_ncat(1, 1, i) = real(sea_ice_thickness(map_sice(1)+i-1), r_um)
    end do

    ! Ocean coupling point
    ! Temporarily set to true for all points in coupled models.
    ! This will need to be fed from ocn_cpl_point lfric variable
    ! as part of LFRIC#3250
    if ( l_esm_couple ) then
      ainfo%ocn_cpl_point(1, 1) = .true.
    else
      ainfo%ocn_cpl_point(1, 1) = .false.
    end if

    do n = 1, npft
      ! Leaf area index
      lai_pft(1, n) = real(leaf_area_index(map_pft(1)+n-1), r_um)
      ! Canopy height
      canht_pft(1, n) = real(canopy_height(map_pft(1)+n-1), r_um)
    end do

    do i = 1, sm_levels
      ! Soil temperature (t_soil_soilt)
      t_soil_soilt(1, i) = real(soil_temperature(map_soil(1)+i-1), r_um)
    end do

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

    ! Net SW on sea-ice
    i_sice = 0
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      fluxes%sw_sicat(1, i_sice) = real(sw_down_surf(map_2d(1)) - &
                                 sw_up_tile(map_tile(1)+i-1), r_um)
    end do

    ! Carbon dioxide
    co2 = co2_mmr

    ! Half of peak-to-trough height over root(2) of orography (ho2r2_orog_gb)
    ho2r2_orog_gb = real(peak_to_trough_orog(map_2d(1)), r_um)

    sil_orog_land_gb = real(silhouette_area_orog(map_2d(1)), r_um)

    ! Canopy water on each tile (canopy_surft)
    do i = 1, n_land_tile
      canopy_surft(1, i) = real(canopy_water(map_tile(1)+i-1), r_um)
    end do

    do i = 1, n_land_tile
      ! Lying snow mass on land tiles
      snow_surft(1, i) = real(tile_snow_mass(map_tile(1)+i-1), r_um)
      ! Number of snow layers on tiles (nsnow_surft)
      progs%nsnow_surft(1, i) = n_snow_layers(map_tile(1)+i-1)
      ! Equivalent snowdepth for surface calculations.
      ! code copied from jules_land_sf_explicit
      ! 4 is a magic number inherited from Jules, meaning radiative canopy
      ! with heat capacity and snow beneath
      if ( (can_model == 4) .and. cansnowtile(i) .and. l_snowdep_surf) then
        jules_vars%snowdep_surft(1, i) = snow_surft(1, i) / rho_snow_const
      else
        jules_vars%snowdep_surft(1, i) = real(snow_depth(map_tile(1)+i-1), r_um)
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

    if ( leonard_term ) then
      do k = 1, nlayers
        thetal_inc_leonard(1,1,k) = thetal_inc_leonard_wth(map_wth(1) + k)
        qw_inc_leonard(1,1,k) = mt_inc_leonard_wth(map_wth(1) + k)
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
    ! near surface moisture fields
    q(1,1,0) = m_v_n(map_wth(1) + 0)
    qcl(1,1,0) = m_cl_n(map_wth(1) + 0)
    qcf(1,1,0) = m_ci_n(map_wth(1) + 0)
    ! surface height
    r_theta_levels(1,1,0) = height_wth(map_wth(1) + 0) + planet_radius
    ! height of levels above surface
    z_rho = r_rho_levels-r_theta_levels(1,1,0)
    z_theta(1,1,:) = r_theta_levels(1,1,1:nlayers)-r_theta_levels(1,1,0)

    !-----------------------------------------------------------------------
    ! Things passed from other parametrization schemes on this timestep
    !-----------------------------------------------------------------------
    cumulus(1,1) = (cumulus_2d(map_2d(1)) == 1_i_def)
    ntml(1,1) = ntml_2d(map_2d(1))

    do k = 1, bl_levels
      rhokh(1,1,k) = rhokh_bl(map_w3(1) + k-1)
      bq_gb(1,1,k) = bq_bl(map_wth(1) + k-1)
      bt_gb(1,1,k) = bt_bl(map_wth(1) + k-1)
      fqw(1,1,k) = moist_flux_bl(map_w3(1) + k-1)
      ftl(1,1,k) = heat_flux_bl(map_w3(1) + k-1)
      dtrdz_charney_grid(1,1,k) = dtrdz_tq_bl(map_wth(1) + k)
      rdz_charney_grid(1,1,k) = rdz_tq_bl(map_w3(1) + k-1)
    end do
    if (fric_heating) then
      do k = 1, bl_levels
        dissip_u(1,1,k) = diss_u(map_w3(1) + k-1)
        dissip_v(1,1,k) = diss_v(map_w3(1) + k-1)
      end do
    end if

    do i = 1, n_land_tile
      alpha1(1, i) = alpha1_tile(map_tile(1)+i-1)
      ashtf_prime_surft(1, i) = ashtf_prime_tile(map_tile(1)+i-1)
      dtstar_surft(1, i) = dtstar_tile(map_tile(1)+i-1)
      fraca(1, i) = fraca_tile(map_tile(1)+i-1)
      fluxes%z0h_surft(1, i) = z0h_tile(map_tile(1)+i-1)
      fluxes%z0m_surft(1, i) = z0m_tile(map_tile(1)+i-1)
      rhokh_surft(1, i) = rhokh_tile(map_tile(1)+i-1)
      chr1p5m(1, i) = chr1p5m_tile(map_tile(1)+i-1)
      resfs(1, i) = resfs_tile(map_tile(1)+i-1)
      canhc_surft(1, i) = canhc_tile(map_tile(1)+i-1)
    end do
    if (flandg(1,1) > 0.0_r_um) then
      ! recalculate the total resistance factor
      resft = fraca + (1.0 - fraca) * resfs
      resft(1,lake) = 1.0
      ! recalculate the total lake fraction
      flake = 0.0
      flake(1,lake) = 1.0
      ! recalculate the surface emissivity
      do i = 1, n_land_tile
        fluxes%emis_surft(1, i) = 1.0_r_um - real(tile_lw_grey_albedo(map_tile(1)+i-1), r_um)
      end do
    end if

    if (emis_method_soil /= emis_method_soil_fixed) fluxes%l_emis_surft_set(soil)=.TRUE.

    i_tile = 0
    do i = 1, n_land_tile
      do n = 1, sm_levels
        wt_ext_surft(1,n,i) = tile_water_extract(map_smtile(1)+i_tile)
        i_tile = i_tile + 1
      end do
    end do

    alpha1_sea(1,1) = alpha1_tile(map_tile(1)+first_sea_tile-1)
    ashtf_prime_sea(1,1) = ashtf_prime_tile(map_tile(1)+first_sea_tile-1)
    dtstar_sea(1,1) = dtstar_tile(map_tile(1)+first_sea_tile-1)
    rhokh_sea(1,1) = rhokh_tile(map_tile(1)+first_sea_tile-1)

    i_sice = 0
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      alpha1_sice(1,1,i_sice) = alpha1_tile(map_tile(1)+i-1)
      ashtf_prime(1,1,i_sice) = ashtf_prime_tile(map_tile(1)+i-1)
      rhokh_sice(1,1,i_sice) = rhokh_tile(map_tile(1)+i-1)
      dtstar_sice(1,1,i_sice) = dtstar_tile(map_tile(1)+i-1)
    end do
    z0hssi(1,1) = z0h_tile(map_tile(1)+first_sea_ice_tile-1)
    z0mssi(1,1) = z0m_tile(map_tile(1)+first_sea_ice_tile-1)
    chr1p5m_sice(1,1) = chr1p5m_tile(map_tile(1)+first_sea_ice_tile-1)

    k_blend_tq(1,1) = blend_height_tq(map_2d(1))
    u_s(1,1) = ustar(map_2d(1))
    smc_soilt(1) = soil_moist_avail(map_2d(1))
    zhnl(1,1) = zh_nonloc(map_2d(1))
    zh(1,1) = zh_2d(map_2d(1))
    zhsc(1,1) = zhsc_2d(map_2d(1))
    zlcl(1,1) = real(z_lcl(map_2d(1)), r_um)
    dzh(1,1) = real(inv_depth(map_2d(1)), r_um)
    qcl_inv_top(1,1) = real(qcl_at_inv_top(map_2d(1)), r_um)
    zlcl_mix = 0.0_r_um

    bl_type_1(1,1) = bl_type_ind(map_bl(1)+0)
    bl_type_2(1,1) = bl_type_ind(map_bl(1)+1)
    bl_type_3(1,1) = bl_type_ind(map_bl(1)+2)
    bl_type_4(1,1) = bl_type_ind(map_bl(1)+3)
    bl_type_5(1,1) = bl_type_ind(map_bl(1)+4)
    bl_type_6(1,1) = bl_type_ind(map_bl(1)+5)
    bl_type_7(1,1) = bl_type_ind(map_bl(1)+6)

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
      q_latest(1,1,k)   = m_v(map_wth(1) + k)
      qcl_latest(1,1,k) = m_cl(map_wth(1) + k)
      qcf_latest(1,1,k) = m_ci(map_wth(1) + k)
      ! Set qcf2_latest to zero for now, until needed by CASIM coupling.
      qcf2_latest(1,1,k) = 0.0_r_um
    end do

    !-----------------------------------------------------------------------
    ! fields for bimodal cloud scheme
    !-----------------------------------------------------------------------
    do k = 1, nlayers
      tgrad_in(1,1,k)    = dsldzm(map_wth(1) + k)
      tau_dec_in(1,1,k)  = tau_dec_bm(map_wth(1) + k)
      tau_hom_in(1,1,k)  = tau_hom_bm(map_wth(1) + k)
      tau_mph_in(1,1,k)  = tau_mph_bm(map_wth(1) + k)
      wvar_in(1,1,k)     = wvar(map_wth(1) + k )
    end do

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
          , u_0, v_0, land_index, cca_3d, lcbase, lcbase0, ccb0, cct0   &
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
          , ashtf_prime_surft,gc_surft,aresist_surft                    &
          , resft,rhokh_sice,rhokh_sea                                  &
          , chr1p5m_sice                                                &
          , fland, flandg, flandg_u,flandg_v                            &
          , t_soil_soilt, snow_surft, sstfrz                            &
    ! IN JULES variables for STASH
          , gs_gb, npp_gb                                               &
          , resp_s_tot_soilt                                            &
          , catch_surft                                                 &
          , co2_emits, co2flux                                          &
    ! INOUT diagnostic info
          , STASHwork3, STASHwork9                                      &
    ! SCM Diagnostics (dummy in full UM) & BL diags
          , nSCMDpkgs, L_SCMDiags, bl_diag, sf_diag                     &
    ! INOUT (Note ti_sice_ncat and ti_sice are IN if l_sice_multilayers=T)
          , TScrnDcl_SSI, TScrnDcl_surft, tStbTrans                     &
          , cca0, ccw0, fqw, ftl, taux, tauy, rhokh                     &
          , dtstar_surft,dtstar_sea,dtstar_sice,ti_sice_ncat &
          , area_cloud_fraction, bulk_cloud_fraction                    &
          , t_latest, q_latest, qcl_latest, qcf_latest, qcf2_latest     &
          , cf_latest, cfl_latest, cff_latest                           &
          , R_u, R_v, R_w, cloud_fraction_liquid, cloud_fraction_frozen &
          , sum_eng_fluxes,sum_moist_flux, rhcpt, dissip_u, dissip_v    &
    ! IN arguments for bimodal scheme (dzh is a duplicate argument for
    ! UM compatibility)
          , tgrad_in, wvar_in, tau_dec_in, tau_hom_in, tau_mph_in, dzh  &
    ! INOUT tracer fields
          , aerosol, free_tracers,  resist_b,  resist_b_surft           &
          , dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6 &
          , drydep2, so2, dms, so4_aitken, so4_accu, so4_diss, nh3      &
          , soot_new, soot_aged, soot_cld, bmass_new, bmass_aged        &
          , bmass_cld, ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss &
          , co2, ozone_tracer                                           &
    ! INOUT additional variables for JULES
          , tstar_surft, epot_surft                                     &
          , radnet_sice,olr,tstar_sice_ncat,tstar_ssi                   &
          , tstar_sea,taux_land,taux_ssi,tauy_land,tauy_ssi,Error_code  &
    ! JULES TYPES (IN OUT)
          , crop_vars, ainfo, aerotype, progs, coast, jules_vars        &
          , fluxes &
          , lake_vars &
          , forcing &
          , trifctltype                                                  &
          !rivers, &
          !veg3_parm, &
          !veg3_field, &
          !chemvars
    ! OUT fields
          , surf_ht_flux_land, zlcl_mix                                 &
          , theta_star_surf, qv_star_surf                               &
    ! OUT additional variables for JULES
          , tstar, ti_sice, tstar_land,tstar_sice                       &
    ! OUT fields for coupling to the wave model
          , uwind_wav, vwind_wav                                        &
    ! OUT Fields for use in conservation of atmospheric CO2
          , co2_flux_tot, target_mass                                   &
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
      if ( scheme == scheme_smith .or. scheme == scheme_bimodal ) then
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

      if (.not. associated(snomlt_surf_htf, empty_real_data) ) then
        snomlt_surf_htf(map_2d(1)) = sf_diag%snomlt_surf_htf(1,1)
      end if

      if (.not. associated(soil_evap, empty_real_data) ) then
        soil_evap(map_2d(1)) = fluxes%esoil_ij_soilt(1,1,1)
      end if

      if (.not. associated(soil_surf_ht_flux, empty_real_data) ) then
        soil_surf_ht_flux(map_2d(1)) = coast%surf_ht_flux_land_ij(1,1)
      end if

      do k = 0, bl_levels-1
        heat_flux_bl(map_w3(1)+k) = ftl(1,1,k+1)
        moist_flux_bl(map_w3(1)+k) = fqw(1,1,k+1)
      end do

      tile_heat_flux(map_tile(1)+first_sea_tile-1) = 0.0_r_def
      tile_moisture_flux(map_tile(1)+first_sea_tile-1) = 0.0_r_def
      ! Update land tiles
      do i = 1, n_land_tile
        tile_temperature(map_tile(1)+i-1) = real(tstar_surft(1, i), r_def)
        screen_temperature(map_tile(1)+i-1) = real(tscrndcl_surft(1, i), r_def)
        tile_heat_flux(map_tile(1)+i-1) = real(fluxes%ftl_surft(1, i), r_def)
        tile_moisture_flux(map_tile(1)+i-1) = real(fluxes%fqw_surft(1, i), r_def)
        ! Sum the fluxes over the land for use in sea point calculation
        if (tile_fraction(map_tile(1)+i-1) > 0.0_r_def) then
        tile_heat_flux(map_tile(1)+first_sea_tile-1) =                         &
             tile_heat_flux(map_tile(1)+first_sea_tile-1)                      &
             + fluxes%ftl_surft(1,i) * tile_fraction(map_tile(1)+i-1)
        tile_moisture_flux(map_tile(1)+first_sea_tile-1) =                     &
             tile_moisture_flux(map_tile(1)+first_sea_tile-1)                  &
             + fluxes%fqw_surft(1,i) * tile_fraction(map_tile(1)+i-1)
        end if
        snowice_sublimation(map_tile(1)+i-1) = real(fluxes%ei_surft(1, i), r_def)
        ! NB - net surface heat flux
        surf_heat_flux(map_tile(1)+i-1) = real(fluxes%surf_htf_surft(1, i), r_def)
        canopy_evap(map_tile(1)+i-1) = real(fluxes%ecan_surft(1, i), r_def)
        snowice_melt(map_tile(1)+i-1) = real(fluxes%melt_surft(1, i), r_def)
      end do

      ! Update sea tile
      tile_temperature(map_tile(1)+first_sea_tile-1) = real(tstar_sea(1,1), r_def)
      screen_temperature(map_tile(1)+first_sea_tile-1) = real(tscrndcl_ssi(1, 1), r_def)
      time_since_transition(map_2d(1)) = tstbtrans(1,1)

      ! Update sea-ice tiles
      i_sice = 0
      do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
        i_sice = i_sice + 1
        tile_temperature(map_tile(1)+i-1) = real(tstar_sice_ncat(1,1,i_sice), r_def)
        tile_heat_flux(map_tile(1)+i-1) = real(fluxes%ftl_sicat(1,1,i_sice), r_def)
        tile_moisture_flux(map_tile(1)+i-1) = real(fluxes%fqw_sicat(1,1,i_sice), r_def)
        snowice_melt(map_tile(1)+i-1) = real(fluxes%sice_melt(1, 1, i_sice), r_def)
        snowice_sublimation(map_tile(1)+i-1) = real(fluxes%ei_sice(1,1,i_sice), r_def)
        canopy_evap(map_tile(1)+i-1) = 0.0_r_def
        ! Sum the fluxes over the sea-ice for use in sea point calculation
        tile_heat_flux(map_tile(1)+first_sea_tile-1) =                         &
             tile_heat_flux(map_tile(1)+first_sea_tile-1)                      &
             + tile_heat_flux(map_tile(1)+i-1) * tile_fraction(map_tile(1)+i-1)
        tile_moisture_flux(map_tile(1)+first_sea_tile-1) =                     &
             tile_moisture_flux(map_tile(1)+first_sea_tile-1)                  &
             + tile_moisture_flux(map_tile(1)+i-1) * tile_fraction(map_tile(1)+i-1)
      end do

      ! Sea-ice bulk temperature
      do i = 1, n_sea_ice_tile
        sea_ice_temperature(map_sice(1)+i-1) = real(ti_sice_ncat(1, 1, i), r_def)
      end do

      ! Sea tile fluxes
      if (tile_fraction(map_tile(1)+first_sea_tile-1) > 0.0_r_def) then
        tile_heat_flux(map_tile(1)+first_sea_tile-1) = ( ftl(1,1,1) -          &
             tile_heat_flux(map_tile(1)+first_sea_tile-1) )                    &
             / tile_fraction(map_tile(1)+first_sea_tile-1)
        tile_moisture_flux(map_tile(1)+first_sea_tile-1) = ( fqw(1,1,1) -      &
             tile_moisture_flux(map_tile(1)+first_sea_tile-1) )                &
             / tile_fraction(map_tile(1)+first_sea_tile-1)
      else
        tile_heat_flux(map_tile(1)+first_sea_tile-1) = 0.0_r_def
        tile_moisture_flux(map_tile(1)+first_sea_tile-1) = 0.0_r_def
      end if
      snowice_sublimation(map_tile(1)+first_sea_tile-1) = 0.0_r_def
      snowice_melt(map_tile(1)+first_sea_tile-1) = 0.0_r_def
      canopy_evap(map_tile(1)+first_sea_tile-1) = 0.0_r_def

      do i = 1, sm_levels
        water_extraction(map_soil(1)+i-1) = real(fluxes%ext_soilt(1, 1, i), r_def)
      end do

      ! diagnostics
      if (.not. associated(t1p5m, empty_real_data) .or.                        &
          .not. associated(q1p5m, empty_real_data) .or.                        &
          .not. associated(rh1p5m, empty_real_data) .or.                       &
          .not. associated(qcl1p5m, empty_real_data) ) then
        call ls_cld(                                                           &
           p_star, rhcpt, 1, bl_levels, 1, 1, ntml, cumulus, .false.,          &
           sf_diag%t1p5m, work_2d_1, sf_diag%q1p5m, qcf_latest, qcl1p5m_loc,   &
           work_2d_2, work_2d_3, error_code )
      end if

      if (.not. associated(t1p5m, empty_real_data) ) then
        t1p5m(map_2d(1)) = sf_diag%t1p5m(1,1)
      end if
      if (.not. associated(q1p5m, empty_real_data) ) then
        q1p5m(map_2d(1)) = sf_diag%q1p5m(1,1)
      end if
      if (.not. associated(qcl1p5m, empty_real_data) ) then
        qcl1p5m(map_2d(1)) = qcl1p5m_loc(1,1)
      end if

      if (.not. associated(rh1p5m, empty_real_data) ) then
        ! qsat needed since q1p5m always a specific humidity
        call qsat(work_2d_1,sf_diag%t1p5m,p_star,pdims%i_end,pdims%j_end)
        rh1p5m(map_2d(1)) = max(0.0, sf_diag%q1p5m(1,1)) * 100.0 / work_2d_1(1,1)
      end if

      if (.not. associated(t1p5m_surft, empty_real_data) ) then
        if (land_field > 0) then
          do i = 1, n_land_tile
            t1p5m_surft(map_tile(1)+i-1) = real(sf_diag%t1p5m_surft(1, i), r_def)
          end do
        else
          do i = 1, n_land_tile
            t1p5m_surft(map_tile(1)+i-1) = 0.0_r_def
          end do
        end if
        do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          t1p5m_surft(map_tile(1)+i-1) = 0.0_r_def
        end do
        t1p5m_surft(map_tile(1)+first_sea_tile-1) = 0.0_r_def
      end if
      if (.not. associated(q1p5m_surft, empty_real_data) ) then
        if (land_field > 0) then
          do i = 1, n_land_tile
            q1p5m_surft(map_tile(1)+i-1) = real(sf_diag%q1p5m_surft(1, i), r_def)
          end do
        else
          do i = 1, n_land_tile
            q1p5m_surft(map_tile(1)+i-1) = 0.0_r_def
          end do
        end if
        do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          q1p5m_surft(map_tile(1)+i-1) = 0.0_r_def
        end do
        q1p5m_surft(map_tile(1)+first_sea_tile-1) = 0.0_r_def
      end if

      if (.not. associated(latent_heat, empty_real_data) ) then
        if (land_field > 0) then
          do i = 1, n_land_tile
            latent_heat(map_tile(1)+i-1) = real(fluxes%le_surft(1, i), r_def)
          end do
        else
          do i = 1, n_land_tile
            latent_heat(map_tile(1)+i-1) = 0.0_r_def
          end do
        end if
        do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          if (tile_fraction(map_tile(1)+i-1) > 0.0_r_def) then
            latent_heat(map_tile(1)+i-1) = (lc + lf) *                    &
               tile_moisture_flux(map_tile(1)+i-1)
          else
            latent_heat(map_tile(1)+i-1) =  0.0_r_def
          end if
        end do
        if (tile_fraction(map_tile(1)+first_sea_tile-1) > 0.0_r_def) then
          latent_heat(map_tile(1)+first_sea_tile-1) = lc *                &
                 tile_moisture_flux(map_tile(1)+first_sea_tile-1)
        else
          latent_heat(map_tile(1)+first_sea_tile-1) =  0.0_r_def
        end if
      end if

      if (.not. associated(surf_ht_flux, empty_real_data) ) then
        if (land_field > 0) then
          do i = 1, n_land_tile
            surf_ht_flux(map_tile(1)+i-1) = real(fluxes%surf_htf_surft(1, i), r_def)
          end do
        else
          do i = 1, n_land_tile
            surf_ht_flux(map_tile(1)+i-1) = 0.0_r_def
          end do
        end if
        i_sice = 0
        do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          i_sice = i_sice + 1
          if (tile_fraction(map_tile(1)+i-1) > 0.0_r_def) then
            surf_ht_flux(map_tile(1)+i-1) =                               &
              real(fluxes%surf_ht_flux_sice(1,1,i_sice), r_def)
          else
            surf_ht_flux(map_tile(1)+i-1) = 0.0_r_def
          end if
        end do
        surf_ht_flux(map_tile(1)+first_sea_tile-1) = 0.0_r_def
      end if

      if (.not. associated(surf_sw_net, empty_real_data) ) then
        if (land_field > 0) then
          do i = 1, n_land_tile
            surf_sw_net(map_tile(1)+i-1) =                                &
                real(sw_down_surf(map_2d(1)), r_um) -                     &
                real(sw_up_tile(map_tile(1)+i-1), r_um)
          end do
        else
          do i = 1, n_land_tile
            surf_sw_net(map_tile(1)+i-1) =  0.0_r_def
          end do
        end if
        do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          if (tile_fraction(map_tile(1)+i-1) > 0.0_r_def) then
            surf_sw_net(map_tile(1)+i-1) =                                &
                real(sw_down_surf(map_2d(1)), r_um) -                     &
                real(sw_up_tile(map_tile(1)+i-1), r_um)
          else
            surf_sw_net(map_tile(1)+i-1) = 0.0_r_def
          end if
        end do
        if (tile_fraction(map_tile(1)+first_sea_tile-1) > 0.0_r_def) then
          surf_sw_net(map_tile(1)+first_sea_tile-1) =                     &
              real(sw_down_surf(map_2d(1)), r_um) -                       &
              real(sw_up_tile(map_tile(1)+first_sea_tile-1), r_um)
        else
          surf_sw_net(map_tile(1)+first_sea_tile-1) = 0.0_r_def
        end if
      end if

      if (.not. associated(surf_radnet, empty_real_data) ) then
        if (land_field > 0) then
          do i = 1, n_land_tile
            surf_radnet(map_tile(1)+i-1) = real(fluxes%radnet_surft(1, i), r_def)
          end do
        else
          do i = 1, n_land_tile
            surf_radnet(map_tile(1)+i-1) = 0.0_r_def
          end do
        end if
        i_sice = 0
        do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          i_sice = i_sice + 1
          if (tile_fraction(map_tile(1)+i-1) > 0.0_r_def) then
            surf_radnet(map_tile(1)+i-1) =                                &
              radnet_sice(1,1,i_sice)/ice_fract_ncat(1,1,i_sice)
          else
            surf_radnet(map_tile(1)+i-1) = 0.0_r_def
          end if
        end do
        if (tile_fraction(map_tile(1)+first_sea_tile-1) > 0.0_r_def) then
          surf_radnet(map_tile(1)+first_sea_tile-1) =                     &
              real(sw_down_surf(map_2d(1)), r_um) -                       &
              real(sw_up_tile(map_tile(1)+first_sea_tile-1), r_um) +      &
              emis_sea * (lw_down(1,1) - sbcon * tstar_sea(1,1) ** 4.0)
        else
          surf_radnet(map_tile(1)+first_sea_tile-1) =  0.0_r_def
        end if
      end if

      if (.not. associated(surf_lw_up, empty_real_data) ) then
        if (land_field > 0) then
          do i = 1, n_land_tile
            surf_lw_up(map_tile(1)+i-1) = real(sf_diag%lw_up_surft(1, i), r_def)
          end do
        else
          do i = 1, n_land_tile
            surf_lw_up(map_tile(1)+i-1) = 0.0_r_def
          end do
        end if
        i_sice = 0
        do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          i_sice = i_sice + 1
          if (tile_fraction(map_tile(1)+i-1) > 0.0_r_def) then
            surf_lw_up(map_tile(1)+i-1) =                                &
              sf_diag%lw_up_sice_weighted_cat(1,1,i_sice)/ice_fract_ncat(1,1,i_sice)
          else
            surf_lw_up(map_tile(1)+i-1) = 0.0_r_def
          end if
        end do
        if (tile_fraction(map_tile(1)+first_sea_tile-1) > 0.0_r_def) then
          surf_lw_up(map_tile(1)+first_sea_tile-1) =                     &
              (1.0 - emis_sea) * lw_down(1,1) +                          &
               emis_sea * sbcon * tstar_sea(1,1) ** 4.0
        else
          surf_lw_up(map_tile(1)+first_sea_tile-1) =  0.0_r_def
        end if
      end if
      if (.not. associated(surf_lw_down, empty_real_data) ) then
        if (land_field > 0) then
          do i = 1, n_land_tile
            surf_lw_down(map_tile(1)+i-1) = real(sf_diag%lw_down_surft(1, i), r_def)
          end do
        else
          do i = 1, n_land_tile
            surf_lw_down(map_tile(1)+i-1) = 0.0_r_def
          end do
        end if
        i_sice = 0
        do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          if (tile_fraction(map_tile(1)+i-1) > 0.0_r_def) then
            surf_lw_down(map_tile(1)+i-1) = lw_down(1,1)
          else
            surf_lw_down(map_tile(1)+i-1) =  0.0_r_def
          end if
        end do
        if (tile_fraction(map_tile(1)+first_sea_tile-1) > 0.0_r_def) then
          surf_lw_down(map_tile(1)+first_sea_tile-1) = lw_down(1,1)
        else
          surf_lw_down(map_tile(1)+first_sea_tile-1) =  0.0_r_def
        end if
      end if

    endif  ! outer = outer_iterations

    ! deallocate diagnostics deallocated in atmos_physics2
    call dealloc_bl_imp(bl_diag)
    call dealloc_sf_imp(sf_diag)
    call dealloc_sf_expl(sf_diag)

    ! set this back to 1 before exit
    land_field = 1

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

    call aero_nullify(aerotype)
    call aero_dealloc(aero_data)

    call jules_vars_dealloc(jules_vars_data)
    call jules_vars_nullify(jules_vars)

    call prognostics_nullify(progs)
    call prognostics_dealloc(progs_data)

    call fluxes_nullify(fluxes)
    call fluxes_dealloc(fluxes_data)

  end subroutine bl_imp_code

end module bl_imp_kernel_mod
