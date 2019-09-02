!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2019.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!> @brief Initialise Jules fields

module initial_jules_kernel_mod

  use argument_mod, only: arg_type, &
    GH_FIELD, GH_INTEGER, GH_WRITE, GH_READ, &
    ANY_SPACE_1, ANY_SPACE_2, ANY_SPACE_3, ANY_SPACE_4, ANY_SPACE_5, &
    CELLS
  use constants_mod, only: r_def, i_def
  use kernel_mod, only: kernel_type
  use surface_config_mod, only: surf_tile_fracs

  implicit none

  private

  !> Kernel metadata for Psyclone 
  type, public, extends(kernel_type) :: initial_jules_kernel_type
      private
      type(arg_type) :: meta_args(40) = (/              &
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1), & ! tile_fraction
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_2), & ! leaf_area_index
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_2), & ! canopy_height
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! sd_orog
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! peak_to_trough_orog
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! silhouette_area_orog
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! soil_albedo
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! soil_roughness
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! albedo_obs_sw
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! albedo_obs_vis
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! albedo_obs_nir
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_4), & ! soil_moist_wilt
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_4), & ! soil_moist_crit
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_4), & ! soil_moist_sat
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_4), & ! soil_cond_sat
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_4), & ! soil_thermal_cap
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! soil_thermal_cond
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_4), & ! soil_suction_sat
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_4), & ! clapp_horn_b
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! soil_carbon_content
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1), & ! tile_temperature
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1), & ! tile_snow_mass
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1), & ! tile_snow_rgrain
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1), & ! n_snow_layers
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1), & ! snow_depth
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! snow_soot
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! surface_conductance
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1), & ! canopy_water
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1), & ! sw_up_tile
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_4), & ! soil_temperature
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_4), & ! soil_moisture
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_4), & ! unfrozen_soil_moisture
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_4), & ! frozen_soil_moisture
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! chloro_sea
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_5), & ! sea_ice_thickness
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_5), & ! sea_ice_pond_frac
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_5), & ! sea_ice_pond_depth
          arg_type(GH_INTEGER, GH_READ               ), & ! n_surf_tile
          arg_type(GH_INTEGER, GH_READ               ), & ! n_sea_ice_tile
          arg_type(GH_INTEGER, GH_READ               )  & ! n_soil_levs
          /)
      integer :: iterates_over = CELLS

  contains
      procedure, nopass :: initial_jules_code
  end type

  ! overload the default structure constructor for function space
  interface initial_jules_kernel_type
      module procedure initial_jules_kernel_constructor
  end interface

  public initial_jules_code
contains

  function initial_jules_kernel_constructor() result(self)
    implicit none
    type(initial_jules_kernel_type) :: self
    return
  end function initial_jules_kernel_constructor

  !> @param[in]  nlayers            The number of layers
  !> @param[out] tile_fraction      Surface tile fractions
  !> @param[out] leaf_area_index    Leaf Area Index
  !> @param[out] canopy_height      Canopy height
  !> @param[out] sd_orog            Standard deviation of orography
  !> @param[out] peak_to_trough_orog  Half of peak-to-trough height over root(2) of orography
  !> @param[out] silhouette_area_orog Silhouette area of orography
  !> @param[out] soil_albedo        Snow-free soil albedo
  !> @param[out] soil_roughness     Bare soil surface roughness length
  !> @param[out] albedo_obs_sw      Observed snow-free shortwave albedo
  !> @param[out] albedo_obs_vis     Observed snow-free visible albedo
  !> @param[out] albedo_obs_nir     Observed snow-free near-IR albedo
  !> @param[out] soil_moist_wilt    Volumetric soil moisture at wilting point
  !> @param[out] soil_moist_crit    Volumetric soil moisture at critical point
  !> @param[out] soil_moist_sat     Volumetric soil moisture at saturation
  !> @param[out] soil_cond_sat      Saturated soil conductivity
  !> @param[out] soil_thermal_cap   Soil thermal capacity
  !> @param[out] soil_thermal_cond  Soil thermal conductivity
  !> @param[out] soil_suction_sat   Saturated soil water suction
  !> @param[out] clapp_horn_b       Clapp and Hornberger b coefficient
  !> @param[out] soil_carbon_content Soil carbon content
  !> @param[out] tile_temperature   Surface tile temperatures
  !> @param[out] tile_snow_mass     Snow mass on tiles (kg/m2)
  !> @param[out] tile_snow_rgrain   Snow grain size on tiles (microns)
  !> @param[out] n_snow_layers      Number of snow layers on tiles
  !> @param[out] snow_depth         Snow depth on tiles
  !> @param[out] snow_soot          Snow soot content (kg/kg)
  !> @param[out] surface_conductance Surface conductance
  !> @param[out] canopy_water       Canopy water on each tile
  !> @param[out] sw_up_tile         Upwelling SW radiation on tiles
  !> @param[out] soil_temperature   Soil temperature
  !> @param[out] soil_moisture      Soil moisture content (kg m-2)
  !> @param[out] unfrozen_soil_moisture Unfrozen soil moisture proportion
  !> @param[out] frozen_soil_moisture Frozen soil moisture proportion
  !> @param[out] chloro_sea         Chlorophyll content of the sea
  !> @param[out] sea_ice_thickness  Sea ice thickness (m)
  !> @param[out] sea_ice_pond_frac  Meltpond fraction on sea ice
  !> @param[out] sea_ice_pond_depth Meltpond depth on sea ice (m)
  !> @param[in]  n_surf_tile        Number of surface tiles
  !> @param[in]  n_sea_ice_tile     Number of sea ice tiles
  !> @param[in]  n_soil_levs        Number of soil levels
  !> @param[in]  ndf_tile           Number of DOFs per cell for tiles
  !> @param[in]  undf_tile          Number of total DOFs for tiles
  !> @param[in]  map_tile           Dofmap for cell for surface tiles
  !> @param[in]  ndf_pft            Number of DOFs per cell for PFTs
  !> @param[in]  undf_pft           Number of total DOFs for PFTs
  !> @param[in]  map_pft            Dofmap for cell for PFTs
  !> @param[in]  ndf_2d             Number of DOFs per cell for 2D fields
  !> @param[in]  undf_2d            Number of total DOFs for 2D fields
  !> @param[in]  map_2d             Dofmap for cell for 2D fields
  !> @param[in]  ndf_soil           Number of DOFs per cell for soil levels
  !> @param[in]  undf_soil          Number of total DOFs for soil levels
  !> @param[in]  map_soil           Dofmap for cell for soil levels
  !> @param[in]  ndf_sice           Number of DOFs per cell for sea ice tiles
  !> @param[in]  undf_sice          Number of total DOFs for sea ice tiles
  !> @param[in]  map_sice           Dofmap for cell for sea ice tiles
  subroutine initial_jules_code(nlayers,                       &
                                tile_fraction,                 &
                                leaf_area_index,               &
                                canopy_height,                 &
                                sd_orog,                       &
                                peak_to_trough_orog,           &
                                silhouette_area_orog,          &
                                soil_albedo,                   &
                                soil_roughness,                &
                                albedo_obs_sw,                 &
                                albedo_obs_vis,                &
                                albedo_obs_nir,                &
                                soil_moist_wilt,               &
                                soil_moist_crit,               &
                                soil_moist_sat,                &
                                soil_cond_sat,                 &
                                soil_thermal_cap,              &
                                soil_thermal_cond,             &
                                soil_suction_sat,              &
                                clapp_horn_b,                  &
                                soil_carbon_content,           &
                                tile_temperature,              &
                                tile_snow_mass,                &
                                tile_snow_rgrain,              &
                                n_snow_layers,                 &
                                snow_depth,                    &
                                snow_soot,                     &
                                surface_conductance,           &
                                canopy_water,                  &
                                sw_up_tile,                    &
                                soil_temperature,              &
                                soil_moisture,                 &
                                unfrozen_soil_moisture,        &
                                frozen_soil_moisture,          &
                                chloro_sea,                    &
                                sea_ice_thickness,             &
                                sea_ice_pond_frac,             &
                                sea_ice_pond_depth,            &
                                n_surf_tile,                   &
                                n_sea_ice_tile,                &
                                n_soil_levs,                   &
                                ndf_tile, undf_tile, map_tile, &
                                ndf_pft, undf_pft, map_pft,    &
                                ndf_2d, undf_2d, map_2d,       &
                                ndf_soil, undf_soil, map_soil, &
                                ndf_sice, undf_sice, map_sice)

      implicit none

      ! Arguments
      integer(kind=i_def), intent(in) :: nlayers, n_surf_tile, n_sea_ice_tile, &
                                         n_soil_levs
      integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
      integer(kind=i_def), intent(in) :: map_tile(ndf_tile)
      integer(kind=i_def), intent(in) :: ndf_pft, undf_pft
      integer(kind=i_def), intent(in) :: map_pft(ndf_pft)
      integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
      integer(kind=i_def), intent(in) :: map_2d(ndf_2d)
      integer(kind=i_def), intent(in) :: ndf_soil, undf_soil
      integer(kind=i_def), intent(in) :: map_soil(ndf_soil)
      integer(kind=i_def), intent(in) :: ndf_sice, undf_sice
      integer(kind=i_def), intent(in) :: map_sice(ndf_sice)

      real(kind=r_def), intent(out) :: tile_fraction(undf_tile)
      real(kind=r_def), intent(out) :: tile_temperature(undf_tile)
      real(kind=r_def), intent(out) :: tile_snow_mass(undf_tile)
      real(kind=r_def), intent(out) :: tile_snow_rgrain(undf_tile)
      real(kind=r_def), intent(out) :: n_snow_layers(undf_tile)
      real(kind=r_def), intent(out) :: snow_depth(undf_tile)
      real(kind=r_def), intent(out) :: canopy_water(undf_tile)
      real(kind=r_def), intent(out) :: sw_up_tile(undf_tile)

      real(kind=r_def), intent(out) :: leaf_area_index(undf_pft)
      real(kind=r_def), intent(out) :: canopy_height(undf_pft)

      real(kind=r_def), intent(out) :: sd_orog(undf_2d)
      real(kind=r_def), intent(out) :: peak_to_trough_orog(undf_2d)
      real(kind=r_def), intent(out) :: silhouette_area_orog(undf_2d)
      real(kind=r_def), intent(out) :: soil_albedo(undf_2d)
      real(kind=r_def), intent(out) :: soil_roughness(undf_2d)
      real(kind=r_def), intent(out) :: albedo_obs_sw(undf_2d)
      real(kind=r_def), intent(out) :: albedo_obs_vis(undf_2d)
      real(kind=r_def), intent(out) :: albedo_obs_nir(undf_2d)
      real(kind=r_def), intent(out) :: soil_thermal_cond(undf_2d)
      real(kind=r_def), intent(out) :: soil_carbon_content(undf_2d)
      real(kind=r_def), intent(out) :: snow_soot(undf_2d)
      real(kind=r_def), intent(out) :: surface_conductance(undf_2d)
      real(kind=r_def), intent(out) :: chloro_sea(undf_2d)

      real(kind=r_def), intent(out) :: soil_moist_wilt(undf_soil)
      real(kind=r_def), intent(out) :: soil_moist_crit(undf_soil)
      real(kind=r_def), intent(out) :: soil_moist_sat(undf_soil)
      real(kind=r_def), intent(out) :: soil_cond_sat(undf_soil)
      real(kind=r_def), intent(out) :: soil_thermal_cap(undf_soil)
      real(kind=r_def), intent(out) :: soil_suction_sat(undf_soil)
      real(kind=r_def), intent(out) :: clapp_horn_b(undf_soil)
      real(kind=r_def), intent(out) :: soil_temperature(undf_soil)
      real(kind=r_def), intent(out) :: soil_moisture(undf_soil)
      real(kind=r_def), intent(out) :: unfrozen_soil_moisture(undf_soil)
      real(kind=r_def), intent(out) :: frozen_soil_moisture(undf_soil)

      real(kind=r_def), intent(out) :: sea_ice_thickness(undf_sice)
      real(kind=r_def), intent(out) :: sea_ice_pond_frac(undf_sice)
      real(kind=r_def), intent(out) :: sea_ice_pond_depth(undf_sice)

      !Internal variables
      integer(kind=i_def) :: i


      ! Tile fraction
      do i = 1, n_surf_tile
          tile_fraction(map_tile(i)) = surf_tile_fracs(i)
      end do

      ! Leaf area index (UKV values: lai_pft)
      leaf_area_index(map_pft(1)) = 5.0_r_def
      leaf_area_index(map_pft(2)) = 4.0_r_def
      leaf_area_index(map_pft(3)) = 2.0_r_def
      leaf_area_index(map_pft(4)) = 4.0_r_def
      leaf_area_index(map_pft(5)) = 1.0_r_def

      ! Canopy height (UKV values: canht_pft)
      canopy_height(map_pft(1)) = 19.01_r_def
      canopy_height(map_pft(2)) = 16.38_r_def
      canopy_height(map_pft(3)) =  1.46_r_def
      canopy_height(map_pft(4)) =  1.26_r_def
      canopy_height(map_pft(5)) =  1.59_r_def

      ! Standard deviation of orography
      sd_orog(map_2d(1)) = 0.0_r_def

      ! Half of peak-to-trough height over root(2) of orography (ho2r2_orog_gb)
      peak_to_trough_orog(map_2d(1)) = 0.0_r_def

      ! Silhouette area of orography (sil_orog_land_gb)
      silhouette_area_orog(map_2d(1)) = 0.0_r_def
      
      ! Snow-free soil albedo (Cardington UKV value: albsoil_soilt)
      soil_albedo(map_2d(1)) = 0.1065_r_def

      ! Bare soil surface roughness length
      soil_roughness(map_2d(1)) = 0.0_r_def
      
      ! Observed snow-free shortwave albedo
      albedo_obs_sw(map_2d(1)) = 0.0_r_def

      ! Observed snow-free visible albedo
      albedo_obs_vis(map_2d(1)) = 0.0_r_def

      ! Observed snow-free near infra-red albedo
      albedo_obs_nir(map_2d(1)) = 0.0_r_def

      ! Volumetric soil moisture at wilting point (smvcwt_soilt)
      do i = 1, n_soil_levs
        soil_moist_wilt(map_soil(i)) = 0.226216_r_def
      end do

      ! Volumetric soil moisture at critical point (smvccl_soilt)
      do i = 1, n_soil_levs
        soil_moist_crit(map_soil(i)) = 0.340808_r_def
      end do

      ! Volumetric soil moisture at saturation (smvcst_soilt)
      do i = 1, n_soil_levs
        soil_moist_sat(map_soil(i)) = 0.4489_r_def
      end do

      ! Saturated soil conductivity (satcon_soilt)
      do i = 1, n_soil_levs
        soil_cond_sat(map_soil(i)) = 0.00286_r_def
      end do

      ! Soil thermal capacity (hcap_soilt)
      do i = 1, n_soil_levs
        soil_thermal_cap(map_soil(i)) = 1228100.0_r_def
      end do

      ! Soil thermal conductivity (hcon_soilt)
      soil_thermal_cond(map_2d(1)) = 0.230211_r_def

      ! Saturated soil water suction (sathh_soilt)
      do i = 1, n_soil_levs
        soil_suction_sat(map_soil(i)) = 0.31282_r_def
      end do

      ! Clapp and Hornberger b coefficient (bexp_soilt)
      do i = 1, n_soil_levs
        clapp_horn_b(map_soil(i)) = 9.3080_r_def
      end do

      ! Soil carbon content (cs_pool_gb_um)
      soil_carbon_content(map_2d(1)) = 12.1_r_def



      ! Tile temperature
      do i = 1, n_surf_tile
        tile_temperature(map_tile(i)) = 295.0_r_def
      end do

      ! Snow mass on tiles (kg/m2, snow_surft)
      do i = 1, n_surf_tile
        tile_snow_mass(map_tile(i)) = 0.0_r_def
      end do

      ! Snow grain size on tiles (microns)
      do i = 1, n_surf_tile
        tile_snow_rgrain(map_tile(i)) = 300.0_r_def
      end do

      ! Number of snow layers on tiles (nsnow_surft)
      do i = 1, n_surf_tile
        n_snow_layers(map_tile(i)) = 0.0_r_def
      end do

      ! Snow depth on tiles (snowdepth_surft)
      do i = 1, n_surf_tile
        snow_depth(map_tile(i)) = 0.0_r_def
      end do      

      ! Snow soot content (kg/kg)
      snow_soot(map_2d(1)) = 0.0_r_def

      ! Surface conductance (gs_gb)
      surface_conductance(map_2d(1)) = 1.0_r_def

      ! Canopy water on each tile (canopy_surft)
      do i = 1, n_surf_tile
        canopy_water(map_tile(i)) = 0.0_r_def
      end do

      ! Upwelling SW on tiles - set to 0 for testing without radiation
      do i = 1, n_surf_tile
        sw_up_tile(map_tile(i)) = 0.0_r_def
      end do

      ! Soil temperature (t_soil_soilt)
      soil_temperature(map_soil(1)) = 285.0_r_def
      soil_temperature(map_soil(2)) = 280.0_r_def
      soil_temperature(map_soil(3)) = 275.0_r_def
      soil_temperature(map_soil(4)) = 275.0_r_def

      ! Soil moisture content (kg m-2, soil_layer_moisture)
      soil_moisture(map_soil(1)) = 40.0_r_def
      soil_moisture(map_soil(2)) = 90.0_r_def
      soil_moisture(map_soil(3)) = 150.0_r_def
      soil_moisture(map_soil(4)) = 460.0_r_def

      ! Unfrozen soil moisture proportion (sthu_soilt)
      unfrozen_soil_moisture(map_soil(1)) = 0.89107_r_def
      unfrozen_soil_moisture(map_soil(2)) = 0.80196_r_def
      unfrozen_soil_moisture(map_soil(3)) = 0.51408_r_def
      unfrozen_soil_moisture(map_soil(4)) = 0.51236_r_def

      ! Frozen soil moisture proportion (sthf_soilt)
      frozen_soil_moisture(map_soil(1)) = 0.0_r_def
      frozen_soil_moisture(map_soil(2)) = 0.0_r_def
      frozen_soil_moisture(map_soil(3)) = 0.0_r_def
      frozen_soil_moisture(map_soil(4)) = 0.0_r_def

      ! Chlorophyll content of the sea near-surface (kg/m3)
      chloro_sea(map_2d(1)) = 5.0e-7_r_def

      ! Variables on sea ice tiles
      do i = 1, n_sea_ice_tile
        ! Sea ice thickness (m)
        sea_ice_thickness(map_sice(i)) = 0.0_r_def
        ! Meltpond fraction
        sea_ice_pond_frac(map_sice(i)) = 0.0_r_def
        ! Meltpond depth (m)
        sea_ice_pond_depth(map_sice(i)) = 0.0_r_def
      end do

  end subroutine initial_jules_code

end module initial_jules_kernel_mod
