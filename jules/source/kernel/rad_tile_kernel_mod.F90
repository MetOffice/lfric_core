!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Interface to Jules for surface tile radiative properties

module rad_tile_kernel_mod

use argument_mod,      only : arg_type, func_type,       &
                              GH_FIELD, GH_READ,         &
                              GH_WRITE, GH_INC, CELLS,   &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_3, &
                              ANY_DISCONTINUOUS_SPACE_4, &
                              ANY_DISCONTINUOUS_SPACE_5, &
                              ANY_DISCONTINUOUS_SPACE_6, &
                              ANY_DISCONTINUOUS_SPACE_7
use fs_continuity_mod, only:  W3, WTheta
use constants_mod,     only : r_def, i_def, r_um, i_um
use kernel_mod,        only : kernel_type

implicit none

private

public :: rad_tile_kernel_type
public :: rad_tile_code

!------------------------------------------------------------------------------
! Public types
!------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, extends(kernel_type) :: rad_tile_kernel_type
  private
  type(arg_type) :: meta_args(26) = (/                          &
       arg_type(GH_FIELD, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! tile_sw_direct_albedo
       arg_type(GH_FIELD, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! tile_sw_diffuse_albedo
       arg_type(GH_FIELD, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! tile_lw_albedo
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! tile_fraction
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! leaf_area_index
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! canopy_height
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! sd_orog
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! soil_albedo
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! soil_roughness
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! albedo_obs_vis
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! albedo_obs_nir
       arg_type(GH_FIELD, GH_WRITE, ANY_DISCONTINUOUS_SPACE_6), & ! albedo_obs_scaling
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! tile_temperature
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! tile_snow_mass
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! tile_snow_rgrain
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! snow_depth
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! snowpack_density
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! snow_soot
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! chloro_sea
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_7), & ! sea_ice_thickness
       arg_type(GH_FIELD, GH_READ,  W3),                        & ! u1_in_w3
       arg_type(GH_FIELD, GH_READ,  W3),                        & ! u2_in_w3
       arg_type(GH_FIELD, GH_READ,  W3),                        & ! height_w3
       arg_type(GH_FIELD, GH_READ,  WTHETA),                    & ! height_wth
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_5), & ! z0msea
       arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_5)  & ! cos_zenith_angle
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: rad_tile_code
end type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

! @param[in]  nlayers                Number of layers
! @param[out] tile_sw_direct_albedo  SW direct tile albedos
! @param[out] tile_sw_diffuse_albedo SW diffuse tile albedos
! @param[out] tile_lw_albedo         LW tile albedos
! @param[in]  tile_fraction          Surface tile fractions
! @param[in]  leaf_area_index        Leaf Area Index
! @param[in]  canopy_height          Canopy height
! @param[in]  sd_orog                Standard deviation of orography
! @param[in]  soil_albedo            Snow-free soil albedo
! @param[in]  soil_roughness         Bare soil surface roughness length
! @param[in]  albedo_obs_vis         Observed snow-free visible albedo
! @param[in]  albedo_obs_nir         Observed snow-free near-IR albedo
! @param[out] albedo_obs_scaling     Scaling factor to adjust albedos by
! @param[in]  tile_temperature       Surface tile temperatures
! @param[in]  tile_snow_mass         Snow mass on tiles (kg/m2)
! @param[in]  tile_snow_rgrain       Snow grain size on tiles (microns)
! @param[in]  snow_depth             Snow depth on tiles (m)
! @param[in]  snowpack_density       Density of snow on ground (kg m-3)
! @param[in]  snow_soot              Snow soot content (kg/kg)
! @param[in]  chloro_sea             Chlorophyll content of the sea
! @param[in]  sea_ice_thickness      Sea ice thickness (m)
! @param[in]  u1_in_w3               'Zonal' wind in density space
! @param[in]  u2_in_w3               'Meridional' wind in density space
! @param[in]  height_w3              Height of w3 levels above mean sea level
! @param[in]  height_wth             Height of wth levels above mean sea level
! @param[in]  z0msea                 Roughness length of sea
! @param[in]  cos_zenith_angle       Cosine of the stellar zenith angle
! @param[in]  ndf_sw_tile            DOFs per cell for tiles and sw bands
! @param[in]  undf_sw_tile           total DOFs for tiles and sw bands
! @param[in]  map_sw_tile            Dofmap for cell at the base of the column
! @param[in]  ndf_lw_tile            DOFs per cell for tiles and lw bands
! @param[in]  undf_lw_tile           total DOFs for tiles and lw bands
! @param[in]  map_lw_tile            Dofmap for cell at the base of the column
! @param[in]  ndf_tile               Number of DOFs per cell for tiles
! @param[in]  undf_tile              Number of total DOFs for tiles
! @param[in]  map_tile               Dofmap for cell at the base of the column
! @param[in]  ndf_pft                Number of DOFs per cell for PFTs
! @param[in]  undf_pft               Number of total DOFs for PFTs
! @param[in]  map_pft                Dofmap for cell at the base of the column
! @param[in]  ndf_2d                 Number of DOFs per cell for 2D fields
! @param[in]  undf_2d                Number of total DOFs for 2D fields
! @param[in]  map_2d                 Dofmap for cell at the base of the column
! @param[in]  ndf_scal               Number of DOFs per cell for albedo scaling
! @param[in]  undf_scal              Number of total DOFs for albedo scaling
! @param[in]  map_scal               Dofmap for cell at the base of the column
! @param[in]  ndf_sice               Number of DOFs per cell for sea ice tiles
! @param[in]  undf_sice              Number of total DOFs for sea ice tiles
! @param[in]  map_sice               Dofmap for cell at the base of the column
! @param[in]  ndf_w3                 Number of DOFs per cell for density space
! @param[in]  undf_w3                Number of unique DOFs for density space
! @param[in]  map_w3                 Dofmap for cell at the base of the column
! @param[in]  ndf_wth                Number of DOFs per cell for theta space
! @param[in]  undf_wth               Number of unique DOFs for theta space
! @param[in]  map_wth                Dofmap for cell at the base of the column
subroutine rad_tile_code(nlayers,                                &
                         tile_sw_direct_albedo,                  &
                         tile_sw_diffuse_albedo,                 &
                         tile_lw_albedo,                         &
                         tile_fraction,                          &
                         leaf_area_index,                        &
                         canopy_height,                          &
                         sd_orog,                                &
                         soil_albedo,                            &
                         soil_roughness,                         &
                         albedo_obs_vis,                         &
                         albedo_obs_nir,                         &
                         albedo_obs_scaling,                     &
                         tile_temperature,                       &
                         tile_snow_mass,                         &
                         tile_snow_rgrain,                       &
                         snow_depth,                             &
                         snowpack_density,                       &
                         snow_soot,                              &
                         chloro_sea,                             &
                         sea_ice_thickness,                      &
                         u1_in_w3,                               &
                         u2_in_w3,                               &
                         height_w3,                              &
                         height_wth,                             &
                         z0msea,                                 &
                         cos_zenith_angle,                       &
                         ndf_sw_tile, undf_sw_tile, map_sw_tile, &
                         ndf_lw_tile, undf_lw_tile, map_lw_tile, &
                         ndf_tile, undf_tile, map_tile,          &
                         ndf_pft, undf_pft, map_pft,             &
                         ndf_2d, undf_2d, map_2d,                &
                         ndf_scal, undf_scal, map_scal,          &
                         ndf_sice, undf_sice, map_sice,          &
                         ndf_w3, undf_w3, map_w3,                &
                         ndf_wth, undf_wth, map_wth)

  use socrates_init_mod, only: &
    n_sw_band, sw_wavelength_short, sw_wavelength_long, sw_weight_blue, &
    n_lw_band
  use jules_control_init_mod, only: &
    n_surf_tile, n_land_tile, n_sea_tile, n_sea_ice_tile, &
    first_sea_tile, first_sea_ice_tile
  use surface_config_mod, only: albedo_obs
  use nlsizes_namelist_mod, only: row_length, rows, land_field, ntiles
  use jules_surface_types_mod, only: ntype, npft, ice
  use nvegparm, only: emis_nvg
  use pftparm, only: emis_pft
  use jules_sea_seaice_mod, only: nice, nice_use, emis_sea, emis_sice
  use ancil_info, only: sea_pts, sea_index, ssi_index, sice_pts_ncat,       &
                        sice_index_ncat, sice_frac_ncat, rad_nband,         &
                        l_lice_point
  use tilepts_mod, only: tilepts
  use sparm_mod, only: sparm
  use surf_couple_radiation_mod, only: surf_couple_radiation

  ! Horizontally varying information used from modules in Jules
  use jules_mod, only: albobs_scaling_surft
  use prognostics, only: snowdepth_surft, rho_snow_grnd_surft

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_sw_tile, undf_sw_tile
  integer(i_def), intent(in) :: map_sw_tile(ndf_sw_tile)
  integer(i_def), intent(in) :: ndf_lw_tile, undf_lw_tile
  integer(i_def), intent(in) :: map_lw_tile(ndf_lw_tile)
  integer(i_def), intent(in) :: ndf_tile, undf_tile
  integer(i_def), intent(in) :: map_tile(ndf_tile)
  integer(i_def), intent(in) :: ndf_pft, undf_pft
  integer(i_def), intent(in) :: map_pft(ndf_pft)
  integer(i_def), intent(in) :: ndf_2d, undf_2d
  integer(i_def), intent(in) :: map_2d(ndf_2d)
  integer(i_def), intent(in) :: ndf_scal, undf_scal
  integer(i_def), intent(in) :: map_scal(ndf_scal)
  integer(i_def), intent(in) :: ndf_sice, undf_sice
  integer(i_def), intent(in) :: map_sice(ndf_sice)
  integer(i_def), intent(in) :: ndf_w3, undf_w3
  integer(i_def), intent(in) :: map_w3(ndf_w3)
  integer(i_def), intent(in) :: ndf_wth, undf_wth
  integer(i_def), intent(in) :: map_wth(ndf_wth)

  real(r_def), intent(out) :: tile_sw_direct_albedo(undf_sw_tile)
  real(r_def), intent(out) :: tile_sw_diffuse_albedo(undf_sw_tile)
  real(r_def), intent(out) :: tile_lw_albedo(undf_lw_tile)

  real(r_def), intent(in) :: tile_fraction(undf_tile)
  real(r_def), intent(in) :: tile_temperature(undf_tile)
  real(r_def), intent(in) :: tile_snow_mass(undf_tile)
  real(r_def), intent(in) :: tile_snow_rgrain(undf_tile)
  real(r_def), intent(in) :: snow_depth(undf_tile)
  real(r_def), intent(in) :: snowpack_density(undf_tile)

  real(r_def), intent(in) :: leaf_area_index(undf_pft)
  real(r_def), intent(in) :: canopy_height(undf_pft)

  real(r_def), intent(in) :: sd_orog(undf_2d)
  real(r_def), intent(in) :: soil_albedo(undf_2d)
  real(r_def), intent(in) :: soil_roughness(undf_2d)
  real(r_def), intent(in) :: albedo_obs_vis(undf_2d)
  real(r_def), intent(in) :: albedo_obs_nir(undf_2d)
  real(r_def), intent(out):: albedo_obs_scaling(undf_scal)
  real(r_def), intent(in) :: snow_soot(undf_2d)
  real(r_def), intent(in) :: chloro_sea(undf_2d)
  real(r_def), intent(in) :: z0msea(undf_2d)
  real(r_def), intent(in) :: cos_zenith_angle(undf_2d)

  real(r_def), intent(in) :: sea_ice_thickness(undf_sice)

  real(r_def), intent(in) :: u1_in_w3(undf_w3)
  real(r_def), intent(in) :: u2_in_w3(undf_w3)
  real(r_def), intent(in) :: height_w3(undf_w3)
  real(r_def), intent(in) :: height_wth(undf_wth)

  ! Local variables for the kernel
  integer(i_def) :: i, i_tile, i_pft, i_sice, i_band, n
  integer(i_def) :: df_rtile

  ! Inputs to surf_couple_radiation
  real(r_um), dimension(row_length, rows) :: &
    tstar_sea, ws_10m_sea, chloro, ice_fract, cos_zen_rts, flandg, &
    soot
  real(r_um), dimension(row_length, rows, nice_use) :: &
    ice_fract_cat, ice_thick_cat, tstar_sice_cat, snow_sice_cat
  real(r_um), dimension(row_length, rows, nice) :: &
    pond_frac_cat, pond_depth_cat
  real(r_um), dimension(land_field) :: &
    albobs_sw, albobs_vis, albobs_nir, albsoil, sd_orog_land, z0m_soil
  real(r_um), dimension(land_field, ntype) :: &
    frac_tile
  real(r_um), dimension(land_field, ntiles) :: &
    z0_tile, rgrain, snow_tile, tstar_tile, &
    z0h_bare_tile, catch_snow_tile, catch_tile
  real(r_um), dimension(land_field, npft) :: &
    lai, canht
  integer, dimension(land_field) :: &
    land_index
  integer, dimension(ntype) :: &
    type_pts
  integer, dimension(land_field, ntype) :: &
    type_index

  ! Outputs from surf_couple_radiation
  real(r_um), dimension(row_length, rows, 4) :: &
    sea_ice_albedo, land_albedo
  real(r_um), dimension(land_field, ntiles, 4) :: &
    alb_tile
  real(r_um), dimension(row_length, rows, ntiles, 2) :: &
    albobs_sc
  real(r_um), dimension(row_length, rows, 2, n_sw_band) :: &
    open_sea_albedo


  ! ---------------------------------------------------------------------------
  ! SW tile albedos
  ! ---------------------------------------------------------------------------

  ! Land tile fractions
  flandg = 0.0_r_um
  do i = 1, n_land_tile
    flandg = flandg + real(tile_fraction(map_tile(i)), r_um)
    frac_tile(1, i) = real(tile_fraction(map_tile(i)), r_um)
  end do

  ! Jules requires fractions with respect to the land area
  if (flandg(1, 1) > 0.0_r_um) then
    land_field = 1
    land_index = 1
    frac_tile(1, 1:n_land_tile) = frac_tile(1, 1:n_land_tile) / flandg(1, 1)
  else
    land_field = 0
    land_index = 0
  end if

  if (tile_fraction(map_tile(ice)) > 0.5_r_def) then
    l_lice_point = .true.
  else
    l_lice_point = .false.
  end if

  ! Set type_pts and type_index
  call tilepts(land_field, frac_tile, type_pts, type_index)

  ! Land tile temperatures
  do i = 1, n_land_tile
    tstar_tile(1, i) = real(tile_temperature(map_tile(i)), r_um)
  end do

  ! Sea temperature
  tstar_sea = real(tile_temperature(map_tile(first_sea_tile)), r_um)

  ! Sea-ice temperatures
  i_sice = 0
  do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
    i_sice = i_sice + 1
    tstar_sice_cat(1, 1, i_sice) = real(tile_temperature(map_tile(i)), r_um)
  end do

  ! Sea-ice fraction
  i_sice = 0
  ice_fract = 0.0_r_um
  do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
    i_sice = i_sice + 1
    ice_fract = ice_fract + real(tile_fraction(map_tile(i)), r_um)
    ice_fract_cat(1, 1, i_sice) = real(tile_fraction(map_tile(i)), r_um)
  end do

  ! Jules requires sea-ice fractions with respect to the sea area
  if (ice_fract(1, 1) > 0.0_r_um) then
    ice_fract(1, 1) = ice_fract(1, 1) / (1.0_r_um - flandg(1, 1))
    ice_fract_cat(1, 1, 1:n_sea_ice_tile) &
      = ice_fract_cat(1, 1, 1:n_sea_ice_tile) / (1.0_r_um - flandg(1, 1))
  end if

  ! Combined sea and sea-ice index
  if (flandg(1, 1) < 1.0_r_um) then
    ssi_index = 1
  else
    ssi_index = 0
  end if

  ! Individual sea and sea-ice indices
  if (ssi_index(1) > 0) then
    if (ice_fract(1, 1) < 1.0_r_um) then
      sea_pts = 1
      sea_index = 1
    else
      sea_pts = 0
      sea_index = 0
    end if
  end if

  ! Multi-category sea-ice index
  do n = 1, nice_use
    if (ssi_index(1) > 0 .and. ice_fract_cat(1, 1, n) > 0.0_r_um) then
      sice_pts_ncat(n) = 1
      sice_index_ncat(1, n) = 1
      sice_frac_ncat(1, n) = ice_fract_cat(1, 1, n)
    else
      sice_pts_ncat(n) = 0
      sice_index_ncat(1, n) = 0
      sice_frac_ncat(1, n) = 0.0_r_um
    end if
  end do

  do i_sice = 1, n_sea_ice_tile
    ! Sea-ice thickness
    ice_thick_cat(1,1,i_sice) = real(sea_ice_thickness(map_sice(i_sice)),r_um)
  end do

  ! Leaf area index
  do i_pft = 1, npft
    lai(1, i_pft) = real(leaf_area_index(map_pft(i_pft)), r_um)
  end do

  ! Canopy height
  do i_pft = 1, npft
    canht(1, i_pft) = real(canopy_height(map_pft(i_pft)), r_um)
  end do

  ! Roughness length (z0_tile)
  z0m_soil = real(soil_roughness(map_2d(1)), r_um)
  call sparm(land_field, n_land_tile, type_pts, type_index, &
             frac_tile, canht, lai, z0m_soil, &
             catch_snow_tile, catch_tile, z0_tile, z0h_bare_tile)

  ! Snow-free soil albedo
  albsoil = real(soil_albedo(map_2d(1)), r_um)

  ! Cosine of the solar zenith angle
  cos_zen_rts = real(cos_zenith_angle(map_2d(1)), r_um)

  ! Standard deviation of orography
  sd_orog_land = real(sd_orog(map_2d(1)), r_um)

  ! 10m wind speed over the sea
  ws_10m_sea = sqrt(u1_in_w3(map_w3(1))**2 + u2_in_w3(map_w3(1))**2) &
    * log(10.0_r_def / z0msea(map_2d(1))) &
    / log((height_w3(map_w3(1))-height_wth(map_wth(1))) / z0msea(map_2d(1)))

  ! Chlorophyll content of the sea
  chloro = real(chloro_sea(map_2d(1)), r_um)

  ! Observed albedo
  albobs_vis = real(albedo_obs_vis(map_2d(1)), r_um)
  albobs_nir = real(albedo_obs_nir(map_2d(1)), r_um)

  ! Lying snow mass on land tiles
  do i = 1, n_land_tile
    snow_tile(1, i) = real(tile_snow_mass(map_tile(i)), r_um)
    rgrain(1, i) = real(tile_snow_rgrain(map_tile(i)), r_um)
    snowdepth_surft(1, i) = real(snow_depth(map_tile(i)), r_um)
    rho_snow_grnd_surft(1, i) = real(snowpack_density(map_tile(i)), r_um)
  end do

  ! Lying snow mass on sea ice categories
  i_sice = 0
  do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
    i_sice = i_sice + 1
    snow_sice_cat(1, 1, i_sice) = real(tile_snow_mass(map_tile(i)), r_um)
  end do

  ! Snow soot content
  soot = real(snow_soot(map_2d(1)), r_um)

  CALL surf_couple_radiation(                                   &
  ! Fluxes INTENT(IN)
    tstar_sea,                                                  &
  ! Misc INTENT(IN)
    ws_10m_sea, chloro,                                         &
    n_sw_band, n_sw_band,                                       &
    sw_wavelength_short, sw_wavelength_long,                    &
  ! Misc INTENT(OUT)
    sea_ice_albedo,                                             &
  ! Fluxes INTENT(OUT)
    alb_tile, land_albedo,                                      &
  ! UM-only args: INTENT(IN)
    pond_frac_cat, pond_depth_cat,                              &
  ! (ancil_info mod)
    ntiles, land_field, land_index, type_pts, type_index,       &
    row_length, rows, ice_fract, frac_tile,                     &
  ! (p_s_parms mod)
    cos_zen_rts, albobs_sw, albobs_vis, albobs_nir,             &
    z0_tile, albsoil,                                           &
  ! (coastal mod)
    flandg, tstar_sice_cat,                                     &
  ! (prognostics mod)
    snow_sice_cat, ice_thick_cat,                               &
    lai, canht, rgrain,                                         &
    snow_tile, soot, tstar_tile, sd_orog_land,                  &
  ! UM-only args: INTENT(OUT)
    albobs_sc, open_sea_albedo)

  df_rtile = 0
  do i_band = 1, n_sw_band
    ! Land tile albedos
    df_rtile = n_surf_tile*(i_band-1)
    do i_tile = 1, n_land_tile
      df_rtile = df_rtile + 1
      tile_sw_direct_albedo(map_sw_tile(1)+df_rtile-1) &
        = sw_weight_blue(i_band) &
        * real(alb_tile(1, i_tile, 1), r_def)  & ! visible direct albedo
        + (1.0_r_def - sw_weight_blue(i_band)) &
        * real(alb_tile(1, i_tile, 3), r_def)    ! near-ir direct albedo
      tile_sw_diffuse_albedo(map_sw_tile(1)+df_rtile-1) &
        = sw_weight_blue(i_band) &
        * real(alb_tile(1, i_tile, 2), r_def)  & ! visible diffuse albedo
        + (1.0_r_def - sw_weight_blue(i_band)) &
        * real(alb_tile(1, i_tile, 4), r_def)    ! near-ir diffuse albedo
    end do

    ! Sea tile albedos
    df_rtile = first_sea_tile-1 + n_surf_tile*(i_band-1)
    do i_tile = 1, n_sea_tile
      df_rtile = df_rtile + 1
      tile_sw_diffuse_albedo(map_sw_tile(1)+df_rtile-1) &
        = real(open_sea_albedo(1, 1, 1, i_band), r_def)
      tile_sw_direct_albedo(map_sw_tile(1)+df_rtile-1) &
        = real(open_sea_albedo(1, 1, 2, i_band), r_def)
    end do

    ! Sea-ice tile albedos
    df_rtile = first_sea_ice_tile-1 + n_surf_tile*(i_band-1)
    do i_tile = 1, n_sea_ice_tile
      df_rtile = df_rtile + 1
      tile_sw_direct_albedo(map_sw_tile(1)+df_rtile-1) &
        = sw_weight_blue(i_band) &
        * real(sea_ice_albedo(1, 1, 1), r_def) &
        + (1.0_r_def - sw_weight_blue(i_band)) &
        * real(sea_ice_albedo(1, 1, 3), r_def)
      tile_sw_diffuse_albedo(map_sw_tile(1)+df_rtile-1) &
        = sw_weight_blue(i_band) &
        * real(sea_ice_albedo(1, 1, 2), r_def) &
        + (1.0_r_def - sw_weight_blue(i_band)) &
        * real(sea_ice_albedo(1, 1, 4), r_def)
    end do
  end do

  ! ---------------------------------------------------------------------------
  ! LW tile albedos
  ! ---------------------------------------------------------------------------

  do i_band = 1, n_lw_band
    ! Land tile albedos
    df_rtile = n_surf_tile*(i_band-1)
    do i_tile = 1, n_land_tile
      df_rtile = df_rtile + 1
      if (i_tile <= npft) then
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
          = 1.0_r_def - real(emis_pft(i_tile), r_def)
      else
        tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
          = 1.0_r_def - real(emis_nvg(i_tile-npft), r_def)
      end if
    end do

    ! Sea tile albedos
    df_rtile = first_sea_tile-1 + n_surf_tile*(i_band-1)
    do i_tile = 1, n_sea_tile
      df_rtile = df_rtile + 1
      tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
        = 1.0_r_def - real(emis_sea, r_def)
    end do

    ! Sea-ice tile albedos
    df_rtile = first_sea_ice_tile-1 + n_surf_tile*(i_band-1)
    do i_tile = 1, n_sea_ice_tile
      df_rtile = df_rtile + 1
      tile_lw_albedo(map_lw_tile(1)+df_rtile-1) &
        = 1.0_r_def - real(emis_sice, r_def)
    end do
  end do

  ! Scaling factors needed for use in surface exchange code
  if (albedo_obs .and. flandg(1, 1) > 0.0_r_um) then
    df_rtile = 0
    do i_band = 1, rad_nband
      do i_tile = 1, n_land_tile
        albedo_obs_scaling(map_scal(1)+df_rtile) = &
             albobs_scaling_surft(1,i_tile,i_band)
        ! Counting from 0 so increment index here
        df_rtile = df_rtile + 1
      end do
    end do
  end if

  ! set this back to 1 before exit
  land_field = 1

end subroutine rad_tile_code

end module rad_tile_kernel_mod
