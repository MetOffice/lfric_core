!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Interface to Socrates for shortwave fluxes (external illumination)

module sw_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_SCALAR,       &
                              GH_REAL, GH_INTEGER,       &
                              GH_READ, GH_WRITE,         &
                              GH_READWRITE, CELL_COLUMN, &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_3, &
                              ANY_DISCONTINUOUS_SPACE_4, &
                              ANY_DISCONTINUOUS_SPACE_5, &
                              ANY_DISCONTINUOUS_SPACE_6
use fs_continuity_mod, only : W3, Wtheta
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type

implicit none

private

public :: sw_kernel_type
public :: sw_code

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, extends(kernel_type) :: sw_kernel_type
  private
  type(arg_type) :: meta_args(66) = (/ &
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     Wtheta),                    & ! sw_heating_rate
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_surf
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_surf
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_blue_surf
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_blue_surf
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_surf
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_toa
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_toa
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2), & ! sw_up_tile
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2), & ! sw_up_blue_tile
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, Wtheta),                    & ! sw_heating_rate_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_blue_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_blue_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_toa_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! sw_direct_toa_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! sw_up_tile_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! sw_up_blue_tile_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_6), & ! sw_direct_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_6), & ! sw_down_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_6), & ! sw_up_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! theta
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                        & ! exner
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! exner_in_wth
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! rho_in_wth
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! dz_in_wth
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! cos_zenith_angle
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! lit_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! cos_zenith_angle_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! lit_fraction_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! stellar_irradiance_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! ozone
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! mv
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! mcl
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! mci
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! area_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! liquid_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! ice_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! sigma_qcw
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! cca
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! ccw
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! cloud_drop_no_conc
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! tile_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! tile_sw_direct_albedo
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! tile_sw_diffuse_albedo
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      Wtheta),                    & ! sulphuric
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_4), & ! aer_mix_ratio
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_5), & ! aer_sw_absorption
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_5), & ! aer_sw_scattering
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_5), & ! aer_sw_asymmetry
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! latitude
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! longitude
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ                                ), & ! timestep
    ! Diagnostics (section radiation__)
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! cloud_top_re_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! cloud_top_weight_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! warm_cloud_top_re_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! warm_cloud_top_weight_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     Wtheta),                    & ! sw_aer_optical_depth_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_6), & ! sw_direct_clear_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_6), & ! sw_down_clear_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_6), & ! sw_up_clear_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_down_clear_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! sw_up_clear_surf_rts
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1)  & ! sw_up_clear_toa_rts
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: sw_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

! @param[in]     nlayers                   Number of layers
! @param[in,out] sw_heating_rate           SW heating rate
! @param[in,out] sw_down_surf              SW downward surface flux
! @param[in,out] sw_direct_surf            SW unscattered surface flux
! @param[in,out] sw_down_blue_surf         SW blue downward surface flux
! @param[in,out] sw_direct_blue_surf       SW blue unscattered surface flux
! @param[in,out] sw_up_surf                SW upward surface flux
! @param[in,out] sw_up_toa                 SW upward top-of-atmosphere flux
! @param[in,out] sw_direct_toa             SW unscattered top-of-atmosphere flux
! @param[in,out] sw_up_tile                SW upward tiled surface flux
! @param[in,out] sw_up_blue_tile           SW blue upward tiled surface flux
! @param[in,out] sw_heating_rate_rts       SW heating rate
! @param[in,out] sw_down_surf_rts          SW downward surface flux
! @param[in,out] sw_direct_surf_rts        SW unscattered surface flux
! @param[in,out] sw_down_blue_surf_rts     SW blue downward surface flux
! @param[in,out] sw_direct_blue_surf_rts   SW blue unscattered surface flux
! @param[in,out] sw_up_surf_rts            SW upward surface flux
! @param[in,out] sw_up_toa_rts             SW upward top-of-atmosphere flux
! @param[in,out] sw_direct_toa_rts         SW unscattered top-of-atmosphere flux
! @param[in,out] sw_up_tile_rts            SW upward tiled surface flux
! @param[in,out] sw_up_blue_tile_rts       SW blue upward tiled surface flux
! @param[in,out] sw_direct_rts             SW direct flux on radiation levels
! @param[in,out] sw_down_rts               SW downwards flux on radiation levels
! @param[in,out] sw_up_rts                 SW upwards flux on radiation levels
! @param[in]     theta                     Potential temperature field
! @param[in]     exner                     Exner pressure in density space
! @param[in]     exner_in_wth              Exner pressure in wth space
! @param[in]     rho_in_wth                Density in wth space
! @param[in]     dz_in_wth                 Depth of wth levels
! @param[in]     cos_zenith_angle          Cosine of the stellar zenith angle
! @param[in]     lit_fraction              Lit fraction of the timestep
! @param[in]     cos_zenith_angle_rts      Cosine of the stellar zenith angle
! @param[in]     lit_fraction_rts          Lit fraction of the timestep
! @param[in]     stellar_irradiance_rts    Stellar irradaince at the planet
! @param[in]     ozone                     Ozone field
! @param[in]     mv                        Water vapour field
! @param[in]     mcl                       Cloud liquid field
! @param[in]     mci                       Cloud ice field
! @param[in]     area_fraction             Total cloud area fraction field
! @param[in]     liquid_fraction           Liquid cloud fraction field
! @param[in]     ice_fraction              Ice cloud fraction field
! @param[in]     sigma_qcw                 Fractional standard deviation of condensate
! @param[in]     cca                       Convective cloud amount (fraction)
! @param[in]     ccw                       Convective cloud water (kg/kg) (can be ice or liquid)
! @param[in]     cloud_drop_no_conc        Cloud Droplet Number Concentration
! @param[in]     tile_fraction             Surface tile fractions
! @param[in]     tile_sw_direct_albedo     SW direct tile albedos
! @param[in]     tile_sw_diffuse_albedo    SW diffuse tile albedos
! @param[in]     sulphuric                 Sulphuric acid aerosol
! @param[in]     aer_mix_ratio             MODE aerosol mixing ratios
! @param[in]     aer_sw_absorption         MODE aerosol SW absorption
! @param[in]     aer_sw_scattering         MODE aerosol SW scattering
! @param[in]     aer_sw_asymmetry          MODE aerosol SW asymmetry
! @param[in]     latitude                  Latitude field
! @param[in]     longitude                 Longitude field
! @param[in]     timestep                  Timestep number
! @param[in,out] cloud_top_re_rts          Diagnostic: Cloud top effective radius
! @param[in,out] cloud_top_weight_rts      Diagnostic: Cloud top weight
! @param[in,out] warm_cloud_top_re_rts     Diagnostic: Warm cloud top effective radius
! @param[in,out] warm_cloud_top_weight_rts Diagnostic: Warm cloud top weight
! @param[in,out] sw_aer_optical_depth_rts  Diagnostic: Aerosol optical depth in the visible
! @param[in,out] sw_direct_clear_rts       Diagnostic: Clear-sky SW direct flux on radiation levels
! @param[in,out] sw_down_clear_rts         Diagnostic: Clear-sky SW downwards flux on radiation levels
! @param[in,out] sw_up_clear_rts           Diagnostic: Clear-sky SW upwards flux on radiation levels
! @param[in,out] sw_down_clear_surf_rts    Diagnostic: Clear-sky SW downwards surface flux
! @param[in,out] sw_up_clear_surf_rts      Diagnostic: Clear-sky SW upwards surface flux
! @param[in,out] sw_up_clear_toa_rts       Diagnostic: Clear-sky SW upwards top-of-atmosphere flux
! @param[in]     ndf_wth                   No. DOFs per cell for wth space
! @param[in]     undf_wth                  No. unique of DOFs for wth space
! @param[in]     map_wth                   Dofmap for wth space column base cell
! @param[in]     ndf_2d                    No. of DOFs per cell for 2D space
! @param[in]     undf_2d                   No. unique of DOFs for 2D space
! @param[in]     map_2d                    Dofmap for 2D space column base cell
! @param[in]     ndf_tile                  Number of DOFs per cell for tiles
! @param[in]     undf_tile                 Number of total DOFs for tiles
! @param[in]     map_tile                  Dofmap for tile space column base cell
! @param[in]     ndf_flux                  No. of DOFs per cell for flux space
! @param[in]     undf_flux                 No. unique of DOFs for flux space
! @param[in]     map_flux                  Dofmap for flux space column base cell
! @param[in]     ndf_w3                    No. of DOFs per cell for w3 space
! @param[in]     undf_w3                   No. unique of DOFs for w3 space
! @param[in]     map_w3                    Dofmap for w3 space column base cell
! @param[in]     ndf_rtile                 No. of DOFs per cell for rtile space
! @param[in]     undf_rtile                No. unique of DOFs for rtile space
! @param[in]     map_rtile                 Dofmap for rtile space column base cell
! @param[in]     ndf_mode                  No. of DOFs per cell for mode space
! @param[in]     undf_mode                 No. unique of DOFs for mode space
! @param[in]     map_mode                  Dofmap for mode space column base cell
! @param[in]     ndf_rmode                 No. of DOFs per cell for rmode space
! @param[in]     undf_rmode                No. unique of DOFs for rmode space
! @param[in]     map_rmode                 Dofmap for rmode space column base cell
subroutine sw_code(nlayers,                                                    &
                   sw_heating_rate, sw_down_surf, sw_direct_surf,              &
                   sw_down_blue_surf, sw_direct_blue_surf,                     &
                   sw_up_surf, sw_up_toa, sw_direct_toa,                       &
                   sw_up_tile, sw_up_blue_tile,                                &
                   sw_heating_rate_rts, sw_down_surf_rts, sw_direct_surf_rts,  &
                   sw_down_blue_surf_rts, sw_direct_blue_surf_rts,             &
                   sw_up_surf_rts, sw_up_toa_rts, sw_direct_toa_rts,           &
                   sw_up_tile_rts, sw_up_blue_tile_rts,                        &
                   sw_direct_rts, sw_down_rts, sw_up_rts,                      &
                   theta, exner, exner_in_wth, rho_in_wth, dz_in_wth,          &
                   cos_zenith_angle, lit_fraction,                             &
                   cos_zenith_angle_rts, lit_fraction_rts,                     &
                   stellar_irradiance_rts,                                     &
                   ozone, mv, mcl, mci,                                        &
                   area_fraction, liquid_fraction, ice_fraction,               &
                   sigma_qcw, cca, ccw, cloud_drop_no_conc,                    &
                   tile_fraction,                                              &
                   tile_sw_direct_albedo, tile_sw_diffuse_albedo,              &
                   sulphuric, aer_mix_ratio,                                   &
                   aer_sw_absorption, aer_sw_scattering, aer_sw_asymmetry,     &
                   latitude, longitude, timestep,                              &
                   cloud_top_re_rts, cloud_top_weight_rts,                     &
                   warm_cloud_top_re_rts, warm_cloud_top_weight_rts,           &
                   sw_aer_optical_depth_rts,                                   &
                   sw_direct_clear_rts, sw_down_clear_rts, sw_up_clear_rts,    &
                   sw_down_clear_surf_rts, sw_up_clear_surf_rts,               &
                   sw_up_clear_toa_rts,                                        &
                   ndf_wth, undf_wth, map_wth,                                 &
                   ndf_2d, undf_2d, map_2d,                                    &
                   ndf_tile, undf_tile, map_tile,                              &
                   ndf_flux, undf_flux, map_flux,                              &
                   ndf_w3, undf_w3, map_w3,                                    &
                   ndf_rtile, undf_rtile, map_rtile,                           &
                   ndf_mode, undf_mode, map_mode,                              &
                   ndf_rmode, undf_rmode, map_rmode)

  use well_mixed_gases_config_mod, only: &
    co2_mix_ratio, n2o_mix_ratio, ch4_mix_ratio, o2_mix_ratio
  use radiation_config_mod, only: n_radstep,  &
    l_rayleigh_sw, l_trans_zen_correction,    &
    i_cloud_ice_type_sw, i_cloud_liq_type_sw, &
    cloud_vertical_decorr
  use aerosol_config_mod, only: l_radaer, sulphuric_strat_climatology
  use set_thermodynamic_mod, only: set_thermodynamic
  use set_cloud_field_mod, only: set_cloud_field
  use jules_control_init_mod, only: n_surf_tile, sw_band_tile
  use um_physics_init_mod, only: n_aer_mode, mode_dimen, sw_band_mode
  use socrates_runes, only: runes, ip_source_illuminate, StrDiag
  use socrates_bones, only: bones
  use empty_data_mod, only: empty_real_data

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, timestep
  integer(i_def), intent(in) :: ndf_wth, ndf_w3, ndf_2d, ndf_flux
  integer(i_def), intent(in) :: ndf_tile, ndf_rtile, ndf_mode, ndf_rmode
  integer(i_def), intent(in) :: undf_wth, undf_w3, undf_2d, undf_flux
  integer(i_def), intent(in) :: undf_tile, undf_rtile, undf_mode, undf_rmode

  integer(i_def), dimension(ndf_wth),   intent(in) :: map_wth
  integer(i_def), dimension(ndf_w3),    intent(in) :: map_w3
  integer(i_def), dimension(ndf_2d),    intent(in) :: map_2d
  integer(i_def), dimension(ndf_tile),  intent(in) :: map_tile
  integer(i_def), dimension(ndf_rtile), intent(in) :: map_rtile
  integer(i_def), dimension(ndf_mode),  intent(in) :: map_mode
  integer(i_def), dimension(ndf_rmode), intent(in) :: map_rmode
  integer(i_def), dimension(ndf_flux),  intent(in) :: map_flux

  real(r_def), dimension(undf_wth),  intent(inout) :: sw_heating_rate
  real(r_def), dimension(undf_2d),   intent(inout) :: sw_down_surf, &
    sw_direct_surf, sw_down_blue_surf, sw_direct_blue_surf, &
    sw_up_surf, sw_up_toa, sw_direct_toa
  real(r_def), dimension(undf_tile), intent(inout) :: sw_up_tile, &
    sw_up_blue_tile

  real(r_def), dimension(undf_wth),  intent(inout), target :: &
    sw_heating_rate_rts
  real(r_def), dimension(undf_2d),   intent(inout) :: &
    sw_down_surf_rts, sw_direct_surf_rts, &
    sw_up_surf_rts, sw_up_toa_rts, sw_direct_toa_rts
  real(r_def), dimension(undf_2d),   intent(inout), target :: &
    sw_down_blue_surf_rts, sw_direct_blue_surf_rts
  real(r_def), dimension(undf_tile), intent(inout), target :: &
    sw_up_tile_rts, sw_up_blue_tile_rts
  real(r_def), dimension(undf_flux), intent(inout), target :: &
    sw_direct_rts, sw_down_rts, sw_up_rts

  real(r_def), dimension(undf_w3),   intent(in) :: exner
  real(r_def), dimension(undf_wth),  intent(in) :: theta, exner_in_wth, &
    rho_in_wth, dz_in_wth, ozone, mv, mcl, mci, &
    area_fraction, liquid_fraction, ice_fraction, sigma_qcw, &
    cca, ccw, cloud_drop_no_conc
  real(r_def), dimension(undf_2d), intent(in) :: &
    cos_zenith_angle, lit_fraction, &
    cos_zenith_angle_rts, lit_fraction_rts, stellar_irradiance_rts
  real(r_def), dimension(undf_tile),  intent(in) :: tile_fraction
  real(r_def), dimension(undf_rtile), intent(in) :: &
    tile_sw_direct_albedo, tile_sw_diffuse_albedo
  real(r_def), dimension(undf_wth),   intent(in) :: sulphuric
  real(r_def), dimension(undf_mode),  intent(in) :: aer_mix_ratio
  real(r_def), dimension(undf_rmode), intent(in) :: &
    aer_sw_absorption, aer_sw_scattering, aer_sw_asymmetry
  real(r_def), dimension(undf_2d), intent(in) :: latitude, longitude

  ! Conditional Diagnostics
  real(r_def), pointer, dimension(:), intent(inout) :: & ! 2d
    cloud_top_re_rts, cloud_top_weight_rts, &
    warm_cloud_top_re_rts, warm_cloud_top_weight_rts, &
    sw_down_clear_surf_rts, sw_up_clear_surf_rts, sw_up_clear_toa_rts
  real(r_def), pointer, dimension(:), intent(inout) :: & ! wth
    sw_aer_optical_depth_rts
  real(r_def), pointer, dimension(:), intent(inout) :: & ! flux
    sw_direct_clear_rts, sw_down_clear_rts, sw_up_clear_rts

  ! Local variables for the kernel
  integer(i_def), parameter :: n_profile = 1
  integer(i_def) :: i_cloud_representation, i_overlap, i_inhom, i_drop_re
  integer(i_def) :: rand_seed(n_profile)
  integer(i_def) :: n_cloud_layer
  integer(i_def) :: wth_0, wth_1, wth_nlayers, w3_1, w3_nlayers
  integer(i_def) :: tile_1, tile_last, rtile_1, rtile_last
  integer(i_def) :: mode_1, mode_last, rmode_1, rmode_last
  integer(i_def) :: flux_0, flux_nlayers
  real(r_def), dimension(nlayers) :: layer_heat_capacity
    ! Heat capacity for each layer
  real(r_def), dimension(nlayers) :: p_layer, t_layer
    ! Layer pressure and temperature
  real(r_def), dimension(nlayers) :: d_mass
    ! Mass of layer per square metre
  real(r_def), dimension(nlayers) :: cloud_frac, liq_dim
    ! Large-scale cloud fields
  real(r_def), dimension(nlayers) :: conv_frac, &
    liq_conv_frac, ice_conv_frac, liq_conv_mmr, ice_conv_mmr, liq_conv_dim
    ! Convective cloud fields
  type(StrDiag) :: sw_diag


  ! Set indexing
  wth_0 = map_wth(1)
  wth_1 = map_wth(1)+1
  wth_nlayers = map_wth(1)+nlayers
  w3_1 = map_w3(1)
  w3_nlayers = map_w3(1)+nlayers-1
  tile_1 = map_tile(1)
  tile_last = map_tile(1)+n_surf_tile-1
  rtile_1 = map_rtile(1)
  rtile_last = map_rtile(1)+sw_band_tile-1
  mode_1 = map_mode(1)+1
  mode_last = map_mode(1)+(nlayers+1)*mode_dimen-1
  rmode_1 = map_rmode(1)+1
  rmode_last = map_rmode(1)+(nlayers+1)*sw_band_mode-1
  flux_0 = map_flux(1)
  flux_nlayers = map_flux(1)+nlayers

  if (mod(timestep-1_i_def, n_radstep) == 0) then
    ! Radiation time-step: full calculation of radiative fluxes

    ! Set up pressures, temperatures, masses and heat capacities
    call set_thermodynamic(nlayers, &
      exner(w3_1:w3_nlayers), exner_in_wth(wth_0:wth_nlayers), &
      theta(wth_0:wth_nlayers), rho_in_wth(wth_0:wth_nlayers), &
      dz_in_wth(wth_0:wth_nlayers), &
      p_layer, t_layer, d_mass, layer_heat_capacity)

    ! Set up cloud fields for radiation
    call set_cloud_field(nlayers, n_profile, &
      area_fraction(wth_1:wth_nlayers), &
      cca(wth_1:wth_nlayers), ccw(wth_1:wth_nlayers), t_layer, &
      latitude(map_2d(1):map_2d(1)), longitude(map_2d(1):map_2d(1)), &
      i_cloud_representation, i_overlap, i_inhom, i_drop_re, &
      rand_seed, n_cloud_layer, cloud_frac, liq_dim, conv_frac, &
      liq_conv_frac, ice_conv_frac, liq_conv_mmr, ice_conv_mmr, liq_conv_dim)

    ! Set pointers for diagnostic output (using pointer bound remapping)
    ! Always required:
    sw_diag%heating_rate(1:1,1:nlayers) => &
      sw_heating_rate_rts(wth_1:wth_nlayers)

    sw_diag%flux_direct(1:1,0:nlayers) => sw_direct_rts(flux_0:flux_nlayers)
    sw_diag%flux_down(1:1,0:nlayers) => sw_down_rts(flux_0:flux_nlayers)
    sw_diag%flux_up(1:1,0:nlayers) => sw_up_rts(flux_0:flux_nlayers)

    sw_diag%flux_up_tile(1:1,1:n_surf_tile) => sw_up_tile_rts(tile_1:tile_last)
    sw_diag%flux_up_blue_tile(1:1,1:n_surf_tile) => &
      sw_up_blue_tile_rts(tile_1:tile_last)

    sw_diag%flux_direct_blue_surf(1:1) => &
      sw_direct_blue_surf_rts(map_2d(1):map_2d(1))
    sw_diag%flux_down_blue_surf(1:1) => &
      sw_down_blue_surf_rts(map_2d(1):map_2d(1))

    ! Diagnosed on request:
    if (.not. associated(cloud_top_re_rts, empty_real_data) .and. &
        .not. associated(cloud_top_weight_rts, empty_real_data)) then
      sw_diag%cloud_top_liq_dim(1:1) => &
        cloud_top_re_rts(map_2d(1):map_2d(1))
      sw_diag%cloud_top_liq_weight(1:1) => &
        cloud_top_weight_rts(map_2d(1):map_2d(1))
    end if
    if (.not. associated(warm_cloud_top_re_rts, empty_real_data) .and. &
        .not. associated(warm_cloud_top_weight_rts, empty_real_data)) then
      sw_diag%cloud_top_warm_liq_dim(1:1) => &
        warm_cloud_top_re_rts(map_2d(1):map_2d(1))
      sw_diag%cloud_top_warm_liq_weight(1:1) => &
        warm_cloud_top_weight_rts(map_2d(1):map_2d(1))
    end if

    ! Aerosol optical depth for SW band 3 is output (505-690nm). Once the
    ! diagnostic infrastructure is in place the band(s) may be user defined.
    if (.not. associated(sw_aer_optical_depth_rts, empty_real_data)) &
      sw_diag%aerosol_optical_depth(1:1,1:nlayers,3:3) &
                      => sw_aer_optical_depth_rts(wth_1:wth_nlayers)

    if (.not. associated(sw_direct_clear_rts, empty_real_data)) &
      sw_diag%flux_direct_clear(1:1,0:nlayers) &
                      => sw_direct_clear_rts(flux_0:flux_nlayers)
    if (.not. associated(sw_down_clear_rts, empty_real_data)) &
      sw_diag%flux_down_clear(1:1,0:nlayers) &
                      => sw_down_clear_rts(flux_0:flux_nlayers)
    if (.not. associated(sw_up_clear_rts, empty_real_data)) &
      sw_diag%flux_up_clear(1:1,0:nlayers) &
                      => sw_up_clear_rts(flux_0:flux_nlayers)

    ! Calculate the SW fluxes (RUN the Edwards-Slingo two-stream solver)
    call runes(n_profile, nlayers, sw_diag,                                    &
      spectrum_name          = 'sw',                                           &
      i_source               = ip_source_illuminate,                           &
      n_cloud_layer          = n_cloud_layer,                                  &
      p_layer_1d             = p_layer,                                        &
      t_layer_1d             = t_layer,                                        &
      mass_1d                = d_mass,                                         &
      density_1d             = rho_in_wth(wth_1:wth_nlayers),                  &
      h2o_1d                 = mv(wth_1:wth_nlayers),                          &
      o3_1d                  = ozone(wth_1:wth_nlayers),                       &
      co2_mix_ratio          = co2_mix_ratio,                                  &
      n2o_mix_ratio          = n2o_mix_ratio,                                  &
      ch4_mix_ratio          = ch4_mix_ratio,                                  &
      o2_mix_ratio           = o2_mix_ratio,                                   &
      cos_zenith_angle       = cos_zenith_angle_rts(map_2d(1):map_2d(1)),      &
      solar_irrad            = stellar_irradiance_rts(map_2d(1):map_2d(1)),    &
      n_tile                 = n_surf_tile,                                    &
      frac_tile_1d           = tile_fraction(tile_1:tile_last),                &
      albedo_diff_tile_1d    = tile_sw_diffuse_albedo(rtile_1:rtile_last),     &
      albedo_dir_tile_1d     = tile_sw_direct_albedo(rtile_1:rtile_last),      &
      cloud_frac_1d          = cloud_frac,                                     &
      liq_frac_1d            = liquid_fraction(wth_1:wth_nlayers),             &
      ice_frac_1d            = ice_fraction(wth_1:wth_nlayers),                &
      liq_mmr_1d             = mcl(wth_1:wth_nlayers),                         &
      ice_mmr_1d             = mci(wth_1:wth_nlayers),                         &
      liq_dim_1d             = liq_dim,                                        &
      liq_nc_1d              = cloud_drop_no_conc(wth_1:wth_nlayers),          &
      conv_frac_1d           = conv_frac,                                      &
      liq_conv_frac_1d       = liq_conv_frac,                                  &
      ice_conv_frac_1d       = ice_conv_frac,                                  &
      liq_conv_mmr_1d        = liq_conv_mmr,                                   &
      ice_conv_mmr_1d        = ice_conv_mmr,                                   &
      liq_conv_dim_1d        = liq_conv_dim,                                   &
      liq_conv_nc_1d         = cloud_drop_no_conc(wth_1:wth_nlayers),          &
      liq_rsd_1d             = sigma_qcw(wth_1:wth_nlayers),                   &
      ice_rsd_1d             = sigma_qcw(wth_1:wth_nlayers),                   &
      cloud_vertical_decorr  = cloud_vertical_decorr,                          &
      conv_vertical_decorr   = cloud_vertical_decorr,                          &
      rand_seed              = rand_seed,                                      &
      layer_heat_capacity_1d = layer_heat_capacity,                            &
      l_rayleigh             = l_rayleigh_sw,                                  &
      l_mixing_ratio         = .true.,                                         &
      i_cloud_representation = i_cloud_representation,                         &
      i_overlap              = i_overlap,                                      &
      i_inhom                = i_inhom,                                        &
      i_drop_re              = i_drop_re,                                      &
      i_st_water             = i_cloud_liq_type_sw,                            &
      i_st_ice               = i_cloud_ice_type_sw,                            &
      i_cnv_water            = i_cloud_liq_type_sw,                            &
      i_cnv_ice              = i_cloud_ice_type_sw,                            &
      l_sulphuric            = sulphuric_strat_climatology,                    &
      sulphuric_1d           = sulphuric(wth_1:wth_nlayers),                   &
      l_aerosol_mode         = l_radaer,                                       &
      n_aer_mode             = n_aer_mode,                                     &
      n_aer_layer            = nlayers+1,                                      &
      aer_mix_ratio_1d       = aer_mix_ratio(mode_1:mode_last),                &
      aer_absorption_1d      = aer_sw_absorption(rmode_1:rmode_last),          &
      aer_scattering_1d      = aer_sw_scattering(rmode_1:rmode_last),          &
      aer_asymmetry_1d       = aer_sw_asymmetry(rmode_1:rmode_last),           &
      l_invert               = .true.)

    ! Set level 0 increment such that theta increment will equal level 1
    sw_heating_rate_rts(wth_0) = sw_heating_rate_rts(wth_1) &
                               * exner_in_wth(wth_0) / exner_in_wth(wth_1)

    ! Set surface and top-of-atmosphere fluxes
    sw_down_surf_rts(map_2d(1))   = sw_down_rts(flux_0)
    sw_up_surf_rts(map_2d(1))     = sw_up_rts(flux_0)
    sw_direct_surf_rts(map_2d(1)) = sw_direct_rts(flux_0)
    sw_up_toa_rts(map_2d(1))      = sw_up_rts(flux_nlayers)
    sw_direct_toa_rts(map_2d(1))  = sw_direct_rts(flux_nlayers)

    if (.not. associated(sw_down_clear_surf_rts, empty_real_data)) &
      sw_down_clear_surf_rts(map_2d(1)) = sw_down_clear_rts(flux_0)
    if (.not. associated(sw_up_clear_surf_rts, empty_real_data)) &
      sw_up_clear_surf_rts(map_2d(1)) = sw_up_clear_rts(flux_0)
    if (.not. associated(sw_up_clear_toa_rts, empty_real_data)) &
      sw_up_clear_toa_rts(map_2d(1)) = sw_up_clear_rts(flux_nlayers)

  end if

  if (n_radstep == 1) then
    ! Radiation timestep = model timestep
    sw_heating_rate(wth_0:wth_nlayers) = sw_heating_rate_rts(wth_0:wth_nlayers)
    sw_down_surf(map_2d(1))            = sw_down_surf_rts(map_2d(1))
    sw_up_surf(map_2d(1))              = sw_up_surf_rts(map_2d(1))
    sw_direct_surf(map_2d(1))          = sw_direct_surf_rts(map_2d(1))
    sw_up_toa(map_2d(1))               = sw_up_toa_rts(map_2d(1))
    sw_direct_toa(map_2d(1))           = sw_direct_toa_rts(map_2d(1))
    sw_down_blue_surf(map_2d(1))       = sw_down_blue_surf_rts(map_2d(1))
    sw_direct_blue_surf(map_2d(1))     = sw_direct_blue_surf_rts(map_2d(1))
    sw_up_tile(tile_1:tile_last)       = sw_up_tile_rts(tile_1:tile_last)
    sw_up_blue_tile(tile_1:tile_last)  = sw_up_blue_tile_rts(tile_1:tile_last)
  else
    ! Corrections to model timestep. The radiative fluxes have been calculated
    ! for a mean sun angle over the radiation timestep and must be converted
    ! to fluxes appropriate for the mean sun angle over the shorter model
    ! timestep.

    ! The bare "bones" of a simple radiative transfer calculation.
    ! Apply the solar zenith angle correction.
    call bones(n_profile, nlayers,                                             &
      n_tile                    = n_surf_tile,                                 &
      l_cos_zen_correction      = .true.,                                      &
      cos_zen_rts               = cos_zenith_angle_rts(map_2d(1):map_2d(1)),   &
      lit_frac_rts              = lit_fraction_rts(map_2d(1):map_2d(1)),       &
      cos_zen_mts               = cos_zenith_angle(map_2d(1):map_2d(1)),       &
      lit_frac_mts              = lit_fraction(map_2d(1):map_2d(1)),           &
      l_trans_zen_correction    = l_trans_zen_correction,                      &
      heating_rate_1d_rts       = sw_heating_rate_rts(wth_1:wth_nlayers),      &
      flux_up_tile_1d_rts       = sw_up_tile_rts(tile_1:tile_last),            &
      flux_up_blue_tile_1d_rts  = sw_up_blue_tile_rts(tile_1:tile_last),       &
      flux_direct_toa_rts       = sw_direct_toa_rts(map_2d(1):map_2d(1)),      &
      flux_up_toa_rts           = sw_up_toa_rts(map_2d(1):map_2d(1)),          &
      flux_direct_surf_rts      = sw_direct_surf_rts(map_2d(1):map_2d(1)),     &
      flux_down_surf_rts        = sw_down_surf_rts(map_2d(1):map_2d(1)),       &
      flux_up_surf_rts          = sw_up_surf_rts(map_2d(1):map_2d(1)),         &
      flux_direct_blue_surf_rts = sw_direct_blue_surf_rts(map_2d(1):map_2d(1)),&
      flux_down_blue_surf_rts   = sw_down_blue_surf_rts(map_2d(1):map_2d(1)),  &
      heating_rate_1d_mts       = sw_heating_rate(wth_1:wth_nlayers),          &
      flux_up_tile_1d_mts       = sw_up_tile(tile_1:tile_last),                &
      flux_up_blue_tile_1d_mts  = sw_up_blue_tile(tile_1:tile_last),           &
      flux_direct_toa_mts       = sw_direct_toa(map_2d(1):map_2d(1)),          &
      flux_up_toa_mts           = sw_up_toa(map_2d(1):map_2d(1)),              &
      flux_direct_surf_mts      = sw_direct_surf(map_2d(1):map_2d(1)),         &
      flux_down_surf_mts        = sw_down_surf(map_2d(1):map_2d(1)),           &
      flux_up_surf_mts          = sw_up_surf(map_2d(1):map_2d(1)),             &
      flux_direct_blue_surf_mts = sw_direct_blue_surf(map_2d(1):map_2d(1)),    &
      flux_down_blue_surf_mts   = sw_down_blue_surf(map_2d(1):map_2d(1)))

    ! Set level 0 increment such that theta increment will equal level 1
    sw_heating_rate(wth_0) = sw_heating_rate(wth_1) &
                           * exner_in_wth(wth_0) / exner_in_wth(wth_1)
  end if

end subroutine sw_code
end module sw_kernel_mod
