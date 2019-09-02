!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to Socrates for longwave (thermal) fluxes

module lw_kernel_mod

use argument_mod,      only : arg_type, func_type,                       &
                              GH_FIELD, GH_INTEGER, GH_READ, GH_WRITE,   &
                              GH_READWRITE, GH_INC, CELLS,               &
                              ANY_SPACE_1, ANY_SPACE_2, ANY_SPACE_3
use fs_continuity_mod, only:  W3, Wtheta
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type

implicit none
private
public :: lw_kernel_type
public :: lw_code

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the Psy layer.
type, extends(kernel_type) :: lw_kernel_type
  private
  type(arg_type) :: meta_args(23) = (/               &
    arg_type(GH_FIELD,   GH_WRITE,     Wtheta),      & ! lw_heating_rate
    arg_type(GH_FIELD,   GH_WRITE,     ANY_SPACE_1), & ! lw_down_surf
    arg_type(GH_FIELD,   GH_WRITE,     ANY_SPACE_2), & ! lw_up_tile
    arg_type(GH_FIELD,   GH_READWRITE, Wtheta),      & ! lw_heating_rate_rts
    arg_type(GH_FIELD,   GH_INC,       ANY_SPACE_1), & ! lw_down_surf_rts
    arg_type(GH_FIELD,   GH_INC,       ANY_SPACE_2), & ! lw_up_tile_rts
    arg_type(GH_FIELD,   GH_READ,      Wtheta),      & ! theta
    arg_type(GH_FIELD,   GH_READ,      W3),          & ! theta_in_w3
    arg_type(GH_FIELD,   GH_READ,      W3),          & ! exner
    arg_type(GH_FIELD,   GH_READ,      Wtheta),      & ! exner_in_wth
    arg_type(GH_FIELD,   GH_READ,      Wtheta),      & ! rho_in_wth
    arg_type(GH_FIELD,   GH_READ,      W3),          & ! height_w3
    arg_type(GH_FIELD,   GH_READ,      Wtheta),      & ! height_wth
    arg_type(GH_FIELD,   GH_READ,      Wtheta),      & ! mv
    arg_type(GH_FIELD,   GH_READ,      Wtheta),      & ! mcl
    arg_type(GH_FIELD,   GH_READ,      Wtheta),      & ! mci
    arg_type(GH_FIELD,   GH_READ,      Wtheta),      & ! area_fraction
    arg_type(GH_FIELD,   GH_READ,      Wtheta),      & ! liquid_fraction
    arg_type(GH_FIELD,   GH_READ,      Wtheta),      & ! ice_fraction
    arg_type(GH_FIELD,   GH_READ,      ANY_SPACE_2), & ! tile_fraction
    arg_type(GH_FIELD,   GH_READ,      ANY_SPACE_2), & ! tile_temperature
    arg_type(GH_FIELD,   GH_READ,      ANY_SPACE_3), & ! tile_lw_albedo
    arg_type(GH_INTEGER, GH_READ                  )  & ! timestep
    /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: lw_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface lw_kernel_type
  module procedure lw_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

function lw_kernel_constructor() result(self)
  implicit none
  type(lw_kernel_type) :: self
  return
end function lw_kernel_constructor

!> @param[in]     nlayers                 Number of layers
!> @param[out]    lw_heating_rate         LW heating rate
!> @param[out]    lw_down_surf            LW downward flux at the surface
!> @param[out]    lw_up_tile              LW upward tiled surface flux
!> @param[in,out] lw_heating_rate_rts     LW heating rate
!> @param[in,out] lw_down_surf_rts        LW downward flux at the surface
!> @param[in,out] lw_up_tile_rts          LW upward tiled surface flux
!> @param[in]     theta                   Potential temperature
!> @param[in]     theta_in_w3             Potential temperature in density space
!> @param[in]     exner                   exner pressure in density space
!> @param[in]     exner_in_wth            exner pressure in wth space
!> @param[in]     rho_in_wth              density in potential temperature space
!> @param[in]     height_w3               Height of w3 levels above surface
!> @param[in]     height_wth              Height of wth levels above surface
!> @param[in]     mv                      Water vapour field
!> @param[in]     mcl                     Cloud liquid field
!> @param[in]     mci                     Cloud ice field
!> @param[in]     area_fraction           Total cloud area fraction field
!> @param[in]     liquid_fraction         Liquid cloud fraction field
!> @param[in]     ice_fraction            Ice cloud fraction field
!> @param[in]     tile_fraction           Surface tile fractions
!> @param[in]     tile_temperature        Surface tile temperature
!> @param[in]     tile_lw_albedo          LW tile albedos
!> @param[in]     timestep                Timestep number
!> @param[in]     ndf_wth                 No. DOFs per cell for wth space
!> @param[in]     undf_wth                No. unique of DOFs for wth space
!> @param[in]     map_wth                 Dofmap for wth space column base cell
!> @param[in]     ndf_2d                  No. of DOFs per cell for 2d space
!> @param[in]     undf_2d                 No. unique of DOFs for 2d space
!> @param[in]     map_2d                  Dofmap for 2d space column base cell
!> @param[in]     ndf_tile                Number of DOFs per cell for tiles
!> @param[in]     undf_tile               Number of total DOFs for tiles
!> @param[in]     map_tile                Dofmap for tile space column base cell
!> @param[in]     ndf_w3                  No. of DOFs per cell for w3 space
!> @param[in]     undf_w3                 No. unique of DOFs for w3 space
!> @param[in]     map_w3                  Dofmap for w3 space column base cell
!> @param[in]     ndf_rtile               No. of DOFs per cell for rtile space
!> @param[in]     undf_rtile              No. unique of DOFs for rtile space
!> @param[in]     map_rtile               Dofmap for rtile space column base cell
subroutine lw_code(nlayers,                          & 
                   lw_heating_rate,                  &
                   lw_down_surf,                     &
                   lw_up_tile,                       &
                   lw_heating_rate_rts,              &
                   lw_down_surf_rts,                 &
                   lw_up_tile_rts,                   &
                   theta,                            &
                   theta_in_w3,                      &
                   exner,                            &
                   exner_in_wth,                     &
                   rho_in_wth,                       &
                   height_w3,                        &
                   height_wth,                       &
                   mv,                               &
                   mcl,                              &
                   mci,                              &
                   area_fraction,                    &
                   liquid_fraction,                  &
                   ice_fraction,                     &
                   tile_fraction,                    &
                   tile_temperature,                 &
                   tile_lw_albedo,                   &
                   timestep,                         &
                   ndf_wth, undf_wth, map_wth,       &
                   ndf_2d, undf_2d, map_2d,          &
                   ndf_tile, undf_tile, map_tile,    &
                   ndf_w3, undf_w3, map_w3,          &
                   ndf_rtile, undf_rtile, map_rtile)

  use well_mixed_gases_config_mod, only:         &
    co2_mix_ratio, n2o_mix_ratio, ch4_mix_ratio, &
    cfc11_mix_ratio, cfc12_mix_ratio,            &
    cfc113_mix_ratio, hcfc22_mix_ratio, hfc134a_mix_ratio
  use radiation_config_mod, only: n_radstep,  &
    l_planet_grey_surface, planet_emissivity, &
    i_cloud_ice_type_lw, i_cloud_liq_type_lw, &
    cloud_representation, cloud_overlap,      &
    cloud_inhomogeneity,                      &
    cloud_representation_no_cloud,            &
    cloud_representation_liquid_and_ice,      &
    cloud_representation_conv_strat_liq_ice,  &
    cloud_overlap_maximum_random,             &
    cloud_overlap_random,                     &
    cloud_overlap_exponential_random,         &
    cloud_inhomogeneity_homogeneous,          &
    cloud_inhomogeneity_scaling,              &
    cloud_inhomogeneity_mcica,                &
    cloud_inhomogeneity_cairns
  use set_thermodynamic_mod, only: set_thermodynamic
  use set_cloud_top_mod, only: set_cloud_top
  use init_jules_alg_mod, only: n_surf_tile, first_sea_tile
  use socrates_init_mod, only: n_lw_band
  use socrates_runes, only: runes, ip_source_thermal,                       &
    ip_cloud_representation_off, ip_cloud_representation_ice_water,         &
    ip_cloud_representation_csiw, ip_overlap_max_random, ip_overlap_random, &
    ip_overlap_exponential_random, ip_inhom_homogeneous, ip_inhom_scaling,  &
    ip_inhom_mcica, ip_inhom_cairns
  use socrates_bones, only: bones

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, timestep
  integer(i_def), intent(in) :: ndf_wth, ndf_w3, ndf_2d
  integer(i_def), intent(in) :: ndf_tile, ndf_rtile
  integer(i_def), intent(in) :: undf_wth, undf_w3, undf_2d
  integer(i_def), intent(in) :: undf_tile, undf_rtile

  integer(i_def), dimension(ndf_wth),   intent(in) :: map_wth
  integer(i_def), dimension(ndf_w3),    intent(in) :: map_w3
  integer(i_def), dimension(ndf_2d),    intent(in) :: map_2d
  integer(i_def), dimension(ndf_tile),  intent(in) :: map_tile
  integer(i_def), dimension(ndf_rtile), intent(in) :: map_rtile

  real(r_def), dimension(undf_wth),  intent(out) :: lw_heating_rate
  real(r_def), dimension(undf_2d),   intent(out) :: lw_down_surf
  real(r_def), dimension(undf_tile), intent(out) :: lw_up_tile

  real(r_def), dimension(undf_wth),  intent(inout) :: lw_heating_rate_rts
  real(r_def), dimension(undf_2d),   intent(inout) :: lw_down_surf_rts
  real(r_def), dimension(undf_tile), intent(inout) :: lw_up_tile_rts

  real(r_def), dimension(undf_w3),  intent(in) :: theta_in_w3, exner, height_w3
  real(r_def), dimension(undf_wth), intent(in) :: theta, exner_in_wth, &
    rho_in_wth, height_wth, mv, mcl, mci, &
    area_fraction, liquid_fraction, ice_fraction
  real(r_def), dimension(undf_tile),  intent(in) :: tile_fraction
  real(r_def), dimension(undf_tile),  intent(in) :: tile_temperature
  real(r_def), dimension(undf_rtile), intent(in) :: tile_lw_albedo

  ! Local variables for the kernel
  integer, parameter :: n_profile = 1
  integer :: i_cloud_representation, i_overlap, i_inhom
  integer(i_def) :: k, n_cloud_layer
  integer(i_def) :: df_rtile, i_tile, i_band
  integer(i_def) :: wth_1, wth_nlayers
  real(r_def), dimension(nlayers) :: &
    ! Heat capacity for each layer
    layer_heat_capacity,             &
    ! Layer pressure and temperature
    p_layer, t_layer,                &
    ! Mass of layer per square metre
    d_mass,                          &
    ! Effective radius of droplets
    liq_dim
  real(r_def), dimension(0:nlayers) :: t_layer_boundaries
  real(r_def), dimension(0:nlayers) :: lw_down, lw_up

  ! Tiled surface fields
  real(r_def) :: frac_tile(n_profile, n_surf_tile)
  real(r_def) :: t_tile(n_profile, n_surf_tile)
  real(r_def) :: albedo_diff_tile(n_profile, n_surf_tile, n_lw_band)
  real(r_def) :: flux_up_tile_rts(n_surf_tile)
  real(r_def) :: flux_up_tile(n_surf_tile)


  ! Set wth indexing
  wth_1 = map_wth(1)+1
  wth_nlayers = map_wth(1)+nlayers

  ! Tile fractions & temperatures
  do i_tile = 1, n_surf_tile
    frac_tile(1, i_tile) = tile_fraction(map_tile(i_tile))
    t_tile(1, i_tile) = tile_temperature(map_tile(i_tile))
  end do

  ! Tile albedos
  df_rtile = 0
  do i_band = 1, n_lw_band
    do i_tile = 1, n_surf_tile
      df_rtile = df_rtile + 1
      albedo_diff_tile(1, i_tile, i_band) &
        = tile_lw_albedo(map_rtile(df_rtile))
    end do
  end do

  if (mod(timestep-1_i_def, n_radstep) == 0) then
    ! Radiation time-step: full calculation of radiative fluxes

    ! Properties of clouds
    select case (cloud_representation)
    case (cloud_representation_no_cloud)
      i_cloud_representation = ip_cloud_representation_off
    case (cloud_representation_liquid_and_ice)
      i_cloud_representation = ip_cloud_representation_ice_water
    case (cloud_representation_conv_strat_liq_ice)
      i_cloud_representation = ip_cloud_representation_csiw
    case default
      i_cloud_representation = ip_cloud_representation_off
    end select
    select case (cloud_overlap)
    case (cloud_overlap_maximum_random)
      i_overlap = ip_overlap_max_random
    case (cloud_overlap_random)
      i_overlap = ip_overlap_random
    case (cloud_overlap_exponential_random)
      i_overlap = ip_overlap_exponential_random
    case default
      i_overlap = ip_overlap_max_random
    end select
    select case (cloud_inhomogeneity)
    case (cloud_inhomogeneity_homogeneous)
      i_inhom = ip_inhom_homogeneous
    case (cloud_inhomogeneity_scaling)
      i_inhom = ip_inhom_scaling
    case (cloud_inhomogeneity_mcica)
      i_inhom = ip_inhom_mcica
    case (cloud_inhomogeneity_cairns)
      i_inhom = ip_inhom_cairns
    case default
      i_inhom = ip_inhom_homogeneous
    end select

    ! Hardwire droplet effective radius for now
    do k=1, nlayers
      liq_dim(k) = 7.0e-6_r_def
    end do

    ! Use the highest cloud layer for the number of cloud layers
    call set_cloud_top(nlayers, &
      area_fraction(map_wth(1):map_wth(1)+nlayers), n_cloud_layer)

    ! Set up pressures, temperatures, masses and heat capacities
    call set_thermodynamic(nlayers,                &
      exner(map_w3(1):map_w3(1)+nlayers-1),        &
      exner_in_wth(map_wth(1):map_wth(1)+nlayers), &
      theta(map_wth(1):map_wth(1)+nlayers),        &
      rho_in_wth(map_wth(1):map_wth(1)+nlayers),   &
      height_w3(map_w3(1):map_w3(1)+nlayers-1),    &
      height_wth(map_wth(1):map_wth(1)+nlayers),   &
      p_layer, t_layer, d_mass, layer_heat_capacity)

    ! Calculate temperature at layer boundaries
    t_layer_boundaries(0) = theta(map_wth(1)) * exner_in_wth(map_wth(1))
    do k=1,nlayers-1
      t_layer_boundaries(k) = theta_in_w3(map_w3(1)+k) * exner(map_w3(1)+k)
    end do
    t_layer_boundaries(nlayers) = t_layer(nlayers)

    ! Calculate the LW fluxes (RUN the Edwards-Slingo two-stream solver)
    call runes(n_profile, nlayers,                                             &
      spectrum_name          = 'lw',                                           &
      i_source               = ip_source_thermal,                              &
      n_cloud_layer          = n_cloud_layer,                                  &
      p_layer_1d             = p_layer,                                        &
      t_layer_1d             = t_layer,                                        &
      mass_1d                = d_mass,                                         &
      density_1d             = rho_in_wth(wth_1:wth_nlayers),                  &
      t_level_1d             = t_layer_boundaries,                             &
      h2o_1d                 = mv(wth_1:wth_nlayers),                          &
      co2_mix_ratio          = co2_mix_ratio,                                  &
      n2o_mix_ratio          = n2o_mix_ratio,                                  &
      ch4_mix_ratio          = ch4_mix_ratio,                                  &
      cfc11_mix_ratio        = cfc11_mix_ratio,                                &
      cfc12_mix_ratio        = cfc12_mix_ratio,                                &
      cfc113_mix_ratio       = cfc113_mix_ratio,                               &
      hcfc22_mix_ratio       = hcfc22_mix_ratio,                               &
      hfc134a_mix_ratio      = hfc134a_mix_ratio,                              &
      t_ground               = t_tile(1, first_sea_tile),                      &
      l_grey_albedo          = l_planet_grey_surface,                          &
      grey_albedo            = 1.0_r_def - planet_emissivity,                  &
      l_tile                 = .not.l_planet_grey_surface,                     &
      n_tile                 = n_surf_tile,                                    &
      frac_tile              = frac_tile,                                      &
      t_tile                 = t_tile,                                         &
      albedo_diff_tile       = albedo_diff_tile,                               &
      cloud_frac_1d          = area_fraction(wth_1:wth_nlayers),               &
      liq_frac_1d            = liquid_fraction(wth_1:wth_nlayers),             &
      ice_frac_1d            = ice_fraction(wth_1:wth_nlayers),                &
      liq_mmr_1d             = mcl(wth_1:wth_nlayers),                         &
      ice_mmr_1d             = mci(wth_1:wth_nlayers),                         &
      liq_dim_1d             = liq_dim,                                        &
      layer_heat_capacity_1d = layer_heat_capacity,                            &
      l_mixing_ratio         = .true.,                                         &
      i_cloud_representation = i_cloud_representation,                         &
      i_overlap              = i_overlap,                                      &
      i_inhom                = i_inhom,                                        &
      i_st_water             = i_cloud_liq_type_lw,                            &
      i_st_ice               = i_cloud_ice_type_lw,                            &
      l_invert               = .true.,                                         &
      flux_down_1d           = lw_down,                                        &
      flux_up_1d             = lw_up,                                          &
      flux_up_tile_1d        = flux_up_tile_rts,                               &
      heating_rate_1d        = lw_heating_rate_rts(wth_1:wth_nlayers))

    ! Copy lowest level to surface
    lw_heating_rate_rts(map_wth(1)) = lw_heating_rate_rts(map_wth(1) + 1)

    ! Set surface flux
    lw_down_surf_rts(map_2d(1)) = lw_down(0)

    ! Set tiled fluxes
    do i_tile = 1, n_surf_tile
      lw_up_tile_rts(map_tile(i_tile)) = flux_up_tile_rts(i_tile)
    end do

    ! No corrections needed for model timestep
    lw_heating_rate(map_wth(1):wth_nlayers) &
      = lw_heating_rate_rts(map_wth(1):wth_nlayers)

    lw_down_surf(map_2d(1)) = lw_down_surf_rts(map_2d(1))

    do i_tile = 1, n_surf_tile
      lw_up_tile(map_tile(i_tile)) = lw_up_tile_rts(map_tile(i_tile))
    end do

  else
    ! Not a radiation time-step: apply corrections to the model timestep

    ! The bare "bones" of a simple radiative transfer calculation.
    ! Update the upward fluxes for the change in surface temperature.
    call bones(n_profile, nlayers,                                             &
      n_tile                 = n_surf_tile,                                    &
      l_grey_emis_correction = .true.,                                         &
      grey_albedo_tile       = albedo_diff_tile(:, :, 1),                      &
      t_tile                 = t_tile,                                         &
      heating_rate_1d_rts    = lw_heating_rate_rts(wth_1:wth_nlayers),         &
      flux_down_surf_rts     = lw_down_surf_rts(map_2d(1):map_2d(1)),          &
      heating_rate_1d_mts    = lw_heating_rate(wth_1:wth_nlayers),             &
      flux_up_tile_1d_mts    = flux_up_tile,                                   &
      flux_down_surf_mts     = lw_down_surf(map_2d(1):map_2d(1)))

    ! Copy lowest level to surface
    lw_heating_rate(map_wth(1)) = lw_heating_rate(map_wth(1) + 1)

    ! Set tiled fluxes
    do i_tile = 1, n_surf_tile
      lw_up_tile(map_tile(i_tile)) = flux_up_tile(i_tile)
    end do

  end if

end subroutine lw_code

end module lw_kernel_mod
