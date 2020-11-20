!-------------------------------------------------------------------------------
!(c) Crown copyright 2020 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Process Jules soil ancillaries
module process_soil_kernel_mod

  use argument_mod,  only: arg_type, GH_FIELD, GH_WRITE, CELLS, &
       ANY_DISCONTINUOUS_SPACE_1, ANY_DISCONTINUOUS_SPACE_2, GH_READ
  use constants_mod, only: r_def, i_def, r_um, i_um
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: process_soil_kernel_type
    private
    type(arg_type) :: meta_args(13) = (/                          &
         arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2)  &
         /)
    integer :: iterates_over = CELLS

  contains
    procedure, nopass :: process_soil_code
  end type process_soil_kernel_type

  public process_soil_code
contains

  !> @param[in]  nlayers            The number of layers
  !> @param[in]  soil_moist_sat     Volumetric soil moist at saturation
  !> @param[in]  mean_topog_index   Mean topographic index
  !> @param[in]  stdev_topog_index  Standard dev in topog index
  !> @param[out] a_sat_frac         saturated fraction fitting parameter a
  !> @param[out] c_sat_frac         saturated fraction fitting parameter c
  !> @param[out] a_wet_frac         wet fraction fitting parameter a
  !> @param[out] c_wet_frac         wet fraction fitting parameter c
  !> @param[in]  clapp_horn_b       Clapp and Hornberger b coefficient
  !> @param[in]  soil_suction_sat   Saturated soil water suction (m)
  !> @param[in]  soil_temperature       Soil temperature (K)
  !> @param[in]  soil_moisture          Soil moisture content (kg m-2)
  !> @param[out] unfrozen_soil_moisture Unfrozen soil moisture proportion
  !> @param[out] frozen_soil_moisture   Frozen soil moisture proportion
  !> @param[in]  ndf_2d             Number of DOFs per cell for 2d fields
  !> @param[in]  undf_2d            Number of total DOFs for 2d fields
  !> @param[in]  map_2d             Dofmap for cell for surface 2d fields
  !> @param[in]  ndf_soil           Number of DOFs per cell for soil levels
  !> @param[in]  undf_soil          Number of total DOFs for soil levels
  !> @param[in]  map_soil           Dofmap for cell for soil levels
  subroutine process_soil_code(nlayers,                       &
                               soil_moist_sat,                &
                               mean_topog_index,              &
                               stdev_topog_index,             &
                               a_sat_frac,                    &
                               c_sat_frac,                    &
                               a_wet_frac,                    &
                               c_wet_frac,                    &
                               clapp_horn_b,                  &
                               soil_suction_sat,              &
                               soil_temperature,              &
                               soil_moisture,                 &
                               unfrozen_soil_moisture,        &
                               frozen_soil_moisture,          &
                               ndf_2d, undf_2d, map_2d,       &
                               ndf_soil, undf_soil, map_soil)

    use ancil_info, only: soil_pts
    use calc_fsat_mod, only: calc_fsat
    use jules_soil_mod, only: dzsoil
    use nlsizes_namelist_mod, only: sm_levels

    use jules_physics_init_mod, only: decrease_sath_cond

    use freeze_soil_mod, only: freeze_soil

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)
    integer(kind=i_def), intent(in) :: ndf_soil, undf_soil
    integer(kind=i_def), intent(in) :: map_soil(ndf_soil)

    real(kind=r_def), intent(in)    :: soil_moist_sat(undf_2d)
    real(kind=r_def), intent(in)    :: mean_topog_index(undf_2d)
    real(kind=r_def), intent(in)    :: stdev_topog_index(undf_2d)
    real(kind=r_def), intent(out)   :: a_sat_frac(undf_2d)
    real(kind=r_def), intent(out)   :: c_sat_frac(undf_2d)
    real(kind=r_def), intent(out)   :: a_wet_frac(undf_2d)
    real(kind=r_def), intent(out)   :: c_wet_frac(undf_2d)

    real(kind=r_def), intent(in)  :: soil_suction_sat(undf_2d)
    real(kind=r_def), intent(in)  :: clapp_horn_b(undf_2d)
    real(kind=r_def), intent(in)  :: soil_temperature(undf_soil)
    real(kind=r_def), intent(in)  :: soil_moisture(undf_soil)
    real(kind=r_def), intent(out) :: unfrozen_soil_moisture(undf_soil)
    real(kind=r_def), intent(out) :: frozen_soil_moisture(undf_soil)

    real(r_um), dimension(1) :: fexp, ti_mean, ti_sig,  &
         gamtot, a_fsat, c_fsat, a_fwet, c_fwet, dummy

    real(r_um), dimension(sm_levels) :: sthu, sthf

    integer(i_um), dimension(1) :: soil_index

    real(r_um) :: zdepth

    integer(i_um) :: n

    ! Only process the ancils if this is a soil point
    if ( soil_moist_sat(map_2d(1)) > 0.0_r_def ) then

      soil_pts = 1
      soil_index = 1

      ! Calculate total soil depth
      zdepth = 0.0_r_um
      do n = 1, sm_levels
        zdepth = zdepth + dzsoil(n)
      end do

      ! Map LFRic variables into Jules variables
      ti_mean = real(mean_topog_index(map_2d(1)), r_um)
      ti_sig = real(stdev_topog_index(map_2d(1)), r_um)

      ! First calculate the integrated gamma function
      gamtot = 0.0_r_um
      call calc_fsat(.true., soil_pts, soil_index, soil_pts, ti_mean, &
                     ti_sig, dummy, dummy, gamtot, dummy, dummy)

      ! Then calculate the fitting parameters (a_fsat, a_fwet, c_fsat, c_fwet)
      fexp = decrease_sath_cond
      ! DEPENDS ON: calc_fit_fsat
      call calc_fit_fsat(soil_pts, soil_index, soil_pts, fexp, ti_mean, &
                         ti_sig, gamtot, zdepth, a_fsat, c_fsat,        &
                         a_fwet, c_fwet)

      ! Map the calculated values back to LFRic variables
      a_sat_frac(map_2d(1)) = real(a_fsat(1), r_def)
      c_sat_frac(map_2d(1)) = real(c_fsat(1), r_def)
      a_wet_frac(map_2d(1)) = real(a_fwet(1), r_def)
      c_wet_frac(map_2d(1)) = real(c_fwet(1), r_def)

      ! Calculate frozen and unfrozen soil fractions
      call freeze_soil(soil_pts, sm_levels, clapp_horn_b(map_2d(1)),  &
                       dzsoil, soil_suction_sat(map_2d(1)),           &
                       soil_moisture(map_soil(1:sm_levels)),          &
                       soil_temperature(map_soil(1:sm_levels)),       &
                       soil_moist_sat(map_2d(1)), sthu, sthf)

      unfrozen_soil_moisture(map_soil(1:sm_levels)) = sthu
      frozen_soil_moisture(map_soil(1:sm_levels)) = sthf

    else

      ! Set to 0 if this isn't a soil point
      a_sat_frac(map_2d(1)) = 0.0_r_def
      c_sat_frac(map_2d(1)) = 0.0_r_def
      a_wet_frac(map_2d(1)) = 0.0_r_def
      c_wet_frac(map_2d(1)) = 0.0_r_def
      unfrozen_soil_moisture(map_soil(1:sm_levels)) = 0.0_r_def
      frozen_soil_moisture(map_soil(1:sm_levels)) = 0.0_r_def

    end if

  end subroutine process_soil_code

end module process_soil_kernel_mod
