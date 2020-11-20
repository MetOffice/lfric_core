!-------------------------------------------------------------------------------
!(c) Crown copyright 2020 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Process Jules surface ancillaries
!> @details Kerenel used without propoer psyclone support of multi-data fields
!>          see https://github.com/stfc/PSyclone/issues/868
module process_surface_kernel_mod

  use argument_mod,  only: arg_type, GH_FIELD, GH_WRITE, CELLS, &
       ANY_DISCONTINUOUS_SPACE_1, ANY_DISCONTINUOUS_SPACE_2,    &
       ANY_DISCONTINUOUS_SPACE_3, ANY_DISCONTINUOUS_SPACE_4,    &
       GH_READ, GH_READWRITE
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  use jules_control_init_mod, only: n_land_tile, n_sea_ice_tile, &
       first_sea_tile, first_sea_ice_tile, n_surf_tile

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: process_surface_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                               &
         arg_type(GH_FIELD, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_READ,      ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_FIELD, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_4)  &
         /)
    integer :: iterates_over = CELLS

  contains
    procedure, nopass :: process_surface_code
  end type process_surface_kernel_type

  public process_surface_code
contains

  !> @param[in]     nlayers            The number of layers
  !> @param[in,out] land_area_fraction Fraction of land in grid-box
  !> @param[in]     land_tile_fraction Land tile fractions from ancil
  !> @param[in]     sea_ice_fraction   Fraction of sea-ice in grid-box
  !> @param[out]    tile_fraction      Surface tile fractions
  !> @param[in]     ndf_2d             Number of DOFs per cell for 2d fields
  !> @param[in]     undf_2d            Number of total DOFs for 2d fields
  !> @param[in]     map_2d             Dofmap for cell for surface 2d fields
  !> @param[in]     ndf_land           Number of DOFs per cell for land
  !> @param[in]     undf_land          Number of total DOFs for land
  !> @param[in]     map_land           Dofmap for cell for surface land
  !> @param[in]     ndf_sice           Number of DOFs per cell for sea ice
  !> @param[in]     undf_sice          Number of total DOFs for sea ice
  !> @param[in]     map_sice           Dofmap for cell for surface sea ice
  !> @param[in]     ndf_tile           Number of DOFs per cell for tiles
  !> @param[in]     undf_tile          Number of total DOFs for tiles
  !> @param[in]     map_tile           Dofmap for cell for surface tiles
  subroutine process_surface_code(nlayers,                       &
                                  land_area_fraction,            &
                                  land_tile_fraction,            &
                                  sea_ice_fraction,              &
                                  tile_fraction,                 &
                                  ndf_2d, undf_2d, map_2d,       &
                                  ndf_land, undf_land, map_land, &
                                  ndf_sice, undf_sice, map_sice, &
                                  ndf_tile, undf_tile, map_tile)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile)
    integer(kind=i_def), intent(in) :: ndf_land, undf_land
    integer(kind=i_def), intent(in) :: map_land(ndf_land)
    integer(kind=i_def), intent(in) :: ndf_sice, undf_sice
    integer(kind=i_def), intent(in) :: map_sice(ndf_sice)
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)

    real(kind=r_def), intent(inout)  :: land_area_fraction(undf_2d)
    real(kind=r_def), intent(in)     :: land_tile_fraction(undf_land)
    real(kind=r_def), intent(in)     :: sea_ice_fraction(undf_sice)
    real(kind=r_def), intent(out)    :: tile_fraction(undf_tile)

    !Internal variables
    integer(kind=i_def) :: i, i_sice
    real(kind=r_def) :: tot_ice, tot_land

    ! Set land tile fractions from ancillary
    tot_land = sum(land_tile_fraction(map_land(1):map_land(1)+n_land_tile-1))
    do i = 1, n_land_tile
      tile_fraction(map_tile(i)) = land_tile_fraction(map_land(1)+i-1) &
                                 * land_area_fraction(map_2d(1))       &
                                 / tot_land
    end do

    ! Set the sea ice fraction from an ancillary
    tot_ice = 0.0_r_def
    i_sice = 0
    do i = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = i_sice + 1
      ! Only use where field contains valid data
      if (sea_ice_fraction(map_sice(i_sice)) > 0.1_r_def .and. &
          land_area_fraction(map_2d(1)) < 1.0_r_def) then
        tile_fraction(map_tile(i)) = sea_ice_fraction(map_sice(i_sice)) &
                                   * (1.0_r_def - land_area_fraction(map_2d(1)))
      else
        tile_fraction(map_tile(i)) = 0.0_r_def
      end if
      tot_ice = tot_ice + tile_fraction(map_tile(i))
    end do

    ! Now set the sea fraction
    tile_fraction(map_tile(first_sea_tile)) = max(1.0_r_def &
                                            - land_area_fraction(map_2d(1)) &
                                            - tot_ice, 0.0_r_def)

    ! Calculate the total fraction as a test
    land_area_fraction(map_2d(1)) = sum(tile_fraction(map_tile(1:n_surf_tile)))

  end subroutine process_surface_code

end module process_surface_kernel_mod
