!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

module sci_extend_chi_field_kernel_mod

use kernel_mod,            only: kernel_type
use argument_mod,          only: arg_type,                                     &
                                 GH_FIELD, GH_SCALAR,                          &
                                 GH_REAL, GH_INTEGER,                          &
                                 GH_READ, GH_WRITE,                            &
                                 ANY_DISCONTINUOUS_SPACE_3,                    &
                                 ANY_DISCONTINUOUS_SPACE_5,                    &
                                 OWNED_AND_HALO_CELL_COLUMN
use constants_mod,         only: r_def, i_def, l_def
use log_mod,               only: log_event, LOG_LEVEL_ERROR
implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: extend_chi_field_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                          &
       arg_type(GH_FIELD*3, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_5),  &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3)   &
  /)
  integer :: operates_on = OWNED_AND_HALO_CELL_COLUMN
contains
  procedure, nopass :: extend_chi_field_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: extend_chi_field_code

contains

!> @brief Extend the equiangular coordinate fields from panel boundaries into
!!        the halos
!> @param[in]     nlayers    Number of layers
!> @param[in,out] alpha_ext  Extension of the alpha (chi(1)) coordinate field
!> @param[in,out] beta_ext   Extension of the beta (chi(2)) coordinate field
!> @param[in,out] radius_ext Extension of the radius (chi(3)) coordinate field
!> @param[in]     chi1       Original alpha (chi(1)) coordinate field
!> @param[in]     chi2       Original beta (chi(2)) coordinate field
!> @param[in]     chi3       Original radius (chi(3)) coordinate field
!> @param[in]     panel_id   Indicator of the panel for each cell column
!> @param[in]     ndf_wx     Number of degrees of freedom per cell for the coordinate fields
!> @param[in]     undf_wx    Number of unique degrees of freedom for the coordinate space
!> @param[in]     map_wx     Dofmap for the cell at the base of the column for the coordinate space
!> @param[in]     ndf_pid    Number of degrees of freedom per cell for the panel_id field
!> @param[in]     undf_pid   Number of unique degrees of freedom for the panel_id field
!> @param[in]     map_pid    Dofmap for the cell at the base of the column for the panel_id field
subroutine extend_chi_field_code(nlayers,                         &
                                 alpha_ext, beta_ext, radius_ext, &
                                 chi1, chi2, chi3,                &
                                 panel_id,                        &
                                 ndf_wx, undf_wx, map_wx,         &
                                 ndf_pid, undf_pid, map_pid       &
                                )

  use coord_transform_mod, only: alphabetar2xyz, &
                                 xyz2alphabetar
  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in) :: nlayers
  integer(kind=i_def),                     intent(in) :: ndf_wx, ndf_pid, &
                                                         undf_wx, undf_pid
  integer(kind=i_def), dimension(ndf_wx),  intent(in) :: map_wx
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  real(kind=r_def), dimension(undf_wx),  intent(inout) :: alpha_ext, beta_ext
  real(kind=r_def), dimension(undf_wx),  intent(inout) :: radius_ext
  real(kind=r_def), dimension(undf_wx),  intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_pid), intent(in)    :: panel_id

  integer(kind=i_def) :: owned_panel, halo_panel, panel_edge, df, k, idx

  real(kind=r_def)                    :: x, y, z
  real(kind=r_def), parameter         :: unit_radius = 1.0_r_def
  real(kind=r_def), dimension(ndf_wx) :: alpha_halo, beta_halo
  real(kind=r_def), dimension(ndf_wx) :: alpha_owned, beta_owned
  real(kind=r_def), dimension(ndf_wx) :: alpha_extended, beta_extended

  ! This dummy variable is needed to be an output for the xyz2alphabetar
  ! subroutine. The transformation to change (alpha,beta) is only in the
  ! horizontal so the height value doesn't matter
  real(kind=r_def) :: h_dummy

  ! First, set all of the coordinates to match the standard coordinates
  do df = 1, ndf_wx
    idx = map_wx(df)
    alpha_ext(idx:idx+nlayers-1) = chi1(idx:idx+nlayers-1)
    beta_ext(idx:idx+nlayers-1)  = chi2(idx:idx+nlayers-1)
    radius_ext(idx:idx+nlayers-1) = chi3(idx:idx+nlayers-1)
  end do

  ! Assume the first entry in panel id corresponds to an owned (not halo cell)
  owned_panel = int(panel_id(1),i_def)
  ! Panel id for this column
  halo_panel = int(panel_id(map_pid(1)),i_def)

  ! Only need to extend panels if the halo is on a different panel to the
  ! owned cell
  if ( halo_panel /= owned_panel ) then
    do df = 1, ndf_wx
      alpha_halo(df) = alpha_ext(map_wx(df))
      beta_halo(df)  = beta_ext(map_wx(df))

      ! Convert (alpha,beta,h) on panel halo_panel to xyz
      ! The h value is not important as the transformation is only horizontal
      call alphabetar2xyz(alpha_halo(df), beta_halo(df), unit_radius, &
                          halo_panel, x, y, z)
      ! Now convert back to (alpha, beta, h) on owned panel
      call xyz2alphabetar(x, y, z, owned_panel, &
                          alpha_owned(df), beta_owned(df), h_dummy )
    end do

    ! alpha_halo now contains the (alpha,beta) coordinates from the halo panel
    ! alpha_owned contains the same physical points but in coordinates of the
    ! owned panel

    ! We now want to write the extended coordinates (alpha_extended) this
    ! consists of one component of alpha_owned and one component of alpha_halo
    ! e.g between panels 1 & 2 we have
    ! ------------------------------
    ! |             |              |
    ! | beta        | beta         |
    ! | ^     1     | ^     2      |
    ! | |           | |            |
    ! | --> alpha   | --> alpha    |
    ! |             |              |
    ! ------------------------------
    ! and so when extending panel 1 we set alpha_extended to be the alpha
    ! coordinate from panel 1 (alpha_extended = alpha_owned)
    ! and the beta coordinate is from the beta coordinate on panel 2
    ! (beta_extended = beta_halo)

    ! Now extend coordinates in the halo regions
    panel_edge = 10*owned_panel + halo_panel
    select case ( panel_edge )
      case(12, 21, 36, 45, 54, 63)
        ! 12: EAST edge of panel 1 to 2
        ! 21: WEST edge of panel 2 to 1
        ! 36: SOUTH edge of panel 3 to 6
        ! 45: NORTH edge of panel 4 to 5
        ! 54: WEST edge of panel 5 to 4
        ! 63: EAST edge of panel 6 to 3
        alpha_extended = alpha_owned
        beta_extended = beta_halo
      case(15, 26, 34, 43, 51, 62)
        ! 15: NORTH edge of panel 1 to 5
        ! 26: SOUTH edge of panel 2 to 6
        ! 34: WEST edge of panel 3 to 4
        ! 43: WEST edge of panel 4 to 3
        ! 51: SOUTH edge of panel 5 to 1
        ! 62: NORTH edge of panel 6 to 2
        alpha_extended = alpha_halo
        beta_extended = beta_owned
      case(14, 23, 35, 46, 52, 61)
        ! 14: WEST edge of panel 1 to 4
        ! 23: EAST edge of panel 2 to 3
        ! 35: NORTH edge of panel 3 to 5
        ! 46: SOUTH edge of panel 4 to 6
        ! 52: EAST edge of panel 5 to 2
        ! 61: WEST edge of panel 6 to 1
        beta_extended = alpha_halo
        alpha_extended = alpha_owned
      case(16, 25, 32, 41, 53, 64)
        ! 16: SOUTH edge of panel 1 to 6
        ! 25: NORTH edge of panel 2 to 5
        ! 32: EAST edge of panel 3 to 2
        ! 41: EAST edge of panel 4 to 1
        ! 53: NORTH edge of panel 5 to 3
        ! 64: SOUTH edge of panel 6 to 4
        alpha_extended = beta_halo
        beta_extended = beta_owned
      case default
        call log_event('Invalid panel edge',LOG_LEVEL_ERROR)
    end select

    ! Write back to coordinate fields keeping the same value in the whole column
    do df = 1, ndf_wx
      do k = 0, nlayers - 1
        alpha_ext(map_wx(df)+k) = alpha_extended(df)
        beta_ext(map_wx(df)+k)  = beta_extended(df)
      end do
    end do
  end if

end subroutine extend_chi_field_code

end module sci_extend_chi_field_kernel_mod
