!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the difference between cell centre layer heights,
!!        with the result evaluated at cell face layer heights.
module calc_dz_face_kernel_mod

use argument_mod,  only: arg_type,                              &
                         GH_FIELD, GH_SCALAR, GH_READ, GH_INC,  &
                         CELL_COLUMN, ANY_SPACE_1, ANY_SPACE_2, &
                         GH_INTEGER, GH_REAL, GH_LOGICAL
use constants_mod, only: r_def, i_def
use kernel_mod,    only: kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: calc_dz_face_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                         &
       arg_type(GH_FIELD,  GH_REAL,    GH_INC,  ANY_SPACE_1), &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ, ANY_SPACE_2), &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ, ANY_SPACE_1), &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),              &
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)               &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: calc_dz_face_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: calc_dz_face_code

contains

!> @details For lowest-order fields only, calculates the finite-difference
!!          vertical distance between adjacent cell centres. As such,
!!          the result is evaluated at cell top/bottom face. The code
!!          works for any horizontally located position within the cell
!!          (centre or horizontal faces), by use of the input variable
!!          "df_to_do" which describes how many horizontal locations to
!!          calculate dz for.
!> @param[in]     nlayers     Number of layers
!> @param[in,out] dz          Distance between centres (evaluated at cell face)
!> @param[in]     height      Height field of cell centres
!> @param[in]     height_face Height field of cell faces
!> @param[in]     df_to_do    Number of dofs to calculate dz for
!> @param[in]     include_surface Flag for whether to include level 0 in level 1
!> @param[in]     ndf_dz      Number of degrees of freedom per cell for dz
!> @param[in]     undf_dz     Number of unique degrees of freedom for dz
!> @param[in]     map_dz      Dofmap for the cell at the base of the column for dz
!> @param[in]     ndf_height  Number of degrees of freedom per cell for height
!> @param[in]     undf_height Number of unique degrees of freedom for height
!> @param[in]     map_height  Dofmap for the cell at the base of the column for height
subroutine calc_dz_face_code(nlayers, dz, height, height_face, df_to_do, &
                             include_surface, ndf_dz, undf_dz, map_dz,   &
                             ndf_height, undf_height, map_height)
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, df_to_do
  integer(kind=i_def), intent(in) :: ndf_height, ndf_dz, undf_height, undf_dz

  integer(kind=i_def), dimension(ndf_dz), intent(in)     :: map_dz
  integer(kind=i_def), dimension(ndf_height), intent(in) :: map_height

  real(kind=r_def), dimension(undf_dz), intent(inout)  :: dz
  real(kind=r_def), dimension(undf_height), intent(in) :: height
  real(kind=r_def), dimension(undf_dz), intent(in)     :: height_face

  logical, intent(in) :: include_surface

  ! Internal variables
  integer(kind=i_def) :: df, k

  do df = 1, df_to_do
    ! Lowest level is distance between surface and first cell centre
    dz(map_dz(df)) = height(map_height(df))-height_face(map_dz(df))
    ! Difference between cell centres
    do k = 1, nlayers-1
      dz(map_dz(df)+k) = height(map_height(df)+k) &
                       - height(map_height(df)+k-1)
    end do

    if (include_surface) then
      ! Include the dz(0) value in dz(1)
      dz(map_dz(df)+1) = dz(map_dz(df)+1) + dz(map_dz(df))
    end if

    ! Top-most level is double the distance between top and last cell centre
    dz(map_dz(df)+nlayers) = 2.0_r_def * ( height_face(map_dz(df)+nlayers) - &
                                           height(map_height(df)+nlayers-1) )
  end do

end subroutine calc_dz_face_code

end module calc_dz_face_kernel_mod
