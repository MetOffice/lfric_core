!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Extracts the horizontal dof components from a W2 wind field

!> @details Extracts horizontal dof components from a 3D wind field
!!          on W2 and places them in a W2H field
!!          This kernel only works for the lowest-order elements

module sci_extract_uv_kernel_mod

use kernel_mod,                    only: kernel_type
use argument_mod,                  only: arg_type, GH_INTEGER,      &
                                         GH_FIELD, GH_REAL,         &
                                         GH_READ, GH_WRITE,         &
                                         ANY_DISCONTINUOUS_SPACE_2, &
                                         ANY_DISCONTINUOUS_SPACE_3, &
                                         CELL_COLUMN
use constants_mod,                 only: r_single, r_double, i_def
use fs_continuity_mod,             only: W2
use kernel_mod,                    only: kernel_type
use sci_face_selector_support_mod, only: face_from_face_selector

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: extract_uv_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                       &
       ! NB: This is to be used to write to a continuous W2H field, but using
       ! a discontinuous data pattern, so use discontinuous metadata
       arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
       arg_type(GH_FIELD, GH_REAL,    GH_READ,  W2),                        &
       arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
       arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)  &
       /)
  integer :: operates_on = CELL_COLUMN
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: extract_uv_code

! Generic interface for real32 and real64 types
interface extract_uv_code
  module procedure  &
    extract_uv_code_r_single, &
    extract_uv_code_r_double
end interface

contains

!> @brief Extracts the horizontal dof components from a W2 wind field
!> @param[in] nlayers Integer the number of layers
!> @param[in,out] h_wind Real array, horizontal components of wind
!> @param[in] u_wind Real array, 3d wind field
!> @param[in] face_selector_ew  2D field indicating which W/E faces
!!            to loop over in this column
!> @param[in] face_selector_ns  2D field indicating which N/S faces
!!            to loop over in this column
!> @param[in] ndf_w2h The number of degrees of freedom per cell for w2h
!> @param[in] undf_w2h The number of unique degrees of freedom for w2h
!> @param[in] map_w2h Integer array holding the dofmap for the cell at the
!!            base of the column for w2h
!> @param[in] ndf_w2 The number of degrees of freedom per cell for w2
!> @param[in] undf_w2 The number of unique degrees of freedom for w2
!> @param[in] map_w2 Integer array holding the dofmap for the cell at the
!!            base of the column for w2
!> @param[in] ndf_w3_2d  Num of DoFs for 2D W3 per cell
!> @param[in] undf_w3_2d Num of DoFs for this partition for 2D W3
!> @param[in] map_w3_2d  Map for 2D W3
subroutine extract_uv_code_r_single( nlayers,                         &
                                     h_wind,                          &
                                     u_wind,                          &
                                     face_selector_ew,                &
                                     face_selector_ns,                &
                                     ndf_w2h, undf_w2h, map_w2h,      &
                                     ndf_w2, undf_w2, map_w2,         &
                                     ndf_w3_2d, undf_w3_2d, map_w3_2d &
                                    )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: ndf_w3_2d, undf_w3_2d

  real(kind=r_single), dimension(undf_w2h),   intent(inout) :: h_wind
  real(kind=r_single), dimension(undf_w2),    intent(in)    :: u_wind
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ew
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ns
  integer(kind=i_def), dimension(ndf_w2h),    intent(in)    :: map_w2h
  integer(kind=i_def), dimension(ndf_w2),     intent(in)    :: map_w2
  integer(kind=i_def), dimension(ndf_w3_2d),  intent(in)    :: map_w3_2d

  ! Internal variables
  integer(kind=i_def) :: df, k, j

  ! Loop over horizontal W2 DoFs
  do j = 1, ABS(face_selector_ew(map_w3_2d(1))) + ABS(face_selector_ns(map_w3_2d(1)))
    df = face_from_face_selector(j, face_selector_ew(map_w3_2d(1)), face_selector_ns(map_w3_2d(1)))

    ! Loop over layers
    do k = 0, nlayers-1
      h_wind(map_w2h(df) + k) = u_wind(map_w2(df) + k)
    end do
  end do

end subroutine extract_uv_code_r_single

subroutine extract_uv_code_r_double( nlayers,                         &
                                     h_wind,                          &
                                     u_wind,                          &
                                     face_selector_ew,                &
                                     face_selector_ns,                &
                                     ndf_w2h, undf_w2h, map_w2h,      &
                                     ndf_w2, undf_w2, map_w2,         &
                                     ndf_w3_2d, undf_w3_2d, map_w3_2d &
                                    )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_w2h, undf_w2h
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: ndf_w3_2d, undf_w3_2d

  real(kind=r_double), dimension(undf_w2h),   intent(inout) :: h_wind
  real(kind=r_double), dimension(undf_w2),    intent(in)    :: u_wind
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ew
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ns
  integer(kind=i_def), dimension(ndf_w2h),    intent(in)    :: map_w2h
  integer(kind=i_def), dimension(ndf_w2),     intent(in)    :: map_w2
  integer(kind=i_def), dimension(ndf_w3_2d),  intent(in)    :: map_w3_2d

  ! Internal variables
  integer(kind=i_def) :: df, k, j

  ! Loop over horizontal W2 DoFs
  do j = 1, ABS(face_selector_ew(map_w3_2d(1))) + ABS(face_selector_ns(map_w3_2d(1)))
    df = face_from_face_selector(j, face_selector_ew(map_w3_2d(1)), face_selector_ns(map_w3_2d(1)))

    ! Loop over layers
    do k = 0, nlayers-1
      h_wind(map_w2h(df) + k) = u_wind(map_w2(df) + k)
    end do
  end do

end subroutine extract_uv_code_r_double

end module sci_extract_uv_kernel_mod
