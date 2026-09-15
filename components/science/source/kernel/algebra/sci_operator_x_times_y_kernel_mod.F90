!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Kernel to compute Z = X*Y where Z, X and Y are all operators
module sci_operator_x_times_y_kernel_mod

use, intrinsic :: iso_fortran_env, only: real32, real64

use constants_mod, only: i_def
use kernel_mod,    only: kernel_type
use argument_mod,  only: arg_type, func_type,        &
                         GH_OPERATOR, GH_SCALAR,     &
                         GH_REAL, GH_READ, GH_WRITE, &
                         ANY_SPACE_1, ANY_SPACE_2,   &
                         ANY_SPACE_3,                &
                         CELL_COLUMN
implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: operator_x_times_y_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                      &
       arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, ANY_SPACE_1, ANY_SPACE_2), &
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  ANY_SPACE_1, ANY_SPACE_3), &
       arg_type(GH_OPERATOR, GH_REAL, GH_READ,  ANY_SPACE_3, ANY_SPACE_2)  &
       /)
  integer :: operates_on = CELL_COLUMN
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: operator_x_times_y_kernel_code

  ! Generic interface for real32 and real64 types
  interface operator_x_times_y_kernel_code
    module procedure  &
      operator_x_times_y_kernel_code_real32, &
      operator_x_times_y_kernel_code_real64
  end interface

contains

!> @brief Computes x*y, where x & y are operators.
!! @param[in] cell Cell number
!! @param[in] nlayers Number of layers
!! @param[in] ncell_3d_1 Number of 3D cells for x_plus_ay operator
!! @param[in,out] x_times_y Operator to compute
!! @param[in] ncell_3d_2 Number of 3D cells for x operator
!! @param[in] x First operator
!! @param[in] ncell_3d_3 Number of 3D cells for y operator
!! @param[in] y Second operator
!! @param[in] ndf1 Number of dofs per cell for space 1
!! @param[in] ndf2 Number of dofs per cell for space 2
!! @param[in] ndf3 Number of dofs per cell for space 3
subroutine operator_x_times_y_kernel_code_real32(  cell, nlayers,         &
                                                   ncell_3d_1, x_times_y, &
                                                   ncell_3d_2, x,         &
                                                   ncell_3d_3, y,         &
                                                   ndf1, ndf2, ndf3)
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: cell
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ncell_3d_1, ncell_3d_2, ncell_3d_3
  integer(kind=i_def), intent(in) :: ndf1, ndf2, ndf3

  real(kind=real32),   dimension(ncell_3d_1, ndf1, ndf2), intent(inout) :: x_times_y
  real(kind=real32),   dimension(ncell_3d_2, ndf1, ndf3), intent(in)    :: x
  real(kind=real32),   dimension(ncell_3d_3, ndf3, ndf2), intent(in)    :: y

  ! Internal variables
  integer(kind=i_def) :: k, ik

  do k = 0, nlayers - 1
    ik = k + 1 + (cell-1)*nlayers
    x_times_y(ik,:,:) = matmul(x(ik,:,:), y(ik,:,:))
  end do

end subroutine operator_x_times_y_kernel_code_real32

subroutine operator_x_times_y_kernel_code_real64(  cell, nlayers,         &
                                                   ncell_3d_1, x_times_y, &
                                                   ncell_3d_2, x,         &
                                                   ncell_3d_3, y,         &
                                                   ndf1, ndf2, ndf3)
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: cell
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ncell_3d_1, ncell_3d_2, ncell_3d_3
  integer(kind=i_def), intent(in) :: ndf1, ndf2, ndf3

  real(kind=real64),   dimension(ncell_3d_1, ndf1, ndf2), intent(inout) :: x_times_y
  real(kind=real64),   dimension(ncell_3d_2, ndf1, ndf3), intent(in)    :: x
  real(kind=real64),   dimension(ncell_3d_3, ndf3, ndf2), intent(in)    :: y

  ! Internal variables
  integer(kind=i_def) :: k, ik

  do k = 0, nlayers - 1
    ik = k + 1 + (cell-1)*nlayers
    x_times_y(ik,:,:) = matmul(x(ik,:,:), y(ik,:,:))
  end do

end subroutine operator_x_times_y_kernel_code_real64

end module sci_operator_x_times_y_kernel_mod
