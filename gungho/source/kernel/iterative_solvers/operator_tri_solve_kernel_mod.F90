!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Tridiagonal solver using Thomas algorithm for Wtheta and W2v fields.
!> @details Solve the triadiagonal system of equations
!!          \f[ A^+ y^+ + A^0 y + A^- y^- = x \]f
!!          for known x and matrix A using the Thomas algorithm.
!!          Code is only valid for lowest order elements on Wtheta
!!          and W2v spaces (2 dofs per cell). The tridiagonal operator
!!          is extracted from the LMA operators.
module operator_tri_solve_kernel_mod

  use argument_mod,  only: arg_type,              &
                           GH_FIELD, GH_OPERATOR, &
                           GH_REAL,               &
                           GH_READ, GH_WRITE,     &
                           CELL_COLUMN,           &
                           ANY_DISCONTINUOUS_SPACE_1
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  type, public, extends(kernel_type) :: operator_tri_solve_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                       &
         arg_type(GH_FIELD,    GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                  ANY_DISCONTINUOUS_SPACE_1)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: operator_tri_solve_code
  end type operator_tri_solve_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: operator_tri_solve_code

contains

!> @brief Tridiagonal solver using Thomas algorithm where the tridiagonal
!!        arrays are extracted from an operator.
!> @param[in]     cell        Horizontal cell index
!> @param[in]     nlayers     Number of levels to solve over
!> @param[in,out] y           LHS field to solve for
!> @param[in]     x           RHS field
!> @param[in]     ncell       Total number of cells in 3d mesh
!> @param[in]     mass_matrix Matrix to extract tridiagonal operator from
!> @param[in]     ndf         Number of dofs per cell for all fields
!> @param[in]     undf        Size of all field arrays
!> @param[in]     map         Array containing the address of the first dof in the column
subroutine operator_tri_solve_code(cell, nlayers, &
                                   y, x,          &
                                   ncell, matrix, &
                                   ndf, undf, map)

  implicit none

  integer(kind=i_def), intent(in) :: ndf, undf
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: cell, ncell

  integer(kind=i_def), dimension(ndf),  intent(in) :: map

  real(kind=r_def), dimension(undf),            intent(inout) :: y
  real(kind=r_def), dimension(undf),            intent(in)    :: x
  real(kind=r_def), dimension(ndf, ndf, ncell), intent(in)    :: matrix

  integer(kind=i_def)                    :: ik, k, ij
  real(kind=r_def), dimension(0:nlayers) :: tri_p, tri_0, tri_m, x_new, tri_p_new
  real(kind=r_def)                       :: denom

  ! Extract tridiagonal arrays from operator
  tri_0(0) = 0.0_r_def
  tri_m(0) = 0.0_r_def
  do k = 0, nlayers - 1
    ik = 1 + k + (cell - 1) * nlayers
    tri_0(k)   = matrix(1,1,ik) + tri_0(k)
    tri_p(k)   = matrix(1,2,ik)
    tri_m(k+1) = matrix(2,1,ik)
    tri_0(k+1) = matrix(2,2,ik)
  end do
  tri_p(nlayers) = 0.0_r_def

  k  = 0
  ij = map(1)
  denom = 1.0_r_def/tri_0(k)
  tri_p_new(k) = tri_p(k)*denom
  x_new(k)     = x(ij+k) *denom

  do k = 1, nlayers
    denom = 1.0_r_def/(tri_0(k) - tri_m(k)*tri_p_new(k-1))
    tri_p_new(k) = tri_p(k)*denom
    x_new(k)     = (x(ij+k) - tri_m(k)*x_new(k-1))*denom
  end do

  k = nlayers
  y(ij+k) = x_new(k)
  do k = nlayers-1, 0, -1
    y(ij+k) = x_new(k) - tri_p_new(k)*y(ij+k+1)
  end do

end subroutine operator_tri_solve_code

end module operator_tri_solve_kernel_mod
