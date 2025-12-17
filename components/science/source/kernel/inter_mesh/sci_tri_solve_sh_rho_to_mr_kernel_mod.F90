!-----------------------------------------------------------------------------
! Copyright (c) 2020,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Tridiagonal solver for mapping between Wtheta and shifted W3 spaces.
!>
!> Solve the triadiagonal system of equations
!> \f$A^+ y^+ + A^0 y + A^- y^- = x\f$ for \f$y\f$ in Wtheta given known
!> \f$x\f$ in shifted W3 and matrix \f$A\f$ using the Thomas algorithm.
!>
!> Code is only valid for lowest order elements.
!>
module sci_tri_solve_sh_rho_to_mr_kernel_mod

  use argument_mod,      only: arg_type,          &
                               GH_FIELD, GH_REAL, &
                               GH_READ, GH_WRITE, &
                               CELL_COLUMN
  use constants_mod,     only: r_single, r_double, i_def
  use fs_continuity_mod, only: W3, Wtheta
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  type, public, extends(kernel_type) :: tri_solve_sh_rho_to_mr_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                   &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, Wtheta), & ! field_wt
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),     & ! field_sh_w3
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  W3)      & ! tri_below/diag/above
         /)
    integer :: operates_on = CELL_COLUMN
  end type tri_solve_sh_rho_to_mr_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: tri_solve_sh_rho_to_mr_code

  ! Generic interface for real32 and real64 types
  interface tri_solve_sh_rho_to_mr_code
    module procedure  &
      tri_solve_sh_rho_to_mr_code_r_single, &
      tri_solve_sh_rho_to_mr_code_r_double
  end interface

contains

!> @brief Tridiagonal solver using Thomas algorithm
!> @param[in]  nlayers Number of levels to solve over
!> @param[in,out] field_wt LHS field in Wtheta to solve for.
!> @param[in]  field_sh_w3 RHS field in shifted W3.
!> @param[in]  tri_below Below-diagonal elements of tridiagonal matrix. A field
!> in shifted W3 space.
!> @param[in]  tri_diag Centre diagonal elements of tridiagonal matrix. A field
!> in shifted W3 space.
!> @param[in]  tri_above Above-diagonal elements of tridiagonal matrix. A field
!> in shifted W3 space.
!> @param[in]  ndf_wt Number of dofs per cell for Wtheta
!> @param[in]  undf_wt Size of Wtheta field array
!> @param[in]  map_wt Dofmap for Wtheta
!> @param[in]  ndf_sh_w3 Number of dofs per cell for shifted W3. Should be 1.
!> @param[in]  undf_sh_w3 Size of shifted W3 field arrays.
!> @param[in]  map_sh_w3 Dofmap for shifted W3 space

! R_SINGLE PRECISION
! ==================
subroutine tri_solve_sh_rho_to_mr_code_r_single(                                  &
                                                 nlayers,                         &
                                                 field_wt, field_sh_w3,           &
                                                 tri_below, tri_diag, tri_above,  &
                                                 ndf_wt, undf_wt, map_wt,         &
                                                 ndf_sh_w3, undf_sh_w3, map_sh_w3 &
                                                 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: ndf_sh_w3, undf_sh_w3
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), dimension(ndf_wt),    intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_sh_w3), intent(in) :: map_sh_w3

  real(kind=r_single), dimension(undf_wt),     intent(inout) :: field_wt
  real(kind=r_single), dimension(undf_sh_w3),  intent(in)    :: field_sh_w3
  real(kind=r_single), dimension(undf_sh_w3),  intent(in)    :: tri_above
  real(kind=r_single), dimension(undf_sh_w3),  intent(in)    :: tri_diag
  real(kind=r_single), dimension(undf_sh_w3),  intent(in)    :: tri_below

  ! Internal variables
  integer(kind=i_def)                     :: k, ij, nlayers_shifted
  real(kind=r_single), dimension(nlayers+1)  :: rhs_new, tri_above_new
  real(kind=r_single)                        :: denom

  nlayers_shifted = nlayers + 1

  k  = 0
  ij = map_sh_w3(1)
  denom = 1.0_r_single / tri_diag(ij+k)
  tri_above_new(1) = tri_above(ij+k) * denom
  rhs_new(1) = field_sh_w3(ij+k) * denom

  do k = 1, nlayers_shifted-1
    denom = 1.0_r_single / (tri_diag(ij+k) - tri_below(ij+k) * tri_above_new(k))
    tri_above_new(k+1) = tri_above(ij+k) * denom
    rhs_new(k+1) = (field_sh_w3(ij+k) - tri_below(ij+k) * rhs_new(k)) * denom
  end do

  k = nlayers_shifted-1
  field_wt(ij+k) = rhs_new(k+1)
  do k = nlayers_shifted-2, 0, -1
    field_wt(ij+k) = rhs_new(k+1) - tri_above_new(k+1)*field_wt(ij+k+1)
  end do

end subroutine tri_solve_sh_rho_to_mr_code_r_single

! R_DOUBLE PRECISION
! ==================
subroutine tri_solve_sh_rho_to_mr_code_r_double(                                  &
                                                 nlayers,                         &
                                                 field_wt, field_sh_w3,           &
                                                 tri_below, tri_diag, tri_above,  &
                                                 ndf_wt, undf_wt, map_wt,         &
                                                 ndf_sh_w3, undf_sh_w3, map_sh_w3 &
                                                 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: ndf_sh_w3, undf_sh_w3
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), dimension(ndf_wt),    intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_sh_w3), intent(in) :: map_sh_w3

  real(kind=r_double), dimension(undf_wt),     intent(inout) :: field_wt
  real(kind=r_double), dimension(undf_sh_w3),  intent(in)    :: field_sh_w3
  real(kind=r_double), dimension(undf_sh_w3),  intent(in)    :: tri_above
  real(kind=r_double), dimension(undf_sh_w3),  intent(in)    :: tri_diag
  real(kind=r_double), dimension(undf_sh_w3),  intent(in)    :: tri_below

  ! Internal variables
  integer(kind=i_def)                     :: k, ij, nlayers_shifted
  real(kind=r_double), dimension(nlayers+1)  :: rhs_new, tri_above_new
  real(kind=r_double)                        :: denom

  nlayers_shifted = nlayers + 1

  k  = 0
  ij = map_sh_w3(1)
  denom = 1.0_r_double / tri_diag(ij+k)
  tri_above_new(1) = tri_above(ij+k) * denom
  rhs_new(1) = field_sh_w3(ij+k) * denom

  do k = 1, nlayers_shifted-1
    denom = 1.0_r_double / (tri_diag(ij+k) - tri_below(ij+k) * tri_above_new(k))
    tri_above_new(k+1) = tri_above(ij+k) * denom
    rhs_new(k+1) = (field_sh_w3(ij+k) - tri_below(ij+k) * rhs_new(k)) * denom
  end do

  k = nlayers_shifted-1
  field_wt(ij+k) = rhs_new(k+1)
  do k = nlayers_shifted-2, 0, -1
    field_wt(ij+k) = rhs_new(k+1) - tri_above_new(k+1)*field_wt(ij+k+1)
  end do

end subroutine tri_solve_sh_rho_to_mr_code_r_double

end module sci_tri_solve_sh_rho_to_mr_kernel_mod
