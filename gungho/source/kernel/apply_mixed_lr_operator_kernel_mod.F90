!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Apply the semi-implicit mixed operator to the continuity equation.

module apply_mixed_lr_operator_kernel_mod

use argument_mod,            only : arg_type,                       &
                                    GH_FIELD, GH_OPERATOR, GH_REAL, &
                                    GH_READ, GH_WRITE,              &
                                    CELLS
use constants_mod,           only : r_def, i_def
use kernel_mod,              only : kernel_type
use fs_continuity_mod,       only : W2, W3

implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: apply_mixed_lr_operator_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/           &
       arg_type(GH_FIELD,    GH_WRITE, W3),     & ! lhs_rho
       arg_type(GH_FIELD,    GH_READ,  W3),     & ! rho'
       arg_type(GH_OPERATOR, GH_READ,  W3, W3), & ! M3^-1
       arg_type(GH_OPERATOR, GH_READ,  W3, W2), & ! div
       arg_type(GH_REAL,     GH_READ),          & ! tau_r*dt
       arg_type(GH_FIELD,    GH_READ,  W2)      & ! u'*rho^ref
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: apply_mixed_lr_operator_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public apply_mixed_lr_operator_code

contains

!> @brief Compute the LHS of the continuity equation:
!!        lhs_rho = rho + tau_dt*inv_m3*div*f_star.
!! @param[in] cell Horizontal cell index
!! @param[in] nlayers Number of layers
!! @param[in,out] lhs_rho Mixed operator applied to the continuity equation
!! @param[in] rho Density field
!! @param[in] ncell1 Total number of cells for the inv_m3 operator
!! @param[in] inv_m3 Inverse mass matrix for the continuity equation
!! @param[in] ncell2 Total number of cells for the divergence operator
!! @param[in] div Divergence matrix from W2 to W3
!! @param[in] tau_dt tau_r*dt where tau_r is relaxation parameter and dt is
!!                   the timestep
!! @param[in] f_star Linearised mass flux: rho^ref * u'
!! @param[in] ndf_w3 Number of degrees of freedom per cell for the density space
!! @param[in] undf_w3 Unique number of degrees of freedom for the density space
!! @param[in] map_w3 Dofmap for the cell at the base of the column for the
!!                   density space
!! @param[in] ndf_w2 Number of degrees of freedom per cell for the wind space
!! @param[in] undf_w2 Unique number of degrees of freedom for the wind space
!! @param[in] map_w2 Dofmap for the cell at the base of the column for the wind space
subroutine apply_mixed_lr_operator_code(cell,                    &
                                        nlayers,                 &
                                        lhs_rho,                 &
                                        rho,                     &
                                        ncell1, inv_m3,          &
                                        ncell2, div,             &
                                        tau_dt,                  &
                                        f_star,                  &
                                        ndf_w3, undf_w3, map_w3, &
                                        ndf_w2, undf_w2, map_w2)

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in) :: cell, nlayers
  integer(kind=i_def),                    intent(in) :: ncell1, ncell2
  integer(kind=i_def),                    intent(in) :: undf_w2, ndf_w2
  integer(kind=i_def),                    intent(in) :: undf_w3, ndf_w3
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3

  ! Fields
  real(kind=r_def), dimension(undf_w3), intent(inout) :: lhs_rho
  real(kind=r_def), dimension(undf_w2), intent(in)    :: f_star
  real(kind=r_def), dimension(undf_w3), intent(in)    :: rho

  ! Operators
  real(kind=r_def), dimension(ndf_w3, ndf_w3, ncell1), intent(in) :: inv_m3
  real(kind=r_def), dimension(ndf_w3, ndf_w2, ncell2), intent(in) :: div

  ! Constants
  real(kind=r_def), intent(in) :: tau_dt

  ! Internal variables
  integer(kind=i_def)                 :: df, k, ik
  real(kind=r_def), dimension(ndf_w2) :: f_e
  real(kind=r_def), dimension(ndf_w3) :: r_e, lhs_e

  do k = 0, nlayers-1
    do df = 1, ndf_w2
      f_e(df) = f_star(map_w2(df)+k)
    end do
    do df = 1, ndf_w3
      r_e(df) = rho(map_w3(df)+k)
    end do

    ik = (cell-1)*nlayers + k + 1

    ! lhs_rho for this element
    lhs_e = r_e + tau_dt*matmul(inv_m3(:,:,ik), matmul(div(:,:,ik), f_e))
    do df = 1, ndf_w3
      lhs_rho(map_w3(df)+k) = lhs_e(df)
    end do
  end do

end subroutine apply_mixed_lr_operator_code

end module apply_mixed_lr_operator_kernel_mod
