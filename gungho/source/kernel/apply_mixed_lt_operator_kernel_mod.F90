!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Apply the semi-implicit mixed operator to the potential temperature equation.

module apply_mixed_lt_operator_kernel_mod

use argument_mod,      only : arg_type,              &
                              GH_FIELD, GH_OPERATOR, &
                              GH_READ, GH_READWRITE, &
                              CELLS
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type
use fs_continuity_mod, only : W2, Wtheta

implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: apply_mixed_lt_operator_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                       &
       arg_type(GH_FIELD,    GH_READWRITE, Wtheta),         & ! lhs_theta
       arg_type(GH_FIELD,    GH_READ,      W2),             & ! u'
       arg_type(GH_FIELD,    GH_READ,      Wtheta),         & ! theta'
       arg_type(GH_OPERATOR, GH_READ,      Wtheta, Wtheta), & ! mtheta
       arg_type(GH_OPERATOR, GH_READ,      Wtheta, W2),     & ! ptheta2
       arg_type(GH_FIELD,    GH_READ,      Wtheta)          & ! norm_theta

       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: apply_mixed_lt_operator_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public apply_mixed_lt_operator_code

contains

!> @brief Compute the LHS of the potential temperature equation
!!        lhs_theta = norm_theta*(mtheta*theta + Ptheta2*wind).
!! @param[in] cell Horizontal cell index
!! @param[in] nlayers Number of layers
!! @param[inout] lhs_theta Mixed operator applied to the potential temperature equation
!! @param[in] wind Wind field
!! @param[in] theta Potential temperature field
!! @param[in] ncell1 Total number of cells for the mtheta operator
!! @param[in] mtheta Generalised mass matrix for the potential temperature equation
!! @param[in] ncell2 Total number of cells for the ptheta2 operator
!! @param[in] ptheta2 Generalised projection matrix from W2 to Wtheta
!! @param[in] norm_theta Normalisation field for the potential temperature equation
!! @param[in] ndf_wt Number of degrees of freedom per cell for the potential
!!                   temperature space
!! @param[in] undf_wt Unique number of degrees of freedom for the potential
!!                    temperature space
!! @param[in] map_wt Dofmap for the cell at the base of the column for the
!!                   potential temperature space
!! @param[in] ndf_w2 Number of degrees of freedom per cell for the wind space
!! @param[in] undf_w2 Unique number of degrees of freedom for the wind space
!! @param[in] map_w2 Dofmap for the cell at the base of the column for the wind space
subroutine apply_mixed_lt_operator_code(cell,                    &
                                        nlayers,                 &
                                        lhs_theta,               &
                                        wind, theta,             &
                                        ncell1, mtheta,          &
                                        ncell2, ptheta2,         &
                                        norm_theta,              &
                                        ndf_wt, undf_wt, map_wt, &
                                        ndf_w2, undf_w2, map_w2)

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in) :: cell, nlayers
  integer(kind=i_def),                    intent(in) :: ncell1, ncell2
  integer(kind=i_def),                    intent(in) :: undf_w2, ndf_w2
  integer(kind=i_def),                    intent(in) :: undf_wt, ndf_wt
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt

  ! Fields
  real(kind=r_def), dimension(undf_wt), intent(inout) :: lhs_theta
  real(kind=r_def), dimension(undf_w2), intent(in)    :: wind
  real(kind=r_def), dimension(undf_wt), intent(in)    :: theta, norm_theta

  ! Operators
  real(kind=r_def), dimension(ndf_wt, ndf_wt, ncell1), intent(in) :: mtheta
  real(kind=r_def), dimension(ndf_wt, ndf_w2, ncell2), intent(in) :: ptheta2

  ! Internal variables
  integer(kind=i_def)                 :: df, k, ik
  real(kind=r_def), dimension(ndf_w2) :: u_e
  real(kind=r_def), dimension(ndf_wt) :: t_e, lhs_e

  do k = 0, nlayers-1
    do df = 1, ndf_w2
      u_e(df) = wind(map_w2(df)+k)
    end do
    do df = 1, ndf_wt
      t_e(df) = theta(map_wt(df)+k)
    end do

    ik = (cell-1)*nlayers + k + 1

    ! lhs_theta for this element
    lhs_e = matmul(mtheta(:,:,ik),t_e) + matmul(ptheta2(:,:,ik),u_e)
    do df = 1, ndf_wt
      lhs_theta(map_wt(df)+k) = lhs_theta(map_wt(df)+k) &
                              + norm_theta(map_wt(df)+k)*lhs_e(df)
    end do
  end do

end subroutine apply_mixed_lt_operator_code

end module apply_mixed_lt_operator_kernel_mod
