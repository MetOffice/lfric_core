!-------------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Compute the projection in of the pressure field into the same space as
!>        density
module project_pressure_kernel_mod

use argument_mod,      only : arg_type, func_type,                 &
                              GH_FIELD, GH_OPERATOR,               &
                              GH_READ, GH_WRITE,                   &
                              ANY_SPACE_1, ANY_SPACE_2, W3,        &
                              GH_BASIS, GH_DIFF_BASIS, CELLS,      &                    
                              GH_QUADRATURE_XYoZ, QUADRATURE_XYoZ
use constants_mod,     only : r_def
use kernel_mod,        only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: project_pressure_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE,  W3),                            &
       arg_type(GH_FIELD,   GH_READ,   W3),                            &
       arg_type(GH_FIELD,   GH_READ,   ANY_SPACE_1),                   &
       arg_type(GH_FIELD*3, GH_READ,   ANY_SPACE_2),                   &
       arg_type(GH_OPERATOR, GH_READ,  W3, W3)                         &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(W3,          GH_BASIS),                               &
       func_type(ANY_SPACE_1, GH_BASIS),                               &
       func_type(ANY_SPACE_2, GH_DIFF_BASIS)                           &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = QUADRATURE_XYoZ
contains
  procedure, nopass ::project_pressure_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface project_pressure_kernel_type
   module procedure project_pressure_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public project_pressure_code
contains

type(project_pressure_kernel_type) function project_pressure_kernel_constructor() result(self)
  return
end function project_pressure_kernel_constructor

!> @brief Compute the pressure gradient component of the momentum equation
!! @param[in] cell Horizontal cell index
!! @param[in] nlayers Number of layers
!! @param[inout] exner Pressure field
!! @param[in] rho Density
!! @param[in] theta Potential temperature
!! @param[in] chi1 coordinate field
!! @param[in] chi2 coordinate field
!! @param[in] chi3 coordinate field
!! @param[in] ncell_3d number of cells
!! @param[in] m3_inv Inverse of W3 mass matrix
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number unique of degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Basis functions evaluated at gaussian quadrature points 
!! @param[in] ndf_wt Number of degrees of freedom per cell for theta space
!! @param[in] undf_wt Number unique of degrees of freedom  for theta space
!! @param[in] map_wt Dofmap for the cell at the base of the column for theta space
!! @param[in] wt_basis Basis functions evaluated at gaussian quadrature points 
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi space
!! @param[in] undf_chi Number unique of degrees of freedom  for chi space
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi space
!! @param[in] chi_diff_basis Basis functions evaluated at gaussian quadrature points 
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h horizontal quadrature weights
!! @param[in] wqp_v vertical quadrature weights
subroutine project_pressure_code(cell, nlayers,                                &
                                 exner, rho, theta, chi1, chi2, chi3,          &
                                 ncell_3d, m3_inv,                             &
                                 ndf_w3, undf_w3, map_w3, w3_basis,            &
                                 ndf_wt, undf_wt, map_wt, wt_basis,            &
                                 ndf_chi, undf_chi, map_chi, chi_diff_basis,   &
                                 nqp_h, nqp_v, wqp_h, wqp_v                    &
                                 )
                           
  use calc_exner_pointwise_mod,only: calc_exner_pointwise
  use coordinate_jacobian_mod, only: coordinate_jacobian

  !Arguments
  integer, intent(in) :: nlayers, nqp_h, nqp_v, ncell_3d, cell
  integer, intent(in) :: ndf_wt, ndf_w3, ndf_chi
  integer, intent(in) :: undf_wt, undf_w3, undf_chi
  integer, dimension(ndf_wt),  intent(in) :: map_wt
  integer, dimension(ndf_chi), intent(in) :: map_chi
  integer, dimension(ndf_w3),  intent(in) :: map_w3
  

  real(kind=r_def), dimension(1,ndf_w3, nqp_h,nqp_v), intent(in) :: w3_basis  
  real(kind=r_def), dimension(1,ndf_wt, nqp_h,nqp_v), intent(in) :: wt_basis 
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis   

  real(kind=r_def), dimension(undf_w3),  intent(out) :: exner
  real(kind=r_def), dimension(undf_w3),  intent(in)  :: rho
  real(kind=r_def), dimension(undf_wt),  intent(in)  :: theta
  real(kind=r_def), dimension(undf_chi), intent(in)  :: chi1, chi2, chi3

  real(kind=r_def), dimension(ndf_w3,ndf_w3,ncell_3d), intent(in) :: m3_inv

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer               :: df, k, ik 
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_w3)          :: rho_e
  real(kind=r_def), dimension(ndf_w3)          :: r_exner, exner_e
  real(kind=r_def), dimension(ndf_wt)          :: theta_e
  real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac

  real(kind=r_def) :: exner_at_quad, rho_at_quad, theta_at_quad

  do k = 0, nlayers-1
    do df = 1, ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                             chi_diff_basis, jac, dj)

    do df = 1, ndf_w3
      rho_e(df) = rho( map_w3(df) + k )
      r_exner(df) = 0.0_r_def
    end do    
    do df = 1, ndf_wt
      theta_e(df) = theta( map_wt(df) + k )
    end do   

    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        rho_at_quad = 0.0_r_def 
        do df = 1, ndf_w3
          rho_at_quad  = rho_at_quad + rho_e(df)*w3_basis(1,df,qp1,qp2) 
        end do
        theta_at_quad = 0.0_r_def
        do df = 1, ndf_wt
          theta_at_quad = theta_at_quad + theta_e(df)*wt_basis(1,df,qp1,qp2)
        end do
        exner_at_quad = wqp_h(qp1)*wqp_v(qp2)*dj(qp1,qp2) &
                      *calc_exner_pointwise(rho_at_quad, theta_at_quad)

        do df = 1, ndf_w3
          r_exner(df) = r_exner(df) + w3_basis(1,df,qp1,qp2)*exner_at_quad
        end do
      end do
    end do
    ik = 1 + k + (cell-1)*nlayers
    exner_e = matmul(m3_inv(:,:,ik),r_exner)
    do df = 1, ndf_w3
      exner( map_w3(df) + k ) =  exner_e(df)
    end do 
  end do
  
end subroutine project_pressure_code

end module project_pressure_kernel_mod
