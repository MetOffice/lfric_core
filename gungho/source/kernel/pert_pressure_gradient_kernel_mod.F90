!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the perturbation pressure gradient for lhs of the momentum equation
module pert_pressure_gradient_kernel_mod

use argument_mod,      only : arg_type, func_type,                 &
                              GH_FIELD, GH_READ, GH_INC,           &
                              ANY_SPACE_9, W2, W3,                 &
                              GH_BASIS, GH_DIFF_BASIS,             &
                              CELLS, GH_QUADRATURE_XYoZ,           &
                              QUADRATURE_XYoZ

use constants_mod,     only : r_def
use kernel_mod,        only : kernel_type
use planet_config_mod, only : Cp, kappa

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: pert_pressure_gradient_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W2),                              &
       arg_type(GH_FIELD,   GH_READ, W3),                              &
       arg_type(GH_FIELD,   GH_READ, W3),                              &
       arg_type(GH_FIELD,   GH_READ, ANY_SPACE_9),                     &
       arg_type(GH_FIELD,   GH_READ, ANY_SPACE_9)                      &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(W2, GH_BASIS, GH_DIFF_BASIS),                         &
       func_type(W3, GH_BASIS),                                        &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                 &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = QUADRATURE_XYoZ
contains
  procedure, nopass ::pert_pressure_gradient_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface pert_pressure_gradient_kernel_type
  module procedure pert_pressure_gradient_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public pert_pressure_gradient_code
contains

type(pert_pressure_gradient_kernel_type) function pert_pressure_gradient_kernel_constructor() result(self)
  return
end function pert_pressure_gradient_kernel_constructor
!> @details The kernel computes the perturbation pressure gradient part for the 
!>         lhs of the momentum equation for the nonlinear equations,
!>         pert_pressure_gradient:
!>         \f[  cp\left(\theta^{ref}\nabla\Pi' + \theta'\nabla\Pi^{ref}\right) \f]
!! @param[in] nlayers Number of layers
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2 Number unique of degrees of freedom  for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Basis functions evaluated at quadrature points 
!! @param[in] w2_diff_basis Differntial of the basis functions evaluated at  quadrature points
!! @param[inout] r_u Momentum equation residual
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number unique of degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Basis functions evaluated at gaussian quadrature points 
!! @param[in] rho Density
!! @param[in] rho_ref Reference density
!! @param[in] ndf_wtheta Number of degrees of freedom per cell for wtheta
!! @param[in] undf_wtheta Number unique of degrees of freedom  for wtheta
!! @param[in] map_wtheta Dofmap for the cell at the base of the column for wtheta
!! @param[in] wtheta_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] wtheta_diff_basis Differntial of the basis functions evaluated at gaussian quadrature point
!! @param[in] theta Potential temperature
!! @param[in] theta_ref Reference potential temperature
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine pert_pressure_gradient_code(nlayers,                                     &
                                  r_u, rho, rho_ref, theta, theta_ref,              &
                                  ndf_w2, undf_w2, map_w2, w2_basis, w2_diff_basis, &
                                  ndf_w3, undf_w3, map_w3, w3_basis,                &
                                  ndf_wtheta, undf_wtheta, map_wtheta,              &
                                  wtheta_basis, wtheta_diff_basis,                  &
                                  nqp_h, nqp_v, wqp_h, wqp_v                        &
                                  )
                           
  use calc_exner_pointwise_mod, only: calc_exner_pointwise
  
  !Arguments
  integer, intent(in) :: nlayers,nqp_h, nqp_v
  integer, intent(in) :: ndf_wtheta, ndf_w2, ndf_w3
  integer, intent(in) :: undf_wtheta, undf_w2, undf_w3

  integer, dimension(ndf_wtheta), intent(in) :: map_wtheta
  integer, dimension(ndf_w2), intent(in)     :: map_w2
  integer, dimension(ndf_w3), intent(in)     :: map_w3
  
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in)     :: w3_basis
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in)     :: w2_basis
  real(kind=r_def), dimension(1,ndf_wtheta,nqp_h,nqp_v), intent(in) :: wtheta_basis
  real(kind=r_def), dimension(1,ndf_w2,nqp_h,nqp_v), intent(in)     :: w2_diff_basis
  real(kind=r_def), dimension(3,ndf_wtheta,nqp_h,nqp_v), intent(in) :: wtheta_diff_basis

  real(kind=r_def), dimension(undf_w2), intent(inout)  :: r_u
  real(kind=r_def), dimension(undf_w3), intent(in)     :: rho, rho_ref
  real(kind=r_def), dimension(undf_wtheta), intent(in) :: theta, theta_ref

  real(kind=r_def), dimension(nqp_h), intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) ::  wqp_v

  !Internal variables
  integer               :: df, k 
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_w3)     :: rho_e, rho_ref_e
  real(kind=r_def), dimension(ndf_w2)     :: ru_e
  real(kind=r_def), dimension(ndf_wtheta) :: theta_e, theta_ref_e

  real(kind=r_def) :: grad_theta_at_quad(3), v(3), grad_theta_ref_at_quad(3)
  real(kind=r_def) :: exner_at_quad, rho_at_quad, theta_at_quad, &
                      grad_term, dv, rho_ref_at_quad, exner_ref_at_quad, &
                      theta_ref_at_quad
  
  do k = 0, nlayers-1
    do df = 1, ndf_w3
      rho_e(df)     = rho    ( map_w3(df) + k )
      rho_ref_e(df) = rho_ref( map_w3(df) + k )
    end do    
    do df = 1, ndf_wtheta
      theta_e(df)     = theta    ( map_wtheta(df) + k )
      theta_ref_e(df) = theta_ref( map_wtheta(df) + k )
    end do   
    do df = 1, ndf_w2
      ru_e(df) = 0.0_r_def
    end do 
  ! compute the RHS integrated over one cell    
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        rho_at_quad = 0.0_r_def 
        rho_ref_at_quad = 0.0_r_def 
        do df = 1, ndf_w3
          rho_at_quad      = rho_at_quad     + rho_e(df)    *w3_basis(1,df,qp1,qp2) 
          rho_ref_at_quad  = rho_ref_at_quad + rho_ref_e(df)*w3_basis(1,df,qp1,qp2) 
        end do
        theta_at_quad = 0.0_r_def
        grad_theta_at_quad(:) = 0.0_r_def
        theta_ref_at_quad = 0.0_r_def
        grad_theta_ref_at_quad(:) = 0.0_r_def
        do df = 1, ndf_wtheta
          theta_at_quad   = theta_at_quad                                      &
                          + theta_e(df)*wtheta_basis(1,df,qp1,qp2)
          grad_theta_at_quad(:) = grad_theta_at_quad(:) &
                                + theta_e(df)*wtheta_diff_basis(:,df,qp1,qp2)
          theta_ref_at_quad   = theta_ref_at_quad                              &
                          + theta_ref_e(df)*wtheta_basis(1,df,qp1,qp2)
          grad_theta_ref_at_quad(:) = grad_theta_ref_at_quad(:) &
                                + theta_ref_e(df)*wtheta_diff_basis(:,df,qp1,qp2)
        end do

        exner_ref_at_quad = calc_exner_pointwise(rho_ref_at_quad, &
                                                 theta_ref_at_quad)
        exner_at_quad = (kappa/(1.0_r_def - kappa))*exner_ref_at_quad &
                      *(theta_at_quad/theta_ref_at_quad &
                      + rho_at_quad/rho_ref_at_quad)

        do df = 1, ndf_w2
          v  = w2_basis(:,df,qp1,qp2)
          dv = w2_diff_basis(1,df,qp1,qp2)
          ! theta_prime * grad(exner_ref) term
          grad_term = Cp*exner_ref_at_quad * (                       &
                      theta_at_quad * dv                             &
                    + dot_product( grad_theta_at_quad(:),v)          &
                                         )
          ! theta_ref * grad(exner_prime) term
          grad_term = grad_term + Cp*exner_at_quad * (               &
                      theta_ref_at_quad * dv                         &
                    + dot_product( grad_theta_ref_at_quad(:),v)      &
                                         )

          ru_e(df) = ru_e(df) +  wqp_h(qp1)*wqp_v(qp2)*grad_term
        end do
      end do
    end do
    do df = 1, ndf_w2
      r_u( map_w2(df) + k ) =  r_u( map_w2(df) + k ) + ru_e(df)
    end do 
  end do
  
end subroutine pert_pressure_gradient_code

end module pert_pressure_gradient_kernel_mod
