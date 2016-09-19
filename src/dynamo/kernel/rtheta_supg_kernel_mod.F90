!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief The kernel applies supg to the thermodynamic equation for the
!>        nonlinear equations.
!
module rtheta_supg_kernel_mod

use argument_mod,            only : arg_type, func_type,                     &
                                    GH_FIELD, GH_READ, GH_INC,               &
                                    W0, W2, W3, ANY_SPACE_9,                 &
                                    GH_BASIS, GH_DIFF_BASIS,                 &
                                    CELLS
use constants_mod,           only : r_def, EPS
use coordinate_jacobian_mod, only : coordinate_jacobian
use kernel_mod,              only : kernel_type
use timestepping_config_mod, only : dt

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy
!> layer.
!
type, public, extends(kernel_type) :: rtheta_supg_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W0),                              &
       arg_type(GH_FIELD,   GH_READ, W0),                              &
       arg_type(GH_FIELD,   GH_READ, W0),                              &
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD,   GH_READ, W3),                              &
       arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9)                               &
       /)
  type(func_type) :: meta_funcs(4) = (/                                &
       func_type(W0, GH_BASIS, GH_DIFF_BASIS),                         &
       func_type(W2, GH_BASIS),                                        &
       func_type(W3, GH_BASIS),                                        &
       func_type(ANY_SPACE_9, GH_DIFF_BASIS)                          &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass ::rtheta_supg_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface rtheta_supg_kernel_type
   module procedure rtheta_supg_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public rtheta_supg_code
contains

type(rtheta_supg_kernel_type) function rtheta_supg_kernel_constructor() &
                              result(self)
  return
end function rtheta_supg_kernel_constructor
!> @details Kernel to  compute the application of streamline upwind
!>         Petrov-Galerkin (supg) method to the thermodynamic equation. This
!>         replaces the test function
!>         \f[ 
!!         \gamma \longrightarrow \gamma^* \equiv \gamma
!!         + \frac{\alpha}{\Delta}\frac{\mathbf{F}}{\left|\mathbf{F}\right|}
!!         .\nabla\gamma
!!         \f]
!!
!!         where \f$\alpha\f$ is some parameter and \f$\Delta\f$ is a measure
!!         of the grid spacing.
!! @param[in] nlayers Number of layers
!! @param[in] ndf_w0 Number of degrees of freedom per cell for w0
!! @param[in] undf_w0  Number of unique degrees of freedom  for w0
!! @param[in] map_w0 Dofmap for the cell at the base of the column for w0
!! @param[in] w0_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] w0_diff_basis Differential basis functions evaluated at gaussian
!!                          quadrature points
!! @param[inout] r_theta right hand side of the thermodynamic equation
!! @param[in] theta Potential temperature
!! @param[in] theta_n Potential temperature at timelevel n
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2  Number of unique degrees of freedom  for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] f Mass flux
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number of unique degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi  Number of unique degrees of freedom  for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_diff_basis Differential basis functions evaluated at gaussian
!!                          quadrature points
!! @param[in] rho Density
!! @param[inout] chi1 Chi in the 1st dir
!! @param[inout] chi2 Chi in the 2nd dir
!! @param[inout] chi3 Chi in the 3rd dir
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Weights of the horizontal quadrature points
!! @param[in] wqp_v Weights of the vertical quadrature points
subroutine rtheta_supg_code(nlayers,                                          &
                            r_theta, theta, theta_n,                           &
                            f, rho,                                            &
                            chi1, chi2, chi3,                                  &
                            ndf_w0, undf_w0, map_w0, w0_basis, w0_diff_basis,  &
                            ndf_w2, undf_w2, map_w2, w2_basis,                 &
                            ndf_w3, undf_w3, map_w3, w3_basis,                 &
                            ndf_chi, undf_chi, map_chi, chi_diff_basis,  &
                            nqp_h, nqp_v, wqp_h, wqp_v )

  !Arguments
  integer, intent(in) :: nlayers, nqp_h, nqp_v
  integer, intent(in) :: ndf_w0, ndf_w2, ndf_w3, undf_w0, undf_w2, undf_w3, ndf_chi, undf_chi

  integer, dimension(ndf_w0), intent(in)  :: map_w0
  integer, dimension(ndf_w2), intent(in)  :: map_w2
  integer, dimension(ndf_w3), intent(in)  :: map_w3
  integer, dimension(ndf_chi), intent(in) :: map_chi

  real(kind=r_def), dimension(1,ndf_w0,nqp_h,nqp_v), intent(in)  :: w0_basis
  real(kind=r_def), dimension(3,ndf_w0,nqp_h,nqp_v), intent(in)  :: w0_diff_basis
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in)  :: w2_basis
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in)  :: w3_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w0), intent(inout) :: r_theta
  real(kind=r_def), dimension(undf_w0), intent(in)    :: theta, theta_n
  real(kind=r_def), dimension(undf_chi), intent(in)   :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_w2), intent(in)    :: f
  real(kind=r_def), dimension(undf_w3), intent(in)    :: rho

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer               :: df, k 
  integer               :: qp1, qp2

  real(kind=r_def), dimension(ndf_w0)          :: rtheta_e, theta_e, &
                                                  theta_n_e
  real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
  real(kind=r_def), dimension(ndf_w2)          :: f_e
  real(kind=r_def), dimension(ndf_w3)          :: rho_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def) :: f_at_quad(3), grad_theta_at_quad(3) 
  real(kind=r_def) :: rho_at_quad, abs_f, gamma, dthetadt, supg
  real(kind=r_def) :: advective_term     

  supg = 5.0_r_def

  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_chi
      chi1_e(df) = chi1( map_chi(df) + k )
      chi2_e(df) = chi2( map_chi(df) + k )
      chi3_e(df) = chi3( map_chi(df) + k )
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                             chi_diff_basis, jac, dj)
    do df = 1, ndf_w0
      rtheta_e(df) = 0.0_r_def
      theta_e(df)   = theta( map_w0(df) + k )
      theta_n_e(df) = theta_n( map_w0(df) + k )
    end do
    do df = 1, ndf_w2
      f_e(df) = f( map_w2(df) + k )
    end do
    do df = 1, ndf_w3
      rho_e(df) = rho( map_w3(df) + k )
    end do

    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        f_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w2
          f_at_quad(:)  = f_at_quad(:)  + f_e(df)*w2_basis(:,df,qp1,qp2)
        end do
        abs_f = max(EPS, sqrt(f_at_quad(1)**2 &
                              + f_at_quad(2)**2 + f_at_quad(3)**2))

        grad_theta_at_quad(:) = 0.0_r_def
        dthetadt = 0.0_r_def
        do df = 1, ndf_w0
          dthetadt = dthetadt  &
                   + (theta_e(df) - theta_n_e(df))*w0_basis(1,df,qp1,qp2)
          grad_theta_at_quad(:) = grad_theta_at_quad(:) &
                                + 0.5_r_def*(theta_e(df) + theta_n_e(df)) &
                                 *w0_diff_basis(:,df,qp1,qp2)
        end do
        dthetadt = dthetadt*dj(qp1,qp2)
        rho_at_quad = 0.0_r_def
        do df = 1, ndf_w3
          rho_at_quad  = rho_at_quad + rho_e(df)*w3_basis(1,df,qp1,qp2)
        end do

        advective_term = DT * dot_product( f_at_quad,grad_theta_at_quad ) &
                         / rho_at_quad

        do df = 1, ndf_w0
! Compute test function gamma = F/rho . grad(gamma)
          gamma = supg*dot_product(f_at_quad,&
                                   w0_diff_basis(:,df,qp1,qp2))/abs_f
          rtheta_e(df) = rtheta_e(df) &
                       - wqp_h(qp1)*wqp_v(qp2)*gamma*(dthetadt + advective_term)
        end do
      end do
    end do
    do df = 1, ndf_w0
      r_theta( map_w0(df) + k ) =  r_theta( map_w0(df) + k ) + rtheta_e(df)
    end do
  end do

end subroutine rtheta_supg_code

end module rtheta_supg_kernel_mod
