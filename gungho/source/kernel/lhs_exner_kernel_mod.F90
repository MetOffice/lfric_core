!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes lhs of the equation of state for the nonlinear equations 

!> @details The kernel computes the lhs of the equation of state for the nonlinear equations, 
!>          That is: \f[ lhs_{\Pi} = (1-\kappa)/\kappa * \Pi'/\Pi^* - \rho'/\rho^* - \theta'/\theta^* \f]
!>          Where ' are increments to the fields and ^* are the reference states
module lhs_exner_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,                     &
                                    GH_FIELD, GH_READ, GH_WRITE,             &
                                    ANY_SPACE_1, W3, ANY_SPACE_9,            &
                                    GH_BASIS, GH_DIFF_BASIS,                 &
                                    CELLS, GH_QUADRATURE_XYoZ,               &
                                    QUADRATURE_XYoZ
use constants_mod,           only : r_def, i_def

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: lhs_exner_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, W3),                             &
       arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_1),                    &
       arg_type(GH_FIELD,   GH_READ,  W3),                             &
       arg_type(GH_FIELD,   GH_READ,  W3),                             &
       arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_1),                    &
       arg_type(GH_FIELD,   GH_READ,  W3),                             &
       arg_type(GH_FIELD*3, GH_READ,  ANY_SPACE_9)                     &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(W3, GH_BASIS),                                        &
       func_type(ANY_SPACE_9, GH_DIFF_BASIS),                          &
       func_type(ANY_SPACE_1, GH_BASIS)                                &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = QUADRATURE_XYoZ
contains
  procedure, nopass ::lhs_exner_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface lhs_exner_kernel_type
   module procedure lhs_exner_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!------------^*---------------------------------------------------------------
public lhs_exner_code
contains

type(lhs_exner_kernel_type) function lhs_exner_kernel_constructor() result(self)
  return
end function lhs_exner_kernel_constructor

!> @brief Computes lhs of the equation of state for the nonlinear equations 
!! @param[in] nlayers Number of layers
!! @param[inout] l_exner LHS array for the equation of state
!! @param[in] theta Potential temperature increment
!! @param[in] rho Density increment
!! @param[in] exner Exner pressure increment
!! @param[in] theta_ref Potential temperature reference state
!! @param[in] rho_ref Density reference state
!! @param[in] chi1 First coordinate array
!! @param[in] chi2 Second coordinate array
!! @param[in] chi3 Third coordinate array
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number of (local) unique degrees of freedom
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Basis functions evaluated at quadrature points 
!! @param[in] ndf_w0 Number of degrees of freedom per cell for w0
!! @param[in] undf_w0 Number of (local) unique degrees of freedom
!! @param[in] map_w0 Dofmap for the cell at the base of the column for w0
!! @param[in] w0_basis Basis functions evaluated at quadrature points 
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number of (local) unique degrees of freedom for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_diff_basis Differential basis functions evaluated at quadrature points 
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h horizontal quadrature weights
!! @param[in] wqp_v vertical quadrature weights
subroutine lhs_exner_code(nlayers,                                         &
                          l_exner, theta, rho, exner,                      &
                          theta_ref, rho_ref,                              &
                          chi1, chi2, chi3,                                &
                          ndf_w3, undf_w3, map_w3, w3_basis,               &
                          ndf_w0, undf_w0, map_w0, w0_basis,               &
                          ndf_chi, undf_chi, map_chi, chi_diff_basis,      &
                          nqp_h, nqp_v, wqp_h, wqp_v )

  use coordinate_jacobian_mod,  only: coordinate_jacobian
  use planet_config_mod,        only: kappa
  use calc_exner_pointwise_mod, only: calc_exner_pointwise

  implicit none
  !Arguments
  integer(kind=i_def), intent(in) :: nlayers, nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_w0, ndf_w3, ndf_chi
  integer(kind=i_def), intent(in) :: undf_w0, undf_w3, undf_chi
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_w0), intent(in) :: map_w0
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi

  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in) :: w3_basis
  real(kind=r_def), dimension(1,ndf_w0,nqp_h,nqp_v), intent(in) :: w0_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w3), intent(inout) :: l_exner
  real(kind=r_def), dimension(undf_w0), intent(in)    :: theta
  real(kind=r_def), dimension(undf_w3), intent(in)    :: rho, exner
  real(kind=r_def), dimension(undf_w0), intent(in)    :: theta_ref
  real(kind=r_def), dimension(undf_w3), intent(in)    :: rho_ref
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi1, chi2, chi3

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer(kind=i_def) :: df, k 
  integer(kind=i_def) :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_w0)  :: theta_e, theta_ref_e
  real(kind=r_def), dimension(ndf_chi) :: chi1_e, chi2_e, chi3_e
  real(kind=r_def), dimension(ndf_w3)  :: rho_e, exner_e, rho_ref_e
  real(kind=r_def), dimension(ndf_w3)  :: lhs_exner_e
  real(kind=r_def)                     :: rho_quad, theta_quad, exner_quad
  real(kind=r_def)                     :: rho_ref_quad, theta_ref_quad, exner_ref_quad
  real(kind=r_def)                             :: integrand, eos
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac

  do k = 0, nlayers-1
    do df = 1, ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                             chi_diff_basis, jac, dj)
    do df = 1, ndf_w0
      theta_e(df) = theta(map_w0(df) + k)
      theta_ref_e(df) = theta_ref(map_w0(df) + k)
    end do
    do df = 1, ndf_w3
      exner_e(df)  = exner(map_w3(df) + k)
      rho_e(df)    = rho(map_w3(df) + k)
      rho_ref_e(df)    = rho_ref(map_w3(df) + k)
      lhs_exner_e(df) = 0.0_r_def
    end do
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        theta_quad = 0.0_r_def
        theta_ref_quad = 0.0_r_def
        do df = 1, ndf_w0
          theta_quad = theta_quad + theta_e(df)*w0_basis(1,df,qp1,qp2)
          theta_ref_quad = theta_ref_quad + theta_ref_e(df)*w0_basis(1,df,qp1,qp2)
        end do
        exner_quad = 0.0_r_def
        rho_quad = 0.0_r_def
        rho_ref_quad = 0.0_r_def
        do df = 1, ndf_w3
          exner_quad   = exner_quad   + exner_e(df)  *w3_basis(1,df,qp1,qp2)
          rho_quad     = rho_quad     + rho_e(df)    *w3_basis(1,df,qp1,qp2)
          rho_ref_quad = rho_ref_quad + rho_ref_e(df)*w3_basis(1,df,qp1,qp2)
        end do
        ! Recompute exner_ref from rho_ref & theta_ref
        exner_ref_quad = calc_exner_pointwise(rho_ref_quad, theta_ref_quad)
        eos = (1.0_r_def - kappa)/kappa * exner_quad/exner_ref_quad &
            - rho_quad/rho_ref_quad - theta_quad/theta_ref_quad
        do df = 1, ndf_w3          
          integrand = wqp_h(qp1)*wqp_v(qp2)*w3_basis(1,df,qp1,qp2)*eos*dj(qp1,qp2)
          lhs_exner_e(df) = lhs_exner_e(df) + integrand
        end do
      end do
    end do
    do df = 1, ndf_w3
      l_exner( map_w3(df) + k ) =  lhs_exner_e(df)
    end do 
  end do
  
end subroutine lhs_exner_code

end module lhs_exner_kernel_mod
