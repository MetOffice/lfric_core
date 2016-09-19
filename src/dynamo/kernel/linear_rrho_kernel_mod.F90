!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes rhs of the continuity equation for the linear equations 
!!        with no advection

!> @details The kernel computes thr rhs of the continuity equation 
!>         That is: \f[ r_{\rho} = N^2/g * \rho_s * u.k - \rho_s \nabla.(u) \f]

!>@deprecated The Usefulness of the linear model is to be revaluated at 
!>            the end of the Gung-Ho project and removied if possible
module linear_rrho_kernel_mod

use argument_mod,                   only : arg_type, func_type,         &
                                           GH_FIELD, GH_READ, GH_WRITE, &
                                           W0, W2, W3, ANY_SPACE_9,     &
                                           GH_BASIS, GH_DIFF_BASIS,     &
                                           CELLS
use constants_mod,                  only : r_def
use idealised_config_mod,           only : test
use initial_temperature_config_mod, only : bvf_square
use kernel_mod,                     only : kernel_type
use planet_config_mod,              only : gravity
use reference_profile_mod,          only : reference_profile

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: linear_rrho_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, W3),                             &
       arg_type(GH_FIELD,   GH_READ,  W2),                             &
       arg_type(GH_FIELD,   GH_READ,  W0),                             &
       arg_type(GH_FIELD*3, GH_READ,  ANY_SPACE_9)                     &
       /)
  type(func_type) :: meta_funcs(4) = (/                                &
       func_type(W3, GH_BASIS),                                        &
       func_type(W2, GH_BASIS, GH_DIFF_BASIS),                         &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS),                &
       func_type(W0, GH_DIFF_BASIS)                                    &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass ::linear_rrho_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface linear_rrho_kernel_type
   module procedure linear_rrho_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public linear_rrho_code
contains

type(linear_rrho_kernel_type) function linear_rrho_kernel_constructor() result(self)
  return
end function linear_rrho_kernel_constructor

!> @brief Computes rhs of the continuity equation for the linear equations
!! @param[in] nlayers Number of layers
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number of (local) unique degrees of freedom
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Basis functions evaluated at quadrature points 
!! @param[inout] r_rho Right hand side of the continuity equation
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2 Number of (local) unique degrees of freedom
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Basis functions evaluated at quadrature points 
!! @param[in] w2_diff_basis Differential basis functions evaluated at quadrature points 
!! @param[in] u Velocity
!! @param[in] ndf_w0 Number of degrees of freedom per cell for w0
!! @param[in] undf_w0 Number of (local) unique degrees of freedom
!! @param[in] map_w0 Dofmap for the cell at the base of the column for w0
!! @param[in] w0_diff_basis Differential of the basis functions for w0 evaluated at quadrature points
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number of (local) unique degrees of freedom for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Basis functions for chi evaluated at  quadrature points 
!! @param[in] chi_diff_basis Differential of the basis functions for chi evaluated at quadrature points  
!! @param[in] phi Geopotential
!! @param[in] chi_1 Physical x coordinate in chi
!! @param[in] chi_2 Physical y coordinate in chi
!! @param[in] chi_3 Physical z coordinate in chi
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine linear_rrho_code(nlayers,                                                  &
                            r_rho, u, phi, chi_1, chi_2, chi_3,                       &
                            ndf_w3, undf_w3, map_w3, w3_basis,                        &
                            ndf_w2, undf_w2, map_w2, w2_basis, w2_diff_basis,         &
                            ndf_w0, undf_w0, map_w0, w0_diff_basis,                   &
                            ndf_chi, undf_chi, map_chi, chi_basis, chi_diff_basis,    &
                            nqp_h, nqp_v, wqp_h, wqp_v         )

  use coordinate_jacobian_mod, only: coordinate_jacobian

  !Arguments
  integer, intent(in) :: nlayers, nqp_h, nqp_v
  integer, intent(in) :: ndf_w0, ndf_w2, ndf_w3, ndf_chi
  integer, intent(in) :: undf_w0, undf_w2, undf_w3, undf_chi
  integer, dimension(ndf_w3), intent(in) :: map_w3
  integer, dimension(ndf_w2), intent(in) :: map_w2
  integer, dimension(ndf_w0), intent(in) :: map_w0
  integer, dimension(ndf_chi), intent(in) :: map_chi

  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in)  :: w3_basis
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in)  :: w2_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_basis
  real(kind=r_def), dimension(1,ndf_w2,nqp_h,nqp_v), intent(in)  :: w2_diff_basis
  real(kind=r_def), dimension(3,ndf_w0,nqp_h,nqp_v), intent(in)  :: w0_diff_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w3), intent(inout) :: r_rho
  real(kind=r_def), dimension(undf_chi), intent(in)   :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_w0), intent(in)    :: phi
  real(kind=r_def), dimension(undf_w2), intent(in)    :: u

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer               :: df, k, loc 
  integer               :: qp1, qp2

  real(kind=r_def), dimension(ndf_w0)          :: phi_e
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(ndf_w2)          :: u_e
  real(kind=r_def), dimension(nqp_h,nqp_v)        :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v)    :: jac
  real(kind=r_def), dimension(ndf_w3) :: rrho_e
  real(kind=r_def) :: rho_s_at_quad, exner_s_at_quad, &
                      theta_s_at_quad, div_u_at_quad, &
                      div_term, buoy_term, vec_term
  real(kind=r_def) :: u_at_quad(3), x_at_quad(3), grad_phi_at_quad(3)


  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_chi
      loc = map_chi(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do
    do df = 1, ndf_w3
      rrho_e(df) = 0.0_r_def
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             chi_diff_basis, jac, dj)
    do df = 1, ndf_w0
      phi_e(df)   = phi( map_w0(df) + k )
    end do

    do df = 1, ndf_w2
      u_e(df) = u( map_w2(df) + k )
    end do
  ! compute the RHS integrated over one cell
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        x_at_quad(:) = 0.0_r_def
        grad_phi_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w0
          grad_phi_at_quad(:) = grad_phi_at_quad(:) &
                              + phi_e(df)*w0_diff_basis(:,df,qp1,qp2)
        end do
        do df = 1, ndf_chi
          x_at_quad(1) = x_at_quad(1) + chi_1_e(df)*chi_basis(1,df,qp1,qp2)
          x_at_quad(2) = x_at_quad(2) + chi_2_e(df)*chi_basis(1,df,qp1,qp2)
          x_at_quad(3) = x_at_quad(3) + chi_3_e(df)*chi_basis(1,df,qp1,qp2)
        end do
        call reference_profile(exner_s_at_quad, rho_s_at_quad, & 
                               theta_s_at_quad, x_at_quad, test)
        u_at_quad(:) = 0.0_r_def
        div_u_at_quad = 0.0_r_def
        do df = 1, ndf_w2
          u_at_quad(:)  = u_at_quad(:)  + u_e(df)*w2_basis(:,df,qp1,qp2)
          div_u_at_quad = div_u_at_quad + u_e(df)*w2_diff_basis(1,df,qp1,qp2)
        end do

        div_term  =  - rho_s_at_quad*div_u_at_quad 
        vec_term = dot_product(u_at_quad,grad_phi_at_quad)/gravity
        buoy_term = bvf_square/gravity*rho_s_at_quad*vec_term

        do df = 1, ndf_w3
          rrho_e(df) = rrho_e(df) + wqp_h(qp1)*wqp_v(qp2)*w3_basis(1,df,qp1,qp2)*( buoy_term + div_term )
        end do
      end do
    end do
    do df = 1, ndf_w3
      r_rho( map_w3(df) + k ) =  rrho_e(df)
    end do 
  end do

end subroutine linear_rrho_code

end module linear_rrho_kernel_mod
