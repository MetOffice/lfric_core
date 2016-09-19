!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes rhs of the momentum equation for the linear equations with
!>         no advection


!> @details The kernel computes thr rhs of the momentum equation
!>         That is: \f[ r_u = -Cp * \theta_s * \nabla(\Pi) + \theta/\theta_s * g \f]

!>@deprecated The Usefulness of the linear model is to be revaluated at 
!>            the end of the Gung-Ho project and removied if possible
module linear_ru_kernel_mod

use base_mesh_config_mod,           only : f_lat,    &
                                           geometry, &
                                           base_mesh_geometry_spherical
use argument_mod,                   only : arg_type, func_type,              &
                                           GH_FIELD, GH_READ, GH_INC,        &
                                           W0, W2, W3, ANY_SPACE_9,          &
                                           GH_BASIS, GH_DIFF_BASIS, CELLS
use constants_mod,                  only : r_def
use idealised_config_mod,           only : test
use initial_temperature_config_mod, only : bvf_square
use kernel_mod,                     only : kernel_type
use planet_config_mod,              only : gravity, scaled_omega, Cp

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: linear_ru_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W2),                              &
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD,   GH_READ, W3),                              &
       arg_type(GH_FIELD,   GH_READ, W0),                              &
       arg_type(GH_FIELD,   GH_READ, W0),                              &
       arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9)                      &
       /)
  type(func_type) :: meta_funcs(4) = (/                                &
       func_type(W2, GH_BASIS, GH_DIFF_BASIS),                         &
       func_type(W3, GH_BASIS),                                        &
       func_type(W0, GH_BASIS, GH_DIFF_BASIS),                         &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                 &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass ::linear_ru_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface linear_ru_kernel_type
   module procedure linear_ru_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public linear_ru_code
contains

type(linear_ru_kernel_type) function linear_ru_kernel_constructor() result(self)
  return
end function linear_ru_kernel_constructor

!> @brief Compute the right hand side of the linear momentum equation
!! @param[in] nlayers Number of layers
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2 Number unique of degrees of freedom  for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Basis functions evaluated at quadrature points 
!! @param[in] w2_diff_basis Differntial of the basis functions evaluated at  quadrature points
!! @param[inout] r_u Right hand side of the linear momentum equation
!! @param[in] u Velocity
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number unique of degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Basis functions evaluated at gaussian quadrature points 
!! @param[in] rho Density
!! @param[in] ndf_w0 Number of degrees of freedom per cell for w0
!! @param[in] undf_w0 Number unique of degrees of freedom  for w0
!! @param[in] map_w0 Dofmap for the cell at the base of the column for w0
!! @param[in] w0_basis Basis functions evaluated at gaussian quadrature points 
!! @param[in] w0_diff_basis Differential of the basis functions evaluated at gaussian quadrature point
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number unique of degrees of freedom  for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Basis functions evaluated at gaussian quadrature points 
!! @param[in] chi_diff_basis Differntial of the basis functions evaluated at gaussian quadrature point
!! @param[in] theta Potential temperature
!! @param[in] phi Geopotential
!! @param[in] chi_1 Physical x coordinate in chi
!! @param[in] chi_2 Physical y coordinate in chi
!! @param[in] chi_3 Physical z coordinate in chi
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine linear_ru_code(nlayers, r_u, u, rho, theta, phi, chi_1, chi_2, chi_3, &
                          ndf_w2, undf_w2, map_w2, w2_basis, w2_diff_basis,      &
                          ndf_w3, undf_w3, map_w3, w3_basis,                     &
                          ndf_w0, undf_w0, map_w0, w0_basis, w0_diff_basis,      &
                          ndf_chi, undf_chi, map_chi, chi_basis, chi_diff_basis, &
                          nqp_h, nqp_v, wqp_h, wqp_v                             &
                         )

  use coordinate_jacobian_mod,  only: coordinate_jacobian
  use reference_profile_mod,    only: reference_profile 
  use calc_exner_pointwise_mod, only: linear_calc_exner_pointwise
  use coord_transform_mod,      only: xyz2llr, sphere2cart_vector
  use rotation_vector_mod,      only: rotation_vector_fplane,  &
                                      rotation_vector_sphere
  use cross_product_mod,        only: cross_product

  !Arguments
  integer, intent(in) :: nlayers,nqp_h, nqp_v
  integer, intent(in) :: ndf_w0, ndf_w2, ndf_w3, ndf_chi
  integer, intent(in) :: undf_w0, undf_w2, undf_w3, undf_chi
  integer, dimension(ndf_w0), intent(in) :: map_w0
  integer, dimension(ndf_w2), intent(in) :: map_w2
  integer, dimension(ndf_w3), intent(in) :: map_w3
  integer, dimension(ndf_chi), intent(in) :: map_chi
  
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in)  :: w3_basis
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in)  :: w2_basis
  real(kind=r_def), dimension(1,ndf_w0,nqp_h,nqp_v), intent(in)  :: w0_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_basis
  real(kind=r_def), dimension(1,ndf_w2,nqp_h,nqp_v), intent(in)  :: w2_diff_basis
  real(kind=r_def), dimension(3,ndf_w0,nqp_h,nqp_v), intent(in)  :: w0_diff_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w2), intent(inout) :: r_u
  real(kind=r_def), dimension(undf_w2), intent(in)    :: u
  real(kind=r_def), dimension(undf_w3), intent(in)    :: rho
  real(kind=r_def), dimension(undf_chi), intent(in)   :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_w0), intent(in)    :: theta, phi

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer               :: df, k, loc 
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_chi)          :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(3,nqp_h,nqp_v)   :: rotation_vector
  real(kind=r_def), dimension(ndf_w3)          :: rho_e
  real(kind=r_def), dimension(ndf_w0)          :: theta_e, phi_e
  real(kind=r_def), dimension(ndf_w2)          :: ru_e
  real(kind=r_def) :: grad_theta_s_at_quad(3), jac_v(3), x_at_quad(3), &
                      grad_phi_at_quad(3), v(3), omega_cross_u(3), &
                      u_at_quad(3), j_u(3)
  real(kind=r_def) :: exner_at_quad, rho_at_quad, theta_at_quad,       &
                      exner_s_at_quad, rho_s_at_quad, theta_s_at_quad, &
                      grad_term, buoy_term, coriolis_term

  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_w0
      loc = map_w0(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do
! Calculate rotation and Jacobian
    if ( geometry == base_mesh_GEOMETRY_SPHERICAL ) then
      call rotation_vector_sphere(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e, &
                                  chi_basis, rotation_vector)
    else
      call rotation_vector_fplane(nqp_h, nqp_v, scaled_omega, f_lat, rotation_vector)
    end if
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             chi_diff_basis, jac, dj)
    do df = 1, ndf_w0
      phi_e(df)   = phi( map_w0(df) + k )
    end do
    do df = 1, ndf_w3
      rho_e(df) = rho( map_w3(df) + k )
    end do
    do df = 1, ndf_w0
      theta_e(df) = theta( map_w0(df) + k )
    end do
  ! compute the RHS integrated over one cell
    ru_e(:) = 0.0_r_def
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        rho_at_quad = 0.0_r_def 
        do df = 1, ndf_w3
          rho_at_quad  = rho_at_quad + rho_e(df)*w3_basis(1,df,qp1,qp2) 
        end do
        theta_at_quad = 0.0_r_def
        grad_theta_s_at_quad(:) = 0.0_r_def
        grad_phi_at_quad(:) = 0.0
        do df = 1, ndf_w0
          theta_at_quad   = theta_at_quad                                      &
                          + theta_e(df)*w0_basis(1,df,qp1,qp2)
          grad_phi_at_quad(:) = grad_phi_at_quad(:) &
                              + phi_e(df)*w0_diff_basis(:,df,qp1,qp2)
        end do
        x_at_quad(:) = 0.0_r_def
        do df = 1, ndf_chi
          x_at_quad(1) = x_at_quad(1) + chi_1_e(df)*chi_basis(1,df,qp1,qp2)
          x_at_quad(2) = x_at_quad(2) + chi_2_e(df)*chi_basis(1,df,qp1,qp2)
          x_at_quad(3) = x_at_quad(3) + chi_3_e(df)*chi_basis(1,df,qp1,qp2)
        end do
        u_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w2
          u_at_quad(:) = u_at_quad(:) + u(map_w2(df)+k)*w2_basis(:,df,qp1,qp2)
        end do
        j_u = matmul(jac(:,:,qp1,qp2),u_at_quad)
        omega_cross_u = cross_product(rotation_vector(:,qp1,qp2), j_u)

        call reference_profile(exner_s_at_quad, rho_s_at_quad, &
                               theta_s_at_quad, x_at_quad, test)

        exner_at_quad = linear_calc_exner_pointwise( &
                               rho_at_quad, theta_at_quad,       & 
                               exner_s_at_quad, rho_s_at_quad,   & 
                               theta_s_at_quad)

        do df = 1, ndf_w2
          v = w2_basis(:,df,qp1,qp2)
          jac_v = matmul(jac(:,:,qp1,qp2),w2_basis(:,df,qp1,qp2))
          buoy_term = dot_product( v, grad_phi_at_quad ) &
                    *(theta_at_quad/theta_s_at_quad)
          grad_term = Cp*exner_at_quad*theta_s_at_quad*( &
                      bvf_square/gravity*dot_product( v, grad_phi_at_quad ) + &
                      w2_diff_basis(1,df,qp1,qp2) )

          coriolis_term = dot_product(jac_v/dj(qp1,qp2),omega_cross_u)

          ru_e(df) = ru_e(df) +  wqp_h(qp1)*wqp_v(qp2)*( grad_term + buoy_term - coriolis_term )
        end do
      end do
    end do
    do df = 1, ndf_w2
      r_u( map_w2(df) + k ) =  r_u( map_w2(df) + k ) + ru_e(df)
    end do 
  end do 

end subroutine linear_ru_code

end module linear_ru_kernel_mod
