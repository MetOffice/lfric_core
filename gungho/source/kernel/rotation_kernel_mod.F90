!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes rotation component of the rhs of the momentum equation for the nonlinear equations


!> @details The kernel computes the rotation component rhs of the momentum equation as given by
!>         \f[ 2\Omega \times u \f]
!>         Omega is the rotation vector of the domain and u is the fluid velocity
module rotation_kernel_mod

use argument_mod,      only : arg_type, func_type,                     &
                              GH_FIELD, GH_READ, GH_INC,               &
                              ANY_SPACE_9, W2,                         &
                              GH_BASIS, GH_DIFF_BASIS,                 &
                              CELLS, GH_QUADRATURE_XYoZ,               &
                              QUADRATURE_XYoZ
use constants_mod,     only : r_def
use cross_product_mod, only : cross_product
use kernel_mod,        only : kernel_type
use planet_config_mod, only : scaled_omega

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: rotation_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W2),                              &
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9)                      &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(W2, GH_BASIS),                                        &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                 &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = QUADRATURE_XYoZ
contains
  procedure, nopass ::rotation_code
end type

!-------------------------------------------------------------------------------
! Constrotationctors
!-------------------------------------------------------------------------------

! overload the default strotationcture constrotationctor for function space
interface rotation_kernel_type
   module procedure rotation_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public rotation_code
contains

type(rotation_kernel_type) function rotation_kernel_constructor() result(self)
  return
end function rotation_kernel_constructor

!> @brief Compute the rotational part of the momentum equation
!! @param[in] nlayers Number of layers
!! @param[inout] r_u Right hand side of the momentum equation
!! @param[in] u Velocity
!! @param[in] chi_1 Physical x coordinate in w0
!! @param[in] chi_2 Physical y coordinate in w0
!! @param[in] chi_3 Physical z coordinate in w0
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2 Number unique of degrees of freedom  for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Basis functions evaluated at quadrature points
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number unique of degrees of freedom  for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] chi_diff_basis Differential basis functions evaluated at gaussian quadrature points
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine rotation_code(nlayers,                                              &
                         r_u, u,                                               &
                         chi_1, chi_2, chi_3,                                  &
                         ndf_w2, undf_w2, map_w2, w2_basis,                    &
                         ndf_chi, undf_chi, map_chi, chi_basis, chi_diff_basis, &
                         nqp_h, nqp_v, wqp_h, wqp_v                            &
                         )

  use base_mesh_config_mod,    only: geometry,                     &
                                     base_mesh_geometry_spherical, &
                                     f_lat
  use coordinate_jacobian_mod, only: coordinate_jacobian
  use rotation_vector_mod,     only: rotation_vector_fplane,  &
                                     rotation_vector_sphere

  !Arguments
  integer, intent(in) :: nlayers, nqp_h, nqp_v
  integer, intent(in) :: ndf_chi, ndf_w2
  integer, intent(in) :: undf_chi, undf_w2
  integer, dimension(ndf_chi), intent(in) :: map_chi
  integer, dimension(ndf_w2), intent(in) :: map_w2

  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w2), intent(inout) :: r_u
  real(kind=r_def), dimension(undf_w2), intent(in)    :: u
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer               :: df, k, loc
  integer               :: qp1, qp2

  real(kind=r_def), dimension(ndf_chi)          :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)         :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v)     :: jac
  real(kind=r_def), dimension(ndf_w2)              :: ru_e, u_e
  real(kind=r_def), dimension(3,nqp_h,nqp_v)       :: rotation_vector

  real(kind=r_def) :: u_at_quad(3), jac_v(3), v(3), jac_u(3)
  real(kind=r_def) :: omega_cross_u(3)
  real(kind=r_def) :: coriolis_term

  do k = 0, nlayers-1

    ! Extract element arrays of chi
    do df = 1, ndf_chi
      loc = map_chi(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do

    ! Calculate rotation vector Omega = (0, 2*cos(lat), 2*sin(lat)) and Jacobian
    if ( geometry == base_mesh_geometry_spherical ) then
      call rotation_vector_sphere(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e,      &
                              chi_3_e, chi_basis, rotation_vector)
    else
      call rotation_vector_fplane(nqp_h, nqp_v, scaled_omega, f_lat, rotation_vector)
    end if
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             chi_diff_basis, jac, dj)

    do df = 1, ndf_w2
      u_e(df) = u( map_w2(df) + k )
      ru_e(df) = 0.0_r_def
    end do

    ! Compute the rotation component of RHS integrated over one cell
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h

        u_at_quad(:) = 0.0_r_def

        do df = 1, ndf_w2
          u_at_quad(:) = u_at_quad(:) &
                       + u_e(df)*w2_basis(:,df,qp1,qp2)
        end do

        ! Rotation term
        jac_u = matmul(jac(:,:,qp1,qp2),u_at_quad)
        omega_cross_u = cross_product(rotation_vector(:,qp1,qp2), jac_u)
        do df = 1, ndf_w2
          v  = w2_basis(:,df,qp1,qp2)
          jac_v = matmul(jac(:,:,qp1,qp2),v)/dj(qp1,qp2)

          coriolis_term = dot_product(jac_v,omega_cross_u)

          ru_e(df) = ru_e(df) -  wqp_h(qp1)*wqp_v(qp2)*coriolis_term
        end do

      end do
    end do

    do df = 1, ndf_w2
      r_u( map_w2(df) + k ) =  r_u( map_w2(df) + k ) + ru_e(df)
    end do

  end do

end subroutine rotation_code

end module rotation_kernel_mod
