!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Compute the projection operator from the potential temperature space 
!!        to the velocity space weighted by the pressure gradient
!> @details Compute the projection operator \f[<v,\nabla(\Pi)*\gamma>\f]
!!          where v is in W2 and gamma is in the potential temperature space and
!!          exner is computed pointwise from the equation of state.
!!          This is integrated by parts to give
!!          \f[
!!          \left< \nabla.v,\Pi\gamma\right> + \left< v,\Pi\nabla\gamma\right>
!!          \f]
module weighted_proj_2theta_kernel_mod
use constants_mod,           only: r_def, i_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, func_type,            &
                                   GH_OPERATOR, GH_FIELD,          &
                                   GH_READ, GH_WRITE,              &
                                   ANY_SPACE_9, W2, W3,            &
                                   GH_BASIS,GH_DIFF_BASIS,         &
                                   CELLS, GH_QUADRATURE_XYoZ,      &
                                   QUADRATURE_XYoZ

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: weighted_proj_2theta_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_OPERATOR, GH_WRITE, W2, ANY_SPACE_9),               &
       arg_type(GH_FIELD,    GH_READ,  ANY_SPACE_9),                   &
       arg_type(GH_FIELD,    GH_READ,  W3)                             &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(W2, GH_BASIS, GH_DIFF_BASIS),                         &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS),                &
       func_type(W3, GH_BASIS)                                         &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = QUADRATURE_XYoZ
contains
  procedure, nopass :: weighted_proj_2theta_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface weighted_proj_2theta_kernel_type
   module procedure weighted_proj_2theta_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public weighted_proj_2theta_code
contains

type(weighted_proj_2theta_kernel_type) function weighted_proj_2theta_constructor() result(self)
  return
end function weighted_proj_2theta_constructor

!> @brief Compute the weigthed projection from Wtheta to W2
!! @param[in] cell Cell number
!! @param[in] nlayers Number of layers.
!! @param[in] ncell_3d ncell*ndf
!! @param[inout] projection Projection operator to compute
!! @param[in] theta Potential temperature
!! @param[in] rho Density
!! @param[in] ndf_w3 Number of degrees of freedom per cell.
!! @param[in] undf_w3 Total number of degrees.
!! @param[in] map_w3 Dofmap at the base of the column
!! @param[in] basis_w3 Basis functions evaluated at quadrature points.
!! @param[in] ndf_w2 Number of degrees of freedom per cell.
!! @param[in] basis_w2 Basis functions evaluated at quadrature points.
!! @param[in] diff_basis_w2 Differential vector basis functions evaluated at quadrature points.
!! @param[in] ndf_wtheta Number of degrees of freedom per cell.
!! @param[in] undf_wtheta Total number of degrees.
!! @param[in] map_wtheta Dofmap at the base of the column
!! @param[in] basis_wtheta Basis functions evaluated at quadrature points.
!! @param[in] diff_basis_wtheta Differential vector basis functions evaluated at quadrature points.
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine weighted_proj_2theta_code(cell, nlayers, ncell_3d,             &
                                     projection,                          &
                                     theta,                               &
                                     rho,                                 &
                                     ndf_w2, basis_w2, diff_basis_w2,     &
                                     ndf_wtheta, undf_wtheta, map_wtheta, &
                                     basis_wtheta, diff_basis_wtheta,     &
                                     ndf_w3, undf_w3, map_w3, basis_w3,   &
                                     nqp_h, nqp_v, wqp_h, wqp_v )

  use calc_exner_pointwise_mod, only: calc_exner_pointwise

  implicit none
  !Arguments
  integer(kind=i_def),                     intent(in) :: cell, nqp_h, nqp_v
  integer(kind=i_def),                     intent(in) :: nlayers
  integer(kind=i_def),                     intent(in) :: ncell_3d
  integer(kind=i_def),                     intent(in) :: undf_w3, ndf_w3, ndf_w2, ndf_wtheta, undf_wtheta
  integer(kind=i_def), dimension(ndf_wtheta),  intent(in) :: map_wtheta
  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3

  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v),   intent(in) :: basis_w3
  real(kind=r_def), dimension(1,ndf_w2,nqp_h,nqp_v),   intent(in) :: diff_basis_w2
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v),   intent(in) :: basis_w2
  real(kind=r_def), dimension(1,ndf_wtheta,nqp_h,nqp_v),   intent(in) :: basis_wtheta
  real(kind=r_def), dimension(3,ndf_wtheta,nqp_h,nqp_v),   intent(in) :: diff_basis_wtheta

  real(kind=r_def), dimension(ndf_w2,ndf_wtheta,ncell_3d), intent(inout) :: projection
  real(kind=r_def), dimension(undf_wtheta),                intent(in)    :: theta
  real(kind=r_def), dimension(undf_w3),                intent(in)    :: rho
  real(kind=r_def), dimension(nqp_h),                  intent(in)    :: wqp_h
  real(kind=r_def), dimension(nqp_v),                  intent(in)    :: wqp_v

  !Internal variables
  integer(kind=i_def)                  :: df, df0, df2, k, ik
  integer(kind=i_def)                  :: qp1, qp2
  real(kind=r_def), dimension(ndf_wtheta)  :: theta_e
  real(kind=r_def), dimension(ndf_w3)  :: rho_e
  real(kind=r_def)                     :: integrand
  real(kind=r_def)                     :: div_gamma_v
  real(kind=r_def)                     :: theta_quad, rho_quad, exner_quad
  real(kind=r_def)                     :: i1

  do k = 0, nlayers - 1
    ik = k + 1 + (cell-1)*nlayers
    do df = 1,ndf_wtheta
      theta_e(df) = theta(map_wtheta(df) + k)
    end do
    do df = 1,ndf_w3
      rho_e(df) = rho(map_w3(df) + k)
    end do
    projection(:,:,ik) = 0.0_r_def
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        theta_quad = 0.0_r_def            
        do df = 1, ndf_wtheta
          theta_quad = theta_quad                                      &
                     + theta_e(df)*basis_wtheta(1,df,qp1,qp2)
        end do
        rho_quad = 0.0_r_def
        do df = 1, ndf_w3
          rho_quad = rho_quad                                      &
                   + rho_e(df)*basis_w3(1,df,qp1,qp2)
        end do
        exner_quad = calc_exner_pointwise(rho_quad, theta_quad)
        i1 = exner_quad* wqp_h(qp1)*wqp_v(qp2)
        do df0 = 1, ndf_wtheta
          do df2 = 1, ndf_w2
            div_gamma_v = diff_basis_w2(1,df2,qp1,qp2)*basis_wtheta(1,df0,qp1,qp2) &
                        + dot_product(basis_w2(:,df2,qp1,qp2), &
                                      diff_basis_wtheta(:,df0,qp1,qp2))
            integrand = i1*div_gamma_v 
            projection(df2,df0,ik) = projection(df2,df0,ik) + integrand
          end do
        end do
      end do
    end do
  end do 
end subroutine weighted_proj_2theta_code

end module weighted_proj_2theta_kernel_mod
