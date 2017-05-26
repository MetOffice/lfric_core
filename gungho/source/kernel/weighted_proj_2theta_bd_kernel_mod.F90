!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Compute the boundary part of the projection operator from the potential temperature space
!!        to the velocity space weighted by the pressure gradient
!> @details Compute the boundary projection operator \f[<v.n,{\Pi}*\gamma>\f]
!!          where v is in W2, gamma is in the potential temperature space and
!!          exner is computed pointwise from the equation of state.
module weighted_proj_2theta_bd_kernel_mod
use constants_mod,           only: r_def, i_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, func_type,            &
                                   GH_OPERATOR, GH_FIELD,          &
                                   GH_READ, GH_INC,                &
                                   Wtheta, W2, W3,                 &
                                   GH_BASIS,GH_DIFF_BASIS,         &
                                   CELLS, GH_QUADRATURE_XYoZ,      &
                                   QUADRATURE_XYoZ
use reference_element_mod,   only : nfaces_h, out_face_normal


implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: weighted_proj_2theta_bd_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                     &
       arg_type(GH_OPERATOR, GH_INC, W2, Wtheta),                         &
       arg_type(GH_FIELD,    GH_READ,  Wtheta),                           &
       arg_type(GH_FIELD,    GH_READ,  W3)                                &
       /)
  type(func_type) :: meta_funcs(3) = (/                                   &
       func_type(W2, GH_BASIS),                                           &
       func_type(Wtheta, GH_BASIS),                                       &
       func_type(W3, GH_BASIS)                                            &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = QUADRATURE_XYoZ
contains
  procedure, nopass :: weighted_proj_2theta_bd_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface weighted_proj_2theta_bd_kernel_type
   module procedure weighted_proj_2theta_bd_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public weighted_proj_2theta_bd_code
contains

type(weighted_proj_2theta_bd_kernel_type) function weighted_proj_2theta_bd_constructor() result(self)
  return
end function weighted_proj_2theta_bd_constructor

!> @brief Compute the weigthed projection from Wtheta to W2
!! @param[in] cell Cell number
!! @param[in] nlayers Number of layers.
!! @param[in] ncell_3d ncell*ndf
!! @param[inout] projection Projection operator to compute
!! @param[in] theta Potential temperature
!! @param[in] rho Density
!! @param[in] ndf_w2 Number of degrees of freedom per cell.
!! @param[in] ndf_wtheta Number of degrees of freedom per cell.
!! @param[in] undf_wtheta Total number of degrees.
!! @param[in] stencil_wtheta_map W2 dofmaps for the stencil
!! @param[in] stencil_wtheta_size Size of the W2 stencil (number of cells)
!! @param[in] ndf_w3 Number of degrees of freedom per cell.
!! @param[in] undf_w3 Total number of degrees.
!! @param[in] stencil_w3_map W3 dofmaps for the stencil
!! @param[in] stencil_w3_size Size of the W3 stencil (number of cells)
!! @param[in] nqp_h_1d Number of quadrature points in a single horizontal direction
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_v Vertical quadrature weights
!! @param[in] basis_wtheta Basis functions evaluated at quadrature points.
!! @param[in] w3_basis_face Basis functions evaluated at gaussian quadrature points on horizontal faces
!! @param[in] wtheta_basis_face Basis functions evaluated at gaussian quadrature points on horizontal faces

subroutine weighted_proj_2theta_bd_code(cell, nlayers, ncell_3d,          &
                                     projection,                          &
                                     theta,                               &
                                     rho,                                 &
                                     ndf_w2,                              &
                                     ndf_wtheta, undf_wtheta,             &
                                     stencil_wtheta_map,                  &
                                     stencil_wtheta_size,                 &
                                     ndf_w3, undf_w3,                     &
                                     stencil_w3_map,                      &
                                     stencil_w3_size,                     &
                                     nqp_h_1d, nqp_v, wqp_v,              &
                                     w2_basis_face,                       &
                                     w3_basis_face,                       &
                                     wtheta_basis_face, adjacent_face )

  use calc_exner_pointwise_mod, only: calc_exner_pointwise

  implicit none
  !Arguments
  integer(kind=i_def),                     intent(in) :: cell, nqp_h_1d, nqp_v
  integer(kind=i_def),                     intent(in) :: nlayers
  integer(kind=i_def),                     intent(in) :: ncell_3d
  integer(kind=i_def),                     intent(in) :: undf_w3, ndf_w3, ndf_w2, ndf_wtheta, undf_wtheta

  integer(kind=i_def), intent(in) :: stencil_w3_size
  integer(kind=i_def), dimension(ndf_w3, stencil_w3_size), intent(in)  :: stencil_w3_map

  integer(kind=i_def), intent(in) :: stencil_wtheta_size
  integer(kind=i_def), dimension(ndf_wtheta, stencil_wtheta_size), intent(in)  :: stencil_wtheta_map

  real(kind=r_def), dimension(4,3,ndf_w2,nqp_h_1d,nqp_v), intent(in)     :: w2_basis_face
  real(kind=r_def), dimension(4,1,ndf_w3,nqp_h_1d,nqp_v), intent(in)     :: w3_basis_face
  real(kind=r_def), dimension(4,1,ndf_wtheta,nqp_h_1d,nqp_v), intent(in) :: wtheta_basis_face

  real(kind=r_def), dimension(ndf_w2,ndf_wtheta,ncell_3d), intent(inout) :: projection
  real(kind=r_def), dimension(undf_wtheta),                intent(in)    :: theta
  real(kind=r_def), dimension(undf_w3),                    intent(in)    :: rho
  real(kind=r_def), dimension(nqp_v),                      intent(in)    :: wqp_v

  integer(kind=i_def), dimension(nfaces_h), intent(in) :: adjacent_face

  !Internal variables
  integer(kind=i_def)                      :: df, df0, df2, k, ik, face, face_next
  integer(kind=i_def)                      :: qp1, qp2
  real(kind=r_def), dimension(ndf_wtheta)  :: theta_e, theta_next_e
  real(kind=r_def), dimension(ndf_w3)      :: rho_e, rho_next_e

  real(kind=r_def)                         :: v(3), integrand
  real(kind=r_def)                         :: theta_at_fquad, theta_next_at_fquad
  real(kind=r_def)                         :: rho_at_fquad, rho_next_at_fquad
  real(kind=r_def)                         :: exner_at_fquad, exner_next_at_fquad

  do k = 0, nlayers - 1
    ik = k + 1 + (cell-1)*nlayers
    do face = 1, nfaces_h

      face_next = adjacent_face(face)

      do df = 1,ndf_wtheta
        theta_e(df)           = theta(stencil_wtheta_map(df, 1)      + k)
        theta_next_e(df)      = theta(stencil_wtheta_map(df, face+1)      + k)
      end do
      do df = 1,ndf_w3
        rho_e(df)           = rho(stencil_w3_map(df, 1) + k)
        rho_next_e(df)      = rho(stencil_w3_map(df, face+1) + k)
      end do

      do qp2 = 1, nqp_v
        do qp1 = 1, nqp_h_1d
           theta_at_fquad      = 0.0_r_def
           theta_next_at_fquad = 0.0_r_def
           do df = 1, ndf_wtheta
             theta_at_fquad      = theta_at_fquad + theta_e(df)*wtheta_basis_face(face,1,df,qp1,qp2)
             theta_next_at_fquad = theta_next_at_fquad + theta_next_e(df)*wtheta_basis_face(face_next,1,df,qp1,qp2)
           end do
           rho_at_fquad      = 0.0_r_def
           rho_next_at_fquad = 0.0_r_def
           do df = 1, ndf_w3
             rho_at_fquad      = rho_at_fquad + rho_e(df)*w3_basis_face(face,1,df,qp1,qp2)
             rho_next_at_fquad = rho_next_at_fquad + rho_next_e(df)*w3_basis_face(face_next,1,df,qp1,qp2)
           end do
           exner_at_fquad      = calc_exner_pointwise(rho_at_fquad, theta_at_fquad)
           exner_next_at_fquad = calc_exner_pointwise(rho_next_at_fquad, theta_next_at_fquad)


           do df0 = 1, ndf_wtheta
             do df2 = 1, ndf_w2
              v  = w2_basis_face(face,:,df2,qp1,qp2)

              integrand = wqp_v(qp1)*wqp_v(qp2)* &
                            0.5_r_def*(exner_at_fquad + exner_next_at_fquad)* &
                              dot_product(v, out_face_normal(:,face))*wtheta_basis_face(face,1,df0,qp1,qp2)
              projection(df2,df0,ik) = projection(df2,df0,ik) - integrand
            end do
          end do
        end do
      end do
    end do
  end do
  end subroutine weighted_proj_2theta_bd_code

end module weighted_proj_2theta_bd_kernel_mod
