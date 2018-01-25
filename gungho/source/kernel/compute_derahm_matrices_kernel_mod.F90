!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------

!> @brief Computes the mass matrices of the DeRahm complex
!>
!> @details Compute the mass matrices of the DeRahm cochain,
!> these are the mass matrices for the W0, W1, W2 & W3 function spaces
!> along with the gradient, curl and divergence operators and the 
!> Wtheta mass matrix.
!> Since they all depend upon the mesh Jacobian they are computed as one
!> to reduce the cost
module compute_derahm_matrices_kernel_mod
use constants_mod,           only: r_def, i_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, func_type,            &
                                   GH_OPERATOR, GH_FIELD,          &
                                   GH_READ, GH_WRITE,              &
                                   ANY_SPACE_9, W0, W1, W2, W3,    &
                                   Wtheta,                         &
                                   GH_BASIS, GH_DIFF_BASIS,        &
                                   CELLS, GH_QUADRATURE_XYoZ
use coordinate_jacobian_mod, only: pointwise_coordinate_jacobian, &
                                   pointwise_coordinate_jacobian_inverse
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: compute_derahm_matrices_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                                  &
       arg_type(GH_OPERATOR, GH_WRITE, W0, W0),                        &
       arg_type(GH_OPERATOR, GH_WRITE, W1, W1),                        &
       arg_type(GH_OPERATOR, GH_WRITE, W2, W2),                        &
       arg_type(GH_OPERATOR, GH_WRITE, W3, W3),                        &
       arg_type(GH_OPERATOR, GH_WRITE, Wtheta, Wtheta),                &
       arg_type(GH_OPERATOR, GH_WRITE, W1, W0),                        &
       arg_type(GH_OPERATOR, GH_WRITE, W2, W1),                        &
       arg_type(GH_OPERATOR, GH_WRITE, W3, W2),                        &
       arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_9)                    &
       /)
  type(func_type) :: meta_funcs(6) = (/                                &
       func_type(W0,         GH_BASIS, GH_DIFF_BASIS),                 &
       func_type(W1,         GH_BASIS, GH_DIFF_BASIS),                 &
       func_type(W2,         GH_BASIS, GH_DIFF_BASIS),                 &
       func_type(W3,         GH_BASIS),                                &
       func_type(Wtheta,     GH_BASIS),                                &
       func_type(ANY_SPACE_9,          GH_DIFF_BASIS)                  &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: compute_derahm_matrices_code
end type
!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface compute_derahm_matrices_kernel_type
   module procedure compute_derahm_matrices_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public compute_derahm_matrices_code
contains

type(compute_derahm_matrices_kernel_type) &
                        function compute_derahm_matrices_constructor() result(self)
  return
end function compute_derahm_matrices_constructor

!> @brief This subroutine computes the operator matrices for the modified DeRahm complex
!! @param[in] cell Cell number
!! @param[in] nlayers Number of layers.
!! @param[in] ncell_3d0 ncell*nlayers
!! @param[inout] mm0 Local assembly mass matrix for W0 space
!! @param[in] ncell_3d1 ncell*nlayers
!! @param[inout] mm1 Local assembly mass matrix for W1 space
!! @param[in] ncell_3d2 ncell*nlayers
!! @param[inout] mm2 Local assembly mass matrix for W2 space
!! @param[in] ncell_3d3 ncell*nlayers
!! @param[inout] mm3 Local assembly mass matrix for W3 space
!! @param[in] ncell_3dt ncell*nlayers
!! @param[inout] mmt Local assembly mass matrix for Wtheta space
!! @param[in] ncell_3d4 ncell*nlayers
!! @param[inout] grad Local assembly gradient matrix
!! @param[in] ncell_3d5 ncell*nlayers
!! @param[inout] curl Local assembly curl matrix
!! @param[in] ncell_3d6 ncell*nlayers
!! @param[inout] div Local assembly divergence matrix
!! @param[in] chi1 Physical coordinates in the 1st dir
!! @param[in] chi2 Physical coordinates in the 2nd dir
!! @param[in] chi3 Physical coordinates in the 3rd dir
!! @param[in] ndf_w0 Number of degrees of freedom per cell for W0 space.
!! @param[in] basis_w0 Basis functions evaluated at quadrature points for W0 space.
!! @param[in] diff_basis_w0 Differential of basis functions evaluated at quadrature points for W0 space.
!! @param[in] ndf_w1 Number of degrees of freedom per cell for W1 space.
!! @param[in] basis_w1 Basis functions evaluated at quadrature points for W1 space.
!! @param[in] diff_basis_w1 Differential of basis functions evaluated at quadrature points for W1 space.
!! @param[in] ndf_w2 Number of degrees of freedom per cell for W2 space.
!! @param[in] basis_w2 Basis functions evaluated at quadrature points for W2 space.
!! @param[in] diff_basis_w2 Differential of basis functions evaluated at quadrature points for W2 space.
!! @param[in] ndf_w3 Number of degrees of freedom per cell for W3 space.
!! @param[in] basis_w3 Basis functions evaluated at quadrature points for W3 space.
!! @param[in] ndf_wt Number of degrees of freedom per cell for Wtheta space.
!! @param[in] basis_wt Basis functions evaluated at quadrature points for Wtheta space.
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi field
!! @param[in] undf_chi Number of unique degrees of freedom for chi field
!! @param[in] map_chi Dofmap for the cell at the base of the column, for the
!!                    space on which the chi field lives
!! @param[in] diff_basis_chi Vector differential basis functions evaluated at
!!                           quadrature points.
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine compute_derahm_matrices_code(cell, nlayers,                    & 
                                        ncell_3d0, mm0,                   &
                                        ncell_3d1, mm1,                   &
                                        ncell_3d2, mm2,                   &
                                        ncell_3d3, mm3,                   &
                                        ncell_3dt, mmt,                   &
                                        ncell_3d4, grad,                  &
                                        ncell_3d5, curl,                  &
                                        ncell_3d6, div,                   &
                                        chi1, chi2, chi3,                 &
                                        ndf_w0, basis_w0, diff_basis_w0,  &
                                        ndf_w1, basis_w1, diff_basis_w1,  &
                                        ndf_w2, basis_w2, diff_basis_w2,  &
                                        ndf_w3, basis_w3,                 &
                                        ndf_wt, basis_wt,                 &
                                        ndf_chi, undf_chi,                &
                                        map_chi, diff_basis_chi,          &
                                        nqp_h, nqp_v, wqp_h, wqp_v)

  implicit none
  ! Arguments
  integer(kind=i_def),   intent(in)     :: cell, nqp_h, nqp_v
  integer(kind=i_def),   intent(in)     :: nlayers 
  integer(kind=i_def),   intent(in)     :: ndf_w0, ndf_w1, ndf_w2, ndf_w3, ndf_wt
  integer(kind=i_def),   intent(in)     :: ncell_3d0, ncell_3d1, ncell_3d2, &
                                           ncell_3d3, ncell_3d4, ncell_3d5, &
                                           ncell_3d6, ncell_3dt
  integer(kind=i_def),   intent(in)     :: ndf_chi
  integer(kind=i_def),   intent(in)     :: undf_chi
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  real(kind=r_def), intent(inout) :: mm0(ndf_w0,ndf_w0,ncell_3d0)
  real(kind=r_def), intent(inout) :: mm1(ndf_w1,ndf_w1,ncell_3d1)
  real(kind=r_def), intent(inout) :: mm2(ndf_w2,ndf_w2,ncell_3d2)
  real(kind=r_def), intent(inout) :: mm3(ndf_w3,ndf_w3,ncell_3d3)
  real(kind=r_def), intent(inout) :: mmt(ndf_wt,ndf_wt,ncell_3dt)
  real(kind=r_def), intent(inout) :: grad(ndf_w1,ndf_w0,ncell_3d4)
  real(kind=r_def), intent(inout) :: curl(ndf_w2,ndf_w1,ncell_3d5)
  real(kind=r_def), intent(inout) :: div(ndf_w3,ndf_w2,ncell_3d6)
  real(kind=r_def), intent(in)    :: diff_basis_chi(3,ndf_chi,nqp_h,nqp_v)
  real(kind=r_def), intent(in)    :: basis_w0(1,ndf_w0,nqp_h,nqp_v)
  real(kind=r_def), intent(in)    :: basis_w1(3,ndf_w1,nqp_h,nqp_v)
  real(kind=r_def), intent(in)    :: basis_w2(3,ndf_w2,nqp_h,nqp_v)
  real(kind=r_def), intent(in)    :: basis_w3(1,ndf_w3,nqp_h,nqp_v)
  real(kind=r_def), intent(in)    :: basis_wt(1,ndf_wt,nqp_h,nqp_v)
  real(kind=r_def), intent(in)    :: diff_basis_w0(3,ndf_w0,nqp_h,nqp_v)
  real(kind=r_def), intent(in)    :: diff_basis_w1(3,ndf_w1,nqp_h,nqp_v)
  real(kind=r_def), intent(in)    :: diff_basis_w2(1,ndf_w2,nqp_h,nqp_v)
  real(kind=r_def), intent(in)    :: chi1(undf_chi)
  real(kind=r_def), intent(in)    :: chi2(undf_chi)
  real(kind=r_def), intent(in)    :: chi3(undf_chi)
  real(kind=r_def), intent(in)    :: wqp_h(nqp_h)
  real(kind=r_def), intent(in)    :: wqp_v(nqp_v)

  ! Internal variables
  integer(kind=i_def)                          :: df, df2, k, ik
  integer(kind=i_def)                          :: qp1, qp2

  real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
  real(kind=r_def)                             :: integrand
  real(kind=r_def)                             :: dj
  real(kind=r_def), dimension(3,3)             :: jac, jac_inv
  real(kind=r_def), dimension(3,3)             :: jac_t
  real(kind=r_def), dimension(3)               :: jac_v
  real(kind=r_def)                             :: wt

  do k = 0, nlayers-1

    ! Indirect the chi coord field here
    do df = 1, ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do

    ik = 1 + k + (cell-1)*nlayers
    mm0(:,:,ik) = 0.0_r_def
    mm1(:,:,ik) = 0.0_r_def
    mm2(:,:,ik) = 0.0_r_def
    mm3(:,:,ik) = 0.0_r_def
    mmt(:,:,ik) = 0.0_r_def
    grad(:,:,ik) = 0.0_r_def
    curl(:,:,ik) = 0.0_r_def
    div(:,:,ik) = 0.0_r_def
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        ! Precompute some frequently used terms
        call pointwise_coordinate_jacobian(ndf_chi, chi1_e, chi2_e, chi3_e,  &
                                           diff_basis_chi(:,:,qp1,qp2), &
                                           jac, dj)
        jac_inv = pointwise_coordinate_jacobian_inverse(jac,dj)
        jac_t = transpose(jac_inv)
        wt = wqp_h(qp1)*wqp_v(qp2)

        ! W0 mass matrix
        do df2 = 1, ndf_w0
          do df = df2, ndf_w0 ! mass matrix is symmetric
            integrand = wt  &
                       *basis_w0(1,df,qp1,qp2) &
                       *basis_w0(1,df2,qp1,qp2)*dj
            mm0(df,df2,ik) = mm0(df,df2,ik) + integrand
          end do
        end do
        ! W1 mass matrix
        do df2 = 1, ndf_w1
          jac_v = matmul(jac_t,basis_w1(:,df2,qp1,qp2))
          do df = df2, ndf_w1 ! mass matrix is symmetric
            integrand = wt*                                               &
                 dot_product(                                             &
                 matmul(jac_t,basis_w1(:,df,qp1,qp2)),jac_v               &
                            )*dj
            mm1(df,df2,ik) = mm1(df,df2,ik) + integrand
          end do
        end do
        ! W2 mass matrix
        do df2 = 1, ndf_w2
          jac_v = matmul(jac,basis_w2(:,df2,qp1,qp2))
          do df = df2, ndf_w2 ! mass matrix is symmetric
            integrand = wt*                                               &
                 dot_product(                                             &
                 matmul(jac,basis_w2(:,df,qp1,qp2)),jac_v    &
                            )/dj
            mm2(df,df2,ik) = mm2(df,df2,ik) + integrand
          end do
        end do
        ! W3 mass matrix
        do df2 = 1, ndf_w3
          do df = df2, ndf_w3 ! mass matrix is symmetric
            integrand = wt  &
                       *basis_w3(1,df,qp1,qp2) &
                       *basis_w3(1,df2,qp1,qp2)*dj
            mm3(df,df2,ik) = mm3(df,df2,ik) + integrand
          end do
        end do
        ! Wtheta mass matrix
        do df2 = 1, ndf_wt
          do df = df2, ndf_wt ! mass matrix is symmetric
            integrand = wt  &
                       *basis_wt(1,df,qp1,qp2) &
                       *basis_wt(1,df2,qp1,qp2)*dj
            mmt(df,df2,ik) = mmt(df,df2,ik) + integrand
          end do
        end do
        ! Grad matrix
        do df2 = 1, ndf_w0
          jac_v = matmul(jac_t,diff_basis_w0(:,df2,qp1,qp2))
          do df = 1, ndf_w1
            integrand = wt*                                               &
                 dot_product(                                             &
                 matmul(jac_t,basis_w1(:,df,qp1,qp2)),jac_v               &
                            )*dj
            grad(df,df2,ik) = grad(df,df2,ik) + integrand
          end do
        end do
        ! Curl matrix
        do df2 = 1, ndf_w1
          jac_v = matmul(jac,diff_basis_w1(:,df2,qp1,qp2))
          do df = 1, ndf_w2
            integrand = wt*                                               &
                 dot_product(                                             &
                 matmul(jac,basis_w2(:,df,qp1,qp2)),jac_v                 &
                            )/dj
            curl(df,df2,ik) = curl(df,df2,ik) + integrand
          end do
        end do
        ! Div matrix
        do df2 = 1, ndf_w2
          do df = 1, ndf_w3
            integrand = wt  &
                       *basis_w3(1,df,qp1,qp2) &
                       *diff_basis_w2(1,df2,qp1,qp2)!*dj
            div(df,df2,ik) = div(df,df2,ik) + integrand
          end do
        end do
      end do
    end do
    ! Apply symmetry of mass matrices
    do df2 = 1,ndf_w0
      do df = df2, 1, -1
        mm0(df,df2,ik) = mm0(df2,df,ik)
      end do
    end do
    do df2 = 1,ndf_w1
      do df = df2, 1, -1
        mm1(df,df2,ik) = mm1(df2,df,ik)
      end do
    end do
    do df2 = 1,ndf_w2
      do df = df2, 1, -1
        mm2(df,df2,ik) = mm2(df2,df,ik)
      end do
    end do
    do df2 = 1,ndf_w3
      do df = df2, 1, -1
        mm3(df,df2,ik) = mm3(df2,df,ik)
      end do
    end do
    do df2 = 1,ndf_wt
      do df = df2, 1, -1
        mmt(df,df2,ik) = mmt(df2,df,ik)
      end do
    end do
  end do ! end of k loop
end subroutine compute_derahm_matrices_code

end module compute_derahm_matrices_kernel_mod
