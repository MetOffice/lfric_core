!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Computes the mass matrices of the derham complex.
!>
!> @details Compute the mass matrices of the derham cochain,
!> these are the mass matrices for the W0, W1, W2 & W3 function spaces
!> along with the gradient, curl and divergence operators and the
!> W2broken/Wtheta mass matrices.
!> Since they all depend upon the mesh Jacobian they are computed as one
!> to reduce the cost
!>
module sci_compute_derham_matrices_kernel_mod

  use argument_mod,            only: arg_type, func_type,       &
                                     GH_OPERATOR, GH_FIELD,     &
                                     GH_READ, GH_WRITE,         &
                                     GH_REAL, ANY_SPACE_9,      &
                                     ANY_DISCONTINUOUS_SPACE_3, &
                                     GH_BASIS, GH_DIFF_BASIS,   &
                                     CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,           only: r_def, i_def
  use sci_coordinate_jacobian_mod, only: pointwise_coordinate_jacobian, &
                                         pointwise_coordinate_jacobian_inverse
  use fs_continuity_mod,       only: W0, W1, W2, W2broken, W3, Wtheta
  use kernel_mod,              only: kernel_type

  use base_mesh_config_mod,      only: geometry, topology
  use finite_element_config_mod, only: coord_system
  use planet_config_mod,         only: scaled_radius

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: compute_derham_matrices_kernel_type
    private
    type(arg_type) :: meta_args(12) = (/                                     &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W0, W0),                   &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W1, W1),                   &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2, W2),                   &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2broken, W2broken),       &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W3, W3),                   &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, Wtheta, Wtheta),           &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W1, W0),                   &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2, W1),                   &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W3, W2),                   &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W3, W2broken),             &
         arg_type(GH_FIELD*3,  GH_REAL, GH_READ,  ANY_SPACE_9),              &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3) &
         /)
    type(func_type) :: meta_funcs(7) = (/                                    &
         func_type(W0,          GH_BASIS, GH_DIFF_BASIS),                    &
         func_type(W1,          GH_BASIS, GH_DIFF_BASIS),                    &
         func_type(W2,          GH_BASIS, GH_DIFF_BASIS),                    &
         func_type(W2broken,    GH_BASIS, GH_DIFF_BASIS),                    &
         func_type(W3,          GH_BASIS),                                   &
         func_type(Wtheta,      GH_BASIS),                                   &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                     &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: compute_derham_matrices_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_derham_matrices_code

contains

!> @brief This subroutine computes the operator matrices for the modified derham complex
!! @param[in] cell Cell number.
!! @param[in] nlayers Number of layers.
!! @param[in] ncell_3d0 ncell*nlayers
!! @param[in,out] mm0 Local assembly mass matrix for W0 space.
!! @param[in] ncell_3d1 ncell*nlayers
!! @param[in,out] mm1 Local assembly mass matrix for W1 space.
!! @param[in] ncell_3d2 ncell*nlayers
!! @param[in,out] mm2 Local assembly mass matrix for W2 space.
!! @param[in] ncell_3db ncell*nlayers
!! @param[in,out] mm2b Local assembly mass matrix for W2broken space.
!! @param[in] ncell_3d3 ncell*nlayers
!! @param[in,out] mm3 Local assembly mass matrix for W3 space.
!! @param[in] ncell_3dt ncell*nlayers
!! @param[in,out] mmt Local assembly mass matrix for Wtheta space.
!! @param[in] ncell_3d4 ncell*nlayers
!! @param[in,out] grad Local assembly gradient matrix.
!! @param[in] ncell_3d5 ncell*nlayers
!! @param[in,out] curl Local assembly curl matrix.
!! @param[in] ncell_3d6 ncell*nlayers
!! @param[in,out] div Local assembly divergence matrix.
!! @param[in] ncell_3d7 ncell*nlayers
!! @param[in,out] broken_div Local assembly divergence matrix for W2broken.
!! @param[in] chi1 Physical coordinates in the 1st dir.
!! @param[in] chi2 Physical coordinates in the 2nd dir.
!! @param[in] chi3 Physical coordinates in the 3rd dir.
!! @param[in] panel_id Field giving the ID for mesh panels.
!! @param[in] ndf_w0 Number of degrees of freedom per cell for W0 space.
!! @param[in] basis_w0 Basis functions evaluated at quadrature points for W0 space.
!! @param[in] diff_basis_w0 Differential of basis functions evaluated at quadrature points for W0 space.
!! @param[in] ndf_w1 Number of degrees of freedom per cell for W1 space.
!! @param[in] basis_w1 Basis functions evaluated at quadrature points for W1 space.
!! @param[in] diff_basis_w1 Differential of basis functions evaluated at quadrature points for W1 space.
!! @param[in] ndf_w2 Number of degrees of freedom per cell for W2 space.
!! @param[in] basis_w2 Basis functions evaluated at quadrature points for W2 space.
!! @param[in] diff_basis_w2 Differential of basis functions evaluated at quadrature points for W2 space.
!! @param[in] ndf_w2b Number of degrees of freedom per cell for W2broken space.
!! @param[in] basis_w2b Basis functions evaluated at quadrature points for W2broken space.
!! @param[in] diff_basis_w2b Differential of basis functions evaluated at quadrature points for W2broken space.
!! @param[in] ndf_w3 Number of degrees of freedom per cell for W3 space.
!! @param[in] basis_w3 Basis functions evaluated at quadrature points for W3 space.
!! @param[in] ndf_wt Number of degrees of freedom per cell for Wtheta space.
!! @param[in] basis_wt Basis functions evaluated at quadrature points for Wtheta space.
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi field
!! @param[in] undf_chi Number of unique degrees of freedom for chi field
!! @param[in] map_chi Dofmap for the cell at the base of the column, for the
!!                    space on which the chi field lives
!! @param[in] basis_chi Basis functions for Wchi evaluated at quadrature points.
!! @param[in] diff_basis_chi Vector differential basis functions evaluated at
!!                           quadrature points.
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine compute_derham_matrices_code(cell, nlayers,                      &
                                        ncell_3d0, mm0,                     &
                                        ncell_3d1, mm1,                     &
                                        ncell_3d2, mm2,                     &
                                        ncell_3db, mm2b,                    &
                                        ncell_3d3, mm3,                     &
                                        ncell_3dt, mmt,                     &
                                        ncell_3d4, grad,                    &
                                        ncell_3d5, curl,                    &
                                        ncell_3d6, div,                     &
                                        ncell_3d7, broken_div,              &
                                        chi1, chi2, chi3,                   &
                                        panel_id,                           &
                                        ndf_w0, basis_w0, diff_basis_w0,    &
                                        ndf_w1, basis_w1, diff_basis_w1,    &
                                        ndf_w2, basis_w2, diff_basis_w2,    &
                                        ndf_w2b, basis_w2b, diff_basis_w2b, &
                                        ndf_w3, basis_w3,                   &
                                        ndf_wt, basis_wt,                   &
                                        ndf_chi, undf_chi,                  &
                                        map_chi, basis_chi, diff_basis_chi, &
                                        ndf_pid, undf_pid, map_pid,         &
                                        nqp_h, nqp_v, wqp_h, wqp_v)

  implicit none
  ! Arguments
  integer(kind=i_def),   intent(in)     :: cell, nqp_h, nqp_v
  integer(kind=i_def),   intent(in)     :: nlayers
  integer(kind=i_def),   intent(in)     :: ndf_w0, ndf_w1, ndf_w2, ndf_w3, ndf_wt, ndf_w2b
  integer(kind=i_def),   intent(in)     :: ncell_3d0, ncell_3d1, ncell_3d2
  integer(kind=i_def),   intent(in)     :: ncell_3d3, ncell_3d4, ncell_3d5
  integer(kind=i_def),   intent(in)     :: ncell_3d6, ncell_3dt, ncell_3db
  integer(kind=i_def),   intent(in)     :: ncell_3d7
  integer(kind=i_def),   intent(in)     :: ndf_chi
  integer(kind=i_def),   intent(in)     :: undf_chi
  integer(kind=i_def),   intent(in)     :: ndf_pid
  integer(kind=i_def),   intent(in)     :: undf_pid
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  real(kind=r_def), intent(inout) :: mm0(ncell_3d0,ndf_w0,ndf_w0)
  real(kind=r_def), intent(inout) :: mm1(ncell_3d1,ndf_w1,ndf_w1)
  real(kind=r_def), intent(inout) :: mm2(ncell_3d2,ndf_w2,ndf_w2)
  real(kind=r_def), intent(inout) :: mm2b(ncell_3db,ndf_w2b,ndf_w2b)
  real(kind=r_def), intent(inout) :: mm3(ncell_3d3,ndf_w3,ndf_w3)
  real(kind=r_def), intent(inout) :: mmt(ncell_3dt,ndf_wt,ndf_wt)
  real(kind=r_def), intent(inout) :: grad(ncell_3d4,ndf_w1,ndf_w0)
  real(kind=r_def), intent(inout) :: curl(ncell_3d5,ndf_w2,ndf_w1)
  real(kind=r_def), intent(inout) :: div(ncell_3d6,ndf_w3,ndf_w2)
  real(kind=r_def), intent(inout) :: broken_div(ncell_3d7,ndf_w3,ndf_w2b)

  real(kind=r_def), intent(in)  :: basis_chi(1,ndf_chi,nqp_h,nqp_v)
  real(kind=r_def), intent(in)  :: diff_basis_chi(3,ndf_chi,nqp_h,nqp_v)
  real(kind=r_def), intent(in)  :: basis_w0(1,ndf_w0,nqp_h,nqp_v)
  real(kind=r_def), intent(in)  :: basis_w1(3,ndf_w1,nqp_h,nqp_v)
  real(kind=r_def), intent(in)  :: basis_w2(3,ndf_w2,nqp_h,nqp_v)
  real(kind=r_def), intent(in)  :: basis_w2b(3,ndf_w2,nqp_h,nqp_v)
  real(kind=r_def), intent(in)  :: basis_w3(1,ndf_w3,nqp_h,nqp_v)
  real(kind=r_def), intent(in)  :: basis_wt(1,ndf_wt,nqp_h,nqp_v)
  real(kind=r_def), intent(in)  :: diff_basis_w0(3,ndf_w0,nqp_h,nqp_v)
  real(kind=r_def), intent(in)  :: diff_basis_w1(3,ndf_w1,nqp_h,nqp_v)
  real(kind=r_def), intent(in)  :: diff_basis_w2(1,ndf_w2,nqp_h,nqp_v)
  real(kind=r_def), intent(in)  :: diff_basis_w2b(1,ndf_w2,nqp_h,nqp_v)
  real(kind=r_def), intent(in)  :: chi1(undf_chi)
  real(kind=r_def), intent(in)  :: chi2(undf_chi)
  real(kind=r_def), intent(in)  :: chi3(undf_chi)
  real(kind=r_def), intent(in)  :: panel_id(undf_pid)
  real(kind=r_def), intent(in)  :: wqp_h(nqp_h)
  real(kind=r_def), intent(in)  :: wqp_v(nqp_v)

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

  integer(kind=i_def) :: ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1

    ! Indirect the chi coord field here
    do df = 1, ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do

    ik = 1 + k + (cell-1)*nlayers
    mm0(ik,:,:) = 0.0_r_def
    mm1(ik,:,:) = 0.0_r_def
    mm2(ik,:,:) = 0.0_r_def
    mm2b(ik,:,:) = 0.0_r_def
    mm3(ik,:,:) = 0.0_r_def
    mmt(ik,:,:) = 0.0_r_def
    grad(ik,:,:) = 0.0_r_def
    curl(ik,:,:) = 0.0_r_def
    div(ik,:,:) = 0.0_r_def
    broken_div(ik,:,:) = 0.0_r_def

    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        ! Precompute some frequently used terms
        call pointwise_coordinate_jacobian(coord_system, geometry,          &
                                           topology, scaled_radius,         &
                                           ndf_chi, chi1_e, chi2_e, chi3_e, &
                                           ipanel, basis_chi(:,:,qp1,qp2),  &
                                           diff_basis_chi(:,:,qp1,qp2),     &
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
            mm0(ik,df,df2) = mm0(ik,df,df2) + integrand
          end do
        end do
        ! W1 mass matrix
        do df2 = 1, ndf_w1
          jac_v = matmul(jac_t,basis_w1(:,df2,qp1,qp2))
          do df = df2, ndf_w1 ! mass matrix is symmetric
            integrand = wt*                                  &
                 dot_product(                                &
                 matmul(jac_t,basis_w1(:,df,qp1,qp2)),jac_v  &
                            )*dj
            mm1(ik,df,df2) = mm1(ik,df,df2) + integrand
          end do
        end do
        ! W2 mass matrix
        do df2 = 1, ndf_w2
          jac_v = matmul(jac,basis_w2(:,df2,qp1,qp2))
          do df = df2, ndf_w2 ! mass matrix is symmetric
            integrand = wt*                                  &
                 dot_product(                                &
                 matmul(jac,basis_w2(:,df,qp1,qp2)),jac_v    &
                            )/dj
            mm2(ik,df,df2) = mm2(ik,df,df2) + integrand
          end do
        end do
        ! W2broken mass matrix
        do df2 = 1, ndf_w2
          jac_v = matmul(jac,basis_w2b(:,df2,qp1,qp2))
          do df = df2, ndf_w2 ! mass matrix is symmetric
            integrand = wt*                                   &
                 dot_product(                                 &
                 matmul(jac,basis_w2b(:,df,qp1,qp2)),jac_v    &
                            )/dj
            mm2b(ik,df,df2) = mm2b(ik,df,df2) + integrand
          end do
        end do
        ! W3 mass matrix
        do df2 = 1, ndf_w3
          do df = df2, ndf_w3 ! mass matrix is symmetric
            integrand = wt                                    &
                       *basis_w3(1,df,qp1,qp2)                &
                       *basis_w3(1,df2,qp1,qp2)*dj
            mm3(ik,df,df2) = mm3(ik,df,df2) + integrand
          end do
        end do
        ! Wtheta mass matrix
        do df2 = 1, ndf_wt
          do df = df2, ndf_wt ! mass matrix is symmetric
            integrand = wt                                    &
                       *basis_wt(1,df,qp1,qp2)                &
                       *basis_wt(1,df2,qp1,qp2)*dj
            mmt(ik,df,df2) = mmt(ik,df,df2) + integrand
          end do
        end do
        ! Grad matrix
        do df2 = 1, ndf_w0
          jac_v = matmul(jac_t,diff_basis_w0(:,df2,qp1,qp2))
          do df = 1, ndf_w1
            integrand = wt*                                   &
                 dot_product(                                 &
                 matmul(jac_t,basis_w1(:,df,qp1,qp2)),jac_v   &
                            )*dj
            grad(ik,df,df2) = grad(ik,df,df2) + integrand
          end do
        end do
        ! Curl matrix
        do df2 = 1, ndf_w1
          jac_v = matmul(jac,diff_basis_w1(:,df2,qp1,qp2))
          do df = 1, ndf_w2
            integrand = wt*                                   &
                 dot_product(                                 &
                 matmul(jac,basis_w2(:,df,qp1,qp2)),jac_v     &
                            )/dj
            curl(ik,df,df2) = curl(ik,df,df2) + integrand
          end do
        end do
        ! Div matrix
        do df2 = 1, ndf_w2
          do df = 1, ndf_w3
            integrand = wt                                    &
                       *basis_w3(1,df,qp1,qp2)                &
                       *diff_basis_w2(1,df2,qp1,qp2)!*dj
            div(ik,df,df2) = div(ik,df,df2) + integrand
          end do
        end do
        ! Broken div matrix
        do df2 = 1, ndf_w2
          do df = 1, ndf_w3
            integrand = wt  &
                       *basis_w3(1,df,qp1,qp2) &
                       *diff_basis_w2b(1,df2,qp1,qp2)!*dj
            broken_div(ik,df,df2) = broken_div(ik,df,df2) + integrand
          end do
        end do
      end do
    end do
    ! Apply symmetry of mass matrices
    do df2 = 1,ndf_w0
      do df = df2, 1, -1
        mm0(ik,df,df2) = mm0(ik,df2,df)
      end do
    end do
    do df2 = 1,ndf_w1
      do df = df2, 1, -1
        mm1(ik,df,df2) = mm1(ik,df2,df)
      end do
    end do
    do df2 = 1,ndf_w2
      do df = df2, 1, -1
        mm2(ik,df,df2) = mm2(ik,df2,df)
      end do
    end do
    do df2 = 1,ndf_w2
      do df = df2, 1, -1
        mm2b(ik,df,df2) = mm2b(ik,df2,df)
      end do
    end do
    do df2 = 1,ndf_w3
      do df = df2, 1, -1
        mm3(ik,df,df2) = mm3(ik,df2,df)
      end do
    end do
    do df2 = 1,ndf_wt
      do df = df2, 1, -1
        mmt(ik,df,df2) = mmt(ik,df2,df)
      end do
    end do
  end do ! end of k loop
end subroutine compute_derham_matrices_code

end module sci_compute_derham_matrices_kernel_mod
