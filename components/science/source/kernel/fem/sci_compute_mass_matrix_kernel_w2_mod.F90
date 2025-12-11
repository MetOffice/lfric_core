!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Provides access to the members of the w2_kernel class.
!>
!> Accessor functions for the w2_kernel class are defined in this module.
!> Module can be used for W2 and W2broken.
!>
module sci_compute_mass_matrix_kernel_w2_mod

  use argument_mod,            only: arg_type, func_type,       &
                                     GH_OPERATOR, GH_FIELD,     &
                                     GH_READ, GH_WRITE,         &
                                     GH_REAL, ANY_W2,           &
                                     ANY_DISCONTINUOUS_SPACE_3, &
                                     GH_BASIS, GH_DIFF_BASIS,   &
                                     CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,           only: i_def, r_def
  use sci_coordinate_jacobian_mod, only: coordinate_jacobian
  use fs_continuity_mod,       only: Wchi
  use kernel_mod,              only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: compute_mass_matrix_kernel_w2_type
    private
    type(arg_type) :: meta_args(3) = (/                                      &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, ANY_W2, ANY_W2),           &
         arg_type(GH_FIELD*3,  GH_REAL, GH_READ,  Wchi),                     &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3) &
         /)
    type(func_type) :: meta_funcs(2) = (/                                    &
         func_type(ANY_W2, GH_BASIS),                                        &
         func_type(Wchi,   GH_BASIS, GH_DIFF_BASIS)                          &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: compute_mass_matrix_w2_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_mass_matrix_w2_code

contains

!> @brief Computes the mass matrix for the W2 space.
!!
!! @param[in] cell     Identifying number of cell
!! @param[in] nlayers  Number of layers
!! @param[in] ncell_3d ncell*ndf
!! @param[in,out] mm   Local stencil or mass matrix
!! @param[in] chi1     1st coordinate field in Wchi
!! @param[in] chi2     2nd coordinate field in Wchi
!! @param[in] chi3     3rd coordinate field in Wchi
!! @param[in] panel_id Field giving the ID for mesh panels
!! @param[in] ndf_w2   Degrees of freedom per cell
!! @param[in] basis_w2 Vector basis functions evaluated at quadrature points
!! @param[in] ndf_chi  Degrees of freedom per cell for chi field
!! @param[in] undf_chi Unique degrees of freedom for chi field
!! @param[in] map_chi  Dofmap for the cell at the base of the column, for the
!!                     space on which the chi field lives
!! @param[in] basis_chi Wchi basis functions evaluated at quadrature points
!! @param[in] diff_basis_chi Vector differential basis functions evaluated at
!!                           quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h    Number of horizontal quadrature points
!! @param[in] nqp_v    Number of vertical quadrature points
!! @param[in] wqp_h    Horizontal quadrature weights
!! @param[in] wqp_v    Vertical quadrature weights
subroutine compute_mass_matrix_w2_code(cell, nlayers, ncell_3d,     &
                                       mm,                          &
                                       chi1, chi2, chi3,            &
                                       panel_id,                    &
                                       ndf_w2, basis_w2,            &
                                       ndf_chi, undf_chi, map_chi,  &
                                       basis_chi,                   &
                                       diff_basis_chi,              &
                                       ndf_pid, undf_pid,           &
                                       map_pid,                     &
                                       nqp_h, nqp_v, wqp_h, wqp_v )

  implicit none

  ! Arguments
  integer(kind=i_def),   intent(in) :: cell, nqp_h, nqp_v
  integer(kind=i_def),   intent(in) :: nlayers
  integer(kind=i_def),   intent(in) :: ncell_3d
  integer(kind=i_def),   intent(in) :: ndf_w2
  integer(kind=i_def),   intent(in) :: ndf_chi
  integer(kind=i_def),   intent(in) :: undf_chi
  integer(kind=i_def),   intent(in) :: ndf_pid, undf_pid
  integer(kind=i_def),   intent(in) :: map_chi(ndf_chi)
  integer(kind=i_def),   intent(in) :: map_pid(ndf_pid)

  real(kind=r_def),   intent(inout) :: mm(ncell_3d,ndf_w2,ndf_w2)
  real(kind=r_def),      intent(in) :: basis_w2(3,ndf_w2,nqp_h,nqp_v)
  real(kind=r_def),      intent(in) :: diff_basis_chi(3,ndf_chi,nqp_h,nqp_v)
  real(kind=r_def),      intent(in) :: basis_chi(1,ndf_chi,nqp_h,nqp_v)
  real(kind=r_def),      intent(in) :: chi1(undf_chi)
  real(kind=r_def),      intent(in) :: chi2(undf_chi)
  real(kind=r_def),      intent(in) :: chi3(undf_chi)
  real(kind=r_def),      intent(in) :: panel_id(undf_pid)
  real(kind=r_def),      intent(in) :: wqp_h(nqp_h)
  real(kind=r_def),      intent(in) :: wqp_v(nqp_v)

  ! Internal variables
  integer(kind=i_def)                          :: df, df2, k, ik, ipanel
  integer(kind=i_def)                          :: qp1, qp2

  real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
  real(kind=r_def)                             :: integrand
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! loop over layers: Start from 1 as in this loop k is not an offset
  do k = 1, nlayers
     ik = k + (cell-1)*nlayers

     ! indirect the chi coord field here
     do df = 1, ndf_chi
        chi1_e(df) = chi1(map_chi(df) + k - 1)
        chi2_e(df) = chi2(map_chi(df) + k - 1)
        chi3_e(df) = chi3(map_chi(df) + k - 1)
     end do

    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                             ipanel, basis_chi, diff_basis_chi, jac, dj)

    do df2 = 1, ndf_w2
       do df = df2, ndf_w2 ! mass matrix is symmetric
          mm(ik,df,df2) = 0.0_r_def
          do qp2 = 1, nqp_v
             do qp1 = 1, nqp_h
                integrand = wqp_h(qp1) * wqp_v(qp2) *                   &
                     dot_product(                                       &
                     matmul(jac(:,:,qp1,qp2),basis_w2(:,df,qp1,qp2)),   &
                     matmul(jac(:,:,qp1,qp2),basis_w2(:,df2,qp1,qp2)) ) &
                     /dj(qp1,qp2)
                mm(ik,df,df2) = mm(ik,df,df2) + integrand
             end do
          end do
       end do
       do df = df2, 1, -1
          mm(ik,df,df2) = mm(ik,df2,df)
       end do
    end do
  end do ! end of k loop

end subroutine compute_mass_matrix_w2_code

end module sci_compute_mass_matrix_kernel_w2_mod
