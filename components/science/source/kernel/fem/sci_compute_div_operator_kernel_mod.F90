!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
module sci_compute_div_operator_kernel_mod

  use argument_mod,              only: arg_type, func_type,       &
                                       GH_OPERATOR, GH_FIELD,     &
                                       GH_READ, GH_WRITE,         &
                                       GH_REAL, ANY_SPACE_1,      &
                                       GH_BASIS, GH_DIFF_BASIS,   &
                                       ANY_DISCONTINUOUS_SPACE_3, &
                                       GH_QUADRATURE_XYoZ, CELL_COLUMN
  use constants_mod,             only: r_def, i_def
  use sci_coordinate_jacobian_mod, only: coordinate_jacobian
  use fs_continuity_mod,         only: W2, W3
  use finite_element_config_mod, only: rehabilitate
  use kernel_mod,                only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: compute_div_operator_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                      &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W3, W2),                   &
         arg_type(GH_FIELD*3,  GH_REAL, GH_READ,  ANY_SPACE_1),              &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3) &
         /)
    type(func_type) :: meta_funcs(3) = (/                                    &
         func_type(W3,          GH_BASIS),                                   &
         func_type(W2,          GH_DIFF_BASIS),                              &
         func_type(ANY_SPACE_1, GH_BASIS, GH_DIFF_BASIS)                     &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: compute_div_operator_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: compute_div_operator_code

contains

!> @brief Computes the divergence operator
!! @param[in] cell     Cell number
!! @param[in] nlayers  Number of layers.
!! @param[in] ncell_3d ncell*ndf
!! @param[in,out] div  Local stencil of the div operator
!! @param[in] chi_1    1st coordinate field in Wchi
!! @param[in] chi_2    2nd coordinate field in Wchi
!! @param[in] chi_3    3rd coordinate field in Wchi
!! @param[in] panel_id Field giving the ID for mesh panels.
!! @param[in] ndf_w3   Number of degrees of freedom per cell.
!! @param[in] basis_w3 Scalar basis functions
!!                     evaluated at quadrature points.
!! @param[in] ndf_w2   Number of degrees of freedom per cell.
!! @param[in] diff_basis_w2 Differential vector basis
!!                     functions evaluated at quadrature points.
!! @param[in] ndf_chi  Number of degrees of freedom per cell for chi field
!! @param[in] undf_chi Number of unique degrees of freedom for chi field
!! @param[in] map_chi Dofmap for the cell at the base of the column, for the
!!                    space on which the chi field lives
!! @param[in] basis_chi Wchi basis functions evaluated at quadrature points.
!! @param[in] diff_basis_chi Wchi vector differential
!!                    basis functions evaluated at quadrature points.
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine compute_div_operator_code(cell, nlayers, ncell_3d,     &
                                     div,                         &
                                     chi1, chi2, chi3,            &
                                     panel_id,                    &
                                     ndf_w3,                      &
                                     basis_w3,                    &
                                     ndf_w2, diff_basis_w2,       &
                                     ndf_chi, undf_chi,           &
                                     map_chi,                     &
                                     basis_chi,                   &
                                     diff_basis_chi,              &
                                     ndf_pid,                     &
                                     undf_pid,                    &
                                     map_pid,                     &
                                     nqp_h, nqp_v, wqp_h, wqp_v )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: cell, nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ncell_3d
  integer(kind=i_def), intent(in) :: ndf_w3, ndf_w2
  integer(kind=i_def), intent(in) :: ndf_chi, undf_chi
  integer(kind=i_def), intent(in) :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in) :: map_chi(ndf_chi)
  integer(kind=i_def), intent(in) :: map_pid(ndf_pid)

  real(kind=r_def), intent(in) :: diff_basis_chi(3,ndf_chi,nqp_h,nqp_v)
  real(kind=r_def), intent(in) :: basis_chi(1,ndf_chi,nqp_h,nqp_v)
  real(kind=r_def), intent(in) :: basis_w3(1,ndf_w3,nqp_h,nqp_v)
  real(kind=r_def), intent(in) :: diff_basis_w2(1,ndf_w2,nqp_h,nqp_v)

  real(kind=r_def), dimension(ncell_3d,ndf_w3,ndf_w2), intent(inout) :: div
  real(kind=r_def), dimension(undf_pid),               intent(in)    :: panel_id
  real(kind=r_def), dimension(undf_chi),               intent(in)    :: chi1
  real(kind=r_def), dimension(undf_chi),               intent(in)    :: chi2
  real(kind=r_def), dimension(undf_chi),               intent(in)    :: chi3
  real(kind=r_def), dimension(nqp_h),                  intent(in)    :: wqp_h
  real(kind=r_def), dimension(nqp_v),                  intent(in)    :: wqp_v

  ! Internal variables
  integer(kind=i_def)                          :: df, df2, df3, k, ik
  integer(kind=i_def)                          :: qp1, qp2, ipanel
  real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
  real(kind=r_def)                             :: integrand
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers - 1
    ik = k + 1 + (cell-1)*nlayers
    do df = 1, ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                             ipanel, basis_chi, diff_basis_chi, jac, dj)

    do df2 = 1, ndf_w2
      do df3 = 1, ndf_w3
        div(ik,df3,df2) = 0.0_r_def
        do qp2 = 1, nqp_v
          do qp1 = 1, nqp_h
            if ( rehabilitate ) then
              ! With rehabilitation
              !   divergence mapping is div(x) -> ! \hat{div}(\hat{x})
              integrand = wqp_h(qp1)*wqp_v(qp2)                               &
                        *basis_w3(1,df3,qp1,qp2)*diff_basis_w2(1,df2,qp1,qp2)
            else
              ! Without rehabilitation
              !   divergence mapping is div(x) -> ! \hat{div}(\hat{x})/det(J)
              integrand = wqp_h(qp1)*wqp_v(qp2)                               &
                        *basis_w3(1,df3,qp1,qp2)*diff_basis_w2(1,df2,qp1,qp2) &
                        /dj(qp1,qp2)
            end if
            div(ik,df3,df2) = div(ik,df3,df2) + integrand
          end do
        end do
      end do
    end do
  end do
end subroutine compute_div_operator_code

end module sci_compute_div_operator_kernel_mod
