!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

module compute_div_operator_kernel_mod

use argument_mod,              only: arg_type, func_type,            &
                                     GH_OPERATOR, GH_FIELD,          &
                                     GH_READ, GH_WRITE,              &
                                     W2, W3, ANY_SPACE_9,            &
                                     GH_BASIS,GH_DIFF_BASIS,         &
                                     CELLS, GH_QUADRATURE_XYoZ,      &
                                     QUADRATURE_XYoZ
use constants_mod,             only: r_def
use coordinate_jacobian_mod,   only: coordinate_jacobian
use finite_element_config_mod, only: rehabilitate
use kernel_mod,                only: kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: compute_div_operator_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_OPERATOR, GH_WRITE, W3, W2),                        &
       arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_9)                    &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(W3, GH_BASIS),                                        &
       func_type(W2, GH_DIFF_BASIS),                                   &
       func_type(ANY_SPACE_9, GH_DIFF_BASIS)                           &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = QUADRATURE_XYoZ
contains
  procedure, nopass :: compute_div_operator_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface compute_div_operator_kernel_type
   module procedure compute_div_operator_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public compute_div_operator_code
contains

type(compute_div_operator_kernel_type) &
                       function compute_div_operator_constructor() result(self)
  return
end function compute_div_operator_constructor

!> @brief Computes the divergence operator 
!! @param[in] cell Cell number
!! @param[in] nlayers Number of layers.
!! @param[in] ncell_3d ncell*ndf
!! @param[in] ndf_w3 Number of degrees of freedom per cell.
!! @param[in] basis_w3 Scalar basis functions
!!                    evaluated at quadrature points.
!! @param[in] ndf_w2 Number of degrees of freedom per cell.
!! @param[in] diff_basis_w2 Differential vector basis
!!                    functions evaluated at quadrature points.
!! @param[in] div Local stencil of the div operator
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!!                    field
!! @param[in] undf_chi Number of unique degrees of freedom  for chi
!!                    field
!! @param[in] map_chi Dofmap for the cell at the
!!                    base of the column, for the space on which the chi field
!!                    lives
!! @param[in] diff_basis_chi Vector differential
!!                    basis functions evaluated at quadrature points.
!! @param[inout] chi1 Data array for chi in the 1st dir
!! @param[inout] chi2 Data array for chi in the 2nd dir
!! @param[inout] chi3 Data array for chi in the 3rd dir
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine compute_div_operator_code(cell, nlayers, ncell_3d,          &
                                     div,                              &
                                     chi1, chi2, chi3,                 &
                                     ndf_w3, basis_w3,                 &
                                     ndf_w2, diff_basis_w2,            &
                                     ndf_chi, undf_chi,                &
                                     map_chi, diff_basis_chi,          &
                                     nqp_h, nqp_v, wqp_h, wqp_v )

  !Arguments
  integer,                     intent(in) :: cell, nqp_h, nqp_v
  integer,                     intent(in) :: nlayers
  integer,                     intent(in) :: ncell_3d
  integer,                     intent(in) :: ndf_w3, ndf_w2
  integer,                     intent(in) :: ndf_chi, undf_chi
  integer, dimension(ndf_chi), intent(in) :: map_chi

  real(kind=r_def), intent(in) :: diff_basis_chi(3,ndf_chi,nqp_h,nqp_v)
  real(kind=r_def), intent(in) :: basis_w3(1,ndf_w3,nqp_h,nqp_v)
  real(kind=r_def), intent(in) :: diff_basis_w2(1,ndf_w2,nqp_h,nqp_v)

  real(kind=r_def), dimension(ndf_w3,ndf_w2,ncell_3d), intent(inout) :: div
  real(kind=r_def), dimension(undf_chi),               intent(in)    :: chi1
  real(kind=r_def), dimension(undf_chi),               intent(in)    :: chi2
  real(kind=r_def), dimension(undf_chi),               intent(in)    :: chi3
  real(kind=r_def), dimension(nqp_h),                  intent(in)    :: wqp_h
  real(kind=r_def), dimension(nqp_v),                  intent(in)    :: wqp_v

  !Internal variables
  integer                                      :: df, df2, df3, k, ik
  integer                                      :: qp1, qp2
  real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
  real(kind=r_def)                             :: integrand
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac

  do k = 0, nlayers - 1
    ik = k + 1 + (cell-1)*nlayers
    do df = 1, ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                             diff_basis_chi, jac, dj)
    do df2 = 1, ndf_w2
      do df3 = 1, ndf_w3
        div(df3,df2,ik) = 0.0_r_def
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
            div(df3,df2,ik) = div(df3,df2,ik) + integrand
          end do
        end do
      end do
    end do
  end do 
end subroutine compute_div_operator_code

end module compute_div_operator_kernel_mod
