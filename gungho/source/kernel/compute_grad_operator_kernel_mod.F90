!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

module compute_grad_operator_kernel_mod
use constants_mod,           only: r_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, func_type,            &
                                   GH_OPERATOR, GH_FIELD,          &
                                   GH_READ, GH_WRITE,              &
                                   W0, W1, ANY_SPACE_1,            &
                                   GH_BASIS,GH_DIFF_BASIS,         &
                                   CELLS, GH_QUADRATURE_XYoZ,      &
                                   QUADRATURE_XYoZ
use coordinate_jacobian_mod, only: coordinate_jacobian, &
                                   coordinate_jacobian_inverse
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: compute_grad_operator_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_OPERATOR, GH_WRITE, W1, W0),                        &
       arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_1)                    &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(W1, GH_BASIS),                                        &
       func_type(W0, GH_DIFF_BASIS),                                   &
       func_type(ANY_SPACE_1, GH_DIFF_BASIS)                           &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = QUADRATURE_XYoZ
contains
  procedure, nopass :: compute_grad_operator_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface compute_grad_operator_kernel_type
   module procedure compute_grad_operator_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public compute_grad_operator_code
contains

type(compute_grad_operator_kernel_type) &
                      function compute_grad_operator_constructor() result(self)
  return
end function compute_grad_operator_constructor

!> @brief Computes the grad operator 
!! @param[in] cell Cell number
!! @param[in] nlayers Number of layers.
!! @param[in] ncell_3d ncell*ndf
!! @param[in] ndf_w1 Number of degrees of freedom per cell.
!! @param[in] basis_w1 Vector basis functions
!!                    evaluated at quadrature points.
!! @param[in] ndf_w0 Number of degrees of freedom per cell.
!! @param[in] diff_basis_w0 Differential of scalar
!!                    basis functions evaluated at quadrature points.
!! @param[in] grad Local stencil of the grad operator
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!!                    field
!! @param[in] undf_chi Number of unique degrees of freedom  for chi
!!                    field
!! @param[in] map_chi Array holding the dofmap for the cell at the
!!                    base of the column, for the space on which the chi field
!!                    lives
!! @param[in] diff_basis_chi Vector differential
!!                    basis functions evaluated at quadrature points.
!! @param[inout] chi1 Physical coordinates in the 1st dir
!! @param[inout] chi2 Physical coordinates in the 2nd dir
!! @param[inout] chi3 Physical coordinates in the 3rd dir
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine compute_grad_operator_code(cell, nlayers, ncell_3d,          &
                                      grad,                             &
                                      chi1, chi2, chi3,                 &
                                      ndf_w1, basis_w1,                 &
                                      ndf_w0, diff_basis_w0,            &
                                      ndf_chi, undf_chi,                &
                                      map_chi, diff_basis_chi,          &
                                      nqp_h, nqp_v, wqp_h, wqp_v )

  !Arguments
  integer,                     intent(in) :: cell, nqp_h, nqp_v
  integer,                     intent(in) :: nlayers
  integer,                     intent(in) :: ncell_3d
  integer,                     intent(in) :: ndf_w1, ndf_w0
  integer,                     intent(in) :: ndf_chi, undf_chi
  integer, dimension(ndf_chi), intent(in) :: map_chi

  real(kind=r_def), intent(in) :: diff_basis_chi(3,ndf_chi,nqp_h,nqp_v)
  real(kind=r_def), intent(in) :: basis_w1(3,ndf_w1,nqp_h,nqp_v)
  real(kind=r_def), intent(in) :: diff_basis_w0(3,ndf_w0,nqp_h,nqp_v)

  real(kind=r_def), dimension(ndf_w1,ndf_w0,ncell_3d), intent(inout) :: grad
  real(kind=r_def), dimension(undf_chi),               intent(in)    :: chi1
  real(kind=r_def), dimension(undf_chi),               intent(in)    :: chi2
  real(kind=r_def), dimension(undf_chi),               intent(in)    :: chi3
  real(kind=r_def), dimension(nqp_h),                  intent(in)    :: wqp_h
  real(kind=r_def), dimension(nqp_v),                  intent(in)    :: wqp_v

  !Internal variables
  integer                                      :: df, df0, df1, k, ik
  integer                                      :: qp1, qp2
  real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
  real(kind=r_def)                             :: integrand
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac, jac_inv
  real(kind=r_def), dimension(3)               :: c, dgamma


  do k = 0, nlayers - 1
    ik = k + 1 + (cell-1)*nlayers
    do df = 1, ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                             diff_basis_chi, jac, dj)
    call coordinate_jacobian_inverse(nqp_h, nqp_v, jac, dj, jac_inv)
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        jac_inv(:,:,qp1,qp2) = transpose(jac_inv(:,:,qp1,qp2))
      end do
    end do
     
    do df0 = 1, ndf_w0
      do df1 = 1, ndf_w1
        grad(df1,df0,ik) = 0.0_r_def
        do qp2 = 1, nqp_v
          do qp1 = 1, nqp_h
            c      = matmul(jac_inv(:,:,qp1,qp2),     basis_w1(:,df1,qp1,qp2))
            dgamma = matmul(jac_inv(:,:,qp1,qp2),diff_basis_w0(:,df0,qp1,qp2))
            integrand = wqp_h(qp1)*wqp_v(qp2)*dot_product(c,dgamma)*dj(qp1,qp2)
            grad(df1,df0,ik) = grad(df1,df0,ik) + integrand
          end do
        end do
      end do
    end do
  end do 
end subroutine compute_grad_operator_code

end module compute_grad_operator_kernel_mod
