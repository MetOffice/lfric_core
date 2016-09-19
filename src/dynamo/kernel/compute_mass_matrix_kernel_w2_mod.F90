!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides access to the members of the w2_kernel class.
!>
!> @details Accessor functions for the w2_kernel class are defined in this
!>         module.
!
module compute_mass_matrix_kernel_w2_mod

use constants_mod,           only: i_def, r_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, func_type,             &
                                   GH_OPERATOR, GH_FIELD,           &
                                   GH_READ, GH_WRITE,               &
                                   ANY_SPACE_9, W2, GH_BASIS, GH_DIFF_BASIS, &
                                   CELLS

use coordinate_jacobian_mod, only: coordinate_jacobian
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: compute_mass_matrix_kernel_w2_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_OPERATOR, GH_WRITE, W2, W2),                        &
       arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_9)                             &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(ANY_SPACE_9, GH_DIFF_BASIS),                                   &
       func_type(W2, GH_BASIS)                                         &
       /)
  integer :: iterates_over = CELLS

contains
  procedure, nopass :: compute_mass_matrix_w2_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface compute_mass_matrix_kernel_w2_type
   module procedure compute_mass_matrix_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public compute_mass_matrix_w2_code
contains

type(compute_mass_matrix_kernel_w2_type) &
                        function compute_mass_matrix_constructor() result(self)
  return
end function compute_mass_matrix_constructor

!> @brief Computes the mass matrix for the w2 space.
!!
!! @param[in] cell     Identifying number of cell.
!! @param[in] nlayers  Number of layers.
!! @param[in] ndf_w2   Degrees of freedom per cell.
!! @param[in] ncell_3d ncell*ndf
!! @param[in] basis_w2 Vector basis functions evaluated at quadrature points.
!! @param[in] mm       Local stencil or mass matrix.
!! @param[in] ndf_chi  Degrees of freedum per cell for chi field.
!! @param[in] undf_chi Unique degrees of freedum  for chi field.
!! @param[in] map_chi  Dofmap for the cell at the base of the column, for the
!!                     space on which the chi field lives
!! @param[in] diff_basis_chi Vector differential basis functions evaluated at
!!                           quadrature points.
!! @param[inout] chi1  Physical coordinate in the 1st dir.
!! @param[inout] chi2  Physical coordinate in the 2nd dir.
!! @param[inout] chi3  Physical coordinate in the 3rd dir.
!! @param[in] nqp_h    Number of horizontal quadrature points.
!! @param[in] nqp_v    Number of vertical quadrature points.
!! @param[in] wqp_h    Horizontal quadrature weights.
!! @param[in] wqp_v    Vertical quadrature weights.
subroutine compute_mass_matrix_w2_code(cell, nlayers, ncell_3d,     &
                                       mm,                          &
                                       chi1, chi2, chi3,            &
                                       ndf_w2, basis_w2,            &
                                       ndf_chi, undf_chi, map_chi,  &
                                       diff_basis_chi,              &
                                       nqp_h, nqp_v, wqp_h, wqp_v )

  !Arguments
  integer,   intent(in)     :: cell, nqp_h, nqp_v
  integer,   intent(in)     :: nlayers
  integer,   intent(in)     :: ncell_3d
  integer,   intent(in)     :: ndf_w2
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in) :: basis_w2
  integer,   intent(in)     :: ndf_chi
  integer,   intent(in)     :: undf_chi
  integer(i_def), intent(in)    :: map_chi(ndf_chi)
  real(r_def),    intent(inout) :: mm(ndf_w2,ndf_w2,ncell_3d)
  real(r_def),    intent(in)    :: diff_basis_chi(3,ndf_chi,nqp_h,nqp_v)
  real(r_def),    intent(inout) :: chi1(undf_chi)
  real(r_def),    intent(inout) :: chi2(undf_chi)
  real(r_def),    intent(inout) :: chi3(undf_chi)
  real(r_def),    intent(in)    :: wqp_h(nqp_h)
  real(r_def),    intent(in)    :: wqp_v(nqp_v)

  !Internal variables
  integer                                      :: df, df2, k, ik
  integer                                      :: qp1, qp2

  real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
  real(kind=r_def)                             :: integrand
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac

  !loop over layers: Start from 1 as in this loop k is not an offset
  do k = 1, nlayers
     ik = k + (cell-1)*nlayers

     ! indirect the chi coord field here
     do df = 1, ndf_chi
        chi1_e(df) = chi1(map_chi(df) + k - 1)
        chi2_e(df) = chi2(map_chi(df) + k - 1)
        chi3_e(df) = chi3(map_chi(df) + k - 1)
     end do

    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                             diff_basis_chi, jac, dj)
    do df2 = 1, ndf_w2
       do df = df2, ndf_w2 ! mass matrix is symmetric
          mm(df,df2,ik) = 0.0_r_def
          do qp2 = 1, nqp_v
             do qp1 = 1, nqp_h
                integrand = wqp_h(qp1) * wqp_v(qp2) *                   & 
                     dot_product(                                       &
                     matmul(jac(:,:,qp1,qp2),basis_w2(:,df,qp1,qp2)),   &
                     matmul(jac(:,:,qp1,qp2),basis_w2(:,df2,qp1,qp2)) ) &
                     /dj(qp1,qp2) 
                mm(df,df2,ik) = mm(df,df2,ik) + integrand
             end do
          end do
       end do
       do df = df2, 1, -1
          mm(df,df2,ik) = mm(df2,df,ik)
       end do
    end do
  end do ! end of k loop

end subroutine compute_mass_matrix_w2_code

end module compute_mass_matrix_kernel_w2_mod
