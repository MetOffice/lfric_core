!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides access to the members of the w2_kernel class.

!> @details Accessor functions for the w2_kernel class are defined in this module.

!> @param RHS_w2_code              Code to implement the RHS for a w2 field
!> @param gaussian_quadrature      Contains result of gaussian quadrature

module compute_mass_matrix_kernel_w2_mod
use constants_mod,           only: r_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, &           ! the type
                                   gh_inc, w2,fe, cells ! the enums

use coordinate_jacobian_mod, only: coordinate_jacobian
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: compute_mass_matrix_kernel_w2_type
  private
  type(arg_type) :: meta_args(1) = (/ &
       arg_type(gh_inc,w2,fe,.true.,.false.,.false.,.true.) &
       /)
  integer :: iterates_over = cells

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

type(compute_mass_matrix_kernel_w2_type) function compute_mass_matrix_constructor() result(self)
  return
end function compute_mass_matrix_constructor
  
subroutine compute_mass_matrix_w2_code(cell, nlayers, ndf, ncell_3d,          &
                               basis,  mm,                                    &
                               ndf_chi, undf_chi,                             &
                               map_chi, diff_basis_chi, &   ! arrays
                               chi1, chi2, chi3,                & ! dof vectors
                               nqp_h, nqp_v, wqp_h, wqp_v )

!> @brief This subroutine computes the mass matrix for the w2 space
!! @param[in] cell Integer: The cell number
!! @param[in] nlayers Integer: The number of layers.
!! @param[in] ndf Integer: The number of degrees of freedom per cell.
!! @param[in] ncell_3d Integer: ncell*ndf
!! @param[in] basis Real: 4-dim array holding VECTOR basis functions evaluated at quadrature points.
!! @param[in] mm Real array, the local stencil or mass matrix
!! @param[in] ndf_chi Integer: number of degrees of freedum per cell for chi field
!! @param[in] undf_chi Integer: number of unique degrees of freedum  for chi field
!! @param[in] map_chi Integer: Array holding the dofmap for the cell at the base of the column, for the space on which the chi field lives
!! @param[in] diff_basis_chi Real: 4-dim array holding VECTOR differential basis functions evaluated at quadrature points.
!! @param[inout] chi1 Real: The data array for chi in the first dir
!! @param[inout] chi2 Real: The data array for chi in the 2nd dir
!! @param[inout] chi3 Real: The data array for chi in the 3rd dir
!! @param[in] nqp_h Integer number of horizontal quadrature points
!! @param[in] nqp_v Integer number of vertical quadrature points
!! @param[in] wqp_h Real array. Quadrature weights horizontal
!! @param[in] wqp_v Real array. Quadrature weights vertical

  !Arguments
  integer,   intent(in)     :: cell, nqp_h, nqp_v
  integer,   intent(in)     :: nlayers, ndf
  integer,   intent(in)     :: ncell_3d
  integer,   intent(in)     :: ndf_chi
  integer,   intent(in)     :: undf_chi
  integer, dimension(ndf_chi), intent(in)     :: map_chi
  real(kind=r_def), dimension(ndf,ndf,ncell_3d),  intent(inout)  :: mm
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: diff_basis_chi
  real(kind=r_def), dimension(3,ndf,nqp_h,nqp_v), intent(in)     :: basis
  real(kind=r_def), dimension(undf_chi), intent(inout) :: chi1
  real(kind=r_def), dimension(undf_chi), intent(inout) :: chi2
  real(kind=r_def), dimension(undf_chi), intent(inout) :: chi3
  real(kind=r_def), dimension(nqp_h),    intent(in)    :: wqp_h
  real(kind=r_def), dimension(nqp_v),    intent(in)    :: wqp_v

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
    do df2 = 1, ndf
       do df = df2, ndf ! mass matrix is symmetric
          mm(df,df2,ik) = 0.0_r_def
          do qp2 = 1, nqp_v
             do qp1 = 1, nqp_h
                integrand = wqp_h(qp1) * wqp_v(qp2) *                 & 
                     dot_product(                                     &
                     matmul(jac(:,:,qp1,qp2), basis(:,df,qp1,qp2)),   &
                     matmul(jac(:,:,qp1,qp2), basis(:,df2,qp1,qp2)) ) &
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
