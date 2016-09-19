!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides access to the members of the w0_kernel class.

!> @details Accessor functions for the w0_kernel class are defined in this module.

!> @param RHS_w0_code              Code to implement the RHS for a w0 field
!> @param gaussian_quadrature      Contains result of gaussian quadrature

module compute_mass_matrix_kernel_w0_mod
use constants_mod,           only: r_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, func_type,                      &
                                   GH_OPERATOR, GH_FIELD, GH_READ, GH_WRITE, &
                                   W0, ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS, &
                                   CELLS
use coordinate_jacobian_mod, only: coordinate_jacobian

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public, extends(kernel_type) :: compute_mass_matrix_kernel_w0_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_OPERATOR, GH_WRITE, W0, W0),                        &
       arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_9)                    &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(W0, GH_BASIS),                         &
       func_type(ANY_SPACE_9, GH_DIFF_BASIS)                           &
       /)
  integer :: iterates_over = CELLS

contains
  procedure, nopass :: compute_mass_matrix_w0_code
end type compute_mass_matrix_kernel_w0_type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface compute_mass_matrix_kernel_w0
   module procedure compute_mass_matrix_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public compute_mass_matrix_w0_code
contains

type(compute_mass_matrix_kernel_w0_type) function compute_mass_matrix_constructor() result(self)
  return
end function compute_mass_matrix_constructor

!> @brief This subroutine computes the mass matrix for the w0 space
!! @param[in] cell The cell number
!! @param[in] nlayers The number of layers.
!! @param[in] ndf_w0 The number of degrees of freedom per cell for w0.
!! @param[in] ndf_chi The number of degrees of freedom per cell for chi.
!! @param[in] ncell_3d ncell*ndf
!! @param[in] basis_w0 4-dim array holding SCALAR basis functions evaluated at quadrature points for w0.
!! @param[inout] mm The mass matrix data array
!! @param[in] undf_chi Number of unique degrees of freedom for chi
!! @param[in] map_chi Array holding the dofmap for the cell at the base of the column for chi
!! @param[in] diff_basis_chi 4-dim array holding VECTOR differential basis functions evaluated at quadrature points for chi
!! @param[inout] chi1 The data array for chi in the 1st dir
!! @param[inout] chi2 The data array for chi in the 2nd dir
!! @param[inout] chi3 The data array for chi in the 3rd dir
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Quadrature weights horizontal
!! @param[in] wqp_v Quadrature weights vertical
subroutine compute_mass_matrix_w0_code(cell, nlayers, ncell_3d,            &
                                       mm, chi1, chi2, chi3,               &
                                       ndf_w0, basis_w0, ndf_chi,          &
                                       undf_chi, map_chi,  diff_basis_chi, &
                                       nqp_h, nqp_v, wqp_h, wqp_v )

  !Arguments
  integer, intent(in)     :: cell, nqp_h, nqp_v
  integer, intent(in)     :: nlayers, ndf_w0, ndf_chi, undf_chi
  integer, intent(in)     :: ncell_3d

  integer, dimension(ndf_chi), intent(in) :: map_chi

  real(kind=r_def), dimension(ndf_w0,ndf_w0,ncell_3d),  intent(inout)  :: mm

  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: diff_basis_chi
  real(kind=r_def), dimension(1,ndf_w0,nqp_h,nqp_v), intent(in)  :: basis_w0

  real(kind=r_def), dimension(undf_chi), intent(inout)           :: chi1
  real(kind=r_def), dimension(undf_chi), intent(inout)           :: chi2
  real(kind=r_def), dimension(undf_chi), intent(inout)           :: chi3

  real(kind=r_def), dimension(nqp_h), intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) :: wqp_v

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

    do df2 = 1, ndf_w0
       do df = df2, ndf_w0 ! mass matrix is symmetric
          mm(df,df2,ik) = 0.0_r_def
          do qp2 = 1, nqp_v
             do qp1 = 1, nqp_h
                integrand = wqp_h(qp1) * wqp_v(qp2) * & 
                     basis_w0(1,df,qp1,qp2)*basis_w0(1,df2,qp1,qp2)     * &
                     dj(qp1,qp2) 
                mm(df,df2,ik) = mm(df,df2,ik) + integrand
             end do
          end do
       end do
       do df = df2, 1, -1  
          mm(df,df2,ik) = mm(df2,df,ik)
       end do
       
    end do
    
  end do ! end of k loop



end subroutine compute_mass_matrix_w0_code

end module compute_mass_matrix_kernel_w0_mod
