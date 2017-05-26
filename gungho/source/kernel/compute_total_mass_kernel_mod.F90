!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the cell integrated mass

module compute_total_mass_kernel_mod

use argument_mod,      only : arg_type, func_type,             &
                              GH_FIELD, GH_WRITE, GH_READ,     &
                              ANY_SPACE_9, W3,                 &
                              GH_BASIS, GH_DIFF_BASIS,         &
                              CELLS, QUADRATURE_XYoZ,          &
                              GH_QUADRATURE_XYoZ
use constants_mod,     only : r_def
use kernel_mod,        only : kernel_type
use planet_config_mod, only : scaled_radius

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: compute_total_mass_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, W3),                             &
       arg_type(GH_FIELD,   GH_READ,  W3),                             &
       arg_type(GH_FIELD*3, GH_READ,  ANY_SPACE_9)                     &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(W3, GH_BASIS),                                        &
       func_type(ANY_SPACE_9, GH_DIFF_BASIS)                           &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = QUADRATURE_XYoZ
contains
  procedure, nopass ::compute_total_mass_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface compute_total_mass_kernel_type
   module procedure compute_total_mass_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public compute_total_mass_code
contains

type(compute_total_mass_kernel_type) function compute_total_mass_kernel_constructor() result(self)
  return
end function compute_total_mass_kernel_constructor

!> @brief Computes the cell integrated mass for a single model column
!! @param[in] nlayers Number of layers
!! @param[in] ndf_w3 Number of degrees of freedom per cell
!! @param[in] undf_w3 Total number of degrees of freedom for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column
!! @param[in] w3_basis Basis functions evaluated at gaussian quadrature points 
!! @param[out] mass Cell integrated mass
!! @param[in] rho Density
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Total number of degrees of freedom for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_diff_basis Basis functions evaluated at gaussian quadrature points
!! @param[inout] chi_1 X component of the chi coordinate field
!! @param[inout] chi_2 Y component of the chi coordinate field
!! @param[inout] chi_3 Z component of the chi coordinate field
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Weights of horizontal quadrature points
!! @param[in] wqp_v Weights of vertical quadrature points
subroutine compute_total_mass_code(                                            &
                                   nlayers,                                    &
                                   mass, rho, chi_1, chi_2, chi_3,             &
                                   ndf_w3, undf_w3, map_w3, w3_basis,          &
                                   ndf_chi, undf_chi, map_chi, chi_diff_basis, &
                                   nqp_h, nqp_v, wqp_h, wqp_v                  &
                                 )

  use coordinate_jacobian_mod, only : coordinate_jacobian 

  !Arguments
  integer, intent(in) :: nlayers, nqp_h, nqp_v
  integer, intent(in) :: ndf_w3, undf_w3, ndf_chi, undf_chi
  integer, dimension(ndf_w3), intent(in) :: map_w3
  integer, dimension(ndf_chi), intent(in) :: map_chi
  real(kind=r_def), intent(in), dimension(1,ndf_w3,nqp_h,nqp_v) :: w3_basis
  real(kind=r_def), intent(in), dimension(3,ndf_chi,nqp_h,nqp_v) :: chi_diff_basis
  real(kind=r_def), dimension(undf_w3), intent(in)    :: rho
  real(kind=r_def), dimension(undf_w3), intent(out)   :: mass
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v  

  !Internal variables
  integer               :: df1, df2, k
  integer               :: qp1, qp2

  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(ndf_w3) :: rho_e, mass_e
  real(kind=r_def) :: integrand

  ! compute the LHS integrated over one cell and solve
  do k = 0, nlayers-1
    do df1 = 1, ndf_chi
      chi_1_e(df1) = chi_1( map_chi(df1) + k)
      chi_2_e(df1) = chi_2( map_chi(df1) + k)
      chi_3_e(df1) = chi_3( map_chi(df1) + k)
    end do
    do df1 = 1, ndf_w3
      rho_e(df1) = rho( map_w3(df1) + k )
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e, chi_diff_basis, jac, dj)
    do df1 = 1, ndf_w3
       mass_e(df1) = 0.0_r_def
       do df2 = 1, ndf_w3
          do qp2 = 1, nqp_v
             do qp1 = 1, nqp_h
                integrand =  w3_basis(1,df1,qp1,qp2) * &
                             w3_basis(1,df2,qp1,qp2) * rho_e(df2) * dj(qp1,qp2)
                 mass_e(df1) = mass_e(df1) + wqp_h(qp1)*wqp_v(qp2)*integrand/scaled_radius**2
             end do
          end do
       end do
       mass( map_w3(df1) + k ) = mass_e(df1)
    end do
  end do

end subroutine compute_total_mass_code

end module compute_total_mass_kernel_mod
