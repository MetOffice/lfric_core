!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes LHS of Galerkin projection and solves equation in W3 space

module set_rho_kernel_mod

use argument_mod,               only : arg_type, func_type,            &
                                       GH_FIELD, GH_READ, GH_WRITE,    &
                                       ANY_SPACE_9, W3,                &
                                       GH_BASIS, GH_DIFF_BASIS,        &
                                       CELLS, QUADRATURE_XYoZ, GH_REAL
use constants_mod,              only : r_def
use idealised_config_mod,       only : test
use kernel_mod,                 only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: set_rho_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE,  W3),                            &
       arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9),                     &
       arg_type(GH_REAL,    GH_READ)                                   &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(W3, GH_BASIS),                                        &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                 &
       /)
  integer :: iterates_over = CELLS
  integer :: evaluator_shape = QUADRATURE_XYoZ
contains
  procedure, nopass ::set_rho_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface set_rho_kernel_type
   module procedure set_rho_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public set_rho_code
contains

type(set_rho_kernel_type) function set_rho_kernel_constructor() result(self)
  return
end function set_rho_kernel_constructor

!> @brief Computes LHS of Galerkin projection and solves equation in W3 space
!! @param[in] nlayers Number of layers
!! @param[in] ndf_w3 Number of degrees of freedom per cell
!! @param[in] undf_w3 Total number of degrees of freedom
!! @param[in] map_w3 Dofmap for the cell at the base of the column
!! @param[in] w3_basis Basis functions evaluated at gaussian quadrature points 
!! @param[inout] rho Density
!! @param[inout] time Time evaluated as a real value
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number of degrees of freedom for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] chi_diff_basis Basis functions evaluated at gaussian quadrature points
!! @param[inout] chi_1 X component of the chi coordinate field
!! @param[inout] chi_2 Y component of the chi coordinate field
!! @param[inout] chi_3 Z component of the chi coordinate field
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Weights of horizontal quadrature points
!! @param[in] wqp_v Weights of vertical quadrature points
subroutine set_rho_code(nlayers, rho, chi_1, chi_2, chi_3, time, &
                            ndf_w3, undf_w3, map_w3, w3_basis, &
                            ndf_chi, undf_chi, map_chi, chi_basis, chi_diff_basis, &
                            nqp_h, nqp_v, wqp_h, wqp_v )

   use matrix_invert_mod,             only : matrix_invert
   use coordinate_jacobian_mod,       only : coordinate_jacobian
   use analytic_density_profiles_mod, only : analytic_density
  ! needs to compute the integral of rho_df * P
  ! P_analytic over a single column

  !Arguments
  integer, intent(in) :: nlayers, ndf_w3, ndf_chi, undf_w3, undf_chi, nqp_h, nqp_v
  integer, dimension(ndf_w3), intent(in) :: map_w3
  integer, dimension(ndf_chi), intent(in) :: map_chi
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v),  intent(in)    :: w3_basis
  real(kind=r_def), dimension(undf_w3),               intent(inout) :: rho
  real(kind=r_def),                                   intent(in)    :: time
  real(kind=r_def), dimension(undf_chi),              intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in)    :: chi_diff_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in)    :: chi_basis
  real(kind=r_def), dimension(nqp_h),                 intent(in)    :: wqp_h
  real(kind=r_def), dimension(nqp_v),                 intent(in)    :: wqp_v

  !Internal variables
  integer               :: df1, df2, k
  integer               :: qp1, qp2

  real(kind=r_def), dimension(ndf_w3)          :: rho_e, rhs_e
  real(kind=r_def), dimension(ndf_w3,ndf_w3)   :: mass_matrix_w3, inv_mass_matrix_w3
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def)                             :: rho_ref, integrand
  real(kind=r_def)                             :: x(3)

  ! compute the RHS & LHS integrated over one cell and solve
  do k = 0, nlayers-1
    do df1 = 1, ndf_chi
      chi_1_e(df1) = chi_1( map_chi(df1) + k)
      chi_2_e(df1) = chi_2( map_chi(df1) + k)
      chi_3_e(df1) = chi_3( map_chi(df1) + k)
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e, &
                             chi_diff_basis, jac, dj)
! Compute RHS
    do df1 = 1, ndf_w3
      rhs_e(df1) = 0.0_r_def
      do qp2 = 1, nqp_v
        do qp1 = 1, nqp_h
          x(:) = 0.0_r_def
          do df2 = 1, ndf_chi
            x(1) = x(1) + chi_1_e(df2)*chi_basis(1,df2,qp1,qp2)
            x(2) = x(2) + chi_2_e(df2)*chi_basis(1,df2,qp1,qp2)
            x(3) = x(3) + chi_3_e(df2)*chi_basis(1,df2,qp1,qp2)
          end do
          rho_ref = analytic_density(x, test, time)

          integrand =  w3_basis(1,df1,qp1,qp2) * rho_ref * dj(qp1,qp2)
          rhs_e(df1) = rhs_e(df1) + wqp_h(qp1)*wqp_v(qp2)*integrand
        end do
      end do
    end do
! Commpute LHS
    do df1 = 1, ndf_w3
       do df2 = 1, ndf_w3
          mass_matrix_w3(df1,df2) = 0.0_r_def
          do qp2 = 1, nqp_v
             do qp1 = 1, nqp_h
                integrand =  w3_basis(1,df1,qp1,qp2) * &
                             w3_basis(1,df2,qp1,qp2) * dj(qp1,qp2)
                 mass_matrix_w3(df1,df2) = mass_matrix_w3(df1,df2) &
                                         + wqp_h(qp1)*wqp_v(qp2)*integrand
             end do
          end do
       end do
    end do
! Solve
    call matrix_invert(mass_matrix_w3,inv_mass_matrix_w3,ndf_w3)
    rho_e = matmul(inv_mass_matrix_w3,rhs_e)
    do df1 = 1,ndf_w3
      rho(map_w3(df1)+k) = rho_e(df1)
    end do
  end do

end subroutine set_rho_code

end module set_rho_kernel_mod
