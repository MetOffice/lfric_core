!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief Kernel which projects the components of a vector field into a scalar space 

module initial_u_kernel_mod

use argument_mod,            only : arg_type, func_type,           &
                                    GH_FIELD, GH_INC, GH_READ,     &
                                    ANY_SPACE_9, W2,               &
                                    GH_BASIS, GH_DIFF_BASIS,       &
                                    CELLS, GH_QUADRATURE_XYoZ,     &
                                    QUADRATURE_XYoZ
use constants_mod,           only : r_def, PI
use kernel_mod,              only : kernel_type
use initial_wind_config_mod, only : profile, rotation_angle, u0, v0

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: initial_u_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W2),                              &
       ARG_TYPE(GH_FIELD*3, GH_READ, ANY_SPACE_9)                      &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(W2, GH_BASIS),                                        &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                 &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = QUADRATURE_XYoZ
contains
  procedure, public, nopass :: initial_u_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface initial_u_kernel_type
   module procedure initial_u_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

type(initial_u_kernel_type) function initial_u_kernel_constructor() result(self)
  return
end function initial_u_kernel_constructor

!> @brief Computes the right-hand-side of the Galerkin projection for a vector
!> field into a scalar space by decomposing the vector into orthogonal
!> components in cartesian or spherical polar coordinates.
!> @details Computes rhs_i = int (gamma * f_i dx) for a vector field f which  is
!>          decomposed into orthogonal components and a separate right hand side
!>          field is computed for each component, this allows a vector field to
!>          be projected into three separate scalar fields suitable for further
!>          manipulation
!! @param[in] nlayers Number of layers
!! @param[in] ndf Number of degrees of freedom per cell
!! @param[in] undf Total number of degrees of freedom
!! @param[in] map Dofmap for the cell at the base of the column
!! @param[in] basis Basis functions evaluated at gaussian quadrature points
!! @param[inout] rhs Right hand side field to compute
!! @param[in] ndf_chi Number of dofs per cell for the coordinate field
!! @param[in] undf_chi Total number of degrees of freedom
!! @param[in] map_chi Dofmap for the coordinate field
!! @param[in] chi_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] chi_diff_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] chi_1 X component of the coordinate field
!! @param[in] chi_2 Y component of the coordinate field
!! @param[in] chi_3 Z component of the coordinate field
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine initial_u_code(nlayers, &
                          rhs, & 
                          chi_1, chi_2, chi_3, &
                          ndf, undf, &
                          map, basis, &
                          ndf_chi, undf_chi, &
                          map_chi, chi_basis, chi_diff_basis, &
                          nqp_h, nqp_v, wqp_h, wqp_v &
                          )

  use analytic_wind_profiles_mod, only : analytic_wind
  use base_mesh_config_mod,       only : geometry, &
                                         base_mesh_geometry_spherical
  use coordinate_jacobian_mod,    only : coordinate_jacobian
  use coord_transform_mod,        only : sphere2cart_vector, xyz2llr

  !Arguments
  integer, intent(in) :: nlayers, ndf, ndf_chi
  integer, intent(in) :: undf, undf_chi
  integer, intent(in) :: nqp_h, nqp_v

  integer, dimension(ndf),     intent(in) :: map
  integer, dimension(ndf_chi), intent(in) :: map_chi

  real(kind=r_def), intent(in), dimension(3,ndf,    nqp_h,nqp_v) :: basis 
  real(kind=r_def), intent(in), dimension(3,ndf_chi,nqp_h,nqp_v) :: chi_diff_basis
  real(kind=r_def), intent(in), dimension(1,ndf_chi,nqp_h,nqp_v) :: chi_basis

  real(kind=r_def), dimension(undf),     intent(inout) :: rhs
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer               :: df, k, qp1, qp2
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jacobian
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_cell, chi_2_cell, chi_3_cell
  real(kind=r_def), dimension(3)               :: u_physical, u_spherical, xyz, llr
  real(kind=r_def)                             :: integrand
  real(kind=r_def), dimension(2)               :: option

  do k = 0, nlayers-1
    do df = 1, ndf_chi
      chi_1_cell(df) = chi_1( map_chi(df) + k)
      chi_2_cell(df) = chi_2( map_chi(df) + k)
      chi_3_cell(df) = chi_3( map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi, &
                             nqp_h, &
                             nqp_v, &
                             chi_1_cell, &
                             chi_2_cell, &
                             chi_3_cell, &
                             chi_diff_basis, &
                             jacobian, &
                             dj)
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        ! Compute analytical vector wind in physical space
        xyz(:) = 0.0_r_def
        do df = 1, ndf_chi
          xyz(1) = xyz(1) + chi_1_cell(df)*chi_basis(1,df,qp1,qp2)
          xyz(2) = xyz(2) + chi_2_cell(df)*chi_basis(1,df,qp1,qp2)
          xyz(3) = xyz(3) + chi_3_cell(df)*chi_basis(1,df,qp1,qp2)
        end do
        if ( geometry == base_mesh_geometry_spherical ) then
          call xyz2llr(xyz(1), xyz(2), xyz(3), llr(1), llr(2), llr(3))
          option = (/U0, rotation_angle/)
          u_spherical = analytic_wind(llr, profile, 2, option)
          u_physical = sphere2cart_vector(u_spherical,llr) 
        else
          option = (/U0, V0/)
          u_physical = analytic_wind(xyz, profile, 2, option)
        end if
        do df = 1, ndf 
          integrand = dot_product(matmul(jacobian(:,:,qp1,qp2),&
                                         basis(:,df,qp1,qp2)),u_physical)
          rhs(map(df) + k) = rhs(map(df) + k) &
                           + wqp_h(qp1)*wqp_v(qp2)*integrand
        end do
      end do
    end do
  end do

end subroutine initial_u_code

end module initial_u_kernel_mod
