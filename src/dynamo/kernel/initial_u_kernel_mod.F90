!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Kernel which projects the components of a vector field into a scalar space 

module initial_u_kernel_mod
use kernel_mod,              only : kernel_type
use constants_mod,           only : r_def, PI, earth_radius
use argument_mod,            only : arg_type, func_type,           &
                                    GH_FIELD, GH_INC, GH_READ,     &
                                    W0, W2,                        &
                                    GH_BASIS, GH_DIFF_BASIS,       &
                                    GH_ORIENTATION, CELLS
use constants_mod,           only : INITIAL_U_PROFILE, U0, V0, ROTATION_ANGLE
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: initial_u_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W2),                              &
       ARG_TYPE(GH_FIELD*3, GH_READ, W0)                               &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(W2, GH_BASIS, GH_ORIENTATION),                        &
       func_type(W0, GH_BASIS, GH_DIFF_BASIS)                          &
       /)
  integer :: iterates_over = CELLS

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
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf The number of degrees of freedom per cell
!! @param[in] map Integer array holding the dofmap for the cell at the base of the column
!! @param[in] basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points
!! @param[in] orientation The orientation array for the vector field
!! @param[inout] rhs Real array, the ths field to compute
!! @param[in] gq Type, gaussian quadrature rule
!! @param[in] ndf_chi the numbe rof dofs per cell for the coordinate field
!! @param[in] map_chi the dofmap for the coordinate field
!! @param[in] chi_diff_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points
!! @param[in] chi_1 Real array, the x component of the coordinate field
!! @param[in] chi_2 Real array, the y component of the coordinate field
!! @param[in] chi_3 Real array, the z component of the coordinate field
!! @param[in] wqp_h Real array. Quadrature weights horizontal
!! @param[in] wqp_v Real array. Quadrature weights vertical
subroutine initial_u_code(nlayers, &
                          rhs, & 
                          chi_1, chi_2, chi_3, &
                          ndf, undf, &
                          map, basis, &
                          orientation,                        &
                          ndf_chi, undf_chi, &
                          map_chi, chi_basis, chi_diff_basis, &
                          nqp_h, nqp_v, wqp_h, wqp_v &
                          )
                       
  use coordinate_jacobian_mod, only: coordinate_jacobian
  use slush_mod,               only: l_spherical
  use coord_transform_mod,     only: sphere2cart_vector, xyz2llr
  use analytic_profiles_mod,   only: analytic_wind
  
  !Arguments
  integer, intent(in) :: nlayers, ndf, ndf_chi
  integer, intent(in) :: undf, undf_chi
  integer, intent(in) :: nqp_h, nqp_v

  integer, dimension(ndf),     intent(in) :: map
  integer, dimension(ndf_chi), intent(in) :: map_chi
  integer, dimension(ndf),     intent(in) :: orientation

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
        if ( l_spherical ) then
          xyz(:) = 0.0_r_def
          do df = 1, ndf_chi
            xyz(1) = xyz(1) + chi_1_cell(df)*chi_basis(1,df,qp1,qp2)
            xyz(2) = xyz(2) + chi_2_cell(df)*chi_basis(1,df,qp1,qp2)
            xyz(3) = xyz(3) + chi_3_cell(df)*chi_basis(1,df,qp1,qp2)          
          end do
          call xyz2llr(xyz(1), xyz(2), xyz(3), llr(1), llr(2), llr(3))
          u_spherical = analytic_wind(llr, INITIAL_U_PROFILE, &
                                      2, (/U0, ROTATION_ANGLE/))
          u_physical = sphere2cart_vector(u_spherical,llr) 
        else
          u_physical = analytic_wind(xyz, INITIAL_U_PROFILE, &
                                     2, (/U0, V0/))
        end if
        do df = 1, ndf 
          integrand = dot_product(matmul(jacobian(:,:,qp1,qp2),&
                                         basis(:,df,qp1,qp2)), &
                                  u_physical)*real(orientation(df),r_def)
          rhs(map(df) + k) = rhs(map(df) + k) &
                           + wqp_h(qp1)*wqp_v(qp2)*integrand
        end do
      end do
    end do
  end do
  
end subroutine initial_u_code

end module initial_u_kernel_mod
