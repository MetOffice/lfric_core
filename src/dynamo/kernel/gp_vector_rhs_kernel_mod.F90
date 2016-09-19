!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Kernel which projects the components of a vector field into a scalar space 

module gp_vector_rhs_kernel_mod

use argument_mod,            only : arg_type, func_type,           &
                                    GH_FIELD, GH_INC, GH_READ,     &
                                    W0, W2, ANY_SPACE_1,           &
                                    ANY_SPACE_2, ANY_SPACE_9,      &
                                    GH_BASIS, GH_DIFF_BASIS,       &
                                    CELLS
use base_mesh_config_mod,    only : geometry, &
                                    base_mesh_geometry_spherical
use constants_mod,           only : r_def
use coordinate_jacobian_mod, only : coordinate_jacobian, &
                                    coordinate_jacobian_inverse
use coord_transform_mod,     only : cart2sphere_vector
use kernel_mod,              only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: gp_vector_rhs_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                  &
       arg_type(GH_FIELD*3, GH_INC,  ANY_SPACE_1),                     &
       ARG_TYPE(GH_FIELD,   GH_READ, ANY_SPACE_2),                     &
       ARG_TYPE(GH_FIELD*3, GH_READ, ANY_SPACE_9),                     &
       ARG_TYPE(GH_FIELD,   GH_READ, W2)                               &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(ANY_SPACE_1, GH_BASIS),                               &
       func_type(ANY_SPACE_2, GH_BASIS),                               &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                 &
       /)
  integer :: iterates_over = CELLS

contains
  procedure, public, nopass :: gp_vector_rhs_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface gp_vector_rhs_kernel_type
   module procedure gp_vector_rhs_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

type(gp_vector_rhs_kernel_type) function gp_vector_rhs_kernel_constructor() result(self)
  return
end function gp_vector_rhs_kernel_constructor

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
!! @param[in] undf Number of degrees of freedom
!! @param[in] map Dofmap for the cell at the base of the column
!! @param[in] basis Basis functions evaluated at gaussian quadrature points
!! @param[inout] rhs1 Field to compute
!! @param[inout] rhs2 Field to compute
!! @param[inout] rhs3 Field to compute
!! @param[in] ndf_f Number of degrees of freedom per cell for the field to be projected
!! @param[in] undf_f Number of degrees of freedom for the field to be projected
!! @param[in] map_f Dofmap for the cell at the base of the column
!! @param[in] f_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] field Field to be projected 
!! @param[in] w2_field W2_field needed to get function space components
!! @param[in] ndf_w2 Number of degrees of freedom per cell for the w2 field
!! @param[in] undf_w2 Total number of degrees of freedom for the w2 field
!! @param[in] map_w2 Dofmap for the cell at the base of the column
!! @param[in] ndf_chi Number of dofs per cell for the coordinate field
!! @param[in] undf_chi Total number of dofs for the coordinate field
!! @param[in] map_chi Dofmap for the coordinate field
!! @param[in] chi_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] chi_diff_basis Differential basis functions evaluated at gaussian quadrature points
!! @param[in] chi_1 X component of the coordinate field
!! @param[in] chi_2 Y component of the coordinate field
!! @param[in] chi_3 Z component of the coordinate field
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine gp_vector_rhs_code(nlayers, &
                              rhs1, rhs2, rhs3, field, &
                              chi_1, chi_2, chi_3, &
                              w2_field, &
                              ndf, undf, &
                              map, basis, &
                              ndf_f, undf_f, map_f, f_basis, &
                              ndf_chi, undf_chi, &
                              map_chi, chi_basis, chi_diff_basis, &
                              ndf_w2, undf_w2, map_w2, &
                              nqp_h, nqp_v, wqp_h, wqp_v &
                             )

  !Arguments
  integer, intent(in) :: nlayers, ndf, ndf_f, ndf_chi, ndf_w2
  integer, intent(in) :: undf, undf_f, undf_chi, undf_w2
  integer, intent(in) :: nqp_h, nqp_v

  integer, dimension(ndf),     intent(in) :: map
  integer, dimension(ndf_f),   intent(in) :: map_f
  integer, dimension(ndf_chi), intent(in) :: map_chi
  integer, dimension(ndf_w2),  intent(in) :: map_w2


  real(kind=r_def), intent(in), dimension(1,ndf,    nqp_h,nqp_v) :: basis 
  real(kind=r_def), intent(in), dimension(3,ndf_f,  nqp_h,nqp_v) :: f_basis 
  real(kind=r_def), intent(in), dimension(3,ndf_chi,nqp_h,nqp_v) :: chi_diff_basis
  real(kind=r_def), intent(in), dimension(1,ndf_chi,nqp_h,nqp_v) :: chi_basis

  real(kind=r_def), dimension(undf),     intent(inout) :: rhs1, rhs2, rhs3
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_f),   intent(in)    :: field
  real(kind=r_def), dimension(undf_w2),  intent(in)    :: w2_field


  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer                                      :: df, df2, k, qp1, qp2
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jacobian, jacobian_inv
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_cell, chi_2_cell, chi_3_cell
  real(kind=r_def), dimension(3)               :: u_at_quad, x_at_quad, u_physical
  real(kind=r_def)                             :: integrand
  logical                                      :: hdiv

  ! Check if this is hdiv (W2) field or a hcurl (W1) field
  if ( ndf_f == ndf_w2 ) then 
    hdiv = .true.
  else
    hdiv = .false.
  end if

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
    if ( .not. hdiv) call coordinate_jacobian_inverse(nqp_h, nqp_v, &
                                                      jacobian, dj, &
                                                      jacobian_inv)

    do df = 1, ndf
      do qp2 = 1, nqp_v
        do qp1 = 1, nqp_h
          ! Compute vector in computational space
          u_at_quad(:) = 0.0_r_def
          do df2 = 1,ndf_f
            u_at_quad(:) = u_at_quad(:) &
                         + f_basis(:,df2,qp1,qp2)*field(map_f(df2) + k)
          end do
          if ( hdiv ) then
            ! For W2 space
            u_at_quad(:) = matmul(jacobian(:,:,qp1,qp2),u_at_quad(:))/dj(qp1,qp2)
          else
            ! For W1 space
            u_at_quad(:) = matmul(transpose(jacobian_inv(:,:,qp1,qp2)),u_at_quad(:))
          end if
          ! Compute physical coordinate of quadrature point
          if ( geometry == base_mesh_geometry_spherical ) then
            x_at_quad(:) = 0.0_r_def
            do df2 = 1,ndf_chi
              x_at_quad(1) = x_at_quad(1) + chi_1_cell(df2)*chi_basis(1,df2,qp1,qp2)
              x_at_quad(2) = x_at_quad(2) + chi_2_cell(df2)*chi_basis(1,df2,qp1,qp2)
              x_at_quad(3) = x_at_quad(3) + chi_3_cell(df2)*chi_basis(1,df2,qp1,qp2)
            end do
            u_physical(:) = cart2sphere_vector(x_at_quad, u_at_quad)
          else
            u_physical(:) = u_at_quad(:)
          end if
          integrand = wqp_h(qp1)*wqp_v(qp2)*basis(1,df,qp1,qp2)*dj(qp1,qp2)
          rhs1(map(df) + k) = rhs1(map(df) + k) + integrand * u_physical(1)
          rhs2(map(df) + k) = rhs2(map(df) + k) + integrand * u_physical(2)
          rhs3(map(df) + k) = rhs3(map(df) + k) + integrand * u_physical(3)
        end do
      end do
    end do
  end do

end subroutine gp_vector_rhs_code

end module gp_vector_rhs_kernel_mod
