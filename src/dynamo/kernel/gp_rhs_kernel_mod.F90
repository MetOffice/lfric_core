!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Kernel which projects a field into into a given space 

module gp_rhs_kernel_mod
use kernel_mod,              only : kernel_type
use constants_mod,           only : r_def
use quadrature_mod,          only : quadrature_type
use argument_mod,            only : arg_type, func_type,           &
                                    GH_FIELD, GH_INC, GH_READ,     &
                                    ANY_SPACE_9, ANY_SPACE_1, ANY_SPACE_2,  &
                                    GH_BASIS, GH_DIFF_BASIS,       &
                                    CELLS

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: gp_rhs_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  ANY_SPACE_1),                     &
       ARG_TYPE(GH_FIELD,   GH_READ, ANY_SPACE_2),                     &
       ARG_TYPE(GH_FIELD*3, GH_READ, ANY_SPACE_9)                               &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(ANY_SPACE_1, GH_BASIS),                               &
       FUNC_TYPE(ANY_SPACE_2, GH_BASIS),                               &
       FUNC_TYPE(ANY_SPACE_9, GH_DIFF_BASIS)                           &
       /)
  integer :: iterates_over = CELLS

contains
  procedure, public, nopass :: gp_rhs_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface gp_rhs_kernel_type
   module procedure gp_rhs_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

type(gp_rhs_kernel_type) function gp_rhs_kernel_constructor() result(self)
  return
end function gp_rhs_kernel_constructor

!> @brief     Subroutine to compute right hand side of a galerkin projection of
!>            a field from one space to  a different space
!> @details   Computes int( gamma * f  dx) to compute the right hand side of the
!>            galerkin projection of scalar field f into another space of which gamma 
!>            is the test function
!! @param[in] nlayers Number of layers
!! @param[in] ndf Number of degrees of freedom per cell
!! @param[in] undf Number of (local) unique degrees of freedom of the field rhs
!! @param[in] map Dofmap for the cell at the base of the column
!! @param[in] basis Basis functions evaluated at quadrature points
!! @param[inout] rhs Field containing the intergral of test_function * field
!! @param[in] ndf_f Number of degrees of freedom per cell for the field to be projected
!! @param[in] undf_f Number of (local) unique degrees of freedom of the proj. field 
!! @param[in] map_f Dofmap for the cell at the base of the column
!! @param[in] f_basis Basis functions evaluated at quadrature points
!! @param[in] field Field to be projected
!! @param[in] ndf_chi Number of dofs per cell for the coordinate field
!! @param[in] undf_chi Number of (local) unique degrees of freedom of the chi field 
!! @param[in] map_chi Dofmap for the coordinate field
!! @param[in] chi_diff_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] chi_1 X component of the coordinate field
!! @param[in] chi_2 Y component of the coordinate field
!! @param[in] chi_3 Z component of the coordinate field
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine gp_rhs_code(nlayers, &
                       rhs, field, &
                       chi_1, chi_2, chi_3, &
                       ndf, undf, map, basis, &
                       ndf_f, undf_f, map_f, f_basis, &
                       ndf_chi, undf_chi, map_chi, chi_diff_basis, & 
                       nqp_h, nqp_v, wqp_h, wqp_v                    )
                       
  use coordinate_jacobian_mod, only: coordinate_jacobian                       
                         
  !Arguments
  integer,          intent(in) :: nlayers, ndf, undf 
  integer,          intent(in) :: ndf_f, undf_f, ndf_chi, undf_chi
  integer,          intent(in) :: nqp_h, nqp_v

  integer, dimension(ndf),     intent(in) :: map 
  integer, dimension(ndf_f),   intent(in) :: map_f
  integer, dimension(ndf_chi), intent(in) :: map_chi

  real(kind=r_def), intent(in), dimension(1,ndf,    nqp_h,nqp_v) :: basis 
  real(kind=r_def), intent(in), dimension(1,ndf_f,  nqp_h,nqp_v) :: f_basis 
  real(kind=r_def), intent(in), dimension(3,ndf_chi,nqp_h,nqp_v) :: chi_diff_basis 

  real(kind=r_def), dimension(undf),     intent(inout) :: rhs
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_f),   intent(in)    :: field  

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v
  
  !Internal variables
  integer                                      :: df, df2, k, qp1, qp2
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jacobian
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_cell, chi_2_cell, chi_3_cell
  real(kind=r_def)                             :: rhs_cell

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
    do df = 1, ndf
       do qp2 = 1, nqp_v
          do qp1 = 1, nqp_h
            rhs_cell = 0.0_r_def
            do df2 = 1, ndf_f
              rhs_cell = rhs_cell + f_basis(1,df2,qp1,qp2)*field(map_f(df2) + k)
            end do
            rhs(map(df) + k) = rhs(map(df) + k) &
                             + wqp_h(qp1)*wqp_v(qp2)*basis(1,df,qp1,qp2) &
                             * rhs_cell * dj(qp1,qp2)
          end do
       end do      
    end do
  end do
  
end subroutine gp_rhs_code

end module gp_rhs_kernel_mod
