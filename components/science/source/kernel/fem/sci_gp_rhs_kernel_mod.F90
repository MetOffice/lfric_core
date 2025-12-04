!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief Kernel which projects a field into into a given space

module sci_gp_rhs_kernel_mod
use kernel_mod,              only : kernel_type
use constants_mod,           only : r_def, i_def
use argument_mod,            only : arg_type, func_type,       &
                                    GH_FIELD, GH_REAL, GH_INC, &
                                    GH_READ, ANY_SPACE_9,      &
                                    ANY_SPACE_1, ANY_SPACE_2,  &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    GH_BASIS, GH_DIFF_BASIS,   &
                                    CELL_COLUMN, GH_QUADRATURE_XYoZ

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: gp_rhs_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                    &
       arg_type(GH_FIELD,   GH_REAL, GH_INC,  ANY_SPACE_1),              &
       arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_SPACE_2),              &
       arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),              &
       arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3) &
       /)
  type(func_type) :: meta_funcs(3) = (/                                  &
       func_type(ANY_SPACE_1, GH_BASIS),                                 &
       func_type(ANY_SPACE_2, GH_BASIS),                                 &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                   &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: gp_rhs_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: gp_rhs_code

contains

!> @brief     Subroutine to compute right hand side of a galerkin projection of
!>            a field from one space to  a different space
!> @details   Computes int( gamma * f  dx) to compute the right hand side of the
!>            galerkin projection of scalar field f into another space of which gamma
!>            is the test function
!! @param[in] nlayers Number of layers
!! @param[in,out] rhs Field containing the intergral of test_function * field
!! @param[in] field Field to be projected
!! @param[in] chi_1 1st coordinate field in Wchi
!! @param[in] chi_2 2nd coordinate field in Wchi
!! @param[in] chi_3 3rd coordinate field in Wchi
!! @param[in] panel_id Field giving the ID for mesh panels
!! @param[in] ndf Number of degrees of freedom per cell
!! @param[in] undf Number of (local) unique degrees of freedom of the field rhs
!! @param[in] map Dofmap for the cell at the base of the column
!! @param[in] basis Basis functions evaluated at quadrature points
!! @param[in] ndf_f Number of degrees of freedom per cell for the field to be projected
!! @param[in] undf_f Number of (local) unique degrees of freedom of the proj. field
!! @param[in] map_f Dofmap for the cell at the base of the column
!! @param[in] f_basis Basis functions evaluated at quadrature points
!! @param[in] ndf_chi Number of dofs per cell for the coordinate field
!! @param[in] undf_chi Number of (local) unique degrees of freedom of the chi field
!! @param[in] map_chi Dofmap for the coordinate field
!! @param[in] chi_basis Wchi basis functions evaluated at gaussian quadrature points
!! @param[in] chi_diff_basis Derivatives of Wchi basis functions
!!                           evaluated at gaussian quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine gp_rhs_code(nlayers,                       &
                       rhs, field,                    &
                       chi_1, chi_2, chi_3, panel_id, &
                       ndf, undf, map, basis,         &
                       ndf_f, undf_f, map_f, f_basis, &
                       ndf_chi, undf_chi, map_chi,    &
                       chi_basis, chi_diff_basis,     &
                       ndf_pid, undf_pid, map_pid,    &
                       nqp_h, nqp_v, wqp_h, wqp_v     )

  use sci_coordinate_jacobian_mod, only: coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, undf, ndf_pid, undf_pid
  integer(kind=i_def), intent(in) :: ndf_f, undf_f, ndf_chi, undf_chi
  integer(kind=i_def), intent(in) :: nqp_h, nqp_v

  integer(kind=i_def), dimension(ndf),     intent(in) :: map
  integer(kind=i_def), dimension(ndf_f),   intent(in) :: map_f
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  real(kind=r_def), intent(in), dimension(1,ndf,    nqp_h,nqp_v) :: basis
  real(kind=r_def), intent(in), dimension(1,ndf_f,  nqp_h,nqp_v) :: f_basis
  real(kind=r_def), intent(in), dimension(1,ndf_chi,nqp_h,nqp_v) :: chi_basis
  real(kind=r_def), intent(in), dimension(3,ndf_chi,nqp_h,nqp_v) :: chi_diff_basis

  real(kind=r_def), dimension(undf),     intent(inout) :: rhs
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid), intent(in)    :: panel_id
  real(kind=r_def), dimension(undf_f),   intent(in)    :: field
  real(kind=r_def), dimension(nqp_h),    intent(in)    ::  wqp_h
  real(kind=r_def), dimension(nqp_v),    intent(in)    ::  wqp_v

  ! Internal variables
  integer(kind=i_def)                          :: df, df2, k, qp1, qp2, ipanel
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jacobian
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_cell, chi_2_cell, chi_3_cell
  real(kind=r_def)                             :: rhs_cell

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    do df = 1, ndf_chi
      chi_1_cell(df) = chi_1( map_chi(df) + k )
      chi_2_cell(df) = chi_2( map_chi(df) + k )
      chi_3_cell(df) = chi_3( map_chi(df) + k )
    end do
    call coordinate_jacobian(ndf_chi,        &
                             nqp_h,          &
                             nqp_v,          &
                             chi_1_cell,     &
                             chi_2_cell,     &
                             chi_3_cell,     &
                             ipanel,         &
                             chi_basis,      &
                             chi_diff_basis, &
                             jacobian,       &
                             dj              )
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

end module sci_gp_rhs_kernel_mod
