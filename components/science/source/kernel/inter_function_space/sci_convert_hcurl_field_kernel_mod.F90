!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel to compute the apply the curl conforming Piola transform to a
!! computational vector field and return the 3 components of the physical field as
!! separate fields

module sci_convert_hcurl_field_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,       &
                                    GH_FIELD, GH_REAL,         &
                                    GH_READ, GH_INC,           &
                                    ANY_SPACE_9, ANY_SPACE_1,  &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    GH_DIFF_BASIS, GH_BASIS,   &
                                    CELL_COLUMN, GH_EVALUATOR
use constants_mod,           only : r_def, i_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: convert_hcurl_field_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                    &
       arg_type(GH_FIELD*3, GH_REAL, GH_INC,  ANY_SPACE_1),              &
       arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_SPACE_1),              &
       arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),              &
       arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3) &
       /)
  type(func_type) :: meta_funcs(2) = (/                                  &
       func_type(ANY_SPACE_1, GH_BASIS),                                 &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                   &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: convert_hcurl_field_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: convert_hcurl_field_code
contains

!> @param[in] nlayers Number of layers
!> @param[in,out] physical_field1 First component of the output field in physical units
!> @param[in,out] physical_field2 Second component of the  output field in physical units
!> @param[in,out] physical_field3 Third component of the  output field in physical units
!> @param[in] computational_field Input field in computational units
!> @param[in] chi_1 1st coordinate field in Wchi
!> @param[in] chi_2 2nd coordinate field in Wchi
!> @param[in] chi_3 3rd coordinate field in Wchi
!> @param[in] panel_id  Field giving the ID for mesh panels
!> @param[in] ndf Number of degrees of freedom per cell for the output field
!> @param[in] undf Number of unique degrees of freedom for the output field
!> @param[in] map Dofmap for the cell at the base of the column for the output field
!> @param[in] basis Basis functions of the output field evaluated at its nodal points
!> @param[in] ndf_chi Number of degrees of freedom per cell for the coordinate field
!> @param[in] undf_chi Number of unique degrees of freedom for the coordinate field
!> @param[in] map_chi Dofmap for the cell at the base of the column for the coordinate field
!> @param[in] basis_chi Basis functions of the coordinate space evaluated at the nodal points
!> @param[in] diff_basis_chi Differential basis functions of the coordinate space
!!                           evaluated at the nodal points
!> @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!> @param[in] undf_pid Number of unique degrees of freedom for panel_id
!> @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
subroutine convert_hcurl_field_code(nlayers,                                  &
                                    physical_field1,                          &
                                    physical_field2,                          &
                                    physical_field3,                          &
                                    computational_field,                      &
                                    chi1, chi2, chi3, panel_id,               &
                                    ndf, undf, map,                           &
                                    basis,                                    &
                                    ndf_chi, undf_chi, map_chi,               &
                                    basis_chi, diff_basis_chi,                &
                                    ndf_pid, undf_pid, map_pid                &
                                  )

  use sci_coordinate_jacobian_mod, only: coordinate_jacobian, coordinate_jacobian_inverse
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf, undf, ndf_chi, undf_chi, ndf_pid, undf_pid
  integer(kind=i_def), dimension(ndf),     intent(in) :: map
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  real(kind=r_def), dimension(undf),     intent(in)    :: computational_field
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_pid), intent(in)    :: panel_id
  real(kind=r_def), dimension(undf),     intent(inout) :: physical_field1
  real(kind=r_def), dimension(undf),     intent(inout) :: physical_field2
  real(kind=r_def), dimension(undf),     intent(inout) :: physical_field3

  real(kind=r_def), dimension(1,ndf_chi,ndf), intent(in) :: basis_chi
  real(kind=r_def), dimension(3,ndf_chi,ndf), intent(in) :: diff_basis_chi
  real(kind=r_def), dimension(3,ndf,ndf),     intent(in) :: basis

  ! Internal variables
  integer(kind=i_def) :: df, df2, k, ipanel
  real(kind=r_def) :: jacobian(3,3,ndf), jacobian_inv(3,3,ndf), dj(ndf)
  real(kind=r_def) :: vector_in(3), vector_out(3)
  real(kind=r_def), dimension(ndf_chi) :: chi1_e, chi2_e, chi3_e

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    do df = 1,ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi, ndf,chi1_e, chi2_e, chi3_e, &
                             ipanel, basis_chi, diff_basis_chi, jacobian, dj)
    call coordinate_jacobian_inverse(ndf, jacobian, dj, jacobian_inv)
    do df = 1,ndf
      vector_in(:) = 0.0_r_def
      do df2 = 1,ndf
        vector_in(:) = vector_in(:) + computational_field(map(df2)+k)*basis(:,df2,df)
      end do
      vector_out(:) = matmul(transpose(jacobian_inv(:,:,df)),vector_in)
      physical_field1(map(df)+k) = physical_field1(map(df)+k) + vector_out(1)
      physical_field2(map(df)+k) = physical_field2(map(df)+k) + vector_out(2)
      physical_field3(map(df)+k) = physical_field3(map(df)+k) + vector_out(3)
    end do
  end do

end subroutine convert_hcurl_field_code

end module sci_convert_hcurl_field_kernel_mod
