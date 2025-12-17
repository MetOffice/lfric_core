!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the multiplicity of the dofs for a field
!> @details Computes how many times each dof in a field is visited when looping
!>          over cells and all dof's asscociated with that cell
module sci_multiplicity_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type,            &
                                    GH_FIELD, GH_REAL,   &
                                    GH_INC, ANY_SPACE_1, &
                                    CELL_COLUMN
use constants_mod,           only : r_single, r_double, i_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: multiplicity_kernel_type
  private
  type(arg_type) :: meta_args(1) = (/                   &
       arg_type(GH_FIELD, GH_REAL, GH_INC, ANY_SPACE_1) &
       /)
  integer :: operates_on = CELL_COLUMN
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: multiplicity_code

  ! Generic interface for real32 and real64 types
  interface multiplicity_code
    module procedure  &
      multiplicity_code_r_single, &
      multiplicity_code_r_double
  end interface

contains

!> @brief Compute the multiplicity of a field (number of cells each DoF is shared by)
!! @param[in] nlayers Number of layers
!! @param[in,out] field Input/output field
!! @param[in] ndf Number of degrees of freedom per cell for the function space
!! @param[in] undf Number of unique degrees of freedom  for the function space
!! @param[in] map Dofmap for the cell at the base of the column for the function space

! R_SINGLE PRECISION
! ==================
subroutine multiplicity_code_r_single(nlayers, &
                                      field,   &
                                      ndf, undf, map)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf
  integer(kind=i_def), intent(in) :: undf
  integer(kind=i_def), dimension(ndf),  intent(in) :: map
  real(kind=r_single), dimension(undf), intent(inout) :: field

  ! Internal variables
  integer(kind=i_def) :: k, df

  do k = 0, nlayers - 1
    do df = 1,ndf
      field(map(df) + k) = field(map(df) + k) + 1.0_r_single
    end do
  end do

end subroutine multiplicity_code_r_single


! R_DOUBLE PRECISION
! ==================
subroutine multiplicity_code_r_double(nlayers, &
                                      field,   &
                                      ndf, undf, map)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf
  integer(kind=i_def), intent(in) :: undf
  integer(kind=i_def), dimension(ndf),  intent(in) :: map
  real(kind=r_double), dimension(undf), intent(inout) :: field

  ! Internal variables
  integer(kind=i_def) :: k, df

  do k = 0, nlayers - 1
    do df = 1,ndf
      field(map(df) + k) = field(map(df) + k) + 1.0_r_double
    end do
  end do

end subroutine multiplicity_code_r_double

end module sci_multiplicity_kernel_mod
