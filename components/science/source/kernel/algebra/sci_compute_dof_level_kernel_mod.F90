!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel to compute the fractional level a dof lives on. This is given
!! by the layer index (k) + the nodal coordinate (hat{chi})
module sci_compute_dof_level_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type,            &
                                    GH_FIELD, GH_REAL,   &
                                    GH_INC, ANY_SPACE_1, &
                                    CELL_COLUMN
use constants_mod,           only : r_def, i_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: compute_dof_level_kernel_type
  private
  type(arg_type) :: meta_args(1) = (/                   &
       arg_type(GH_FIELD, GH_REAL, GH_INC, ANY_SPACE_1) &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: compute_dof_level_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: compute_dof_level_code
contains

!> @param[in]     nlayers Number of layers
!> @param[in]     ndf Number of degrees of freedom per cell for the output field
!> @param[in]     undf Number of unique degrees of freedom for the output field
!> @param[in]     map Dofmap for the cell at the base of the column for the output field
!> @param[in,out] level Fractional level of the dof's
!> @param[in]     nodes Nodal coordinates of the dofs
subroutine compute_dof_level_code(nlayers,                                  &
                                  level,                                    &
                                  ndf, undf, map,                           &
                                  nodes                                     &
                                  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf
  integer(kind=i_def), intent(in)                    :: undf
  integer(kind=i_def), dimension(ndf), intent(in)    :: map
  real(kind=r_def), dimension(undf),   intent(inout) :: level
  real(kind=r_def), dimension(3,ndf),  intent(in)    :: nodes

  ! Internal variables
  integer(kind=i_def) :: df, k

  do k = 0,nlayers - 1
    do df = 1,ndf
      level(map(df) + k) = real(k,r_def) + nodes(3,df)
    end do
  end do

end subroutine compute_dof_level_code

end module sci_compute_dof_level_kernel_mod
