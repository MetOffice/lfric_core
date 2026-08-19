!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Applies boundary conditions to a lma operator
!> @details Wrapper code for applying boundary conditions to a operator
module sci_enforce_operator_bc_kernel_mod

use, intrinsic :: iso_fortran_env, only: real32, real64

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type,                 &
                                    GH_OPERATOR, GH_REAL,     &
                                    GH_READWRITE,             &
                                    ANY_SPACE_1, ANY_SPACE_2, &
                                    CELL_COLUMN
use constants_mod,           only : i_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: enforce_operator_bc_kernel_type
  private
  type(arg_type) :: meta_args(1) = (/                                         &
       arg_type(GH_OPERATOR, GH_REAL, GH_READWRITE, ANY_SPACE_1, ANY_SPACE_2) &
       /)
  integer :: operates_on = CELL_COLUMN
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: enforce_operator_bc_code

  ! Generic interface for real32 and real64 types
  interface enforce_operator_bc_code
    module procedure  &
      enforce_operator_bc_code_real32, &
      enforce_operator_bc_code_real64
  end interface

contains

!> @brief Applies boundary conditions to an operator
!> @param[in] cell Horizontal cell index
!> @param[in] nlayers Number of layers
!> @param[in] ncell_3d Total number of cells
!> @param[in,out] op Operator data array to map from space 1 to space 2
!> @param[in] ndf1 Number of degrees of freedom per cell for to space
!> @param[in] ndf2 Number of degrees of freedom per cell for from space
!> @param[in] boundary_value Flags (= 0) for dofs that live on the
!!            vertical boundaries of the cell (=1 for other dofs)

! REAL32 PRECISION
! ==================
subroutine enforce_operator_bc_code_real32(  cell, nlayers,                   &
                                             ncell_3d, op,                    &
                                             ndf1, ndf2, boundary_value       &
                                            )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, cell, ncell_3d
  integer(kind=i_def), intent(in) :: ndf1, ndf2
  integer(kind=i_def), dimension(ndf1,2), intent(in) :: boundary_value

  real(kind=real32),   dimension(ncell_3d,ndf1,ndf2), intent(inout) :: op

  ! Local variables
  integer(kind=i_def) :: df, k, ik

  k = 1
  ik = (cell-1)*nlayers + k
  do df = 1,ndf1
    op(ik,df,:) = op(ik,df,:)*real(boundary_value(df,1), real32)
  end do
  k = nlayers
  ik = (cell-1)*nlayers + k
  do df = 1,ndf1
    op(ik,df,:) = op(ik,df,:)*real(boundary_value(df,2), real32)
  end do

end subroutine enforce_operator_bc_code_real32

! REAL64 PRECISION
! ==================
subroutine enforce_operator_bc_code_real64(  cell, nlayers,                   &
                                             ncell_3d, op,                    &
                                             ndf1, ndf2, boundary_value       &
                                            )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, cell, ncell_3d
  integer(kind=i_def), intent(in) :: ndf1, ndf2
  integer(kind=i_def), dimension(ndf1,2), intent(in) :: boundary_value

  real(kind=real64),   dimension(ncell_3d,ndf1,ndf2), intent(inout) :: op

  ! Local variables
  integer(kind=i_def) :: df, k, ik

  k = 1
  ik = (cell-1)*nlayers + k
  do df = 1,ndf1
    op(ik,df,:) = op(ik,df,:)*real(boundary_value(df,1), real64)
  end do
  k = nlayers
  ik = (cell-1)*nlayers + k
  do df = 1,ndf1
    op(ik,df,:) = op(ik,df,:)*real(boundary_value(df,2), real64)
  end do

end subroutine enforce_operator_bc_code_real64
end module sci_enforce_operator_bc_kernel_mod
