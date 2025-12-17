!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Assign field to a value only at a single column
!> @details Sets field to zero everwhere except for the first column, where they are
!>          set to one. The purpose of this kernel is to have a way of setting the field
!>          to zero almost everwhere to test the kernel in
!>          columnwise_op_asm_diag_hmht_kernel_mod.F90
!>          Note that in parallel this might not work, i.e. the values in
!>          several columns will be set to one as a module variable is used to determine
!>          whether the field has already been set.
!>

module sci_assign_field_single_column_kernel_mod

use argument_mod,            only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    GH_READWRITE,              &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    DOMAIN
use constants_mod,           only : r_single, r_double, i_def
use kernel_mod,              only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: assign_field_single_column_kernel_type
  private
  type(arg_type) :: meta_args(1) = (/                                       &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1) &
       /)
  integer :: operates_on = DOMAIN
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: assign_field_single_column_code

  ! Generic interface for real32 and real64 types
  interface assign_field_single_column_code
    module procedure  &
      assign_field_single_column_code_r_single, &
      assign_field_single_column_code_r_double
  end interface

contains

!> @brief Sets the values in the first column.
!> @param[in]     nlayers Number of vertical layers
!> @param[in]     ncols Number of columns
!> @param[in,out] x Output data
!> @param[in]     ndf Number of degrees of freedom per cell for the output field
!> @param[in]     undf Unique number of degrees of freedom  for the output field
!> @param[in]     map Dofmap for the cell at the base of the column for the output field

! R_SINGLE PRECISION
! ==================
subroutine assign_field_single_column_code_r_single(nlayers, ncols, &
                                                    x,              &
                                                    ndf, undf, map)
  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in)    :: nlayers
  integer(kind=i_def),                           intent(in)    :: ncols
  integer(kind=i_def),                           intent(in)    :: undf, ndf
  real   (kind=r_single), dimension(undf),       intent(inout) :: x
  integer(kind=i_def),    dimension(ndf, ncols), intent(in)    :: map

  ! Internal variables
  integer(kind=i_def) :: df, k

  ! Assign values in column number 1
  do k = 0, nlayers-1
    do df = 1, ndf
       x(map(df,1)+k) = 1.0_r_single
    end do
  end do

end subroutine assign_field_single_column_code_r_single

! R_DOUBLE PRECISION
! ==================
subroutine assign_field_single_column_code_r_double(nlayers, ncols, &
                                                    x,              &
                                                    ndf, undf, map)
  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in)    :: nlayers
  integer(kind=i_def),                           intent(in)    :: ncols
  integer(kind=i_def),                           intent(in)    :: undf, ndf
  real   (kind=r_double), dimension(undf),       intent(inout) :: x
  integer(kind=i_def),    dimension(ndf, ncols), intent(in)    :: map

  ! Internal variables
  integer(kind=i_def) :: df, k

  ! Assign values in column number 1
  do k = 0, nlayers-1
    do df = 1, ndf
       x(map(df,1)+k) = 1.0_r_double
    end do
  end do

end subroutine assign_field_single_column_code_r_double


end module sci_assign_field_single_column_kernel_mod
