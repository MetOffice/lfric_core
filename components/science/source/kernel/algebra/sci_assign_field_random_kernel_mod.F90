!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Assign random values to a field
!> @details Sets all field values to uniformly distributed random numbers
!> in the range [0,scale] for some specified scale.

module sci_assign_field_random_kernel_mod

use argument_mod,            only : arg_type,            &
                                    GH_FIELD, GH_REAL,   &
                                    GH_INC, ANY_SPACE_1, &
                                    CELL_COLUMN, GH_READ, &
                                    GH_SCALAR
use constants_mod,           only : r_single, r_double, i_def
use kernel_mod,              only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: assign_field_random_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                  &
    arg_type(GH_FIELD,  GH_REAL, GH_INC, ANY_SPACE_1), &
    arg_type(GH_SCALAR, GH_REAL, GH_READ )             &
  /)
  integer :: operates_on = CELL_COLUMN
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: assign_field_random_code

  ! Generic interface for real32 and real64 types
  interface assign_field_random_code
    module procedure  &
      assign_field_random_code_r_single, &
      assign_field_random_code_r_double
  end interface

contains

!> @brief Sets all field entries to random values
!> @param[in] nlayers Number of layers
!> @param[in,out] x Output data
!> @param[in] scale Output values are in range [0, scale]
!> @param[in] ndf Number of degrees of freedom per cell for the output field
!> @param[in] undf Unique number of degrees of freedom  for the output field
!> @param[in] map Dofmap for the cell at the base of the column for the output field

! R_SINGLE PRECISION
! ==================
subroutine assign_field_random_code_r_single(nlayers,       &
                                             x,             &
                                             scale,         &
                                             ndf, undf, map)
  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in)    :: nlayers
  integer(kind=i_def),                     intent(in)    :: undf, ndf
  integer(kind=i_def), dimension(ndf),     intent(in)    :: map
  real   (kind=r_single), dimension(undf), intent(inout) :: x
  real   (kind=r_single),                  intent(in)    :: scale

  ! Internal variables
  integer(kind=i_def)              :: df, k
  real(kind=r_single), dimension(ndf) :: random_values

  do k = 0, nlayers-1
    call random_number(random_values(:))
    do df = 1, ndf
      x(map(df)+k) = random_values(df) * scale
    end do
  end do

end subroutine assign_field_random_code_r_single

! R_DOUBLE PRECISION
! ==================
subroutine assign_field_random_code_r_double(nlayers,       &
                                             x,             &
                                             scale,         &
                                             ndf, undf, map)
  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in)    :: nlayers
  integer(kind=i_def),                     intent(in)    :: undf, ndf
  integer(kind=i_def),    dimension(ndf),  intent(in)    :: map
  real   (kind=r_double), dimension(undf), intent(inout) :: x
  real   (kind=r_double),                  intent(in)    :: scale

  ! Internal variables
  integer(kind=i_def)              :: df, k
  real(kind=r_double), dimension(ndf) :: random_values

  do k = 0, nlayers-1
    call random_number(random_values(:))
    do df = 1, ndf
      x(map(df)+k) = random_values(df) * scale
    end do
  end do

end subroutine assign_field_random_code_r_double

end module sci_assign_field_random_kernel_mod
