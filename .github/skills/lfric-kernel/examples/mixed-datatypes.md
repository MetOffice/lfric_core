# Kernel with multiple field datatypes

A cell-column kernel can combine real and integer fields in the same invocation.
Each field maintains its own datatype and access mode. Support arguments (ndf, undf, map)
are defined per function space and datatype.

This example has two real fields and one integer field on the same space:

| Code dummy | Access | Type |
| --- | --- | --- |
| `temperature` | `GH_READ` | real field on `ANY_DISCONTINUOUS_SPACE_1` |
| `pressure` | `GH_READ` | real field on `ANY_DISCONTINUOUS_SPACE_1` |
| `quality_flag` | `GH_WRITE` | integer field on `ANY_DISCONTINUOUS_SPACE_1` |

```fortran
!> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
!!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
!>
!> @brief   Module for a mixed-datatype quality-flag kernel.
!> @details Combines real input fields with integer output field.
module quality_flag_kernel_mod

  use argument_mod,  only: arg_type, GH_FIELD, GH_REAL, GH_INTEGER, &
                           GH_READ, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1, &
                           CELL_COLUMN
  use constants_mod, only: i_def, r_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
  !!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
  !>
  !> @brief   Metadata for the quality-flag kernel.
  !> @details Combines real and integer fields with separate support arguments.
  type, public, extends(kernel_type) :: quality_flag_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! temperature
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! pressure
         arg_type(GH_FIELD, GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)  & ! quality_flag
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: quality_flag_code
  end type

  public :: quality_flag_code

contains

!> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
!!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
!>
!> @brief   Assigns integer quality flags based on physical field values.
!> @details Demonstrates a kernel combining real input fields with an integer
!!          output field for quality control or categorization.
!> @param[in]  nlayers      Number of layers
!> @param[in]  temperature  Input temperature field
!> @param[in]  pressure     Input pressure field
!> @param[out] quality_flag Integer quality indicators
!> @param[in]  ndf_r        Degrees of freedom per cell for real fields
!> @param[in]  undf_r       Total degrees of freedom for real fields
!> @param[in]  map_r        Cell dofmap for real fields
!> @param[in]  ndf_i        Degrees of freedom per cell for integer field
!> @param[in]  undf_i       Total degrees of freedom for integer field
!> @param[in]  map_i        Cell dofmap for integer field
subroutine quality_flag_code(nlayers,      &
                             temperature,  &
                             pressure,     &
                             quality_flag, &
                             ndf_r, undf_r, map_r, &
                             ndf_i, undf_i, map_i)

  implicit none

  integer(i_def), intent(in)  :: nlayers
  integer(i_def), intent(in)  :: ndf_r, undf_r
  integer(i_def), intent(in)  :: ndf_i, undf_i
  integer(i_def), intent(in)  :: map_r(ndf_r)
  integer(i_def), intent(in)  :: map_i(ndf_i)
  real(r_def),    intent(in)  :: temperature(undf_r)
  real(r_def),    intent(in)  :: pressure(undf_r)
  integer(i_def), intent(out) :: quality_flag(undf_i)

  ! Local loop counters
  integer(i_def) :: df, k

  do k=0, nlayers-1
    do df=1, ndf_r

      quality_flag(map_i(df) + k) = 0_i_def

      if (temperature(map_r(df) + k) < 200.0_r_def) then
        quality_flag(map_i(df) + k) = 1_i_def
      end if

      if (pressure(map_r(df) + k) < 50000.0_r_def) then
        quality_flag(map_i(df) + k) = quality_flag(map_i(df) + k) + 10_i_def
      end if

    end do
  end do

end subroutine quality_flag_code

end module quality_flag_kernel_mod
```
