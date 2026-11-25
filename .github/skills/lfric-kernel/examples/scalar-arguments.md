# Kernel with scalar arguments

A cell-column kernel can have scalar arguments in addition to fields. Scalars are
passed by value and represent runtime data that does not vary spatially within a
column. This example has one field argument and two scalar arguments.

| Code dummy | Access | Type |
| --- | --- | --- |
| `field` | `GH_READINC` | field on `ANY_DISCONTINUOUS_SPACE_1` |
| `min_threshold` | `GH_READ` | scalar real |
| `max_threshold` | `GH_READ` | scalar real |

```fortran
!> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
!!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
!>
!> @brief   Module for a threshold kernel with scalar arguments.
!> @details Demonstrates scalar argument patterns for runtime threshold values.
module threshold_kernel_mod

  use argument_mod,  only: arg_type, GH_FIELD, GH_REAL, GH_READINC, &
                           GH_SCALAR, ANY_DISCONTINUOUS_SPACE_1, CELL_COLUMN
  use constants_mod, only: i_def, r_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
  !!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
  !>
  !> @brief   Metadata for the threshold kernel.
  !> @details Demonstrates scalar real arguments for threshold control.
  type, public, extends(kernel_type) :: threshold_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                       &
         arg_type(GH_FIELD,  GH_REAL, GH_READINC, ANY_DISCONTINUOUS_SPACE_1), & ! field
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                               & ! min_threshold
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                                & ! max_threshold
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: threshold_code
  end type

  public :: threshold_code

contains

!> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
!!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
!>
!> @brief   Applies min and max thresholds to a field using scalar bounds.
!> @details Scalar arguments allow threshold values to be passed at runtime
!!          without requiring separate single-element fields.
!> @param[in]     nlayers        Number of layers
!> @param[in,out] field          Field values to threshold
!> @param[in]     min_threshold  Minimum allowed value
!> @param[in]     max_threshold  Maximum allowed value
!> @param[in]     ndf            Degrees of freedom per cell
!> @param[in]     undf           Total degrees of freedom
!> @param[in]     map            Cell dofmap
subroutine threshold_code(nlayers,        &
                          field,          &
                          min_threshold,  &
                          max_threshold,  &
                          ndf, undf, map)

  implicit none

  integer(i_def), intent(in)  :: nlayers
  integer(i_def), intent(in)  :: ndf
  integer(i_def), intent(in)  :: undf
  integer(i_def), intent(in)  :: map(ndf)
  real(r_def),    intent(inout) :: field(undf)
  real(r_def),    intent(in)  :: min_threshold
  real(r_def),    intent(in)  :: max_threshold

  ! Local loop counters
  integer(i_def) :: df, k

  do k=0, nlayers-1
    do df=1, ndf
      field(map(df) + k) = max(min_threshold, min(field(map(df) + k), &
                                                   max_threshold))
    end do
  end do

end subroutine threshold_code

end module threshold_kernel_mod
```
