# Kernel with logical arguments

A kernel can have logical scalar arguments to control conditional logic or enable
optional features. Logical scalars are always read-only. This example uses two
logical flags to control processing:

| Code dummy | Access | Type |
| --- | --- | --- |
| `field` | `GH_READINC` | field on `ANY_DISCONTINUOUS_SPACE_1` |
| `enable` | `GH_READ` | scalar logical |
| `strict` | `GH_READ` | scalar logical |

```fortran
!> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
!!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
!>
!> @brief   Module for a logical scalar argument kernel.
!> @details Demonstrates passing logical scalars to control branching logic.
module logical_scalar_kernel_mod

  use argument_mod,  only: arg_type, GH_FIELD, GH_REAL, GH_READINC, &
                           GH_SCALAR, GH_LOGICAL, ANY_DISCONTINUOUS_SPACE_1, &
                           CELL_COLUMN
  use constants_mod, only: i_def, r_def, l_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
  !!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
  !>
  !> @brief   Metadata for the logical-scalar kernel.
  !> @details Shows real field with logical scalar control parameters.
  type, public, extends(kernel_type) :: logical_scalar_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                          &
         arg_type(GH_FIELD,  GH_REAL,    GH_READINC, ANY_DISCONTINUOUS_SPACE_1), & ! field
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                               & ! enable
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                                & ! strict
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: logical_scalar_code
  end type

  public :: logical_scalar_code

contains

!> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
!!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
!>
!> @brief   Applies conditional scaling based on logical flags.
!> @details Demonstrates use of logical scalars to gate computation.
!> @param[in]     nlayers  Number of layers
!> @param[in,out] field    Field to conditionally modify
!> @param[in]     enable   Enable the scaling operation
!> @param[in]     strict   Use strict enforcement rules
!> @param[in]     ndf      Degrees of freedom per cell
!> @param[in]     undf     Total degrees of freedom
!> @param[in]     map      Cell dofmap
subroutine logical_scalar_code(nlayers, &
                               field,   &
                               enable,  &
                               strict,  &
                               ndf, undf, map)

  implicit none

  integer(i_def), intent(in)  :: nlayers
  integer(i_def), intent(in)  :: ndf
  integer(i_def), intent(in)  :: undf
  integer(i_def), intent(in)  :: map(ndf)

  real(r_def),    intent(inout) :: field(undf)
  logical(l_def), intent(in)    :: enable
  logical(l_def), intent(in)    :: strict

  ! Local loop counters
  integer(i_def) :: df, k

  if (enable) then
    do k=0, nlayers-1
      do df=1, ndf

        if (strict) then
          field(map(df) + k) = max(0.0_r_def, field(map(df) + k))
        else
          field(map(df) + k) = field(map(df) + k) * 0.5_r_def
        end if

      end do
    end do
  end if

end subroutine logical_scalar_code

end module logical_scalar_kernel_mod
```
