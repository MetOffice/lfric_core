# Kernel with integer arguments

A kernel can have integer fields or integer scalar arguments, useful for categorical
data, counts, or identifiers. Integer fields support `GH_READ`, `GH_WRITE`, and
`GH_READWRITE` access modes. Scalar integers are read-only.

This example has two integer field arguments and one integer scalar:

| Code dummy | Access | Type |
| --- | --- | --- |
| `source` | `GH_READ` | field on `ANY_DISCONTINUOUS_SPACE_1`, kind `i_def` |
| `result` | `GH_WRITE` | field on `ANY_DISCONTINUOUS_SPACE_1`, kind `i_def` |
| `offset` | `GH_READ` | scalar integer |

```fortran
!> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
!!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
!>
!> @brief   Module for a scalar integer argument kernel.
!> @details Demonstrates passing integer scalars and working with integer fields.
module scalar_int_copy_kernel_mod

  use argument_mod,  only: arg_type, GH_FIELD, GH_INTEGER, GH_READ, &
                           GH_WRITE, GH_SCALAR, ANY_DISCONTINUOUS_SPACE_1, &
                           CELL_COLUMN
  use constants_mod, only: i_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
  !!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
  !>
  !> @brief   Metadata for the integer-scalar kernel.
  !> @details Shows integer field and integer scalar argument patterns.
  type, public, extends(kernel_type) :: scalar_int_copy_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! result
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! source
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                              & ! offset
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: scalar_int_copy_code
  end type

  public :: scalar_int_copy_code

contains

!> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
!!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
!>
!> @brief   Copies source values and applies an integer offset.
!> @details Demonstrates integer field access and scalar integer parameters.
!> @param[in]  nlayers  Number of layers
!> @param[out] result   Destination integer field
!> @param[in]  source   Source integer field
!> @param[in]  offset   Integer scalar offset
!> @param[in]  ndf      Degrees of freedom per cell
!> @param[in]  undf     Total degrees of freedom
!> @param[in]  map      Cell dofmap
subroutine scalar_int_copy_code(nlayers, &
                                result,  &
                                source,  &
                                offset,  &
                                ndf, undf, map)

  implicit none

  integer(i_def), intent(in)  :: nlayers
  integer(i_def), intent(in)  :: ndf
  integer(i_def), intent(in)  :: undf
  integer(i_def), intent(in)  :: map(ndf)
  integer(i_def), intent(out) :: result(undf)
  integer(i_def), intent(in)  :: source(undf)
  integer(i_def), intent(in)  :: offset

  ! Local loop counters
  integer(i_def) :: df, k

  do k=0, nlayers-1
    do df=1, ndf
      result(map(df) + k) = source(map(df) + k) + offset
    end do
  end do

end subroutine scalar_int_copy_code

end module scalar_int_copy_kernel_mod
```
