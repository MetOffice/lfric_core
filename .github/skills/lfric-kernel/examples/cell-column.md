# Cell-column field kernel

This minimal example uses two fields on the same generic discontinuous
function space. It illustrates the common `nlayers`, field-data, and dofmap
interface only for that metadata shape. Replace the function space and access
modes to match the real algorithm and parallel semantics.

The loop assumes `nlayers` vertical values for each dof. That bound is not
universal; adapt it to the selected function space and verify it against a
matching kernel. The counted loops use `do index=start, end`, with no spaces
around `=`, one space after commas, and the simple upper-bound decrement
`nlayers-1` without spaces around `-`. Indent each nested level by two spaces
and align each `end do` with its opener.

```fortran
!> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
!!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
!>
!> @brief   Module for a cell-column field-copy kernel.
!> @details Defines metadata and code for copying cell-column field values.
module copy_field_kernel_mod

  use argument_mod,  only: arg_type, GH_FIELD, GH_REAL, GH_READ, GH_WRITE, &
                           ANY_DISCONTINUOUS_SPACE_1, CELL_COLUMN
  use constants_mod, only: i_def, r_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
  !!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
  !>
  !> @brief   Metadata for the field-copy kernel.
  !> @details Both fields use the same discontinuous function space.
  type, public, extends(kernel_type) :: copy_field_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                          &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! result
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1)  & ! source
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: copy_field_code
  end type

  public :: copy_field_code

contains

!> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
!!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
!>
!> @brief   Copies source values into the result field.
!> @details Uses the supplied map for each dof at every layer.
!> @param[in]  nlayers  Number of layers
!> @param[out] result   Destination field
!> @param[in]  source   Source field
!> @param[in]  ndf      Degrees of freedom per cell
!> @param[in]  undf     Total degrees of freedom
!> @param[in]  map      Cell dofmap
subroutine copy_field_code(nlayers, &
                           result,  &
                           source,  &
                           ndf, undf, map)

  implicit none

  integer(i_def), intent(in)  :: nlayers
  integer(i_def), intent(in)  :: ndf
  integer(i_def), intent(in)  :: undf
  integer(i_def), intent(in)  :: map(ndf)
  real(r_def),    intent(out) :: result(undf)
  real(r_def),    intent(in)  :: source(undf)

  ! Local loop counters
  integer(i_def) :: df, k

  do k=0, nlayers-1
    do df=1, ndf
      result(map(df) + k) = source(map(df) + k)
    end do
  end do

end subroutine copy_field_code

end module copy_field_kernel_mod
```

In this example, both metadata entries use the same function space, so the
support arguments share `ndf`, `undf`, and `map`. Different spaces, stencils,
operators, field vectors, or basis requirements change the generated interface.

