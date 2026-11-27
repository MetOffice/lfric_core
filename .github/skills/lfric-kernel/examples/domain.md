# Domain kernel example

This interface sketch shows a segment-based domain kernel where both fields
share one function space and `operates_on = DOMAIN`. The metadata is:

| Code dummy | Access | Function space |
| --- | --- | --- |
| `result` | `GH_WRITE` | `Wtheta` |
| `source` | `GH_READ` | `Wtheta` |

The corresponding code-procedure argument shape is:

```fortran
!> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
!!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
!>
!> @brief   Code-procedure interface for a segment-based domain kernel.
!> @details Shows support arguments for fields on the same function space.
!> @param[in]  nlayers  Number of vertical layers.
!> @param[in]  seg_len  Number of domain segments.
!> @param[out] result   Result field values.
!> @param[in]  source   Source field values.
!> @param[in]  ndf      Degrees of freedom per cell.
!> @param[in]  undf     Total degrees of freedom.
!> @param[in]  map      Segment dofmap.
subroutine domain_example_code(nlayers, &
                               seg_len, &
                               result,  &
                               source,  &
                               ndf, undf, map)

  implicit none

  integer(i_def), intent(in)  :: nlayers
  integer(i_def), intent(in)  :: seg_len
  integer(i_def), intent(in)  :: ndf
  integer(i_def), intent(in)  :: undf
  integer(i_def), intent(in)  :: map(ndf, seg_len)
  real(r_def),    intent(out) :: result(undf)
  real(r_def),    intent(in)  :: source(undf)

end subroutine domain_example_code
```

The example's `nlayers`, `seg_len`, and two-dimensional map distinguish it from
the cell-column and DOF examples. Domain kernels can have different maps for
multiple function spaces and may process vertical extents differently. Confirm
the configured API's restrictions, generated argument order, and data ownership
before turning this interface sketch into a kernel. The support-argument
section is specific to both fields using the same function space.

