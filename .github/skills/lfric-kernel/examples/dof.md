# DOF kernel

For `operates_on = DOF`, a field value is supplied to the code procedure as a
scalar for the current degree of freedom. Unlike a cell-column kernel, this
minimal example has no `nlayers`, dofmap, or `ndf` support arguments.

```fortran
!> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
!!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
!>
!> @brief   Module for a DOF field-copy kernel.
!> @details Provides a minimal copy kernel with one value per dof.
module copy_dof_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD, GH_REAL, GH_READ, GH_WRITE, DOF
  use constants_mod,     only: r_def
  use fs_continuity_mod, only: Wtheta
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
  !!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
  !>
  !> @brief   Metadata for the DOF field-copy kernel.
  !> @details Both field arguments use the Wtheta function space.
  type, public, extends(kernel_type) :: copy_dof_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                 &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta), & ! result
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta)  & ! source
         /)
    integer :: operates_on = DOF
  contains
    procedure, nopass :: copy_dof_code
  end type

  public :: copy_dof_code

contains

!> @note THESE COMMENTS WERE ADDED BY AI AND MUST BE REVIEWED FOR CORRECTNESS.
!!       REMOVE THIS COMMENT AFTER DEVELOPER REVIEW.
!>
!> @brief   Copies the source value to the result degree of freedom.
!> @details Each invocation receives scalar values for one dof.
!> @param[out] result  Result value
!> @param[in]  source  Source value
subroutine copy_dof_code(result, &
                         source)

  implicit none

  real(r_def), intent(out) :: result
  real(r_def), intent(in)  :: source

  result = source

end subroutine copy_dof_code

end module copy_dof_kernel_mod
```

This is a minimal scalar-field example. Check the configured API's DOF
restrictions before adding scalar arguments, fields, access modes, or function
spaces.

