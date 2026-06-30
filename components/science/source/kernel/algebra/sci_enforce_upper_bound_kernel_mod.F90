!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Enforces a upper bound on a field. If any element of the field array
!>        is above this value, then it is set to the upper bound.
!>        This can be used to remove excessively large values or clip any field
!>        to any desired maximum value.
module sci_enforce_upper_bound_kernel_mod

  use argument_mod,  only : arg_type,                  &
                            GH_FIELD, GH_SCALAR,       &
                            GH_REAL, GH_READ,          &
                            GH_READWRITE,              &
                            ANY_DISCONTINUOUS_SPACE_1, &
                            CELL_COLUMN
  use constants_mod, only : i_def, r_def
  use, intrinsic :: iso_fortran_env, only: real32, real64
  use kernel_mod,    only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  ! The type declaration for the kernel. Contains the metadata needed by the
  ! Psy layer.
  !
  type, public, extends(kernel_type) :: enforce_upper_bound_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                                &
         arg_type(GH_SCALAR, GH_REAL, GH_READ )                                &
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: enforce_upper_bound_code

  ! Generic interface for real32 and real64 types
  interface enforce_upper_bound_code
    module procedure  &
      enforce_upper_bound_code_real32, &
      enforce_upper_bound_code_real64
  end interface

contains

!> @brief Returns field=max(field,upper bound)
!! @param[in] nlayers Number of layers
!! @param[in,out] field Field
!! @param[in] upper_bound The upper bound
!! @param[in] enforce_val The value which will be enforced if above upper_bound
!! @param[in] ndf Number of degrees of freedom per cell
!! @param[in] undf Total number of degrees of freedom
!! @param[in] map Dofmap for the cell at the base of the column


! REAL32 PRECISION
! ==================
subroutine enforce_upper_bound_code_real32(  nlayers, field, upper_bound, enforce_val,&
                                             ndf, undf, map)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, undf
  integer(kind=i_def), dimension(ndf), intent(in) :: map
  real(kind=real32),   dimension(undf), intent(inout) :: field
  real(kind=real32),   intent(in) :: upper_bound
  real(kind=real32),   intent(in) :: enforce_val

  ! Internal variables
  integer(kind=i_def) :: df, k

  do k = 0, nlayers-1
    do df = 1, ndf

      ! Clip field
      if (field(map(df)+k) > upper_bound) field(map(df)+k) = enforce_val

    end do
  end do

end subroutine enforce_upper_bound_code_real32

! REAL64 PRECISION
! ==================
subroutine enforce_upper_bound_code_real64(  nlayers, field, upper_bound, enforce_val,&
                                             ndf, undf, map)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, undf
  integer(kind=i_def), dimension(ndf), intent(in) :: map
  real(kind=real64),   dimension(undf), intent(inout) :: field
  real(kind=real64),   intent(in) :: upper_bound
  real(kind=real64),   intent(in) :: enforce_val

  ! Internal variables
  integer(kind=i_def) :: df, k

  do k = 0, nlayers-1
    do df = 1, ndf

      ! Clip field
      if (field(map(df)+k) > upper_bound) field(map(df)+k) = enforce_val

    end do
  end do

end subroutine enforce_upper_bound_code_real64

end module sci_enforce_upper_bound_kernel_mod
