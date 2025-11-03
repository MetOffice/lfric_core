!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module sci_pointwise_convert_xyz2llr_kernel_mod

use argument_mod,        only : arg_type, GH_FIELD,       &
                                GH_READWRITE, GH_READ,    &
                                GH_REAL, DOF, ANY_SPACE_1
use constants_mod,       only : r_def
use kernel_mod,          only : kernel_type

implicit none

!---------------------------------------------------------------------------
! Public types
!---------------------------------------------------------------------------

type, public, extends(kernel_type) :: pointwise_convert_xyz2llr_kernel_type
  private
  type(arg_type) :: meta_args(1) = (/                            &
        arg_type(GH_FIELD*3, GH_REAL, GH_READWRITE, ANY_SPACE_1) &
        /)
  integer :: operates_on = DOF
contains
  procedure, nopass :: pointwise_convert_xyz2llr_code
end type pointwise_convert_xyz2llr_kernel_type


!---------------------------------------------------------------------------
! Contained functions/subroutines
!---------------------------------------------------------------------------
public pointwise_convert_xyz2llr_code
contains

!--------------------------------------------------------------------
subroutine pointwise_convert_xyz2llr_code(coord_vec_1, coord_vec_2,   &
                                          coord_vec_3)
  use coord_transform_mod, only: xyz2llr
  implicit none

  real(kind=r_def), intent(inout) :: coord_vec_1, coord_vec_2, &
                                     coord_vec_3
  real(kind=r_def)                :: llr_1, llr_2, llr_3

  llr_1 = coord_vec_1
  llr_2 = coord_vec_2
  llr_3 = coord_vec_3

  call xyz2llr(coord_vec_1, coord_vec_2, coord_vec_3, &
               llr_1, llr_2, llr_3)

end subroutine pointwise_convert_xyz2llr_code

end module sci_pointwise_convert_xyz2llr_kernel_mod
