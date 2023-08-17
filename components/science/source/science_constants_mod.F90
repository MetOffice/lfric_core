!----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Universal scientific constants
!----------------------------------------------------------------------------

module science_constants_mod

  use constants_mod, only : r_def

  implicit none

  private
  public :: zero_C_in_K

  ! Celsius to Kelvin offset
  real(r_def), parameter :: zero_C_in_K = 273.15_r_def

end module science_constants_mod
