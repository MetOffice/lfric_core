!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you should have
! received as part of this distribution.
!-----------------------------------------------------------------------------

!> Module containing an interface to the usleep utility from C.
module sleep_mod

  use, intrinsic :: iso_c_binding, only : c_int, c_int32_t
  use constants_mod, only: i_def

  implicit none

  private
  public :: c_sleep

  !> Interface to the C usleep function
  interface
    function bind_c_usleep(useconds) bind(c, name='usleep')
        import :: c_int, c_int32_t
        implicit none
        integer(kind=c_int32_t), value :: useconds
        integer(kind=c_int)            :: bind_c_usleep
    end function bind_c_usleep
  end interface

contains

  !> Sleep for a given number of seconds
  !> @param seconds Number of seconds to sleep
  subroutine c_sleep(seconds)
      integer(i_def), intent(in) :: seconds
      integer(kind=c_int) :: rc
      rc = bind_c_usleep(int(seconds, kind=c_int32_t) * 1000000_c_int32_t )
  end subroutine c_sleep

end module sleep_mod
