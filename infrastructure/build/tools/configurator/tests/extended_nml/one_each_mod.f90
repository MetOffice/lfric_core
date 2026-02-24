!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the one_of_each namelist.
!>
module one_of_each_nml_mod

  use constants_mod, only: i_def, &
                           i_long, &
                           i_short, &
                           l_def, &
                           r_def, &
                           r_double, &
                           r_second, &
                           r_single, &
                           str_def, &
                           str_max_filename

  use namelist_mod, only: namelist_type

  implicit none

  private
  public :: one_of_each_nml_type

  type, extends(namelist_type) :: one_of_each_nml_type
    private
  contains

    procedure :: dint
    procedure :: dlog
    procedure :: dreal
    procedure :: dstr
    procedure :: enum
    procedure :: fstr
    procedure :: lint
    procedure :: lreal
    procedure :: sint
    procedure :: sreal
    procedure :: treal
    procedure :: vint
    procedure :: vreal
    procedure :: vstr

  end type one_of_each_nml_type

contains


  function dint(self) result(answer)

    implicit none

    class(one_of_each_nml_type), intent(in) :: self
    integer(i_def) :: answer

    call self%get_value('dint', answer)

  end function dint


  function dlog(self) result(answer)

    implicit none

    class(one_of_each_nml_type), intent(in) :: self
    logical(l_def) :: answer

    call self%get_value('dlog', answer)

  end function dlog


  function dreal(self) result(answer)

    implicit none

    class(one_of_each_nml_type), intent(in) :: self
    real(r_def) :: answer

    call self%get_value('dreal', answer)

  end function dreal


  function dstr(self) result(answer)

    implicit none

    class(one_of_each_nml_type), intent(in) :: self
    character(str_def) :: answer

    call self%get_value('dstr', answer)

  end function dstr


  function enum(self) result(answer)

    implicit none

    class(one_of_each_nml_type), intent(in) :: self
    integer(i_def) :: answer

    call self%get_value('enum', answer)

  end function enum


  function fstr(self) result(answer)

    implicit none

    class(one_of_each_nml_type), intent(in) :: self
    character(str_max_filename) :: answer

    call self%get_value('fstr', answer)

  end function fstr


  function lint(self) result(answer)

    implicit none

    class(one_of_each_nml_type), intent(in) :: self
    integer(i_long) :: answer

    call self%get_value('lint', answer)

  end function lint


  function lreal(self) result(answer)

    implicit none

    class(one_of_each_nml_type), intent(in) :: self
    real(r_double) :: answer

    call self%get_value('lreal', answer)

  end function lreal


  function sint(self) result(answer)

    implicit none

    class(one_of_each_nml_type), intent(in) :: self
    integer(i_short) :: answer

    call self%get_value('sint', answer)

  end function sint


  function sreal(self) result(answer)

    implicit none

    class(one_of_each_nml_type), intent(in) :: self
    real(r_single) :: answer

    call self%get_value('sreal', answer)

  end function sreal


  function treal(self) result(answer)

    implicit none

    class(one_of_each_nml_type), intent(in) :: self
    real(r_second) :: answer

    call self%get_value('treal', answer)

  end function treal


  function vint(self) result(answer)

    implicit none

    class(one_of_each_nml_type), intent(in) :: self
    integer(i_def) :: answer

    call self%get_value('vint', answer)

  end function vint


  function vreal(self) result(answer)

    implicit none

    class(one_of_each_nml_type), intent(in) :: self
    real(r_def) :: answer

    call self%get_value('vreal', answer)

  end function vreal


  function vstr(self) result(answer)

    implicit none

    class(one_of_each_nml_type), intent(in) :: self
    character(str_def) :: answer

    call self%get_value('vstr', answer)

  end function vstr

end module one_of_each_nml_mod
