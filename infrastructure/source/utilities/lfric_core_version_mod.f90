!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Module containing the lfric_core version string components and
!> character return function (for printing & labelling).

module lfric_core_version_mod

  use constants_mod,                 only: i_def, l_def, str_def

  implicit none

  private
  public lfric_core_version_char

  integer(i_def), public, parameter  :: lfric_core_major_version = 3
  integer(i_def), public, parameter  :: lfric_core_minor_version = 1
  integer(i_def), public, parameter  :: lfric_core_patch_version = 1
  logical(l_def), public, parameter  :: lfric_core_release_version = .false.
  character(2), parameter            :: prefix = 'vn'
  character(4), parameter            :: dev_suffix = '_dev'

contains

  !> Return a character representation of the current version
  !> of lfric_core.
  !> This will only include the patch version if it is greater than zero
  !> and will only skip the '_dev' suffix if the release logical is .true.
  !> e.g. 3.1, 3.1_dev, 3.1.1_dev, 3.2, 3.2_dev
  function lfric_core_version_char(input_major_version, input_minor_version, &
                                   input_patch_version, input_release_version) &
                                   result(core_version_char)

    character(str_def)                   :: core_version_char
    integer(i_def), intent(in), optional :: input_major_version, &
                                            input_minor_version, &
                                            input_patch_version
    logical(l_def), intent(in), optional :: input_release_version
    integer(i_def)                       :: lfric_major_version, &
                                            lfric_minor_version, &
                                            lfric_patch_version
    logical(l_def)                       :: lfric_release_version

    if (present(input_major_version)) then
      lfric_major_version = input_major_version
    else
      lfric_major_version = lfric_core_major_version
    end if
    if (present(input_minor_version)) then
      lfric_minor_version = input_minor_version
    else
      lfric_minor_version = lfric_core_minor_version
    end if
    if (present(input_patch_version)) then
      lfric_patch_version = input_patch_version
    else
      lfric_patch_version = lfric_core_patch_version
    end if
    if (present(input_release_version)) then
      lfric_release_version = input_release_version
    else
      lfric_release_version = lfric_core_release_version
    end if

    if (lfric_patch_version > 0) then
      if (lfric_release_version) then
        write (unit=core_version_char,fmt='(A,I0,A,I0,A,I0)') prefix, &
               lfric_major_version, '.', lfric_minor_version, '.', &
               lfric_patch_version
      else
        write (unit=core_version_char,fmt='(A,I0,A,I0,A,I0,A)') prefix, &
               lfric_major_version, '.', lfric_minor_version, '.', &
               lfric_patch_version, dev_suffix
      end if
    else
      if (lfric_release_version) then
        write (unit=core_version_char,fmt='(A,I0,A,I0)') prefix, &
               lfric_major_version, '.', lfric_minor_version
      else
        write (unit=core_version_char,fmt='(A,I0,A,I0,A)') prefix, &
               lfric_major_version, '.', lfric_minor_version, &
               dev_suffix
      end if
    end if

  end function lfric_core_version_char

end module lfric_core_version_mod
