!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module for the misc_meta_data class
!>
!> @details A key-value pair object to contain miscellaneous information about
!> a field. Extra meta data can be put in here to make the field compliant with
!> another convention, CMIP, for example.

module misc_meta_data_mod

  use constants_mod,          only: str_def

  implicit none

  private

  type, public :: misc_meta_data_type

    !> String containing the name of the type of data being stored
    character(str_def) :: name
    !> Value of the data
    character(str_def) :: value

  end type misc_meta_data_type

  interface misc_meta_data_type
      module procedure misc_meta_data_constructor
  end interface misc_meta_data_type

contains
  !> Construct a <code>meta_data_type</code> object.
  !>
  !> @param [in] name Describes what kind of data is being stored here
  !> @param [in] value The actual data to be stored

  !> @return self the meta_data object
  !> Constructor named shortened for use in meta_mod
  function misc_meta_data_constructor(name, value) result(self)

    implicit none

    character(*),           intent(in) :: name
    character(*),           intent(in) :: value

    type(misc_meta_data_type) :: self

    self%name = name
    self%value = value

  end function misc_meta_data_constructor

end module misc_meta_data_mod