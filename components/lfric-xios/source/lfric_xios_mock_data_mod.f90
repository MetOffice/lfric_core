!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief A mock XIOS data_object used for unit testing containing simple
!>        setter and getter routines which simulate the sending and receiving
!>        of data from XIOS
!>
module lfric_xios_mock_data_mod

  use constants_mod,        only: dp_xios

  implicit none

  private

  !> Container type to hold data for IO unit tests
  type, public :: xios_mock_data_type

    private

    !> The data values of the mock xios data object
    real(dp_xios), allocatable :: data(:,:)
  contains
    !> Setter for the mock XIOS data
    procedure, public :: set_data
    !> Getter for the mock XIOS data
    procedure, public :: get_data
    !> Finaliser for the mock XIOS data type
    final :: finalise
  end type xios_mock_data_type

  contains

  !> Takes input data and stores it as an attribute of the mock data object
  !> @param[in] self      The calling XIOS mock data object
  !> @result    xios_data The input data passed to the XIOS mock data object
  subroutine set_data(self, xios_data)

    implicit none

    class(xios_mock_data_type), intent(inout) :: self
    real(dp_xios),              intent(in)    :: xios_data(:,:)

    self%data = xios_data

  end subroutine set_data

  !> Takes data from mock data object and passes it out of the function
  !> @param[in] self      The calling XIOS mock data object
  !> @result    xios_data The resulting XIOS mock data array
  function get_data(self) result(xios_data)

    implicit none

    class(xios_mock_data_type), intent(in) :: self

    real(dp_xios), allocatable :: xios_data(:,:)

    xios_data = self%data

  end function get_data

  !> Finaliser for XIOS mock data object
  !> @param[in] self The calling XIOS mock data object
  subroutine finalise(self)

    implicit none

    type(xios_mock_data_type), intent(inout) :: self

    if ( allocated(self%data) ) then
      deallocate(self%data)
    end if

  end subroutine finalise

end module lfric_xios_mock_data_mod
