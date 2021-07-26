!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief A mock XIOS module used for unit testing containing dummy routines
!>        for XIOS subroutines used by LFRic's I/O subsystem
!>
module lfric_xios_mock_mod

  use constants_mod,            only: i_def, str_def, dp_xios
  use lfric_xios_mock_data_mod, only: xios_mock_data_type

  implicit none

  private
  public :: xios_recv_field,      &
            xios_send_field,      &
            xios_get_domain_attr, &
            xios_get_axis_attr,   &
            xios_get_field_attr,  &
            get_latest_data

!> Public mock XIOS data type used to hold data for read/write testing
type(xios_mock_data_type), public :: mock_xios_data

  contains

  !> Recieves data from the mock XIOS data object
  !> @param[in] field_id The ID of the field to be tested
  !> @param[in] dat      The data recieved from the XIOS mock data object
  subroutine xios_recv_field(field_id, dat)

    implicit none

    character(len=*),       intent(in)    :: field_id
    real(dp_xios),          intent(inout) :: dat(:)

    integer(i_def) :: i

    dat = real((/(i,i=1,size(dat))/), dp_xios)

  end subroutine xios_recv_field

  !> Sends data to the mock XIOS data object
  !> @param[in] field_id The ID of the field to be tested
  !> @param[in] dat      The data to be sent to the XIOS mock data object
  subroutine xios_send_field(field_id, dat)

    implicit none

    character(len=*),         intent(in) :: field_id
    real(dp_xios), optional,  intent(in) :: dat(:,:)

    call mock_xios_data%set_data(dat)

  end subroutine xios_send_field

  !> Sets hard coded values for domain size based on the domain ID
  !> @param[in] domain_id The ID of the domain to be tested
  !> @param[in] ni        The size of the test domain
  subroutine xios_get_domain_attr(domain_id, ni)

    implicit none

    character(len=*),       intent(in)    :: domain_id
    integer(i_def),         intent(inout) :: ni

    if ( domain_id == "node" ) then
      ni = 9
    else if ( domain_id == "edge" ) then
      ni = 18
    else if ( domain_id == "face" ) then
      ni = 9
    else
      ni = 0
    end if

  end subroutine xios_get_domain_attr

  !> Sets hard coded values for axis size based on the axis ID
  !> @param[in] axis_id The ID of the axis to be tested
  !> @param[in] n_glo   The size of the test axis
  subroutine xios_get_axis_attr(axis_id, n_glo)

    implicit none

    character(len=*),       intent(in)    :: axis_id
    integer(i_def),         intent(inout) :: n_glo

    if ( axis_id == "vert_axis_half_levels" ) then
      n_glo = 3
    else if ( axis_id == "vert_axis_full_levels" ) then
      n_glo = 4
    else
      n_glo = 0
    end if

  end subroutine xios_get_axis_attr

  !> Unused dummy routine required by read_time_data
  !> @param[in] field_id The ID of the field to be tested (unused)
  !> @param[in] unit     The units of the test field (unused)
  !> @param[in] axis_ref The axis reference of the test field (unused)
  subroutine xios_get_field_attr(field_id, unit, axis_ref)

    implicit none

    character(len=*),       intent(in)    :: field_id
    character(str_def),     intent(inout) :: unit
    character(str_def),     intent(inout) :: axis_ref

  end subroutine xios_get_field_attr

  !> Obtains the most recent data added to the mock XIOS data object
  !> @param[in] latest_data Array of the latest data from the XIOS mock data
  !>                        object
  subroutine get_latest_data(latest_data)

    implicit none

    real(dp_xios), intent(inout) :: latest_data(:)

    real(dp_xios), allocatable :: mock_data(:,:)

    mock_data = mock_xios_data%get_data()
    latest_data = mock_data(:,1)

  end subroutine get_latest_data

end module lfric_xios_mock_mod