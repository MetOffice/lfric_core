!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief A mock XIOS module used for unit testing containing dummy routines
!>        for XIOS subroutines used by LFRic's I/O subsystem
!>
module lfric_xios_mock_mod

  use, intrinsic :: iso_fortran_env,  &
                                only: real64
  use constants_mod,            only: i_def, str_def, l_def
  use lfric_xios_constants_mod, only: dp_xios
  use lfric_xios_mock_data_mod, only: xios_mock_data_type
  use log_mod,                  only:             &
                                  log_event,      &
                                  log_level_error

  implicit none

  private

  public :: xios_recv_field,      &
            xios_send_field,      &
            xios_get_domain_attr, &
            xios_get_axis_attr,   &
            xios_set_axis_attr,   &
            xios_is_valid_axis,   &
            xios_get_field_attr,  &
            xios_set_field_attr,  &
            xios_field_is_active, &
            xios_is_valid_file,                   &
            xios_is_defined_file_attr,            &
            xios_get_file_attr,                   &
            xios_is_valid_field,                  &
            xios_is_defined_field_attr,           &
            xios_zoom_axis,                       &
            xios_is_valid_zoom_axis,              &
            xios_get_handle,                      &
            xios_set_attr,                        &
            xios_setvar,                          &
            xios_getvar,                          &
            get_latest_data,                      &
            lfric_xios_mock_pull_in

type :: xios_zoom_axis
  ! intentionally left empty
end type xios_zoom_axis

!> Public mock XIOS data type used to hold data for read/write testing
type(xios_mock_data_type), public :: mock_xios_data

!> Private test data
real(real64), private :: var_var1

  contains


  !> Dummy function for pulling mock implementation into the build.
  function lfric_xios_mock_pull_in() result(pulled_in)
    implicit none
    integer(i_def) :: pulled_in
    pulled_in = 1
  end function lfric_xios_mock_pull_in

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

  !> Sets hard coded value for axis size based on the axis ID
  !> @param[in]     axis_id The ID of the axis to be tested
  !> @param[in,out] n_glo   The size of the test axis
  !> @param[in,out] value   The axis value
  subroutine xios_get_axis_attr(axis_id, n_glo, value)

    implicit none

    character(len=*),       intent(in)    :: axis_id
    integer(i_def), optional, intent(inout) :: n_glo
    real(real64), optional, intent(inout) :: value(:)

    ! This block returns mock sizes of the axis
    if (present(n_glo)) then
      if ( axis_id == "test_axis_99" ) then
        n_glo = 99
      else if ( axis_id == "vert_axis_half_levels" ) then
        n_glo = 3
      else if ( axis_id == "vert_axis_full_levels" ) then
        n_glo = 4
      else
        n_glo = 0
      end if
    end if

    ! This block returns mock values for the axis
    if (present(value)) then
      if ( axis_id == "test_axis_99" ) then
        value = 99
      else if ( axis_id == "vert_axis_half_levels" ) then
        value = 3
      else if ( axis_id == "vert_axis_full_levels" ) then
        value = 4
      else
        value = 0
      end if
    end if

  end subroutine xios_get_axis_attr

  !> Mock implementation for XIOS setter
  !> @param[in] axis_id The ID of the axis to be tested
  !> @param[in] n_glo   The size of the test axis
  subroutine xios_set_axis_attr(axis_id, n_glo)

    implicit none

    character(len=*),       intent(in) :: axis_id
    integer(i_def),         intent(in) :: n_glo

    ! intentionally left blank
  end subroutine xios_set_axis_attr

  !> Mockf implementation for XIOS checker
  !> @param[in] axis_id The ID of the axis to be tested
  function xios_is_valid_axis(axis_id) result(status)

    implicit none

    character(len=*),       intent(in)    :: axis_id
    logical(l_def) :: status
    status = .true.
  end function xios_is_valid_axis

  !> Accept values for field attributes based on the field ID
  !> @param[in]    field_id   The ID of the field to be tested
  !> @param[inout] enabled    The enabled flag of the test field
  subroutine xios_set_field_attr(field_id, enabled)
    implicit none
    character(*), intent(in) :: field_id
    logical(l_def), intent(in) :: enabled

    ! intentionally left bLank
  end subroutine xios_set_field_attr

  !> Sets hard coded values for field attributes based on the field ID
  !> @param[in]    field_id   The ID of the field to be tested
  !> @param[inout] axis_ref   The axis reference of the test field
  !> @param[inout] domain_ref The domain reference of the test field
  !> @param[inout] grid_ref   The grid reference of the test field
  !> @param[inout] enabled    The enabled flag of the test field
  subroutine xios_get_field_attr(field_id, axis_ref, domain_ref, grid_ref, enabled)

    implicit none

    character(len=*),             intent(in)    :: field_id
    character(str_def), optional, intent(inout) :: axis_ref
    character(str_def), optional, intent(inout) :: domain_ref
    character(str_def), optional, intent(inout) :: grid_ref
    logical(l_def), optional,     intent(inout) :: enabled

    if (present(enabled)) then
      enabled = field_id /= 'diag_field_disabled'
    end if

    if (present(grid_ref)) then
      select case(field_id)
      case('diag_field_test_grid')
        grid_ref = 'test_grid'
      case('diag_field_disabled')
        grid_ref = 'half_level_face_grid' ! arbitrary
      case('diag_field_no_flag')
        grid_ref = 'half_level_face_grid' ! arbitrary
      case('diag_field_half_level_face_grid')
        grid_ref = 'half_level_face_grid'
      case('diag_field_full_level_face_grid')
        grid_ref = 'full_level_face_grid'
      case('diag_field_half_level_edge_grid')
        grid_ref = 'half_level_edge_grid'
      case('diag_field_node_grid')
        grid_ref = 'node_grid'
      case default
        call log_event('no grid_ref in field: ' // field_id, log_level_error)
      end select
    end if

    if (present(domain_ref)) then
      select case(field_id)
      case('diag_field_test_domain')
        domain_ref = 'test_domain'
      case ('diag_field_face')
        domain_ref = 'face'
      case ('diag_field_face_tile')
        domain_ref = 'face'
      case('diag_field_face_rad')
        domain_ref = 'face'
      case default
        call log_event('no domain_ref in field: ' // field_id, log_level_error)
      end select
    end if

    if (present(axis_ref)) then
      select case(field_id)
      case('diag_field_test_axis')
        axis_ref = 'test_axis'
      case('diag_field_face_tile')
        axis_ref = 'soil_levels'
      case('diag_field_face_rad')
        axis_ref = 'radiation_levels'
      case default
        call log_event('no axis_ref in field: ' // field_id, log_level_error)
      end select
    end if

  end subroutine xios_get_field_attr

  !> @brief Sets hard coded values for activation status based on the field ID
  !> @param[in]              unique_id       XIOS id of the field
  !> @param[in]              at_current_step Status at current timestep?
  !> @return                    Boolean representing active/inactive status
  function xios_field_is_active(unique_id, at_current_timestep) result(active)
    implicit none
    character(*), intent(in) :: unique_id
    logical(l_def), intent(in) :: at_current_timestep
    logical(l_def) :: active
    if (unique_id == 'diag_field_disabled') then
      active = .false.
    else if (at_current_timestep) then
      active = unique_id /= 'diag_field_inactive_at_current_step'
    else
      active = .true.
    end if
  end function xios_field_is_active

  !> Obtains the most recent data added to the mock XIOS data object
  !> @param[in] latest_data Array of the latest data from the XIOS mock data
  !>                        object
  subroutine get_latest_data(latest_data)

    implicit none

    real(dp_xios), intent(inout) :: latest_data(:)

    real(dp_xios), allocatable :: mock_data(:,:)

    ! The mock_data is in XIOS 2D layered form but it is needed in 1D model
    ! form.
    !
    mock_data = mock_xios_data%get_data()
    latest_data = reshape(mock_data, &
                          (/size(mock_data, 1) * size(mock_data, 2)/))

  end subroutine get_latest_data

  !> Returns hard coded values for validity based on the file ID
  !> @param[in] file_id The ID of the file to be tested
  function xios_is_valid_file(file_id) result(valid)
    implicit none
    character(len=*),       intent(in)    :: file_id
    logical(l_def) :: valid
    select case (file_id)
    case ('lfric_diag')
      valid = .true.
    case ('lfric_diag_no_flag')
      valid = .true.
    case ('lfric_diag_enabled')
      valid = .true.
    case ('lfric_diag_disabled')
      valid = .true.
    case default
      valid = .false.
    end select
  end function xios_is_valid_file

  !> Sets hard coded values for existence of enabled flag based on file ID
  !> @param[in]    file_id The ID of the file to be tested
  !> @param[inout] enabled Existence status of enabled flag of the test file
  !> @param[inout] comment Existence status of comment field of the test file
  !> @param[inout] comment Existence status of name field of the test file
  subroutine xios_is_defined_file_attr(file_id, enabled, comment, name)
    implicit none
    character(len=*),       intent(in)    :: file_id
    logical(l_def),         optional, intent(inout) :: enabled
    logical(l_def),         optional, intent(inout) :: comment
    logical(l_def),         optional, intent(inout) :: name
    if (present(enabled)) enabled = file_id /= 'lfric_diag_no_flag'
    if (present(comment)) comment = .true.
    if (present(name)) name = (file_id == 'lfric_diag')
  end subroutine xios_is_defined_file_attr

  !> Returns hard coded values for validity based on field ID
  !> @param[in] field_id The ID of the field to be tested
  !> @return             Validity status
  function xios_is_valid_field(field_id) result(valid)
    implicit none
    character(len=*),       intent(in)    :: field_id
    logical(l_def) :: valid
    if (field_id(1:6) == 'diag__') then
      valid = .false.
    else
      valid = field_id /= 'diag_no_field'
    end if
  end function xios_is_valid_field

  !> Sets hard coded values for enabled flag based on file ID
  !> @param[in]    file_id The ID of the file to be tested
  !> @param[inout] enabled The enabled flag of the test file
  !> @param[inout] comment The comment describing the test file
  !> @param[inout] name The name of test file
  subroutine xios_get_file_attr(file_id, enabled, comment, name)
    implicit none
    character(len=*),       intent(in)    :: file_id
    logical(l_def),         optional, intent(inout) :: enabled
    character(str_def),     optional, intent(inout) :: comment
    character(len=*),       optional, intent(inout) :: name
    select case (file_id)
    case ('lfric_diag')
      if (present(enabled)) enabled = .true.
      if (present(comment)) comment = '[[diag]]'
      if (present(name)) name = file_id
    case ('lfric_diag_enabled')
      if (present(enabled)) enabled = .true.
      if (present(comment)) comment = ''
    case ('lfric_diag_disabled')
      if (present(enabled)) enabled = .false.
      if (present(comment)) comment = ''
    case ('lfric_diag_no_flag')
      if (present(enabled)) then
        call log_event('no enabled flag in file: ' // file_id, log_level_error)
      end if
      if (present(name)) name = 'lfric_diag'
    case default
      call log_event('unexpected file: ' // file_id, log_level_error)
    end select
  end subroutine xios_get_file_attr

  !> Sets hard coded values for existence of field attributes based on field ID
  !> @param[in]    field_id              ID of test field
  !> @param[inout, optional] axis_ref    Existence of axis ref of field
  !> @param[inout, optional] domain_ref  Existenc of domain ref of field
  !> @param[inout, optional] grid_ref    Existence of grid ref of field
  !> @param[inout, optional] enabled     Existence of enabled flag of field
  subroutine xios_is_defined_field_attr(field_id, axis_ref, domain_ref,       &
      grid_ref, enabled)
    implicit none
    character(len=*),         intent(in)    :: field_id
    logical(l_def), optional, intent(inout) :: axis_ref
    logical(l_def), optional, intent(inout) :: domain_ref
    logical(l_def), optional, intent(inout) :: grid_ref
    logical(l_def), optional, intent(inout) :: enabled
    if (                                                                      &
      field_id /= 'diag_field_half_level_face_grid' .and.                     &
      field_id /= 'diag_field_full_level_face_grid' .and.                     &
      field_id /= 'diag_field_half_level_edge_grid' .and.                     &
      field_id /= 'diag_field_node_grid' .and.                                &
      field_id /= 'diag_field_face' .and.                                     &
      field_id /= 'diag_field_face_tile' .and.                                &
      field_id /= 'diag_field_face_rad' .and.                                 &
      field_id /= 'diag_field_no_flag' .and.                                  &
      field_id /= 'diag_field_enabled' .and.                                  &
      field_id /= 'diag_field_disabled' .and.                                 &
      field_id /= 'diag_field' .and.                                          &
      field_id /= 'diag_field_test_axis' .and.                                &
      field_id /= 'diag_field_test_domain' .and.                              &
      field_id /= 'diag_field_test_grid'  .and.                               &
      field_id /= 'W0_field' .and.                                            &
      field_id /= 'W2H_field' .and.                                           &
      field_id /= 'WTheta_field' .and.                                        &
      field_id /= 'W3_field') then
      call log_event('xios_is_defined_field_attr - unexpected field: '        &
        // field_id, log_level_error)
    end if
    if (present(enabled)) then
      enabled = field_id /= 'diag_field_no_flag'
    end if
    if (present(axis_ref)) then
      axis_ref =                                                              &
        field_id == 'diag_field_test_axis' .or.                               &
        field_id == 'diag_field_face_tile'
    end if
    if (present(domain_ref)) then
      domain_ref =                                                            &
        field_id == 'diag_field_test_domain' .or.                             &
        field_id == 'diag_field_face' .or.                                    &
        field_id == 'diag_field_face_tile' .or.                               &
        field_id == 'diag_field_face_rad'
    end if
    if (present(grid_ref)) then
      grid_ref =                                                              &
        field_id == 'diag_field_test_grid' .or.                               &
        field_id == 'diag_field_half_level_face_grid' .or.                    &
        field_id == 'diag_field_full_level_face_grid' .or.                    &
        field_id == 'diag_field_half_level_edge_grid' .or.                    &
        field_id == 'diag_field_node_grid' .or.                               &
        field_id == 'diag_field_disabled' .or.                                &
        field_id == 'diag_field_no_flag'
    end if
  end subroutine xios_is_defined_field_attr

  !> Mock zoom axis existence test
  !> @param[in]    axis_id  ID of test axis
  !> @return                True iff test axis exists
  function xios_is_valid_zoom_axis(axis_id) result(status)
    implicit none
    character(*), intent(in) :: axis_id

    logical(l_def) :: status

    select case (axis_id)
    case('zoom1')
      status = .true.
    case default
      status = .false.
    end select
  end function xios_is_valid_zoom_axis

  !> Mock XIOS handle creator
  !> @param[in]    unique_id  ID of test object
  !> @param[out]   axis_hdl   Test zoom axis handle
  subroutine xios_get_handle(unique_id, axis_hdl)
    implicit none
    character(*), intent(in) :: unique_id
    type(xios_zoom_axis), optional, intent(inout) :: axis_hdl

    select case (unique_id)
    case ('zoom1')
      ! intentionally left empty
    case default
      call log_event('unknown XIOS id ' // unique_id, log_level_error)
    end select
  end subroutine xios_get_handle

  !> Mock XIOS mutator
  !> @param[in]    unique_id  ID of test object
  !> @param[out]   axis_hdl   Test zoom axis handle
  subroutine xios_set_attr(handle, begin, n)
    implicit none
    type(xios_zoom_axis), intent(inout) :: handle
    integer(i_def), optional, intent(in) :: begin
    integer(i_def), optional, intent(in) :: n
    ! intentionally left empty
  end subroutine xios_set_attr

  !> Sets value of a test variable
  !> @param[in]    var_id                ID of test variable
  !> @param[in]    value                 value to be set
  !> @return                             True iff the variable exists
  function xios_setvar(var_id, value) result(status)
    implicit none
    character(*), intent(in) :: var_id
    real(real64), intent(in) :: value

    logical(l_def) :: status

    select case (var_id)
      case ('var1')
        status = .true.
        var_var1 = value
      case default
        status = .false.
    end select
  end function xios_setvar

  !> Get value of a test variable
  !> @param[in]    var_id                ID of test avariable
  !> @param[in]    value                 Value read
  !> @return                             True iff the variable exists
  function xios_getvar(var_id, value) result(status)
    implicit none
    character(*), intent(in) :: var_id
    real(real64), intent(out) :: value

    logical(l_def) :: status

    select case (var_id)
      case ('var1')
        status = .true.
        value = var_var1
      case default
        status = .false.
    end select
  end function xios_getvar

end module lfric_xios_mock_mod
