!-------------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief XIOS wrapper for use in diagnostics
!> @details XIOS wrapper to be used in the initialisation of diagnostic fields

module lfric_xios_diag_mod

  use, intrinsic :: iso_fortran_env,  only: real64
  use constants_mod,                  only: l_def, i_def, str_def, r_def
#ifdef UNIT_TEST
  use lfric_xios_mock_mod,            only:                                   &
                                        xios_is_valid_file,                   &
                                        xios_is_defined_file_attr,            &
                                        xios_get_file_attr,                   &
                                        xios_is_valid_field,                  &
                                        xios_is_defined_field_attr,           &
                                        xios_get_field_attr,                  &
                                        xios_set_field_attr,                  &
                                        xios_field_is_active,                 &
                                        xios_get_axis_attr,                   &
                                        xios_set_axis_attr,                   &
                                        xios_is_valid_axis,                   &
                                        xios_zoom_axis,                       &
                                        xios_is_valid_zoom_axis,              &
                                        xios_get_handle,                      &
                                        xios_set_attr,                        &
                                        xios_setvar,                          &
                                        xios_getvar,                          &
                                        lfric_xios_mock_pull_in
#else
  use lfric_xios_mock_mod,            only: lfric_xios_mock_pull_in
  use xios,                           only:                                   &
                                        xios_is_valid_file,                   &
                                        xios_is_defined_file_attr,            &
                                        xios_get_file_attr,                   &
                                        xios_is_valid_field,                  &
                                        xios_is_defined_field_attr,           &
                                        xios_get_field_attr,                  &
                                        xios_set_field_attr,                  &
                                        xios_field_is_active,                 &
                                        xios_get_axis_attr,                   &
                                        xios_set_axis_attr,                   &
                                        xios_is_valid_axis,                   &
                                        xios_zoom_axis,                       &
                                        xios_is_valid_zoom_axis,              &
                                        xios_get_handle,                      &
                                        xios_set_attr,                        &
                                        xios_setvar,                          &
                                        xios_getvar
#endif

  use log_mod,                         only:                                  &
                                        log_event,                            &
                                        log_level_error,                      &
                                        log_level_warning,                    &
                                        log_scratch_space

  implicit none

  private

  public ::                                                                    &
    file_is_enabled,                                                           &
    file_is_tagged,                                                            &
    get_file_name,                                                             &
    field_is_valid,                                                            &
    field_is_enabled,                                                          &
    field_is_active,                                                           &
    enable_field,                                                              &
    get_field_order,                                                           &
    get_field_grid_ref,                                                        &
    get_field_domain_ref,                                                      &
    get_field_axis_ref,                                                        &
    get_axis_dimension,                                                        &
    get_axis_values,                                                           &
    set_axis_dimension,                                                        &
    set_zoom_axis_attr,                                                        &
    set_variable,                                                              &
    get_variable

contains

  !> @brief Return true if an only if an XIOS output file is enabled.
  !> @param[in]    unique_id  XIOS id of the file
  !> @return                  Boolean representing enabled/disabled status
  function file_is_enabled(file_id) result(enabled)
    implicit none
    character(*), intent(in) :: file_id
    logical(l_def) :: has_enabled_flag
    logical(l_def) :: enabled
    if (.not. xios_is_valid_file(file_id)) then
      write(log_scratch_space, '(A, A)')                                      &
        'Invalid XIOS file:', file_id
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if
    call xios_is_defined_file_attr(file_id, enabled=has_enabled_flag)
    if (has_enabled_flag) then
       call xios_get_file_attr(file_id, enabled=enabled)
       ! explicitly enabled/disabled
    else
       enabled = .true. ! enabled by default
    end if
  end function file_is_enabled

  !> @brief Scan comment field of a file for (boolean) tag, e.g, [[diag]]
  !> @param[in]    file_id      XIOS id of the field
  !> @param[in]    tag          Tag string, e.g., "diag"
  !> @return                    True if and only if the file is tagged
  function file_is_tagged(file_id, tag) result(tagged)
    implicit none
    character(len=*), intent(in) :: file_id
    character(len=*), intent(in) :: tag
    logical(l_def) :: tagged
    logical(l_def) :: def
    character(str_def) :: comment
    character(str_def) :: pattern
    if (.not. xios_is_valid_file(file_id)) then
      write(log_scratch_space, '(A, A)')                                      &
        'Invalid XIOS file:', file_id
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if
    tagged = .false.
    call xios_is_defined_file_attr(file_id, comment=def)
    if (def) then
      call xios_get_file_attr(file_id, comment=comment)
      pattern = '[[' // tag //']]'
      tagged = (index(comment, pattern) /= 0)
    end if
  end function file_is_tagged

  !> @brief Return the name of an XIOS file
  !> @param[in]    file_id      XIOS id of the file
  !> @return                    The file name
  function get_file_name(file_id) result(name)
    implicit none
    character(len=*), intent(in) :: file_id
    logical(l_def) :: def
    character(str_def) :: name
    if (.not. xios_is_valid_file(file_id)) then
      write(log_scratch_space, '(A, A)')                                      &
        'Invalid XIOS file:', file_id
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if
    call xios_is_defined_file_attr(file_id, name=def)
    if (def) then
      call xios_get_file_attr(file_id, name=name)
    else
      name = file_id
    end if
  end function get_file_name

  !> @brief Return true if an only if an XIOS field id is valid.
  !> @param[in]    unique_id    XIOS id of the field
  !> @return                    Boolean representing validity
  function field_is_valid(unique_id) result(valid)
    implicit none
    character(*), intent(in) :: unique_id
    logical(l_def) :: valid
    valid = xios_is_valid_field(unique_id)
  end function field_is_valid

  !> @brief Return true if an only if an XIOS output field is enabled.
  !> @param[in]              unique_id    XIOS id of the field
  !> @return                    Boolean representing enabled/disabled status
  function field_is_enabled(unique_id) result(enabled)
    implicit none
    character(*), intent(in) :: unique_id
    logical(l_def) :: has_enabled_flag
    logical(l_def) :: enabled
    if (.not. xios_is_valid_field(unique_id)) then
      write(log_scratch_space, '(A, A)')                                      &
      'Invalid XIOS field:', unique_id
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    else
      call xios_is_defined_field_attr(unique_id, enabled=has_enabled_flag)
      if (has_enabled_flag) then
        call xios_get_field_attr(unique_id, enabled=enabled)
        ! explicitly enabled/disabled
      else
        enabled = .true. ! enabled by default
      end if
    end if
  end function field_is_enabled

  !> @brief Return true if and only if an XIOS field is active.
  !> @param[in]              unique_id            XIOS id of the field
  !> @param[in]              at_current_timestep  Active at current timestep?
  !> @return                 Boolean representing active/inactive status
  function field_is_active(unique_id, at_current_timestep) result(active)
    implicit none
    character(*), intent(in) :: unique_id
    logical(l_def), intent(in) :: at_current_timestep
    logical(l_def) :: active
    if (.not. xios_is_valid_field(unique_id)) then
      write(log_scratch_space, '(A, A)')                                      &
      'Invalid XIOS field:', unique_id
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    else
      active = xios_field_is_active(unique_id, at_current_timestep)
    end if
  end function field_is_active

  !> @brief Set the enabled flag of an XIOS field to true
  !> @param[in]              unique_id    XIOS id of the field
  subroutine enable_field(unique_id)
    implicit none
    character(*), intent(in) :: unique_id
    if (.not. xios_is_valid_field(unique_id)) then
      write(log_scratch_space, '(A, A)')                                      &
      'Invalid XIOS field:', unique_id
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    else
      call xios_set_field_attr(unique_id, enabled=.true.)
    end if
  end subroutine enable_field

  !> @brief Return the interpolation order of a field.
  !> @param[in]    unique_id    XIOS id of the field
  !> @param[in]    force_order  Override order with
  !> @return                    The interpolation order
  function get_field_order(unique_id, force_order) result(order)
   implicit none
   character(*), intent(in) :: unique_id
   integer(i_def), optional, intent(in) :: force_order
   integer(i_def) :: order
   order = 0 ! FOR NOW, to be extracted from the field comment?
   if (present(force_order)) order = force_order
  end function get_field_order

  !> @brief Return the XIOS grid reference of a field.
  !> @param[in]    unique_id    XIOS id of the field
  !> @return                    The grid reference
  function get_field_grid_ref(unique_id) result(grid_ref)
    implicit none
    character(*), intent(in) :: unique_id
    character(str_def) :: grid_ref
    logical(l_def) :: has_grid
    call xios_is_defined_field_attr(unique_id, grid_ref=has_grid)
    if (has_grid) then
       call xios_get_field_attr(unique_id, grid_ref=grid_ref)
    else
       grid_ref = ""
    end if
  end function get_field_grid_ref

  !> @brief Return the XIOS domain reference of a field.
  !> @param[in]    unique_id    XIOS id of the field
  !> @return                    The domain reference
  function get_field_domain_ref(unique_id) result(domain_ref)
    implicit none
    character(*), intent(in) :: unique_id
    character(str_def) :: domain_ref
    logical(l_def) :: has_domain
    call xios_is_defined_field_attr(unique_id, domain_ref=has_domain)
    if (has_domain) then
       call xios_get_field_attr(unique_id, domain_ref=domain_ref)
    else
       domain_ref = ""
    end if
  end function get_field_domain_ref

  !> @brief Return the XIOS axis reference of a field.
  !> @param[in]    unique_id    XIOS id of the field
  !> @return                    The axis reference
  function get_field_axis_ref(unique_id) result(axis_ref)
    implicit none
    character(*), intent(in) :: unique_id
    character(str_def) :: axis_ref
    logical(l_def) :: has_axis
    call xios_is_defined_field_attr(unique_id, axis_ref=has_axis)
    if (has_axis) then
       call xios_get_field_attr(unique_id, axis_ref=axis_ref)
    else
       axis_ref = ""
    end if
  end function get_field_axis_ref

  !> @brief Return the dimension of an XIOS axis object.
  !> @param[in]    unique_id    XIOS id of the axis
  !> @return                    The dimension
  function get_axis_dimension(unique_id) result(dim)
    implicit none
    character(*), intent(in) :: unique_id
    integer(i_def) :: dim
    call xios_get_axis_attr(unique_id, n_glo=dim)
  end function get_axis_dimension

  !> @brief Return the values in an XIOS axis object.
  !> @param[in]    unique_id    XIOS id of the axis
  !> @param[in]    nval         Dimension of the axis object
  !> @return                    Axis values
  function get_axis_values(unique_id, nval) result(value)
    implicit none
    character(*), intent(in) :: unique_id
    integer(i_def), intent(in) :: nval
    real(real64) :: value(nval)
    call xios_get_axis_attr(unique_id, value=value)
  end function get_axis_values

  !> @brief Set the dimension of an XIOS axis object.
  !> @param[in]    unique_id    XIOS id of the axis
  !> @param[in]    dim          Dimension value to be set
  !> @param[in]    tolerant     Ignore missing axes?
  subroutine set_axis_dimension(unique_id, dim, tolerant)
    implicit none
    character(*), intent(in) :: unique_id
    integer(i_def), intent(in) :: dim
    logical(l_def), optional, intent(in) :: tolerant
    logical(l_def) :: strict
    strict = .true.
    if (present(tolerant)) strict = .not. tolerant
    if (xios_is_valid_axis(unique_id)) then
      call xios_set_axis_attr(unique_id, n_glo=dim)
    else
      if (strict) then
        write(log_scratch_space, '(A, A)')                                    &
        'Invalid XIOS axis:', unique_id
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if
    end if
  end subroutine set_axis_dimension

  !> @brief Set the dimension of an XIOS zoom_axis object.
  !> @param[in]    unique_id    XIOS id of the zoom_axis
  !> @param[in]    begin        Begin of zoom range
  !> @param]in]    n            Dimension of zoom range
  !> @param[in]    tolerant     Ignore missing axes?
  subroutine set_zoom_axis_attr(unique_id, begin, n, tolerant)
    implicit none
    character(*), intent(in) :: unique_id
    integer(i_def), intent(in) :: begin
    integer(i_def), intent(in) :: n
    logical(l_def), optional, intent(in) :: tolerant
    type(xios_zoom_axis) :: zoom_axis_hdl
    logical(l_def) :: strict
    strict = .true.
    if (present(tolerant)) strict = .not. tolerant
    if (xios_is_valid_zoom_axis(unique_id)) then
      call xios_get_handle(unique_id, zoom_axis_hdl)
      call xios_set_attr(zoom_axis_hdl, begin=begin, n=n)
    else
      if (strict) then
        write(log_scratch_space, '(A, A)')                                    &
        'Invalid XIOS zoom_axis:', unique_id
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if
    end if
  end subroutine set_zoom_axis_attr

  !> @brief Set the value of a real XIOS variable
  !> @param[in]              unique_id    XIOS id of the variable
  !> @param[in]              value        Value to be set
  !> @param[in, optional]    tolerant     Tolerate missing variable?
  subroutine set_variable(unique_id, value, tolerant)
    implicit none
    character(*), intent(in) :: unique_id
    real(real64), intent(in) :: value
    logical(l_def), optional, intent(in) :: tolerant
    ! The XIOS call was ifdef'ed out to avoid having to write a mocking
    ! version of the XIOS routine called. Such a mocking version would
    ! currently add little value. Note that the error handling could not
    ! be tested, because it is not possible to regain control after logging
    ! an error, cf. ticket #3887.
    logical(l_def) :: strict
    integer(i_def) :: level
    strict = .true.
    if (present(tolerant)) strict = .not. tolerant
    if (.not. xios_setvar(unique_id, value)) then
      if (strict) then
        level = log_level_error
      else
        level = log_level_warning
      end if
      call log_event('Failed to set variable ' // trim(unique_id), level)
    end if
  end subroutine set_variable

  !> @brief Get the value of a real XIOS variable
  !> @param[in]              var_id    XIOS id of the variable
  !> @param[in, optional]    dflt      Value to return if no such variable
  !> @return                           Value of the variable
  function get_variable(unique_id, dflt) result(value)
    implicit none
    character(*), intent(in) :: unique_id
    real(real64), optional, intent(in) :: dflt
    real(real64) :: value
    if (.not. xios_getvar(unique_id, value)) then
      if (present(dflt)) then
        value = dflt
      else
        call log_event('Failed to get value of variable ' // &
          trim(unique_id), log_level_error)
      end if
    end if
  end function get_variable

end module lfric_xios_diag_mod
