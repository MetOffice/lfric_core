!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module for the class to hold information about LFRic axes
!>
!> @details Contains information about all of the axes used in an LFRic
!> model run. Pointers to axisspec objects are included in the fieldspec for
!> each dimension for every diagnostic. These axisspecs are used to describe
!> both vertical and any non-spatial axes that a diagnostic may have


module axisspec_mod

  use constants_mod,             only: i_def, str_def, r_def, i_native, l_def
  use linked_list_data_mod,      only: linked_list_data_type
  use non_spatial_dimension_mod, only: NUMERICAL, CATEGORICAL

  implicit none

  private

  type, extends(linked_list_data_type), public :: axisspec_type

    private

    !> Unique id used by the diagnostic system to identify the field
    character(str_def)              :: unique_id
    !> Length of axis definition array
    integer(i_def)                  :: size
    !> Type of axis definition, either NUMERICAL or CATEGORICAL
    integer(i_native)               :: axis_category
    !> Whether axis defines vertical extrusion of lfric model
    logical(l_def)                  :: is_primary_axis
    !> Array containing numerical definition of axis, only for NUMERICAL axes
    real(r_def),        allocatable :: numeric_axis_def(:)
    !> Array containing strings describing categories on axis, only for
    !> CATEGORICAL axes
    character(str_def), allocatable :: label_axis_def(:)


  contains

    !> Getter to return the unique_id
    procedure, public :: get_unique_id

    !> Getter to return the size
    procedure, public :: get_size

    !> Getter to return the numeric_axis_def
    procedure, public :: get_numeric_axis_def

    !> Getter to return the label_axis_def
    procedure, public :: get_label_axis_def

    !> Getter to return the axis_category
    procedure, public :: get_axis_category

    !> Getter to return is_primary_axis
    procedure, public :: get_is_primary_axis

  end type axisspec_type

  interface axisspec_type
    module procedure axisspec_constructor
  end interface


contains


  !> Construct an <code>axisspec_type</code> object.
  !>
  !> Should contain only one axis definition, either:
  !>   numeric_axis_def for NUMERICAL axes (a real array of axis values)
  !>   label_axis_def for CATEGORICAL axes (a character array of labels)
  !>
  !> @param [in] unique_id A unique identifer for the field
  !> @param [in] axis_category Type of the axis definition. This should
  !>                           be either NUMERICAL or CATEGORICAL
  !> @param [in] numeric_axis_def Real array defining NUMERICAL axes
  !> @param [in] label_axis_def Character array defining CATEGORICAL axes
  !> @return self the axisspec object
  function axisspec_constructor(unique_id,          &
                                axis_category,      &
                                numeric_axis_def,   &
                                label_axis_def)     &
                                result(self)

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none

    character(*),                       intent(in)  :: unique_id
    integer(i_native),                  intent(in)  :: axis_category
    real(r_def),              optional, intent(in)  :: numeric_axis_def(:)
    character(str_def),       optional, intent(in)  :: label_axis_def(:)

    type(axisspec_type), target :: self

    self%unique_id = trim(unique_id)
    self%axis_category = axis_category
    if (index(unique_id, 'model_vert_axis') == 1) then
      self%is_primary_axis = .true.
    else
      self%is_primary_axis = .false.
    end if

    if (present(numeric_axis_def) .and. present(label_axis_def)) then
      call log_event("Axis cannot have both a 'numeric' and a 'label'" // &
              "axis_def", LOG_LEVEL_ERROR)
    end if

    if (axis_category == NUMERICAL) then
      if (present(label_axis_def)) then
        call log_event("Numerical axis cannot have a label_axis_def", &
                       LOG_LEVEL_ERROR)
      else
        self%numeric_axis_def = numeric_axis_def
        self%size = size(numeric_axis_def)
      end if

    else if (axis_category == CATEGORICAL) then
      if (present(numeric_axis_def)) then
        call log_event("Categorical axis cannot have a numeric_axis_def", &
                       LOG_LEVEL_ERROR)
      else
        self%label_axis_def = label_axis_def
        self%size = size(label_axis_def)
      end if

    else
      call log_event("Unrecognised axis category", LOG_LEVEL_ERROR)
    end if


  end function axisspec_constructor

  !> Getter for unique_id
  !> @param[in]  self  axisspec_type
  !> @return unique_id
  function get_unique_id(self) result(unique_id)

    implicit none

    class(axisspec_type), intent(in) :: self
    character(str_def)               :: unique_id

    unique_id = trim(self%unique_id)

  end function get_unique_id

  !> Getter for size
  !> @param[in]  self axisspec_type
  !> @return size
  function get_size(self) result(size)

    implicit none

    class(axisspec_type), intent(in) :: self
    integer(i_def) :: size

    size = self%size

  end function get_size

  !> Getter for numeric_axis_def, which is only set for NUMERICAL axes
  !> @param[in]  self  axisspec_type
  !> @return numeric_axis_def
  function get_numeric_axis_def(self) result(numeric_axis_def)

    implicit none

    class(axisspec_type), intent(in) :: self
    real(r_def),         allocatable :: numeric_axis_def(:)

    numeric_axis_def = self%numeric_axis_def

  end function get_numeric_axis_def

  !> Getter for label_axis_def, which is only set for CATEGORICAL axes
  !> @param[in]  self  axisspec_type
  !> @return label_axis_def
  function get_label_axis_def(self) result(label_axis_def)

    implicit none

    class(axisspec_type), intent(in) :: self
    character(str_def), allocatable  :: label_axis_def(:)

    label_axis_def = self%label_axis_def

  end function get_label_axis_def

  !> Getter for axis_category
  !> @param[in]  self  axisspec_type
  !> @return axis_category
  function get_axis_category(self) result(axis_category)

    implicit none

    class(axisspec_type), intent(in) :: self
    integer(i_native)                :: axis_category

    axis_category = self%axis_category

  end function get_axis_category

  !> Getter for is_primary_axis
  !> @param[in]  self  axisspec_type
  !> @return is_primary_axis
  function get_is_primary_axis(self) result(is_primary_axis)

    implicit none

    class(axisspec_type), intent(in) :: self
    logical(l_def)                   :: is_primary_axis

    is_primary_axis = self%is_primary_axis

  end function get_is_primary_axis

end module axisspec_mod
