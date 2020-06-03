!!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module for the vertical dimension types
!>
!> @details This allows for the creation of an object that represents the
!> vertical dimension meta data for a field object.
!> A meta_data_type can be supplied with, as an optional argument, a
!> vertical_dimension_type to define its vertical dimension meta data.


module vertical_dimension_types_mod

  use constants_mod,        only: i_def, str_def, r_def, l_def

  implicit none

  private

  type, abstract, public :: abstract_vertical_meta_data_type

    !> standard_name used by the dimension
    character(str_def) :: standard_name
    !> unit of measure used by the dimension
    character(str_def) :: units
    !> The direction of movement for positive numbers
    integer(i_def)     :: positive
    !> The axis that this dimension represents
    character(str_def) :: axis

  end type abstract_vertical_meta_data_type


  !> Defines the range vertical dimension
  type, public, extends(abstract_vertical_meta_data_type) :: range_vertical_dimension_meta_data_type

    !> The lowest level the field is calculated on
    integer(i_def) :: lower_level

    !> The highest level the field is calculated on
    integer(i_def) :: upper_level

  end type range_vertical_dimension_meta_data_type

  !> Defines the list vertical dimension
  type, public, extends(abstract_vertical_meta_data_type) :: list_vertical_dimension_meta_data_type

    !> Array that defines the vertical level
    real(r_def), allocatable     :: level_definition(:)

  end type list_vertical_dimension_meta_data_type

  interface range_vertical_dimension_meta_data_type
    module procedure range_vertical_dimension_meta_data_constructor
  end interface range_vertical_dimension_meta_data_type

  interface list_vertical_dimension_meta_data_type
    module procedure list_vertical_dimension_meta_data_constructor
   end interface list_vertical_dimension_meta_data_type


  contains

  !> Construct a <code>range_vertical_dimension_meta_data_type</code> object.
  !> @param [in] standard_name The standard name of the verticality
  !> @param [in] units SI Unit of measure for the field
  !> @param [in] positive The direction for positive numbers in the field has a height
  !> @param [in] lower level The lower level that the field will be calculated on
  !> @param [in] upper_level The upper level that the field will be calculated on
  !> @return self the meta_data object

  function range_vertical_dimension_meta_data_constructor(&
                                 standard_name,             &
                                 units,                     &
                                 positive,                  &
                                 upper_level,               &
                                 lower_level)               &
                                 result(self)
    implicit none

    character(*),                 intent(in)  :: standard_name
    character(*),                 intent(in)  :: units
    integer(i_def),               intent(in)  :: positive
    integer(i_def),               intent(in)  :: lower_level
    integer(i_def),               intent(in)  :: upper_level

    type(range_vertical_dimension_meta_data_type)   :: self

    self%standard_name = standard_name
    self%units = units
    self%positive = positive
    self%lower_level = lower_level
    self%upper_level = upper_level
    self%axis = "Z"

  end function range_vertical_dimension_meta_data_constructor

  !> Construct a <code>list_vertical_dimension_meta_data_type</code> object.
  !> @param [in] standard_name The standard name of the verticality
  !> @param [in] units SI Unit of measure for the field
  !> @param [in] positive The direction for positive numbers in the field has a height
  !> @param [in] level_definition The definition of vertical axis. Defined with an array of reals
  !> @return self the meta_data object

  function list_vertical_dimension_meta_data_constructor(&
                                 standard_name,               &
                                 units,                       &
                                 positive,                    &
                                 level_definition)            &
                                 result(self)
    implicit none

    character(*),                         intent(in)  :: standard_name
    character(*),                         intent(in)  :: units
    integer(i_def),                       intent(in)  :: positive
    real(r_def),                          intent(in)  :: level_definition(:)

    type(list_vertical_dimension_meta_data_type)   :: self

    self%standard_name = standard_name
    self%units = units
    self%positive = positive
    self%axis = "Z"

    allocate(self%level_definition(size(level_definition)))
    self%level_definition = level_definition

  end function list_vertical_dimension_meta_data_constructor

end module vertical_dimension_types_mod