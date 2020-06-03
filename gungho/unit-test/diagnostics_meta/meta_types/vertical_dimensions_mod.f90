!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module for creating vertical dimension types
!>
!> @details This module contains functions that create vertical_dimension
!> objects. These functions are used in field meta data definition files to
!> create appropriate vertical dimension objects.


module vertical_dimensions_mod

  use constants_mod,                only : i_def, r_def
  use vertical_dimension_types_mod, only : range_vertical_dimension_meta_data_type, &
                                           list_vertical_dimension_meta_data_type
  use levels_enum_mod,              only : BOTTOM_ATMOSPHERIC_LEVEL, &
                                           TOP_ATMOSPHERIC_LEVEL, &
                                           BOTTOM_SOIL_LEVEL, &
                                           TOP_SOIL_LEVEL
  use positive_enum_mod,            only : POSITIVE_UP, &
                                           POSITIVE_DOWN

  implicit none

  private
  public ::   model_custom_dimension, &
              model_height_dimension, &
              model_depth_dimension,  &
              fixed_custom_dimension, &
              fixed_depth_dimension,  &
              fixed_height_dimension

  contains
  !> @brief Returns a range_vertical_dimension_meta_data_type
  !> @details This is for is the developer wants to specifiy everthing, when
  !> using model levels, instead of using default values
  !> @param [in] name The standard name of the vertical dimension
  !> @param [in] units The unit of measure
  !> @param [in] positive The direction the a positive moves, up or down
  !> @param [in] bottom bottom The bottom level as per levels_enum
  !> @param [in] top The top level as per levels_enum
  function model_custom_dimension(name, units, positive, bottom, top) result (vertical_type)

    implicit none

    character(*),   intent(in)  :: name, units
    integer(i_def), intent(in)  :: positive, bottom, top

    type(range_vertical_dimension_meta_data_type) :: vertical_type

    vertical_type = range_vertical_dimension_meta_data_type(&
                    standard_name = name,                     &
                    units = units,                            &
                    positive = positive,                      &
                    lower_level = bottom,                     &
                    upper_level = top)

  end function model_custom_dimension

  !> @brief Returns a range_vertical_dimension_meta_data_type
  !> @details For creating a depth vertical dimension on model levels using
  !> default settings. Developer can optionally supply highest and lowest model
  !> level. Since the majority of fields will be over all levels, a user can
  !> use the syntax "model_depth_dimension() to obtain a dimension with
  !> default settings. But can provide arguments if a different levels are
  !> required.
  !> @param [in,optional] bottom The bottom level as per levels_enum
  !> @param [in,optional] top The top level as per levels_enum
  function model_depth_dimension(bottom, top) result (vertical_type)

    implicit none

    integer(i_def), intent(in), optional  :: bottom, top

    integer(i_def)                        :: upper_level, lower_level

    type(range_vertical_dimension_meta_data_type) :: vertical_type

    !> Set defaults if arguments not supplied
    if(present(bottom)) then
      lower_level = bottom
    else
      lower_level = BOTTOM_SOIL_LEVEL
    endif

    if(present(top)) then
      upper_level = top
    else
      upper_level = TOP_SOIL_LEVEL
    endif

    vertical_type = range_vertical_dimension_meta_data_type(&
                    standard_name = 'depth',                  &
                    units = 'm',                              &
                    positive = POSITIVE_DOWN,                 &
                    lower_level = lower_level,                &
                    upper_level = upper_level)

  end function model_depth_dimension


  !> @brief Returns a range_vertical_dimension_meta_data_type
  !> @details For creating a depth vertical dimension on model levels using
  !> default settings. Developer can optionally supply highest and lowest model
  !> level. Since the majority of fields will be over all levels, a user can
  !> use the syntax "model_height_dimension() to obtain a dimension with
  !> default settings. But can provide arguments if a different levels are
  !> required.
  !> @param [in,optional] bottom The bottom level as per levels_enum
  !> @param [in,optional] top The top level as per levels_enum
  function model_height_dimension(bottom, top) result (vertical_type)

    implicit none

    integer(i_def), intent(in), optional :: bottom, top

    integer(i_def)                       :: upper_level, lower_level

    type(range_vertical_dimension_meta_data_type)  :: vertical_type

    !> Set defaults if arguments not supplied
    if(present(bottom)) then
      lower_level = bottom
    else
      lower_level = BOTTOM_ATMOSPHERIC_LEVEL
    endif

    if(present(top)) then
      upper_level = top
    else
      upper_level = TOP_ATMOSPHERIC_LEVEL
    endif

    vertical_type = range_vertical_dimension_meta_data_type(&
                    standard_name = 'height',                 &
                    positive = POSITIVE_UP,                   &
                    units = 'm',                              &
                    lower_level = lower_level,                &
                    upper_level = upper_level)

  end function model_height_dimension


  !> @brief Returns a list_vertical_dimension_meta_data_type
  !> @details This is for is the developer wants to specifiy everthing, when
  !> using fixed levels, instead of using default values
  !> @param [in] name The standard name of the vertical dimension
  !> @param [in] units The unit of measure
  !> @param [in] positive The direction the a positive moves, up or down
  !> @param [in] level_definition An array of real numbers that define the
  !> height, or depth of levels
  function fixed_custom_dimension(name, positive, units, level_definition) result (vertical_type)

    implicit none

    character(*),                   intent(in) :: name, units
    integer(i_def),                 intent(in) :: positive
    real(r_def),                    intent(in) :: level_definition(:)

    type(list_vertical_dimension_meta_data_type) :: vertical_type

    vertical_type = list_vertical_dimension_meta_data_type(&
                    standard_name = name,                       &
                    positive = positive,                        &
                    units = units,                              &
                    level_definition = level_definition)

  end function fixed_custom_dimension


  !> @brief Returns a list_vertical_dimension_meta_data_type
  !> @details For creating a depth vertical dimension on fixed levels using
  !> default settings.
  !> @param [in,] level_definition An array of real numbers that define the
  !> height, or depth of levels
  function fixed_depth_dimension(level_definition) result(vertical_type)

    implicit none

    real(r_def), intent(in) :: level_definition(:)

    type(list_vertical_dimension_meta_data_type) :: vertical_type

    vertical_type = list_vertical_dimension_meta_data_type(&
                    standard_name = 'depth',                    &
                    positive = POSITIVE_DOWN,                   &
                    units = 'm',                                &
                    level_definition = level_definition)

  end function fixed_depth_dimension


  !> @brief Returns a list_vertical_dimension_meta_data_type
  !> @details For creating a height vertical dimension on fixed levels using
  !> default settings.
  !> @param [in,] level_definition An array of real numbers that define the
  !> height, or depth of levels
  function fixed_height_dimension(level_definition) result (vertical_type)

    implicit none

    real(r_def), intent(in) :: level_definition(:)

    type(list_vertical_dimension_meta_data_type) :: vertical_type

    vertical_type = list_vertical_dimension_meta_data_type(&
                    standard_name = 'height',                   &
                    positive = POSITIVE_UP,                     &
                    units = 'm',                                &
                    level_definition = level_definition)

  end function fixed_height_dimension

end module vertical_dimensions_mod