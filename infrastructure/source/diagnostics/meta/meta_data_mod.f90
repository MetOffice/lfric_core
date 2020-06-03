!!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!
!> @brief A module for the meta_data type
!>
!> @details This file defines the meta_data_type. This type will hold all meta
!> data information for a given field. Objects of this type are created by the
!> field meta data definition files

module meta_data_mod

 use constants_mod,                only: i_def, str_def, l_def, str_long, &
                                         i_native
 use vertical_dimension_types_mod, only: range_vertical_dimension_meta_data_type, &
                                         list_vertical_dimension_meta_data_type, &
                                         abstract_vertical_meta_data_type
 use misc_meta_data_mod,           only: misc_meta_data_type

 implicit none

 private

 type, public :: meta_data_type

   !> Unique id used by the diagnostic system to identify the field
   character(str_def)  ::  unique_id
   !> Long name of the field
   character(str_def)  ::  long_name
   !> the SI unit of measure for the field
   character(str_def)  ::  units

   !> Function space.
   !> For prognostics, which function space does this need to exist on.
   !> For diagnostics, what function space(s) is it possible to request this on
   integer(i_native)  ::  function_space
   !> The IO driver used by this field
   character(str_def) ::  io_driver
   !> The order of the function space. Cannot be negative
   integer(i_def)     ::  order

   !> Contains Rose triggering syntax.
   character(str_def)  ::  trigger
   !> This string is displayed in the Rose gui under the field
   character(str_long) ::  description

   !> ENUM's are used to represent the data type
   integer(i_def)      ::  data_type
   !> Timestep information.
   !> Specifies which timesteps diagnostics or prognostics are available
   !> This is set to an enumerated value, not the actual number of time steps
   integer(i_def)      ::  time_step
   !> Interpolation method. Yet to be implemented
   integer(i_def)      ::  recommended_interpolation
   !> Level of packing for the field
   integer(i_def)      ::  packing

   !> Vertical attribute
   class(abstract_vertical_meta_data_type), allocatable :: vertical_dimension

   !> Other optional information used to create a field
   !> Standard name for the field, if it exists
   character(str_def)  ::  standard_name

   !> Positive direction of the field (vectors)
   character(str_def)  ::  positive

   !> Any misc data for the field
   !class(misc_meta_data_type), allocatable  :: misc_meta_data
   type(misc_meta_data_type), dimension(:), allocatable  :: misc_meta_data

contains

  !> Getters for attributes of the field
  procedure, public :: get_unique_id

  end type meta_data_type

  interface meta_data_type
    module procedure meta_data_constructor
  end interface

contains

  !> Construct a <code>meta_data_type</code> object.
  !> @brief The constructor for a meta_data_type object
  !> @param [in] unique_id A unique identifier of the field
  !> @param [in] units SI Unit of measure for the field
  !> @param [in] function_space The function_space to create field with
  !> @param [in] order The order of the function space
  !> @param [in] io_driver The IO driver used by the field
  !> @param [in] trigger The triggering syntax for the field in Rose
  !> @param [in] description Description of the field shown in Rose
  !> @param [in] data_type The data type of the field
  !> @param [in] time_step The available time step of the field
  !> @param [in] recommended_interpolation The recommended interpolation method
  !> @param [in] packing Packing setting for the field
  !> @param [in,optional] standard_name The standard name of the field if it exists
  !> @param [in,optional] long_name The long name of the field
  !> @param [in,optional] positive The direction for positive numbers
  !> @param [in,optional] vertical_dimension The vertical of the field
  !> @param [in,optional] misc_meta_data Holds a key/value pair of strings
  !> used for any miscellaneous data that a field might need
  !> @return self the meta_data object
  !>
  function meta_data_constructor(unique_id,                     &
                                 units,                         &
                                 function_space,                &
                                 order,                         &
                                 io_driver,                     &
                                 trigger,                       &
                                 description,                   &
                                 data_type,                     &
                                 time_step,                     &
                                 recommended_interpolation,     &
                                 packing,                       &
                                 standard_name,                 &
                                 long_name,                     &
                                 positive,                      &
                                 vertical_dimension,            &
                                 misc_meta_data)                &
                                 result(self)

    implicit none

    character(*),                       intent(in) :: unique_id
    character(*),                       intent(in) :: units
    integer(i_native),                  intent(in) :: function_space
    integer(i_def),                     intent(in) :: order
    character(*),                       intent(in) :: io_driver
    character(*),                       intent(in) :: trigger
    character(*),                       intent(in) :: description
    integer(i_def),                     intent(in) :: data_type
    integer(i_def),                     intent(in) :: time_step
    integer(i_def),                     intent(in) :: recommended_interpolation
    integer(i_def),                     intent(in) :: packing
    character(*), optional,             intent(in) :: standard_name
    character(*), optional,             intent(in) :: long_name
    character(*), optional,             intent(in) :: positive
    class(abstract_vertical_meta_data_type), optional, &
                                        intent(in) :: vertical_dimension
    type(misc_meta_data_type), dimension(:), optional, &
                                        intent(in) :: misc_meta_data

    type(meta_data_type)   :: self

    self%unique_id                 = trim(unique_id)
    self%units                     = units
    self%function_space            = function_space
    self%order                     = order
    self%io_driver                 = io_driver
    self%trigger                   = trigger
    self%description               = description
    self%data_type                 = data_type
    self%time_step                 = time_step
    self%recommended_interpolation = recommended_interpolation

    !> Handling optional arguments
    if(present(standard_name)) then
      self%standard_name = standard_name
    else
      self%standard_name = ""
    end if

    if(present(long_name)) then
      self%long_name = long_name
    else
      self%long_name = ""
    end if

    if(present(positive)) then
      self%positive = positive
    end if

    if(present(vertical_dimension)) then
      !> This is the way polymorphism is done in fortran 2003. You have to
      !> allocate the derived type instead of using the assignment operator (=)
      !> Fortran 2008 would support self%vertical_dimension = vertical_dimension
      allocate(self%vertical_dimension, source=vertical_dimension)
    end if

    if(present(misc_meta_data)) then
      allocate(self%misc_meta_data, source=misc_meta_data)
    end if

  end function meta_data_constructor

  !> Getter for unique_id
  !> @param[in]  self  meta_data_type
  !> @return unique_id

  function get_unique_id(self) result(unique_id)

    implicit none

    class(meta_data_type), intent(in) :: self
    character(str_def) :: unique_id

    unique_id = trim(self%unique_id)

  end function get_unique_id

end module meta_data_mod