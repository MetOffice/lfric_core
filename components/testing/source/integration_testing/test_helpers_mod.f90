!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise.

!> @brief Helper utilities for creating and manipulating test fields
!> @details
!>   Provides generic interfaces for creating test fields with constant values
!>   and printing field values to stdout in a standardized format. Supports
!>   both 32-bit and 64-bit real precision through overloaded procedures.
module test_helpers_mod

  use, intrinsic :: iso_fortran_env, only : real32, real64, output_unit

  use field_parent_mod,   only : field_parent_type
  use field_real32_mod,   only : field_real32_type
  use field_real64_mod,   only : field_real64_type
  use function_space_mod, only : function_space_type
  use lfric_mpi_mod,      only : global_mpi
  use test_field_mod,     only : test_field_32_type, test_field_64_type

  implicit none

  private
  public :: create_test_field_flood_constant, print_test_field

  !> Generic interface for creating constant-valued test fields
  !> Overloads for both 32-bit and 64-bit real precision
  interface create_test_field_flood_constant
    procedure create_test_field_32_flood_constant
    procedure create_test_field_64_flood_constant
  end interface create_test_field_flood_constant

  !> Generic interface for printing test field values
  !> Overloads for both 32-bit and 64-bit real precision
  interface print_test_field
    procedure print_test_field_32
    procedure print_test_field_64
  end interface print_test_field

contains

  !> @brief Create a 32-bit test field with constant value
  !> @param [in]     fs             Function space on which to create the field
  !> @param [in]     name           Name for the new field
  !> @param [out]    lfric_field    LFRic field object to be initialized
  !> @param [out]    test_field     Test field wrapper object
  !> @param [in]     value          Constant value to fill the field with
  !> @details
  !>   Creates a 32-bit real field, initializes it with a constant value,
  !>   and wraps it in a test_field object for easy data manipulation.
  !>   Allocates memory for both the LFRic field and test wrapper.
  subroutine create_test_field_32_flood_constant( fs, name,                &
                                                  lfric_field, test_field, &
                                                  value )

    implicit none

    class(function_space_type), intent(in), pointer      :: fs
    character(*),               intent(in)               :: name
    class(field_real32_type),   intent(out), &
                                     allocatable, target :: lfric_field
    class(test_field_32_type),  intent(out), allocatable :: test_field
    real(real32),               intent(in)               :: value

    class(field_real32_type), pointer :: lfric_field_ptr
    real(real32),             pointer :: test_data(:)

    ! Allocate and initialize the LFRic field on the given function space
    allocate(lfric_field)
    call lfric_field%initialise( fs, name=name )
    ! Create a pointer to the field for the test wrapper
    lfric_field_ptr => lfric_field
    ! Allocate the test field wrapper around the LFRic field
    allocate(test_field, source=test_field_32_type( lfric_field_ptr ))
    ! Set all field values to the constant value
    test_data => test_field%get_data()
    test_data = value
    ! Synchronize the test data back to the LFRic field
    call test_field%copy_to_lfric()

  end subroutine create_test_field_32_flood_constant


  !> @brief Create a 64-bit test field with constant value
  !> @param [in]     fs             Function space on which to create the field
  !> @param [in]     name           Name for the new field
  !> @param [out]    lfric_field    LFRic field object to be initialized
  !> @param [out]    test_field     Test field wrapper object
  !> @param [in]     value          Constant value to fill the field with
  !> @details
  !>   Creates a 64-bit real field, initializes it with a constant value,
  !>   and wraps it in a test_field object for easy data manipulation.
  !>   Allocates memory for both the LFRic field and test wrapper.
  subroutine create_test_field_64_flood_constant( fs, name,                &
                                                  lfric_field, test_field, &
                                                  value )

    implicit none

    type(function_space_type), intent(in), pointer      :: fs
    character(*),              intent(in)               :: name
    type(field_real64_type),   intent(out), &
                                    allocatable, target :: lfric_field
    type(test_field_64_type),  intent(out), allocatable :: test_field
    real(real64),              intent(in)               :: value

    type(field_real64_type), pointer :: lfric_field_ptr
    real(real64),            pointer :: test_data(:)

    ! Allocate and initialize the LFRic field on the given function space
    allocate(lfric_field)
    call lfric_field%initialise( fs, name=name )
    ! Create a pointer to the field for the test wrapper
    lfric_field_ptr => lfric_field
    ! Allocate the test field wrapper around the LFRic field
    allocate(test_field, source=test_field_64_type( lfric_field_ptr ))
    ! Set all field values to the constant value
    test_data => test_field%get_data()
    test_data = value
    ! Synchronize the test data back to the LFRic field
    call test_field%copy_to_lfric()

  end subroutine create_test_field_64_flood_constant


  !> @brief Print 32-bit test field values to stdout
  !> @param [in,out] test_field  The test field to print
  !> @details
  !>   Copies field data from LFRic field storage to test field storage,
  !>   then outputs the field values formatted as:
  !>   rank_number (undf_count) = value1 value2 ... valueN
  subroutine print_test_field_32( test_field )

    implicit none

    type(test_field_32_type), intent(inout) :: test_field

    class(field_parent_type),  pointer :: lfric_field
    type(function_space_type), pointer :: fs
    character(255)                     :: format_str
    real(real32),              pointer :: test_data(:)

    ! Get references to the field and its function space
    lfric_field => test_field%get_lfric_field_ptr()
    fs => lfric_field%get_function_space()
    ! Construct format string: "rank (ndf) = value value value ..."
    ! Number of values to print: fs%get_undf() (total degrees of freedom)
    write(format_str, '("A, "", ("", I0, "") = "", ", I0, "F7.3")') fs%get_undf()
    ! Synchronize field data from LFRic field to test field
    call test_field%copy_from_lfric()
    ! Get pointer to test field data and print with appropriate format
    test_data => test_field%get_data()
    write(output_unit, format_str) global_mpi%get_comm_rank(), test_data(:fs%get_undf())

  end subroutine print_test_field_32


  !> @brief Print 64-bit test field values to stdout
  !> @param [in,out] test_field  The test field to print
  !> @details
  !>   Copies field data from LFRic field storage to test field storage,
  !>   then outputs the field values formatted as:
  !>   field_name (rank_number) = value1 value2 ... valueN
  subroutine print_test_field_64( test_field )

    implicit none

    type(test_field_64_type), intent(inout) :: test_field

    class(field_parent_type),  pointer :: lfric_field
    type(function_space_type), pointer :: fs
    character(255)                     :: format_str
    real(real64),              pointer :: test_data(:)

    ! Get references to the field and its function space
    lfric_field => test_field%get_lfric_field_ptr()
    fs => lfric_field%get_function_space()
    ! Construct format string: "field_name (rank) = value value value ..."
    ! Number of values to print: fs%get_last_dof_owned() (owned degrees of freedom)
    write(format_str, '("(A, "" ("", I0, "") = "", ", I0, "F7.3)")') &
      fs%get_last_dof_owned()
    ! Synchronize field data from LFRic field to test field
    call test_field%copy_from_lfric()
    ! Get pointer to test field data and print field name, rank, and values
    test_data => test_field%get_data()
    write(output_unit, format_str) trim(lfric_field%get_name()), &
                                   global_mpi%get_comm_rank(),   &
                                   test_data(:fs%get_last_dof_owned())

  end subroutine print_test_field_64

end module test_helpers_mod
