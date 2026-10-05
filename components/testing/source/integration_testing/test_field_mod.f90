!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise.

!> @brief Test field wrapper for easy manipulation of field data in tests
!> @details
!>   Provides wrapper types around LFRic fields that allow direct access to
!>   field data arrays for test setup and verification. Supports both 32-bit
!>   and 64-bit real precision. The wrapper maintains a mirror of the field
!>   data that can be synchronized to/from the underlying LFRic field.
module test_field_mod

  use, intrinsic :: iso_fortran_env, only : real32, real64

  use abstract_external_field_mod, only : abstract_external_field_type
  use constants_mod,               only : i_def
  use field_real32_mod,            only : field_real32_type, &
                                          field_real32_proxy_type
  use field_real64_mod,            only : field_real64_type, &
                                          field_real64_proxy_type
  use log_mod,                     only : log_event, log_level_error
  use field_parent_mod,            only : field_parent_type

  implicit none

  private

  !> @brief 32-bit real test field wrapper
  !> @details
  !>   Wraps a 32-bit LFRic field, mirroring its data for test manipulation.
  !>   Allows direct read/write access to field values without breaking
  !>   encapsulation of the underlying LFRic field object.
  type, public, extends(abstract_external_field_type) :: test_field_32_type
    private
    !> Mirror copy of the field data for test access
    real(real32), pointer :: test_data(:) => null()
  contains
    private
    procedure, public :: copy_from_lfric => copy_32_from_lfric
    procedure, public :: copy_to_lfric => copy_32_to_lfric
    procedure, public :: get_data => get_data_32
    final :: destroy_32
  end type test_field_32_type

  interface test_field_32_type
    procedure test_field_32_constructor
  end interface test_field_32_type


  !> @brief 64-bit real test field wrapper
  !> @details
  !>   Wraps a 64-bit LFRic field, mirroring its data for test manipulation.
  !>   Allows direct read/write access to field values without breaking
  !>   encapsulation of the underlying LFRic field object.
  type, public, extends(abstract_external_field_type) :: test_field_64_type
    private
    !> Mirror copy of the field data for test access
    real(real64), pointer :: test_data(:) => null()
  contains
    private
    procedure, public :: copy_from_lfric => copy_64_from_lfric
    procedure, public :: copy_to_lfric => copy_64_to_lfric
    procedure, public :: get_data => get_data_64
    final :: destroy_64
  end type test_field_64_type

  interface test_field_64_type
    procedure test_field_64_constructor
  end interface test_field_64_type

contains

  !> @brief Constructor for 32-bit test field
  !> @param [in] lfric_field  Pointer to the 32-bit LFRic field to wrap
  !> @return                  New test_field_32_type instance
  !> @details
  !>   Creates a test field wrapper around an LFRic field, mirroring its data
  !>   array for manipulation within the test framework.
  function test_field_32_constructor( lfric_field ) result(new_instance)

    implicit none

    class(field_real32_type), intent(in), pointer :: lfric_field
    type(test_field_32_type) :: new_instance

    type(field_real32_proxy_type)     :: proxy
    class(field_parent_type), pointer :: cast_field

    ! Fortran type system dance: cast child class to parent class
    ! This allows passing the specific field type to an argument expecting
    ! a parent class, which some compilers require for proper type checking
    cast_field => lfric_field
    call new_instance%abstract_external_field_initialiser( cast_field )

    ! Mirror the data array from the LFRic field for test manipulation
    proxy = lfric_field%get_proxy()
    allocate( new_instance%test_data(size(proxy%data)) )

  end function test_field_32_constructor


  !> @brief Constructor for 64-bit test field
  !> @param [in] lfric_field  Pointer to the 64-bit LFRic field to wrap
  !> @return                  New test_field_64_type instance
  !> @details
  !>   Creates a test field wrapper around an LFRic field, mirroring its data
  !>   array for manipulation within the test framework.
  function test_field_64_constructor( lfric_field ) result(new_instance)

    implicit none

    class(field_real64_type), intent(in), pointer :: lfric_field
    type(test_field_64_type) :: new_instance

    type(field_real64_proxy_type)     :: proxy
    class(field_parent_type), pointer :: cast_field

    ! Fortran type system dance: cast child class to parent class
    ! This allows passing the specific field type to an argument expecting
    ! a parent class, which some compilers require for proper type checking
    cast_field => lfric_field
    call new_instance%abstract_external_field_initialiser( cast_field )

    ! Mirror the data array from the LFRic field for test manipulation
    proxy = lfric_field%get_proxy()
    allocate( new_instance%test_data(size(proxy%data)) )

  end function test_field_64_constructor


  !> @brief Destructor for 32-bit test field
  !> @param [in,out] this  The test field to destroy
  !> @details
  !>   Deallocates the internal data array and cleans up resources.
  !>   Called automatically when the test field object goes out of scope.
  subroutine destroy_32( this )

    implicit none

    type(test_field_32_type), intent(inout) :: this

    ! Only deallocate if the pointer was actually allocated
    ! Prevents error if destructor called on uninitialized object
    if (associated(this%test_data)) then
      deallocate( this%test_data )
    end if

  end subroutine destroy_32


  !> @brief Destructor for 64-bit test field
  !> @param [in,out] this  The test field to destroy
  !> @details
  !>   Deallocates the internal data array and cleans up resources.
  !>   Called automatically when the test field object goes out of scope.
  subroutine destroy_64( this )

    implicit none

    type(test_field_64_type), intent(inout) :: this

    ! Only deallocate if the pointer was actually allocated
    ! Prevents error if destructor called on uninitialized object
    if (associated(this%test_data)) then
      deallocate( this%test_data )
    end if

  end subroutine destroy_64


  !> @brief Get pointer to 32-bit test field data array
  !> @param [in] this  The test field object
  !> @return           Pointer to the internal data array
  !> @details
  !>   Returns a pointer to the underlying 32-bit real data array.
  !>   Useful for direct data manipulation in tests.
  function get_data_32( this )

    implicit none

    class(test_field_32_type), intent(in) :: this
    real(real32), pointer :: get_data_32(:)

    get_data_32 => this%test_data

  end function get_data_32


  !> @brief Get pointer to 64-bit test field data array
  !> @param [in] this  The test field object
  !> @return           Pointer to the internal data array
  !> @details
  !>   Returns a pointer to the underlying 64-bit real data array.
  !>   Useful for direct data manipulation in tests.
  function get_data_64( this )

    implicit none

    class(test_field_64_type), intent(in) :: this
    real(real64), pointer :: get_data_64(:)

    get_data_64 => this%test_data

  end function get_data_64


  !> @brief Copy 32-bit field data from LFRic field to test field
  !> @param [in,out] self        The test field object
  !> @param [out]    return_code Optional status code (always 0 on success)
  !> @details
  !>   Synchronizes the test field data array with the underlying LFRic field,
  !>   retrieving the current field values for inspection in tests.
  subroutine copy_32_from_lfric(self, return_code)

    implicit none

    class(test_field_32_type), intent(inout) :: self
    integer(i_def),  optional, intent(out)   :: return_code

    type(field_real32_type), pointer :: lfric_field
    type(field_real32_proxy_type)    :: lfric_proxy

    ! Type guard to ensure the wrapped field is actually a 32-bit field
    select type(lfric_field => self%get_lfric_field_ptr())
    class is (field_real32_type)
      ! Get proxy to access the raw data array from the LFRic field
      lfric_proxy = lfric_field%get_proxy()
      ! Copy data from LFRic field to test field mirror
      self%test_data = lfric_proxy%data
    class default
      ! Log error if field type doesn't match test field type
      call log_event("Inconsistent test_field_32_type", log_level_error )
    end select

    if (present(return_code)) then
      return_code = 0
    end if

  end subroutine copy_32_from_lfric


  !> @brief Copy 64-bit field data from LFRic field to test field
  !> @param [in,out] self        The test field object
  !> @param [out]    return_code Optional status code (always 0 on success)
  !> @details
  !>   Synchronizes the test field data array with the underlying LFRic field,
  !>   retrieving the current field values for inspection in tests.
  subroutine copy_64_from_lfric(self, return_code)

    implicit none

    class(test_field_64_type), intent(inout) :: self
    integer(i_def),  optional, intent(out)   :: return_code

    type(field_real64_type), pointer :: lfric_field
    type(field_real64_proxy_type)    :: lfric_proxy

    ! Type guard to ensure the wrapped field is actually a 64-bit field
    select type(lfric_field => self%get_lfric_field_ptr())
    class is (field_real64_type)
      ! Get proxy to access the raw data array from the LFRic field
      lfric_proxy = lfric_field%get_proxy()
      ! Copy data from LFRic field to test field mirror
      self%test_data = lfric_proxy%data
    class default
      ! Log error if field type doesn't match test field type
      call log_event("Inconsistent test_field_64_type", log_level_error )
    end select

    if (present(return_code)) then
      return_code = 0
    end if

  end subroutine copy_64_from_lfric


  !> @brief Copy 32-bit field data from test field to LFRic field
  !> @param [in,out] self        The test field object
  !> @param [out]    return_code Optional status code (always 0 on success)
  !> @details
  !>   Propagates modified test field data back to the underlying LFRic field,
  !>   allowing tests to set field values that can be used by other code.
  subroutine copy_32_to_lfric(self, return_code)

    implicit none

    class(test_field_32_type), intent(inout) :: self
    integer(i_def),  optional, intent(out)   :: return_code

    type(field_real32_type), pointer :: lfric_field
    type(field_real32_proxy_type)    :: lfric_proxy

    ! Type guard to ensure the wrapped field is actually a 32-bit field
    select type(lfric_field => self%get_lfric_field_ptr())
    class is (field_real32_type)
      ! Get proxy to access the raw data array from the LFRic field
      lfric_proxy = lfric_field%get_proxy()
      ! Copy modified data from test field back to LFRic field
      lfric_proxy%data = self%test_data
    class default
      ! Log error if field type doesn't match test field type
      call log_event("Inconsistent test_field_32_type", log_level_error )
    end select

    if (present(return_code)) then
      return_code = 0
    end if

  end subroutine copy_32_to_lfric


  !> @brief Copy 64-bit field data from test field to LFRic field
  !> @param [in,out] self        The test field object
  !> @param [out]    return_code Optional status code (always 0 on success)
  !> @details
  !>   Propagates modified test field data back to the underlying LFRic field,
  !>   allowing tests to set field values that can be used by other code.
  subroutine copy_64_to_lfric(self, return_code)

    implicit none

    class(test_field_64_type), intent(inout) :: self
    integer(i_def),  optional, intent(out)   :: return_code

    type(field_real64_type), pointer :: lfric_field
    type(field_real64_proxy_type)    :: lfric_proxy

    ! Type guard to ensure the wrapped field is actually a 64-bit field
    select type(lfric_field => self%get_lfric_field_ptr())
    class is (field_real64_type)
      ! Get proxy to access the raw data array from the LFRic field
      lfric_proxy = lfric_field%get_proxy()
      ! Copy modified data from test field back to LFRic field
      lfric_proxy%data = self%test_data
    class default
      ! Log error if field type doesn't match test field type
      call log_event("Inconsistent test_field_32_type", log_level_error )
    end select

    if (present(return_code)) then
      return_code = 0
    end if

  end subroutine copy_64_to_lfric

end module test_field_mod
