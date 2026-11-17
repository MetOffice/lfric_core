!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!
!> @brief   Defines key-value pair objects that can be used w/ the IO system.
!> @details The abstract key-value pair object is extended w/ IO procedures.
!>          The following concrete types are defined:
!>
!> * io_value_type (abstract)
!>     - int32_io_value_type,       int64_io_value_type
!>     - real32_io_value_type,      real64_io_value_type
!>     - int32_arr_io_value_type,   int64_arr_io_value_type
!>     - real32_arr_io_value_type,  real64_arr_io_value_type
!----------------------------------------------------------------------------
module io_value_mod

  use, intrinsic :: iso_fortran_env,     only: int32, int64, real32, real64
  use, intrinsic :: ieee_arithmetic,     only: ieee_is_finite

  use constants_mod,                     only : str_def, r_def, r_double, l_def
  use key_value_mod,                     only : key_value_type
  use log_mod,                           only : log_event, log_scratch_space, &
                                                LOG_LEVEL_INFO, LOG_LEVEL_ERROR, &
                                                LOG_LEVEL_WARNING

  implicit none

  private

  public :: io_value_type
  public :: int32_io_value_type,      int64_io_value_type
  public :: int32_arr_io_value_type,  int64_arr_io_value_type
  public :: real32_io_value_type,     real64_io_value_type
  public :: real32_arr_io_value_type, real64_arr_io_value_type
  public :: io_read_interface, io_write_interface
  public :: to_real32, to_real64, to_int32, to_int64

  !> Numeric copies for direct access to concrete public value components.
  !> Scalars and rank-one arrays are supported within each numeric family.
  !> Array results preserve size; allocatable assignment can resize the target.
  !> Narrowing outside the finite destination range warns once per call, then
  !> converts normally: overflow results are compiler-dependent and may trap.
  !> Ordinary rounding is allowed. Existing NaN/infinity values are propagated
  !> without an overflow warning. Initialise logging before using these helpers.
  !> @{
  interface to_real32
    module procedure real32_from_real32, real32_from_real64
    module procedure real32_arr_from_real32, real32_arr_from_real64
  end interface to_real32

  interface to_real64
    module procedure real64_from_real32, real64_from_real64
    module procedure real64_arr_from_real32, real64_arr_from_real64
  end interface to_real64

  interface to_int32
    module procedure int32_from_int32, int32_from_int64
    module procedure int32_arr_from_int32, int32_arr_from_int64
  end interface to_int32

  interface to_int64
    module procedure int64_from_int32, int64_from_int64
    module procedure int64_arr_from_int32, int64_arr_from_int64
  end interface to_int64
  !> @}

  !> @brief Value with associated I/O methods
  !>        that can be stored in a key-value pair
  type, abstract, extends(key_value_type) :: io_value_type

    procedure(io_write_interface), pointer :: write_method => null()
    procedure(io_read_interface),  pointer :: checkpoint_read_method => null()
    procedure(io_write_interface), pointer :: checkpoint_write_method => null()

    contains

      procedure, public  :: set_write_behaviour
      procedure, public  :: set_checkpoint_write_behaviour
      procedure, public  :: set_checkpoint_read_behaviour
      procedure, public  :: can_write_checkpoint
      procedure, public  :: write_value
      procedure, public  :: write_checkpoint
      procedure, public  :: read_checkpoint
      procedure, private :: get_int32_value
      procedure, private :: get_int64_value
      procedure, private :: get_int32_arr_value
      procedure, private :: get_int64_arr_value
      procedure, private :: get_real32_value
      procedure, private :: get_real64_value
      procedure, private :: get_real32_arr_value
      procedure, private :: get_real64_arr_value
      generic            :: get_value => get_int32_value,       &
                                         get_int64_value,       &
                                         get_int32_arr_value,   &
                                         get_int64_arr_value,   &
                                         get_real32_value,      &
                                         get_real64_value,      &
                                         get_real32_arr_value,  &
                                         get_real64_arr_value



  end type io_value_type

!===============================================
! Concrete Types of Abstract io_value_type
!===============================================

  type, extends(io_value_type) :: int32_io_value_type
    integer(int32) :: value
  end type int32_io_value_type

  type, extends(io_value_type) :: int64_io_value_type
    integer(int64) :: value
  end type int64_io_value_type

  type, extends(io_value_type) :: int32_arr_io_value_type
    integer(int32), allocatable :: value(:)
  end type int32_arr_io_value_type

  type, extends(io_value_type) :: int64_arr_io_value_type
    integer(int64), allocatable :: value(:)
  end type int64_arr_io_value_type

  type, extends(io_value_type) :: real32_io_value_type
    real(real32) :: value
  end type real32_io_value_type

  type, extends(io_value_type) :: real64_io_value_type
    real(real64) :: value
  end type real64_io_value_type

  type, extends(io_value_type) :: real32_arr_io_value_type
    real(real32), allocatable :: value(:)
  end type real32_arr_io_value_type

  type, extends(io_value_type) :: real64_arr_io_value_type
    real(real64), allocatable :: value(:)
  end type real64_arr_io_value_type

!===============================================
! constructor interface
!===============================================

  interface io_value_type
    module procedure initialise_io_value_type_sca
    module procedure initialise_io_value_type_arr
  end interface io_value_type

!===============================================
! I/O interfaces
!===============================================

  abstract interface
    subroutine io_read_interface(self, value_name)
      import io_value_type
      class(io_value_type), intent(inout) :: self
      character(*), optional, intent(in)  :: value_name
    end subroutine io_read_interface
  end interface

  abstract interface
    subroutine io_write_interface(self, value_name)
      import io_value_type
      class(io_value_type), intent(in)   :: self
      character(*), optional, intent(in) :: value_name
    end subroutine io_write_interface
  end interface

contains

!===============================================
! numeric conversion functions
!===============================================

pure function real32_from_real32(value) result(converted)
  real(real32), intent(in) :: value
  real(real32) :: converted

  converted = value
end function real32_from_real32

function real32_from_real64(value) result(converted)
  real(real64), intent(in) :: value
  real(real32) :: converted

  if (overflows_real32(value)) then
    call log_event("real64 to real32 conversion overflow", LOG_LEVEL_WARNING)
  end if
  converted = real(value, real32)
end function real32_from_real64

pure function real64_from_real32(value) result(converted)
  real(real32), intent(in) :: value
  real(real64) :: converted

  converted = real(value, real64)
end function real64_from_real32

pure function real64_from_real64(value) result(converted)
  real(real64), intent(in) :: value
  real(real64) :: converted

  converted = value
end function real64_from_real64

pure function real32_arr_from_real32(value) result(converted)
  real(real32), intent(in) :: value(:)
  real(real32) :: converted(size(value))

  converted = value
end function real32_arr_from_real32

function real32_arr_from_real64(value) result(converted)
  real(real64), intent(in) :: value(:)
  real(real32) :: converted(size(value))

  if (any(overflows_real32(value))) then
    call log_event("real64 to real32 array conversion overflow", LOG_LEVEL_WARNING)
  end if
  converted = real(value, real32)
end function real32_arr_from_real64

pure function real64_arr_from_real32(value) result(converted)
  real(real32), intent(in) :: value(:)
  real(real64) :: converted(size(value))

  converted = real(value, real64)
end function real64_arr_from_real32

pure function real64_arr_from_real64(value) result(converted)
  real(real64), intent(in) :: value(:)
  real(real64) :: converted(size(value))

  converted = value
end function real64_arr_from_real64

pure function int32_from_int32(value) result(converted)
  integer(int32), intent(in) :: value
  integer(int32) :: converted

  converted = value
end function int32_from_int32

function int32_from_int64(value) result(converted)
  integer(int64), intent(in) :: value
  integer(int32) :: converted

  if (overflows_int32(value)) then
    call log_event("int64 to int32 conversion overflow", LOG_LEVEL_WARNING)
  end if
  converted = int(value, int32)
end function int32_from_int64

pure function int64_from_int32(value) result(converted)
  integer(int32), intent(in) :: value
  integer(int64) :: converted

  converted = int(value, int64)
end function int64_from_int32

pure function int64_from_int64(value) result(converted)
  integer(int64), intent(in) :: value
  integer(int64) :: converted

  converted = value
end function int64_from_int64

pure function int32_arr_from_int32(value) result(converted)
  integer(int32), intent(in) :: value(:)
  integer(int32) :: converted(size(value))

  converted = value
end function int32_arr_from_int32

function int32_arr_from_int64(value) result(converted)
  integer(int64), intent(in) :: value(:)
  integer(int32) :: converted(size(value))

  if (any(overflows_int32(value))) then
    call log_event("int64 to int32 array conversion overflow", LOG_LEVEL_WARNING)
  end if
  converted = int(value, int32)
end function int32_arr_from_int64

pure function int64_arr_from_int32(value) result(converted)
  integer(int32), intent(in) :: value(:)
  integer(int64) :: converted(size(value))

  converted = int(value, int64)
end function int64_arr_from_int32

pure function int64_arr_from_int64(value) result(converted)
  integer(int64), intent(in) :: value(:)
  integer(int64) :: converted(size(value))

  converted = value
end function int64_arr_from_int64

elemental function overflows_real32(value) result(overflows)
  real(real64), intent(in) :: value
  logical :: overflows

  overflows = .false.
  if (ieee_is_finite(value)) then
    overflows = abs(value) > real(huge(0.0_real32), real64)
  end if
end function overflows_real32

elemental function overflows_int32(value) result(overflows)
  integer(int64), intent(in) :: value
  logical :: overflows
  integer(int64), parameter :: maximum = int(huge(0_int32), int64)
  integer(int64), parameter :: minimum = -maximum - 1_int64

  overflows = value > maximum .or. value < minimum
end function overflows_int32

!===============================================
! constructors
!===============================================

function initialise_io_value_type_sca(key, value) result(instance)

  implicit none

  character(len=*), intent(in) :: key
  class(*), intent(in) :: value

  class(io_value_type), pointer :: instance

  type(int32_io_value_type)  :: concrete_int32_instance
  type(int64_io_value_type)  :: concrete_int64_instance
  type(real32_io_value_type) :: concrete_real32_instance
  type(real64_io_value_type) :: concrete_real64_instance

  select type (value)

    type is (integer(int32))
      call concrete_int32_instance%key_value_initialise( key )
      concrete_int32_instance%value = value
      allocate( instance, source=concrete_int32_instance )

    type is (integer(int64))
      call concrete_int64_instance%key_value_initialise( key )
      concrete_int64_instance%value = value
      allocate( instance, source=concrete_int64_instance )

    type is (real(real32))
      call concrete_real32_instance%key_value_initialise( key )
      concrete_real32_instance%value = value
      allocate( instance, source=concrete_real32_instance )

    type is (real(real64))
      call concrete_real64_instance%key_value_initialise( key )
      concrete_real64_instance%value = value
      allocate( instance, source=concrete_real64_instance )

    class default
      write( log_scratch_space, &
             '("Unhandled value type for io_value_type ''", A, "''.")' ) key
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )

  end select

end function initialise_io_value_type_sca

function initialise_io_value_type_arr(key, value) result(instance)

  implicit none

  character(len=*), intent(in) :: key
  class(*), intent(in) :: value(:)

  class(io_value_type), pointer :: instance

  type(int32_arr_io_value_type)  :: concrete_int32_arr_instance
  type(int64_arr_io_value_type)  :: concrete_int64_arr_instance
  type(real32_arr_io_value_type) :: concrete_real32_arr_instance
  type(real64_arr_io_value_type) :: concrete_real64_arr_instance

  select type (value)

    type is (integer(int32))
      call concrete_int32_arr_instance%key_value_initialise( key )
      allocate( concrete_int32_arr_instance%value, source=value )
      allocate( instance, source=concrete_int32_arr_instance )

    type is (integer(int64))
      call concrete_int64_arr_instance%key_value_initialise( key )
      allocate( concrete_int64_arr_instance%value, source=value )
      allocate( instance, source=concrete_int64_arr_instance )

    type is (real(real32))
      call concrete_real32_arr_instance%key_value_initialise( key )
      allocate( concrete_real32_arr_instance%value, source=value )
      allocate( instance, source=concrete_real32_arr_instance )

    type is (real(real64))
      call concrete_real64_arr_instance%key_value_initialise( key )
      allocate( concrete_real64_arr_instance%value, source=value )
      allocate( instance, source=concrete_real64_arr_instance )

    class default
      write(6,*) "Unhandled value type for io_value_type '", trim(key), "'."
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )

  end select

end function initialise_io_value_type_arr

!===============================================
! I/O subroutines
!===============================================

!> @brief Sets the diagnostic write behaviour for io_value
!> @param[in] write_behaviour Pointer to procedure implementing the write method
subroutine set_write_behaviour(self, write_behaviour)
  class(io_value_type), intent(inout) :: self
  procedure(io_write_interface), pointer, intent(in) :: write_behaviour

  self%write_method => write_behaviour
end subroutine set_write_behaviour

!> @brief Sets the checkpoint write behaviour for the io_value
!> @param[in] write_behaviour A pointer to the checkpoint write behaviour
subroutine set_checkpoint_write_behaviour(self, write_behaviour)
  class(io_value_type), intent(inout) :: self
  procedure(io_write_interface), pointer, intent(in) :: write_behaviour

  self%checkpoint_write_method => write_behaviour
end subroutine set_checkpoint_write_behaviour

!> @brief Sets the checkpoint read behavoiur for the io_value
!> @param[in] read_behaviour A pointer to the checkpoint read behaviour
subroutine set_checkpoint_read_behaviour(self, read_behaviour)
  class(io_value_type), intent(inout) :: self
  procedure(io_read_interface), pointer, intent(in) :: read_behaviour

  self%checkpoint_read_method => read_behaviour
end subroutine set_checkpoint_read_behaviour

!> @brief Subroutine to write to the diagnostic file with write behaviour
subroutine write_value(self, value_name)
  class(io_value_type), intent(inout) :: self
  character(*), optional, intent(in) :: value_name
  if ( associated(self%write_method) ) then
    call self%write_method(value_name)
  else
    call log_event( 'Error trying to write value ' // self%get_key() // &
                    ', write method not set', LOG_LEVEL_ERROR )
  end if

end subroutine write_value

!> @brief Subroutine to write to a checkpoint file with write behaviour
subroutine write_checkpoint(self, value_name)
  class(io_value_type), intent(inout) :: self
  character(*), optional, intent(in) :: value_name

  if ( associated(self%checkpoint_write_method) ) then
    call self%checkpoint_write_method(value_name)
  else
    call log_event( 'Error trying to write value ' // self%get_key() // &
                    ', checkpoint write method not set', LOG_LEVEL_ERROR )
  end if

end subroutine write_checkpoint

!> @brief Subroutine to read data from the checkpoint file to the value
subroutine read_checkpoint(self, value_name)
  class(io_value_type), intent(inout) :: self
  character(*), optional, intent(in) :: value_name

  if ( associated(self%checkpoint_read_method) ) then
    call self%checkpoint_read_method(value_name)
  else
    call log_event( 'Error trying to read value ' // self%get_key() // &
                    ', checkpoint read method not set', LOG_LEVEL_ERROR )
  end if

end subroutine read_checkpoint

!===============================================
! get_value subroutines
!===============================================

subroutine get_int32_value(self, value)

  class(io_value_type), intent(in) :: self
  integer(kind=int32), intent(out) :: value

  ! 'cast' to the expected value type
  select type(concrete_io_value => self)
    type is (int32_io_value_type)
      value = concrete_io_value%value
    class default
      call log_event("Error trying to get int32 value from " // &
                     concrete_io_value%get_key(), LOG_LEVEL_ERROR)
  end select

end subroutine get_int32_value

subroutine get_int64_value(self, value)

  class(io_value_type), intent(in) :: self
  integer(kind=int64), intent(out) :: value

  ! 'cast' to the expected value type
  select type(concrete_io_value => self)
    type is (int64_io_value_type)
      value = concrete_io_value%value
    class default
      call log_event("Error trying to get int64 value from " // &
                     concrete_io_value%get_key(), LOG_LEVEL_ERROR)
  end select

end subroutine get_int64_value

subroutine get_int32_arr_value(self, value)

  class(io_value_type), intent(in) :: self
  integer(kind=int32), intent(out) :: value(:)

  ! 'cast' to the expected value type
  select type(concrete_io_value => self)
    type is (int32_arr_io_value_type)
      value = concrete_io_value%value
    class default
      call log_event("Error trying to get int32 array value from " // &
                     concrete_io_value%get_key(), LOG_LEVEL_ERROR)
  end select

end subroutine get_int32_arr_value

subroutine get_int64_arr_value(self, value)

  class(io_value_type), intent(in) :: self
  integer(kind=int64), intent(out) :: value(:)

  ! 'cast' to the expected value type
  select type(concrete_io_value => self)
    type is (int64_arr_io_value_type)
      value = concrete_io_value%value
    class default
      call log_event("Error trying to get int64 array value from " // &
                     concrete_io_value%get_key(), LOG_LEVEL_ERROR)
  end select

end subroutine get_int64_arr_value

subroutine get_real32_value(self, value)

  class(io_value_type), intent(in) :: self
  real(kind=real32), intent(out) :: value

  ! 'cast' to the expected value type
  select type(concrete_io_value => self)
    type is (real32_io_value_type)
      value = concrete_io_value%value
    class default
      call log_event("Error trying to get real32 value from " // &
                     concrete_io_value%get_key(), LOG_LEVEL_ERROR)
  end select

end subroutine get_real32_value

subroutine get_real64_value(self, value)

  class(io_value_type), intent(in) :: self
  real(kind=real64), intent(out) :: value

  ! 'cast' to the expected value type
  select type(concrete_io_value => self)
    type is (real64_io_value_type)
      value = concrete_io_value%value
    class default
      call log_event("Error trying to get real64 value from " // &
                     concrete_io_value%get_key(), LOG_LEVEL_ERROR)
  end select

end subroutine get_real64_value

subroutine get_real32_arr_value(self, value)

  class(io_value_type), intent(in) :: self
  real(kind=real32), intent(out) :: value(:)

  ! 'cast' to the expected value type
  select type(concrete_io_value => self)
    type is (real32_arr_io_value_type)
      value = concrete_io_value%value
    class default
      call log_event("Error trying to get real32 array value from " // &
                     concrete_io_value%get_key(), LOG_LEVEL_ERROR)
  end select

end subroutine get_real32_arr_value

subroutine get_real64_arr_value(self, value)

  class(io_value_type), intent(in) :: self
  real(kind=real64), intent(out) :: value(:)

  ! 'cast' to the expected value type
  select type(concrete_io_value => self)
    type is (real64_arr_io_value_type)
      value = concrete_io_value%value
    class default
      call log_event("Error trying to get real64 array value from " // &
                     concrete_io_value%get_key(), LOG_LEVEL_ERROR)
  end select

end subroutine get_real64_arr_value

!===============================================
! helper functions
!===============================================

!> @brief A helper function to determine if the io_value_type
!>        can be written to a checkpoint file
!>
!> @return .true. or .false.
function can_write_checkpoint(self) result(checkpointable)

  implicit none

  class(io_value_type), intent(in) :: self
  logical(l_def) :: checkpointable

  if (associated(self%checkpoint_write_method)) then
    checkpointable = .true.
  else
    checkpointable = .false.
  end if

end function can_write_checkpoint

end module io_value_mod
