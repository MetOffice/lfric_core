!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module io_value_conversion_test_support
  use io_value_mod, only: io_value_type
  implicit none
  integer :: checkpoint_calls = 0
contains
  subroutine checkpoint_probe(self, value_name)
    class(io_value_type), intent(in) :: self
    character(*), optional, intent(in) :: value_name

    if (self%get_key() /= 'conversion_test') error stop 'callback received wrong object'
    checkpoint_calls = checkpoint_calls + 1
  end subroutine checkpoint_probe
end module io_value_conversion_test_support

program io_value_conversion_test

  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64, error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
                                          ieee_positive_inf, ieee_is_nan, &
                                          ieee_is_finite, ieee_support_datatype
  use io_value_mod, only: to_real32, to_real64, to_int32, to_int64, &
                          io_value_type, real32_arr_io_value_type
  use log_mod, only: initialise_logging, finalise_logging, log_set_alert_stream
  use io_value_conversion_test_support, only: checkpoint_calls, checkpoint_probe
#ifndef NO_MPI
  use mpi, only: MPI_COMM_WORLD, MPI_Init, MPI_Finalize
#endif

  implicit none

  integer :: alerts
#ifndef NO_MPI
  integer :: status
  call MPI_Init(status)
  call initialise_logging(MPI_COMM_WORLD, 'io_value_conversion_test')
#else
  call initialise_logging(0, 'io_value_conversion_test')
#endif

  open(newunit=alerts, status='scratch', action='readwrite')
  call log_set_alert_stream(alerts)
  call test_numeric_values()
  call test_arrays()
  call test_direct_assignment()
  call expect_warnings(alerts, 0)
  call test_overflow()
  call test_ieee_values()
  call finalise_logging()
#ifndef NO_MPI
  call MPI_Finalize(status)
#endif
  print '(A)', 'Numeric I/O conversion tests passed'

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (.not. condition) then
      write(error_unit, '(A)') message
      error stop 1
    end if
  end subroutine check

  subroutine test_numeric_values()
    real(real64) :: rounded

    call check(to_real32(-2.5_real32) == -2.5_real32, 'real32 identity')
    call check(to_real32(3.25_real64) == 3.25_real32, 'real64 to real32')
    call check(to_real64(-2.5_real32) == -2.5_real64, 'real32 to real64')
    call check(to_real64(3.25_real64) == 3.25_real64, 'real64 identity')
    call check(to_int32(-2_int32) == -2_int32, 'int32 identity')
    call check(to_int32(3_int64) == 3_int32, 'int64 to int32')
    call check(to_int64(-2_int32) == -2_int64, 'int32 to int64')
    call check(to_int64(3_int64) == 3_int64, 'int64 identity')
    call check(kind(to_real32(1.0_real64)) == real32, 'real32 result kind')
    call check(kind(to_real64(1.0_real32)) == real64, 'real64 result kind')
    call check(kind(to_int32(1_int64)) == int32, 'int32 result kind')
    call check(kind(to_int64(1_int32)) == int64, 'int64 result kind')

    rounded = 1.0_real64 + 2.0_real64**(-25)
    call check(to_real32(rounded) == 1.0_real32, 'ordinary narrowing rounds')
    call check(to_int32(2147483647_int64) == huge(0_int32), 'int32 upper bound')
    call check(to_int32(-2147483648_int64) == -huge(0_int32)-1_int32, 'int32 lower bound')
    call check(to_real32(real(huge(0.0_real32), real64)) == huge(0.0_real32), 'real32 upper bound')
  end subroutine test_numeric_values

  subroutine test_arrays()
    real(real32), allocatable :: r32(:)
    real(real64), allocatable :: r64(:)
    integer(int32), allocatable :: i32(:)
    integer(int64), allocatable :: i64(:)

    r32 = to_real32([-2.5_real32, 0.0_real32, 3.25_real32])
    call check(all(r32 == [-2.5_real32, 0.0_real32, 3.25_real32]), 'real32 array identity')
    r32 = to_real32([-2.5_real64, 0.0_real64, 3.25_real64])
    call check(all(r32 == [-2.5_real32, 0.0_real32, 3.25_real32]), 'real array narrowing')
    r64 = to_real64(r32)
    call check(all(r64 == [-2.5_real64, 0.0_real64, 3.25_real64]), 'real array widening')
    call check(all(to_real64(r64) == r64), 'real64 array identity')
    i32 = to_int32([-2_int32, 0_int32, 3_int32])
    call check(all(i32 == [-2_int32, 0_int32, 3_int32]), 'int32 array identity')
    i32 = to_int32([-2_int64, 0_int64, 3_int64])
    call check(all(i32 == [-2_int32, 0_int32, 3_int32]), 'integer array narrowing')
    i64 = to_int64(i32)
    call check(all(i64 == [-2_int64, 0_int64, 3_int64]), 'integer array widening')
    call check(all(to_int64(i64) == i64), 'int64 array identity')

    r32 = to_real32([real(real64) ::])
    r64 = to_real64([real(real32) ::])
    i32 = to_int32([integer(int64) ::])
    i64 = to_int64([integer(int32) ::])
    call check(size(r32)+size(r64)+size(i32)+size(i64) == 0, 'empty cross-kind arrays')
    call check(size(to_real32(r32))+size(to_real64(r64)) == 0, 'empty real identity arrays')
    call check(size(to_int32(i32))+size(to_int64(i64)) == 0, 'empty integer identity arrays')
    call check(all(to_int64([-huge(0_int32)-1_int32, huge(0_int32)]) == &
                   [-2147483648_int64, 2147483647_int64]), 'signed integer widening bounds')
    call check(all(to_real32([2.0_real64]) == [2.0_real32]), 'single-element conversion')
  end subroutine test_arrays

  subroutine test_direct_assignment()
    type(real32_arr_io_value_type), target :: stored
    class(io_value_type), pointer :: polymorphic
    real(real64), allocatable :: work(:)

    call stored%key_value_initialise('conversion_test')
    stored%checkpoint_write_method => checkpoint_probe
    stored%value = [1.0_real32]
    polymorphic => stored
    select type(concrete => polymorphic)
      type is(real32_arr_io_value_type)
        work = to_real64(concrete%value)
        work = [2.5_real64, 3.25_real64]
        concrete%value = to_real32(work)
      class default
        error stop 'unexpected concrete I/O type'
    end select
    call check(size(stored%value) == 2, 'direct assignment resizes')
    call check(all(stored%value == [2.5_real32, 3.25_real32]), 'direct assignment converts')
    call check(stored%get_key() == 'conversion_test', 'key preserved')
    call check(associated(stored%checkpoint_write_method), 'callback preserved')
    call stored%write_checkpoint()
    call check(checkpoint_calls == 1, 'checkpoint behavior preserved')
  end subroutine test_direct_assignment

  subroutine expect_warnings(unit, expected)
    integer, intent(in) :: unit, expected
    integer :: ios, count
    character(512) :: line

    flush(unit)
    rewind(unit)
    count = 0
    do
      read(unit, '(A)', iostat=ios) line
      if (ios < 0) exit
      call check(ios == 0, 'reading captured log')
      call check(index(line, ':WARN') > 0, 'overflow must warn, not error')
      call check(index(line, 'overflow') > 0, 'warning identifies overflow')
      count = count + 1
    end do
    call check(count == expected, 'unexpected overflow warning count')
    call log_set_alert_stream(error_unit)
    close(unit)
  end subroutine expect_warnings

  subroutine test_overflow()
    integer :: unit
    integer(int32) :: narrowed_integer
    integer(int32), allocatable :: integers(:)
    real(real32) :: narrowed_real
    real(real32), allocatable :: reals(:)
    real(real64) :: large_real

    open(newunit=unit, status='scratch', action='readwrite')
    call log_set_alert_stream(unit)
    narrowed_integer = to_int32(2147483648_int64)
    narrowed_integer = to_int32(-2147483649_int64)
    integers = to_int32([2147483648_int64, 0_int64, -2147483649_int64])
    call check(size(integers) == 3, 'integer overflow conversion continues')
    large_real = 2.0_real64 * real(huge(0.0_real32), real64)
    narrowed_real = to_real32(large_real)
    narrowed_real = to_real32(-large_real)
    reals = to_real32([large_real, 0.0_real64, -large_real])
    call check(size(reals) == 3, 'real overflow conversion continues')
    call expect_warnings(unit, 6)
  end subroutine test_overflow

  subroutine test_ieee_values()
    integer :: unit
    real(real64) :: special
    real(real32) :: converted

    if (.not. ieee_support_datatype(0.0_real64)) return
    if (.not. ieee_support_datatype(0.0_real32)) return
    open(newunit=unit, status='scratch', action='readwrite')
    call log_set_alert_stream(unit)
    special = ieee_value(0.0_real64, ieee_quiet_nan)
    converted = to_real32(special)
    call check(ieee_is_nan(converted), 'NaN propagates')
    special = ieee_value(0.0_real64, ieee_positive_inf)
    converted = to_real32(special)
    call check(.not. ieee_is_finite(converted), 'infinity propagates')
    call expect_warnings(unit, 0)
  end subroutine test_ieee_values

end program io_value_conversion_test
