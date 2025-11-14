!-------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Module providing string utility functions.
module lfric_string_mod
  use constants_mod, only : i_def
  implicit none

  private
  public :: split_string

contains

  !> @brief Splits a string into parts based on a specified delimiter.
  !> @detail This function takes an input string and a delimiter string,
  !> and splits the input string into an array of substrings wherever the
  !> delimiter occurs. The resulting substrings are returned as an allocatable
  !> array of strings. Leading and Trailing delimiters result in empty strings
  !> in the output.
  !> @param[in] input The input string to be split.
  !> @param[in] delimiter The delimiter string used to split the input.
  !> @return An allocatable array of strings containing the split parts.
  function split_string(input, delimiter) result(parts)
    implicit none
    character(*), intent(in) :: input
    character(*), intent(in) :: delimiter
    character(:), allocatable :: parts(:)
    integer(i_def) :: i, start, end_pos, delim_len, part_count, pos
    delim_len = len(delimiter)
    part_count = 1
    pos = 1

    if(delim_len /= 0) then
      ! Count number of parts
      do
        end_pos = index(input(pos:), delimiter)
        if (end_pos == 0) exit
        part_count = part_count + 1
        pos = pos + end_pos + delim_len - 1
        if (pos > len_trim(input)) exit
      end do
      allocate(character(len(input)) :: parts(part_count))
      ! Split string into parts.
      start = 1
      do i = 1, part_count
        end_pos = index(input(start:), delimiter)
        if (end_pos == 0) then
          parts(i) = input(start:)
        else
          parts(i) = input(start:start + end_pos - 2)
          start = start + end_pos + delim_len - 1
        end if
      end do
    else
      allocate(character(len(input)) :: parts(1))
      parts = input
    endif

  end function split_string

end module lfric_string_mod