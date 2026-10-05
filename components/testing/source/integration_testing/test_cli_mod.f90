!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise.

!> @brief Utilities for parsing command-line arguments in tests
!> @details
!>   Provides a convenient interface for retrieving command-line arguments
!>   as allocatable strings with automatic length determination. Useful for
!>   test programs that need to accept runtime configuration parameters.
module test_cli_mod

  use, intrinsic :: iso_fortran_env, only : error_unit

  implicit none

  private
  public get_cli_argument

contains

  !> @brief Retrieve a command-line argument as an allocatable string
  !> @param [in] idx  Index of the command-line argument to retrieve.
  !>                  Defaults to 1 if not provided.
  !> @return          Allocatable string containing the requested argument
  !> @details
  !>   Retrieves a command-line argument and returns it as an allocatable string
  !>   with exactly the length needed to store the argument value.
  !>   Terminates with error if argument retrieval fails.
  function get_cli_argument( idx ) result(argument)

    implicit none

    integer, intent(in), optional :: idx
    character(:), allocatable :: argument

    integer :: argument_index
    integer :: argument_length
    integer :: status

    if (present(idx)) then
      argument_index = idx
    else
      argument_index = 1
    end if

    call get_command_argument( argument_index,         &
                               length=argument_length, &
                               status=status )
    ! Check for errors retrieving argument length
    if (status /= 0) then
        write( error_unit, &
               '("Unable to get command line argument (", I0, ")")' ) status
        ! Terminate test with error code
        error stop 1
    end if

    ! Allocate string with exact length needed for the argument
    allocate( character(argument_length) :: argument )
    ! Retrieve the actual argument value
    call get_command_argument( argument_index, value=argument )

  end function get_cli_argument

end module test_cli_mod
