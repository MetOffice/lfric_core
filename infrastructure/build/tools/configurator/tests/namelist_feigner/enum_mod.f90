!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module enumeration_mod

  use constants_mod, only : i_def
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR
  use mpi_mod,       only : global_mpi

  implicit none

  private
  public :: feign_enum_config

  integer(i_def) :: local_rank = -1
  integer(i_def), parameter :: temporary_unit = 3

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_enum_config( thing )

    use enum_config_mod, only : read_enum_namelist, &
                                postprocess_enum_namelist, &
                                key_from_thing, &
                                thing_from_key

    implicit none

    integer(i_def), intent(in) :: thing

    character(*), parameter :: temp_close_message &
      = "feign_enum_config: Unable to close temporary file"

    integer(i_def)  :: condition

    if (local_rank == -1) then
      local_rank = global_mpi%get_comm_rank()
    end if

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )

    if (condition /= 0) then
      write( 6, '("feign_enum_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&enum")' )
    write( temporary_unit, '("thing = ''", A, "''")' ) key_from_thing( thing )
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_enum_namelist( temporary_unit, local_rank, scan=.false. )
    call postprocess_enum_namelist()
    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop temp_close_message

  end subroutine feign_enum_config

end module enumeration_mod
