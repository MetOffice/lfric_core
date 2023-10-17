!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module computed_mod

  use constants_mod, only : i_def
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR
  use mpi_mod,       only : global_mpi

  implicit none

  private
  public :: feign_computed_config

  integer(i_def) :: local_rank = -1
  integer(i_def), parameter :: temporary_unit = 3

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_computed_config( teapot, &
                                    cheese )

    use computed_config_mod, only : read_computed_namelist, &
                                    postprocess_computed_namelist

    implicit none

    integer(i_def), intent(in) :: teapot
    integer(i_def), intent(in) :: cheese

    character(*), parameter :: temp_close_message &
      = "feign_computed_config: Unable to close temporary file"

    integer(i_def)  :: condition

    if (local_rank == -1) then
      local_rank = global_mpi%get_comm_rank()
    end if

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )

    if (condition /= 0) then
      write( 6, '("feign_computed_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&computed")' )
    write( temporary_unit, '("teapot = ", I0)' ) teapot
    write( temporary_unit, '("cheese = ", I0)' ) cheese
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_computed_namelist( temporary_unit, local_rank, scan=.false. )
    call postprocess_computed_namelist()
    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop temp_close_message

  end subroutine feign_computed_config

end module computed_mod
