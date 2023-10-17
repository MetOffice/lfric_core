!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module simple_mod

  use constants_mod, only : i_def, l_def, r_double, str_def
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR
  use mpi_mod,       only : global_mpi

  implicit none

  private
  public :: feign_simple_config

  integer(i_def) :: local_rank = -1
  integer(i_def), parameter :: temporary_unit = 3

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_simple_config( foo, &
                                  bar, &
                                  baz, &
                                  fred )

    use simple_config_mod, only : read_simple_namelist, &
                                  postprocess_simple_namelist

    implicit none

    integer(i_def), intent(in) :: foo
    real(r_double), intent(in) :: bar
    character(*), intent(in) :: baz
    logical(l_def), intent(in) :: fred

    character(*), parameter :: temp_close_message &
      = "feign_simple_config: Unable to close temporary file"

    integer(i_def)  :: condition

    if (local_rank == -1) then
      local_rank = global_mpi%get_comm_rank()
    end if

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )

    if (condition /= 0) then
      write( 6, '("feign_simple_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&simple")' )
    write( temporary_unit, '("foo = ", I0)' ) foo
    write( temporary_unit, '("bar = ", E14.7)' ) bar
    write( temporary_unit, '("baz = ''", A, "''")' ) baz
    write( temporary_unit, '("fred = ", L2)' ) fred
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_simple_namelist( temporary_unit, local_rank, scan=.false. )
    call postprocess_simple_namelist()
    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop temp_close_message

  end subroutine feign_simple_config

end module simple_mod
