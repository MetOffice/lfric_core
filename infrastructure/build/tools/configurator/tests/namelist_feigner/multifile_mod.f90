!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module multifile_mod

  use constants_mod, only : i_def, l_def, r_def, str_max_filename
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR
  use mpi_mod,       only : global_mpi

  implicit none

  private
  public :: feign_first_config, &
            feign_second_config

  integer(i_def) :: local_rank = -1
  integer(i_def), parameter :: temporary_unit = 3

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_first_config( cake, &
                                 teapot, &
                                 cheese )

    use first_config_mod, only : read_first_namelist, &
                                 postprocess_first_namelist, &
                                 key_from_teapot, &
                                 teapot_from_key

    implicit none

    character(*), intent(in) :: cake
    integer(i_def), intent(in) :: teapot
    logical(l_def), intent(in) :: cheese

    character(*), parameter :: temp_close_message &
      = "feign_first_config: Unable to close temporary file"

    integer(i_def)  :: condition

    if (local_rank == -1) then
      local_rank = global_mpi%get_comm_rank()
    end if

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )

    if (condition /= 0) then
      write( 6, '("feign_first_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&first")' )
    write( temporary_unit, '("cake = ''", A, "''")' ) cake
    write( temporary_unit, '("teapot = ''", A, "''")' ) key_from_teapot( teapot )
    write( temporary_unit, '("cheese = ", L2)' ) cheese
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_first_namelist( temporary_unit, local_rank, scan=.false. )
    call postprocess_first_namelist()
    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop temp_close_message

  end subroutine feign_first_config

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_second_config( fish, &
                                  yarn, &
                                  tail )

    use second_config_mod, only : read_second_namelist, &
                                  postprocess_second_namelist, &
                                  key_from_yarn, &
                                  yarn_from_key

    implicit none

    real(r_def), intent(in) :: fish
    integer(i_def), intent(in) :: yarn
    integer(i_def), intent(in) :: tail

    character(*), parameter :: temp_close_message &
      = "feign_second_config: Unable to close temporary file"

    integer(i_def)  :: condition

    if (local_rank == -1) then
      local_rank = global_mpi%get_comm_rank()
    end if

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )

    if (condition /= 0) then
      write( 6, '("feign_second_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&second")' )
    write( temporary_unit, '("fish = ", E14.7)' ) fish
    write( temporary_unit, '("yarn = ''", A, "''")' ) key_from_yarn( yarn )
    write( temporary_unit, '("tail = ", I0)' ) tail
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_second_namelist( temporary_unit, local_rank, scan=.false. )
    call postprocess_second_namelist()
    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop temp_close_message

  end subroutine feign_second_config

end module multifile_mod
