!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module everything_mod

  use constants_mod, only : i_def, l_def, r_def, str_max_filename
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR
  use mpi_mod,       only : global_mpi

  implicit none

  private
  public :: feign_everything_config

  integer(i_def) :: local_rank = -1
  integer(i_def), parameter :: temporary_unit = 3

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_everything_config( cake, &
                                      teapot, &
                                      cheese, &
                                      fish, &
                                      tail, &
                                      school, &
                                      hanger, &
                                      knife )

    use everything_config_mod, only : read_everything_namelist, &
                                      postprocess_everything_namelist, &
                                      key_from_teapot, &
                                      teapot_from_key

    implicit none

    character(*), intent(in) :: cake
    integer(i_def), intent(in) :: teapot
    logical(l_def), intent(in) :: cheese
    real(r_def), intent(in) :: fish
    integer(i_def), intent(in) :: tail
    integer(i_def), intent(in) :: school(:)
    integer(i_def), intent(in) :: hanger(:)
    character(*), intent(in) :: knife(:)

    character(*), parameter :: temp_close_message &
      = "feign_everything_config: Unable to close temporary file"

    integer(i_def)  :: condition
    integer(i_def)  :: i
    character(str_max_filename) :: tmp_str
    character(str_def) :: fmt_str

    if (local_rank == -1) then
      local_rank = global_mpi%get_comm_rank()
    end if

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )

    if (condition /= 0) then
      write( 6, '("feign_everything_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&everything")' )
    write( temporary_unit, '("cake = ''", A, "''")' ) cake
    write( temporary_unit, '("teapot = ''", A, "''")' ) key_from_teapot( teapot )
    write( temporary_unit, '("cheese = ", L2)' ) cheese
    write( temporary_unit, '("fish = ", E14.7)' ) fish
    write( temporary_unit, '("tail = ", I0)' ) tail
    if (size(school) > 1) then
      write( fmt_str,'(A,I0,A)')  "(A,", size(school)-1, "(I0,','),I0)"
      write( temporary_unit, fmt_str ) 'school = ', school
    else
      write( temporary_unit, '("school = ", I0)' ) school
    end if
    if (size(hanger) > 1) then
      write( fmt_str,'(A,I0,A)')  "(A,", size(hanger)-1, "(I0,','),I0)"
      write( temporary_unit, fmt_str ) 'hanger = ', hanger
    else
      write( temporary_unit, '("hanger = ", I0)' ) hanger
    end if
    write( tmp_str,'(A)') "'"//trim(knife(1))//"'"
    if (size(knife) > 1) then
      do i=2, size(knife)
        write( tmp_str,'(A)') trim(tmp_str)//",'"//trim(knife(i))//"'"
      end do
    end if
    write( temporary_unit, '(A)' ) 'knife = '// trim(tmp_str)
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_everything_namelist( temporary_unit, local_rank, scan=.false. )
    call postprocess_everything_namelist()
    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop temp_close_message

  end subroutine feign_everything_config

end module everything_mod
