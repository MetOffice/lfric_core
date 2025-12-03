!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

program cli_mod_test

  use, intrinsic :: iso_fortran_env, only : output_unit
  use cli_mod,         only : get_initial_filename

  implicit none

  character(:), allocatable :: filename

  call get_initial_filename( filename )
  write( output_unit, '(A)' ) filename

end program cli_mod_test
