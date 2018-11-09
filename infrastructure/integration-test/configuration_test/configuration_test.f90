!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

program configuration_test

  use iso_fortran_env,             only : error_unit
  use mpi_mod,                     only : initialise_comm, store_comm,       &
                                          finalise_comm,                     &
                                          get_comm_rank
  use one_of_each_test_config_mod, only : key_from_an_enum,                  &
                                      postprocess_one_of_each_test_namelist, &
                                          read_one_of_each_test_namelist,    &
                                          a_dim,                             &
                                          angle_deg,                         &
                                          angle_rad,                         &
                                          an_enum,                           &
                                          bounded_array_local_dim,           &
                                          bounded_array1_namelist_dim,       &
                                          bounded_array2_namelist_dim,       &
                                          bounded_array_source_dim,          &
                                          closed_array,                      &
                                          max_array_size,                    &
                                          open_array,                        &
                                          some_string,                       &
                                          whole_number

   use another_list_config_mod, only : postprocess_another_list_namelist,    &
                                       read_another_list_namelist,           &
                                       some_other_dim

  implicit none

  integer,      parameter :: file_unit = 13
  character(*), parameter :: filename = 'configuration_test.nml'

  integer       :: comm
  integer       :: rank
  integer       :: condition
  character(20) :: result_filename
  character(40) :: format_string

  ! Initialse mpi and create the default communicator: mpi_comm_world
  call initialise_comm(comm)

  ! Save lfric's part of the split communicator for later use
  call store_comm(comm)

  rank = get_comm_rank()

  open( file_unit, file=filename, iostat=condition )
  if (condition /= 0) then
    write( error_unit, '("Failed to open file: ",A)' ) filename
    stop 3
  end if

  call read_another_list_namelist( file_unit ,rank )
  REWIND(file_unit)
  call read_one_of_each_test_namelist( file_unit ,rank )

  close( file_unit, iostat=condition )
  if (condition /= 0) then
    write( error_unit, '("Failed to close file: ", A)' ) filename
    stop 4
  end if

  call postprocess_another_list_namelist
  call postprocess_one_of_each_test_namelist

  write( result_filename, '("result.", I0, ".txt")' ) rank
  open( file_unit, file=result_filename, iostat=condition )
  if (condition /= 0) then
    write( error_unit, '("Failed to open file: ",A)' ) result_filename
    stop 5
  end if

  write( file_unit, '(I0, A, I0)' ) rank, ' a_dim: ', a_dim
  write( file_unit, '(I0, " angle_deg: ", E14.7)' ) rank, angle_deg
  write( file_unit, '(I0, " angle_rad: ", E14.7)' ) rank, angle_rad
  write( file_unit, '(I0, " an_enum: ", A)' ) rank, &
                                              trim(key_from_an_enum( an_enum ))

  write( format_string, '(A,I0,A)') '(I0, A, ', max_array_size,  'E14.7)'
  write( file_unit, format_string ) &
    rank, ' bounded_array_local_dim:', bounded_array_local_dim
  write( file_unit, format_string ) &
    rank, ' bounded_array1_namelist_dim:', bounded_array1_namelist_dim
  write( file_unit, format_string ) &
    rank, ' bounded_array2_namelist_dim:', bounded_array2_namelist_dim
  write( file_unit, format_string ) &
    rank, ' bounded_array_source_dim:', bounded_array_source_dim
  write( file_unit, format_string ) rank, ' closed_array: ', closed_array

  write( format_string, &
         '("(I0, "" open_array: "", ", I0, "I3)")' ) max_array_size
  write( file_unit, format_string ) rank, open_array
  write( file_unit, '(I0, " some_string: ''", A, "''")' ) rank, &
                                                           trim(some_string)
  write( file_unit, '(I0, " whole_number: ", I0)' ) rank, whole_number

  close( file_unit, iostat=condition )
  if (condition /= 0) then
    write( error_unit, '("Failed to close file: ", A)' ) result_filename
    stop 6
  end if

  call finalise_comm()

end program configuration_test
