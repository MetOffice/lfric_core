!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp skeleton program

!> @brief Main program used to illustrate how to write LFRic applications.

!> @details Calls init, step and finalise routines from a driver module

program skeleton

  use cli_mod,                 only: parse_command_line
  use constants_mod,           only: precision_real
  use driver_collections_mod,  only: init_collections, final_collections
  use driver_comm_mod,         only: init_comm, final_comm
  use driver_config_mod,       only: init_config, final_config
  use driver_log_mod,          only: init_logger, final_logger
  use driver_modeldb_mod,      only: modeldb_type
  use driver_time_mod,         only: init_time, final_time
  use lfric_mpi_mod,           only: global_mpi
  use log_mod,                 only: log_event,       &
                                     log_level_trace, &
                                     log_scratch_space

  use skeleton_mod,        only: skeleton_required_namelists
  use skeleton_driver_mod, only: initialise, step, finalise

  implicit none

  ! The technical and scientific state
  type(modeldb_type) :: modeldb

  character(*), parameter   :: program_name = "skeleton"
  character(:), allocatable :: filename

  call parse_command_line( filename )
  call modeldb%configuration%initialise( program_name, table_len=10 )
  call modeldb%config%initialise(program_name)

  write(log_scratch_space,'(A)')                          &
      'Application built with '// trim(precision_real) // &
      '-bit real numbers.'
  call log_event( log_scratch_space, log_level_trace )

  modeldb%mpi => global_mpi

  call init_comm( "skeleton", modeldb )
  call init_config( filename, skeleton_required_namelists, &
                    configuration=modeldb%configuration,   &
                    config=modeldb%config )

  call init_logger( modeldb%mpi%get_comm(), program_name )
  call init_collections()
  call init_time( modeldb )
  deallocate( filename )

  ! Create the depository field collection and place it in modeldb
  call modeldb%fields%add_empty_field_collection("depository")

  call modeldb%io_contexts%initialise(program_name, 100)

  call log_event( 'Initialising ' // program_name // ' ...', log_level_trace )
  call initialise( program_name, modeldb )

  do while (modeldb%clock%tick())
    call step( program_name, modeldb )
  end do

  call log_event( 'Finalising ' // program_name // ' ...', log_level_trace )
  call finalise( program_name, modeldb )

  call final_time( modeldb )
  call final_collections()
  call final_logger( program_name )
  call final_config()
  call final_comm( modeldb )

end program skeleton
