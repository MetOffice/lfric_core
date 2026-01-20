!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page Miniapp io_demo program
!> @brief Main program used to calculate diffusion of randomly initialised theta field
!> @details Calls init, run and finalise routines from io_demo driver module
program io_demo

  use cli_mod,                     only : parse_command_line
  use driver_collections_mod,      only : init_collections, final_collections
  use constants_mod,               only : precision_real
  use driver_comm_mod,             only : init_comm, final_comm
  use driver_config_mod,           only : init_config, final_config
  use driver_log_mod,              only : init_logger, final_logger
  use driver_modeldb_mod,          only : modeldb_type
  use driver_time_mod,             only : init_time, final_time
  use lfric_mpi_mod,               only : global_mpi
  use log_mod,                     only : log_event,       &
                                          log_level_trace, &
                                          log_scratch_space
  use random_number_generator_mod, only : random_number_generator_type
  use io_demo_mod,                 only : io_demo_required_namelists
  use io_demo_driver_mod,          only : initialise, step, finalise
  use timing_mod,                  only : init_timing, final_timing
  use io_config_mod,               only : timer_output_path
  use namelist_mod,                only : namelist_type

  implicit none

  ! The technical and scientific state
  type(modeldb_type)           :: modeldb
  character(*), parameter      :: program_name = "io_demo"
  character(:), allocatable    :: filename
  type(namelist_type), pointer :: io_nml
  logical                      :: lsubroutine_timers
  integer, parameter           :: default_seed = 123456789
  type(random_number_generator_type), pointer :: rng

  call parse_command_line( filename )
  call modeldb%values%initialise()
  call modeldb%configuration%initialise( program_name, table_len=10 )

  write(log_scratch_space,&
        '("Application built with ", A, "-bit real numbers")') &
        trim(precision_real)
  call log_event( log_scratch_space, log_level_trace )
  modeldb%mpi => global_mpi
  call init_comm(program_name, modeldb)
  call init_config( filename,                            &
                    io_demo_required_namelists, &
                    modeldb%configuration )
  deallocate( filename )

  call init_logger( modeldb%mpi%get_comm(), program_name )
  io_nml => modeldb%configuration%get_namelist('io')
  call io_nml%get_value('subroutine_timers', lsubroutine_timers)
  call init_timing( modeldb%mpi%get_comm(), lsubroutine_timers, program_name, timer_output_path )
  nullify( io_nml )
  call init_collections()
  call init_time(modeldb)

  allocate(rng, source=random_number_generator_type(default_seed))
  call modeldb%values%add_key_value("rng", rng)

  ! Create the depository field collection and place it in modeldb
  call modeldb%fields%add_empty_field_collection("depository")

  call modeldb%io_contexts%initialise(program_name, 100)
  call log_event( 'Initialising ' // program_name // ' ...', log_level_trace )
  call initialise( program_name, modeldb)

  do while (modeldb%clock%tick())
    call step( program_name, modeldb )
  end do

  call log_event( 'Finalising ' // program_name // ' ...', log_level_trace )
  call finalise( program_name, modeldb )
  call final_time(modeldb)
  call final_collections()
  call final_timing( program_name )
  call final_logger( program_name )
  call final_config()
  call final_comm( modeldb )

end program io_demo
