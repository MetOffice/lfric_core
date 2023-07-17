!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp skeleton program

!> @brief Main program used to illustrate how to write LFRic miniapps.

!> @details Calls init, run and finalise routines from a driver module

program skeleton

  use cli_mod,                only : get_initial_filename
  use driver_collections_mod, only : init_collections, final_collections
  use constants_mod,          only : precision_real
  use driver_comm_mod,        only : init_comm, final_comm
  use driver_config_mod,      only : init_config, final_config
  use driver_log_mod,         only : init_logger, final_logger
  use driver_time_mod,        only : init_time, get_calendar
  use log_mod,                only : log_event,       &
                                     log_level_trace, &
                                     log_scratch_space
  use model_clock_mod,        only : model_clock_type
  use mpi_mod,                only : global_mpi
  use skeleton_mod,           only : skeleton_required_namelists
  use skeleton_driver_mod,    only : initialise, step, finalise

  implicit none

  character(*), parameter :: program_name = "skeleton"

  character(:), allocatable :: filename

  type(model_clock_type), allocatable :: model_clock

  write(log_scratch_space,&
        '("Application built with ", A, "-bit real numbers")') &
        trim(precision_real)
  call log_event( log_scratch_space, log_level_trace )

  call init_comm("skeleton", global_mpi)
  call get_initial_filename( filename )
  call init_config( filename, skeleton_required_namelists )
  deallocate( filename )
  call init_logger( global_mpi%get_comm(), program_name )
  call init_collections()
  call init_time( model_clock )

  call log_event( 'Initialising ' // program_name // ' ...', log_level_trace )
  call initialise( global_mpi, model_clock, program_name, get_calendar() )

  do while (model_clock%tick())
    call step( program_name )
  end do

  call log_event( 'Finalising ' // program_name // ' ...', log_level_trace )
  call finalise( program_name )

  call final_collections()
  call final_logger( program_name )
  call final_config()
  call final_comm( global_mpi )

end program skeleton
