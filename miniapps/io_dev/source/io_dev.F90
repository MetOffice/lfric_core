!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @page io_dev io_dev Miniapp
!> Test program for the XIOS IO implementation.

!> @brief Main program used to test XIOS setup and output of a field.

program io_dev

  use cli_mod,                only: get_initial_filename
  use constants_mod,          only: precision_real
  use driver_collections_mod, only: init_collections, final_collections
  use driver_comm_mod,        only: init_comm, final_comm
  use driver_config_mod,      only: init_config, final_config
  use driver_log_mod,         only: init_logger, final_logger
  use driver_time_mod,        only: init_time, get_calendar
  use driver_timer_mod,       only: init_timers, final_timers
  use io_dev_mod,             only: io_dev_required_namelists
  use io_dev_driver_mod,      only: initialise, step, finalise
  use io_dev_data_mod,        only: io_dev_data_type
  use log_mod,                only: log_event,       &
                                    log_level_trace, &
                                    log_scratch_space
  use model_clock_mod,        only: model_clock_type
  use mpi_mod,                only: global_mpi

  implicit none

  character(*), parameter :: program_name = "io_dev"

  type (io_dev_data_type)             :: model_data
  type(model_clock_type), allocatable :: model_clock

  character(:), allocatable :: filename

  write( log_scratch_space,                                       &
         '("Application built with ", A, "-bit real numbers")' ) &
       precision_real
  call log_event( log_scratch_space, log_level_trace )

  call init_comm( "io_dev", global_mpi )
  call get_initial_filename( filename )
  call init_config( filename, io_dev_required_namelists )
  deallocate( filename )
  call init_logger( global_mpi%get_comm(), program_name )
  call init_timers( program_name )
  call init_collections()
  call init_time( model_clock )

  call log_event( 'Initialising '//program_name//' ...', log_level_trace )
  call initialise( model_data, model_clock, global_mpi, &
                   program_name, get_calendar() )

  write(log_scratch_space,'("Running ", A, " ...")') program_name
  call log_event( log_scratch_space, log_level_trace )
  do while( model_clock%tick() )
    call step( model_data, model_clock, program_name )
  end do

  call log_event( 'Finalising '//program_name//' ...', log_level_trace )
  call finalise( model_data )

  call final_collections()
  call final_timers( program_name )
  call final_logger( program_name )
  call final_config()
  call final_comm( global_mpi )

end program io_dev
