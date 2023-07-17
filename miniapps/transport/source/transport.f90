!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page transport Transport miniapp
!> Program file for running transport miniapp. Subroutine calls include initialise_transport(),
!> run_transport() and finalise_transport().
program transport

  use cli_mod,                only: get_initial_filename
  use constants_mod,          only: i_def, r_def
  use driver_collections_mod, only: init_collections, final_collections
  use driver_comm_mod,        only: init_comm, final_comm
  use driver_config_mod,      only: init_config, final_config
  use driver_log_mod,         only: init_logger, final_logger
  use driver_time_mod,        only: init_time, get_calendar
  use driver_timer_mod,       only: init_timers, final_timers
  use log_mod,                only: log_event,       &
                                    log_level_trace, &
                                    log_scratch_space
  use model_clock_mod,        only: model_clock_type
  use mpi_mod,                only: global_mpi
  use transport_mod,          only: transport_required_namelists
  use transport_driver_mod,   only: initialise_transport, &
                                    step_transport,        &
                                    finalise_transport

  implicit none

  character(*), parameter :: program_name = "transport"

  character(:), allocatable :: filename

  type(model_clock_type), allocatable :: model_clock

  call log_event( 'Miniapp will run with default precision set as:', &
                  log_level_trace )
  write(log_scratch_space, '("        r_def kind = ", I0)') kind(1.0_r_def)
  call log_event( log_scratch_space , log_level_trace )
  write(log_scratch_space, '("        i_def kind = ", I0)') kind(1_i_def)
  call log_event( log_scratch_space , log_level_trace )

  call init_comm( program_name, global_mpi )
  call get_initial_filename( filename )
  call init_config( filename, transport_required_namelists )
  deallocate( filename )
  call init_logger( global_mpi%get_comm(), program_name )
  call init_timers( program_name )
  call init_collections()
  call init_time( model_clock )

  call log_event( 'Initialising ' // program_name // ' ...', log_level_trace )
  call initialise_transport( global_mpi, model_clock, &
                             program_name, get_calendar() )

  call log_event( 'Running ' // program_name // ' ...', log_level_trace )
  do while (model_clock%tick())
    call step_transport( model_clock )
  end do

  call log_event( 'Finalising ' // program_name // ' ...', log_level_trace )
  call finalise_transport( program_name )

  call final_collections()
  call final_timers( program_name )
  call final_logger( program_name )
  call final_config()
  call final_comm( global_mpi )

end program transport
