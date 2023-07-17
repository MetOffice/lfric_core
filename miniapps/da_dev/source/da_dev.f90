!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp da_dev program

!> @brief Tests running the linear model in the miniapp

!> @details Calls init, step and finalise routines from a
!>          DA specific linear driver module in the
!>          lfric-da component

program da_dev

  use cli_mod,                only: get_initial_filename
  use gungho_mod,             only: gungho_required_namelists
  use lfric_da_linear_driver_mod,                     &
                              only: initialise, step, finalise
  use driver_collections_mod, only: init_collections, final_collections
  use driver_comm_mod,        only: init_comm, final_comm
  use driver_config_mod,      only: init_config, final_config
  use driver_log_mod,         only: init_logger, final_logger
  use gungho_modeldb_mod,     only: modeldb_type
  use driver_time_mod,        only: init_time
  use constants_mod,          only: i_native
  use log_mod,                only: log_event, log_level_trace
  use model_clock_mod,        only: model_clock_type
  use mpi_mod,                only: global_mpi

  implicit none

  type(modeldb_type)                  :: modeldb
  type(model_clock_type), allocatable :: model_clock

  character(*), parameter :: program_name = "da_dev"

  character(:), allocatable :: filename

  modeldb%mpi => global_mpi

  call init_comm( program_name, modeldb%mpi )
  call get_initial_filename( filename )
  call init_config( filename, gungho_required_namelists )
  call init_logger( global_mpi%get_comm(), program_name )
  call init_collections()
  deallocate( filename )

  ! Create the depository, prognostics and diagnostics field collections
  call modeldb%model_data%depository%initialise(name='depository', table_len=100)
  call modeldb%model_data%prognostic_fields%initialise(name="prognostics", table_len=100)
  call modeldb%model_data%diagnostic_fields%initialise(name="diagnostics", table_len=100)

  call log_event( 'Initialising ' // program_name // ' ...', log_level_trace )
  call init_time( model_clock )
  call initialise( program_name, modeldb, model_clock)

  call log_event( 'Running ' // program_name // ' ...', log_level_trace )
  do while ( model_clock%tick() )
    call step( modeldb, model_clock )
  end do

  call log_event( 'Finalising ' // program_name // ' ...', log_level_trace )
  call finalise( program_name, modeldb, model_clock )
  deallocate( model_clock )

  call final_collections()
  call final_logger( program_name )
  call final_config()
  call final_comm( modeldb%mpi )

end program da_dev
