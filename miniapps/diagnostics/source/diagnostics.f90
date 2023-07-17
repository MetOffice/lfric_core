!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp diagnostics program

!> @brief Program used as proof of concept for testing core infrastructure. Simple miniapp for easy cannibalising.

!> @details Calls init, run and finalise routines from a driver module

program diagnostics

    use cli_mod,                       only : get_initial_filename
    use driver_collections_mod,        only : init_collections, final_collections
    use driver_comm_mod,               only : init_comm, final_comm
    use driver_config_mod,             only : init_config, final_config
    use driver_log_mod,                only : init_logger, final_logger
    use driver_model_data_mod,         only : model_data_type
    use driver_time_mod,               only : init_time, get_calendar
    use diagnostics_configuration_mod, only : required_namelists
    use diagnostics_driver_mod,        only : initialise, step, finalise
    use log_mod,                       only : log_event, log_level_trace
    use model_clock_mod,               only : model_clock_type
    use mpi_mod,                       only : global_mpi

    implicit none

    character(*), parameter :: program_name = "diagnostics"

    ! Model run working data set
    type(model_data_type)               :: model_data
    type(model_clock_type), allocatable :: model_clock

    character(:), allocatable :: filename

    call init_comm( program_name, global_mpi )
    call get_initial_filename( filename )
    call init_config( filename, required_namelists )
    deallocate( filename )
    call init_logger( global_mpi%get_comm(), program_name )
    call init_collections()
    call init_time( model_clock )

    ! Create the depository
    call model_data%depository%initialise(name = 'depository', table_len=100)

    call log_event( 'Initialising ' // program_name // ' ...', &
                    log_level_trace )
    call initialise( model_data, model_clock, global_mpi, &
                     program_name, get_calendar() )

    call log_event( 'Running ' // program_name // ' ...', &
                    log_level_trace )
    do while (model_clock%tick())
      call step( model_data, model_clock )
    end do

    call log_event( 'Finalising ' // program_name // ' ...', &
                    log_level_trace )
    call finalise( model_data )
    call log_event( program_name // ' completed.', log_level_trace )

    call final_collections()
    call final_logger( program_name )
    call final_config()
    call final_comm( global_mpi )

end program diagnostics
