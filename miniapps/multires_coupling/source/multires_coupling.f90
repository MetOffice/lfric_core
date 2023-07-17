!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp multires_coupling program

!> @brief Main program used for multires_coupling miniapp

!> @details Calls init, run and finalise routines from a driver module

program multires_coupling

  use cli_mod,                      only : get_initial_filename
  use coupling_test_alg_mod,        only : coupling_test_alg
  use driver_collections_mod,       only : init_collections, final_collections
  use driver_comm_mod,              only : init_comm, final_comm
  use driver_config_mod,            only : init_config, final_config
  use driver_log_mod,               only : init_logger, final_logger
  use driver_time_mod,              only : init_time, get_calendar
  use driver_timer_mod,             only : init_timers, final_timers
  use gungho_modeldb_mod,           only : modeldb_type
  use log_mod,                      only : log_event,       &
                                           log_level_error, &
                                           log_level_trace
  use mpi_mod,                      only : global_mpi
  use multires_coupling_mod,        only : multires_required_namelists
  use multires_coupling_config_mod, only : multires_coupling_mode, &
                                           multires_coupling_mode_test
  use multires_coupling_driver_mod, only : initialise, step, finalise

  implicit none

  character(*), parameter :: program_name = "multires_coupling"

  ! Model run working data set
  type (modeldb_type) :: dynamics_mesh_modeldb
  type (modeldb_type) :: physics_mesh_modeldb

  character(:), allocatable :: filename

  logical :: physics_running

  dynamics_mesh_modeldb%mpi => global_mpi
  physics_mesh_modeldb%mpi => global_mpi

  call init_comm( program_name, dynamics_mesh_modeldb%mpi )
  call get_initial_filename( filename )
  call init_config( filename, multires_required_namelists )
  deallocate( filename )
  call init_logger( dynamics_mesh_modeldb%mpi%get_comm(), program_name )
  call init_timers( program_name )
  call init_collections()
  !
  ! Running two clocks like this is fairely unpleasent but until we have
  ! clock proxies it seems the best way to progress.
  !
  call init_time( dynamics_mesh_modeldb%clock )
  call init_time( physics_mesh_modeldb%clock )

  ! Create the dynamics depository, prognostics and diagnostics field collections
  call dynamics_mesh_modeldb%model_data%depository%initialise(name='depository', table_len=100)
  call dynamics_mesh_modeldb%model_data%prognostic_fields%initialise(name="prognostics", table_len=100)
  call dynamics_mesh_modeldb%model_data%diagnostic_fields%initialise(name="diagnostics", table_len=100)

  ! Create the physics depository, prognostics and diagnostics field collections
  call physics_mesh_modeldb%model_data%depository%initialise(name='depository', table_len=100)
  call physics_mesh_modeldb%model_data%prognostic_fields%initialise(name="prognostics", table_len=100)
  call physics_mesh_modeldb%model_data%diagnostic_fields%initialise(name="diagnostics", table_len=100)

  call log_event( 'Initialising' // program_name // ' ...', log_level_trace )
  call initialise( dynamics_mesh_modeldb, &
                   physics_mesh_modeldb,  &
                   program_name,          &
                   get_calendar() )

  if ( multires_coupling_mode == multires_coupling_mode_test ) then
    ! Call coupling test algorithm
    call coupling_test_alg()
  else
    ! Do timestep
    do while (dynamics_mesh_modeldb%clock%tick())
      !
      ! Running two clocks like this is fairely unpleasent but until we have
      ! clock proxies it seems the best way to progress.
      !
      physics_running = physics_mesh_modeldb%clock%tick()
      if (.not. physics_running) then
        call log_event( "Physics clock ended before dynamics clock", &
                        log_level_error )
      end if
      call step( dynamics_mesh_modeldb, physics_mesh_modeldb )
    end do
  end if

  call log_event( 'Finalising ' // program_name // ' ...', log_level_trace )
  call finalise( dynamics_mesh_modeldb, &
                 physics_mesh_modeldb,  &
                 program_name )

  call final_collections()
  call final_timers( program_name )
  call final_logger( program_name )
  call final_config()
  call final_comm( dynamics_mesh_modeldb%mpi )

end program multires_coupling
