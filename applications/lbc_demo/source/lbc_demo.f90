!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page lbc_demo program
!> @brief
!> @details Calls init, run and finalise routines from lbc_demo driver module
program lbc_demo

  use cli_mod,                only: parse_command_line
  use driver_collections_mod, only: init_collections, final_collections

  use constants_mod,       only: precision_real
  use driver_comm_mod,     only: init_comm, final_comm
  use driver_config_mod,   only: init_config, final_config
  use driver_log_mod,      only: init_logger, final_logger
  use driver_modeldb_mod,  only: modeldb_type
  use driver_time_mod,     only: init_time, final_time
  use log_mod,             only: log_event,       &
                                 log_level_trace, &
                                 log_level_error, &
                                 log_scratch_space
  use lfric_mpi_mod,       only: global_mpi
  use lbc_demo_mod,        only: required_namelists
  use lbc_demo_driver_mod, only: initialise, step, finalise

  use base_mesh_config_mod, only: geometry_spherical, &
                                  topology_fully_periodic

  implicit none

  type(modeldb_type) :: modeldb
  character(*), parameter   :: program_name = "lbc_demo"
  character(:), allocatable :: filename

  integer :: geometry, topology

  call parse_command_line( filename )

  write(log_scratch_space, '(A)') &
      'Application built with ' // trim(precision_real) // '-bit real numbers'
  call log_event( log_scratch_space, log_level_trace )

  ! The technical and scientific state
  modeldb%mpi => global_mpi
  call modeldb%configuration%initialise( program_name, table_len=10 )
  call modeldb%config%initialise( program_name )

  call init_comm(program_name, modeldb)

  call init_config(filename, required_namelists,        &
                   configuration=modeldb%configuration, &
                   config=modeldb%config)

  ! Before anything else, test that the mesh provided was a regional domain.
  ! This application is not intended for cubed-sphere meshes.
  geometry = modeldb%config%base_mesh%geometry()
  topology = modeldb%config%base_mesh%topology()

  if ( geometry == geometry_spherical .and. &
       topology == topology_fully_periodic ) then
    ! This is assumed to be a cubedsphere and unsupported by this application
    call log_event( 'Cubed-Sphere mesh is not supported.', log_level_error)
  end if

  call init_logger( modeldb%mpi%get_comm(), program_name )
  call init_collections()
  call init_time(modeldb)

  deallocate( filename )

  ! Create the depository field collection and place it in modeldb
  call modeldb%fields%add_empty_field_collection("depository")
  call modeldb%fields%add_empty_field_collection("lbc")

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
  call final_logger( program_name )
  call final_config()
  call final_comm( modeldb )

end program lbc_demo
