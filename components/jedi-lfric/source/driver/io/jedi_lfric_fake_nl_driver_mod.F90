!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the fake nl DA model
!>
module jedi_lfric_fake_nl_driver_mod

  use base_mesh_config_mod,     only: prime_mesh_name
  use checksum_alg_mod,         only: checksum_alg
  use constants_mod,            only: str_def
  use driver_model_data_mod,    only: model_data_type
  use driver_time_mod,          only: init_time, get_calendar
  use driver_mesh_mod,          only: init_mesh
  use driver_fem_mod,           only: init_fem, final_fem
  use inventory_by_mesh_mod,    only: inventory_by_mesh_type
  use field_mod,                only: field_type
  use jedi_lfric_fake_nl_init_mod, &
                                only: create_da_model_data, &
                                      initialise_da_model_data
  use log_mod,                  only: log_event,          &
                                      log_scratch_space,  &
                                      LOG_LEVEL_ALWAYS,   &
                                      LOG_LEVEL_INFO
  use mesh_mod,                 only: mesh_type
  use mesh_collection_mod,      only: mesh_collection
  use extrusion_mod,            only: extrusion_type, TWOD
  use model_clock_mod,          only: model_clock_type
  use mpi_mod,                  only: mpi_type
  use jedi_lfric_increment_alg_mod, &
                                only: jedi_lfric_increment_alg
  !> @todo: Test code should not appear in the component
  !> @{
  use jedi_lfric_tests_config_mod, &
                                only: write_data, test_field
  !> @}
  use jedi_lfric_fake_nl_extrusion_mod, &
                                only: create_extrusion
  use jedi_lfric_fake_nl_init_files_mod, &
                                only: init_jedi_lfric_files

#ifdef USE_XIOS
  use driver_io_mod,            only: init_io, final_io, filelist_populator, &
                                      get_io_context
  use io_config_mod,            only: use_xios_io
  use io_context_mod,           only: io_context_type
  use lfric_xios_context_mod,   only: lfric_xios_context_type, advance
  use lfric_xios_write_mod,     only: write_state
#endif

  implicit none

  private
  ! To run local LFRic mini-app
  public initialise, step, finalise

  ! To be moved at a later date
  type( model_clock_type ), public, allocatable :: model_clock

  ! Coordinate field
  type(field_type), pointer, public :: chi(:)    => null()
  type(field_type), pointer, public :: panel_id  => null()
  type(mesh_type),  pointer, public :: mesh      => null()
  type(mesh_type),  pointer, public :: twod_mesh => null()
  type(inventory_by_mesh_type)      :: chi_inventory
  type(inventory_by_mesh_type)      :: panel_id_inventory

contains


  !> @brief Initialise the model mesh, fem, clock, and IO
  !>
  !> @param [inout] program_name The program name
  !> @param [inout] mpi          The mpi communicator
  subroutine initialise( program_name, mpi )

    implicit none

    character(len=*), intent(inout) :: program_name
    class(mpi_type),  intent(inout) :: mpi

    class(extrusion_type), allocatable :: extrusion
    character(str_def)                 :: base_mesh_names(1)

    ! Create the mesh
    allocate( extrusion, source=create_extrusion() )
    base_mesh_names(1) = prime_mesh_name
    call init_mesh( mpi%get_comm_rank(), mpi%get_comm_size(), &
                    base_mesh_names, input_extrusion = extrusion )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )

    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)
    call chi_inventory%get_field_array(mesh, chi)
    call panel_id_inventory%get_field(mesh, panel_id)

    call init_time( model_clock )

#ifdef USE_XIOS
    if ( use_xios_io ) then
      call initialise_io( program_name, mpi, model_clock )
    end if
#endif

  end subroutine initialise


#ifdef USE_XIOS
  !> @brief Initialise the model IO
  !>
  !> @param [inout] program_name The program name
  !> @param [inout] mpi          The mpi communicator
  !> @param [inout] model_clock  The model clock
  subroutine initialise_io( program_name, mpi, model_clock )

    implicit none

    character(len=*),        intent(inout) :: program_name
    class(mpi_type),         intent(inout) :: mpi
    type(model_clock_type),  intent(inout) :: model_clock

    procedure(filelist_populator), pointer :: fl_populator => null()
    class(io_context_type),        pointer :: model_io_context => null()

    ! Initialise I/O context
    fl_populator => init_jedi_lfric_files
    call init_io( program_name, mpi%get_comm(), chi_inventory, panel_id_inventory, &
    model_clock, get_calendar(), populate_filelist=fl_populator )

    ! Do initial step
    model_io_context => get_io_context()
    if (model_clock%is_initialisation()) then
      select type (model_io_context)
      type is (lfric_xios_context_type)
        call advance(model_io_context, model_clock)
      end select
    end if

  end subroutine initialise_io
#endif


  !> @brief Performs a single timestep of the fake nl model
  !>
  !> @param [inout] model_data  The model data instance
  subroutine step( model_data )

    implicit none

    type(model_data_type), intent(inout) :: model_data

    type(field_type),       pointer :: working_field    => null()
#ifdef USE_XIOS
    class(io_context_type), pointer :: model_io_context => null()
#endif

    call model_data%depository%get_field( test_field, working_field )

#ifdef USE_XIOS
    if ( use_xios_io ) then
      ! Switch to the main model I/O context
      model_io_context => get_io_context()
      call model_io_context%set_current()
    end if
#endif

    call jedi_lfric_increment_alg( working_field )

#ifdef USE_XIOS
    if ( write_data ) then
      call write_state( model_data%depository, prefix="write_" )
    end if
#endif

  end subroutine step


  !> @brief Initialise the model IO
  !>
  !> @param [in]    program_name     The program name
  subroutine finalise( program_name )

    implicit none

    character(len=*),            intent(in)    :: program_name

#ifdef USE_XIOS
    if ( use_xios_io ) then
      call final_io()
    end if
#endif

    call final_fem()

  end subroutine finalise


end module jedi_lfric_fake_nl_driver_mod
