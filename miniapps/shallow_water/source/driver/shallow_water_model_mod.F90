!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Handles initialisation and finalisation of infrastructure, constants
!!        and the shallow_water model simulation.
module shallow_water_model_mod

  use assign_orography_field_mod,     only: assign_orography_field
  use base_mesh_config_mod,           only: prime_mesh_name
  use calendar_mod,                   only: calendar_type
  use check_configuration_mod,        only: get_required_stencil_depth
  use checksum_alg_mod,               only: checksum_alg
  use conservation_algorithm_mod,     only: conservation_algorithm
  use constants_mod,                  only: i_def, str_def, &
                                            PRECISION_REAL, l_def
  use convert_to_upper_mod,           only: convert_to_upper
  use derived_config_mod,             only: set_derived_config
  use driver_fem_mod,                 only: init_fem, final_fem
  use driver_io_mod,                  only: init_io, final_io, &
                                            filelist_populator
  use driver_mesh_mod,                only: init_mesh
  use field_mod,                      only: field_type
  use field_parent_mod,               only: write_interface
  use field_collection_mod,           only: field_collection_type
  use geometric_constants_mod,        only: get_chi_inventory, &
                                            get_panel_id_inventory
  use inventory_by_mesh_mod,          only: inventory_by_mesh_type
  use io_config_mod,                  only: use_xios_io,             &
                                            write_conservation_diag, &
                                            write_dump,              &
                                            write_minmax_tseries
  use lfric_xios_file_mod,            only: lfric_xios_file_type
  use linked_list_mod,                only: linked_list_type
  use log_mod,                        only: log_event,          &
                                            log_set_level,      &
                                            log_scratch_space,  &
                                            LOG_LEVEL_INFO
  use mesh_collection_mod,            only: mesh_collection
  use mesh_mod,                       only: mesh_type
  use minmax_tseries_mod,             only: minmax_tseries,      &
                                            minmax_tseries_init, &
                                            minmax_tseries_final
  use model_clock_mod,                only: model_clock_type
  use mpi_mod,                        only: mpi_type
  use runtime_constants_mod,          only: create_runtime_constants
  use shallow_water_model_data_mod,   only: model_data_type
  use shallow_water_setup_io_mod,     only: init_shallow_water_files
  use timestepping_config_mod,        only: dt, &
                                            spinup_period
  use xios,                           only: xios_update_calendar

  implicit none

  private
  public :: initialise_infrastructure, &
            initialise_model,          &
            finalise_infrastructure,   &
            finalise_model

  contains

  !=============================================================================
  !> @brief Initialises the infrastructure and sets up constants used by the model.
  !> @param[in]     program_name An identifier given to the model begin run
  !> @param [out]   model_clock  Time within the model
  !> @param [in]    mpi          Communication object
  subroutine initialise_infrastructure(program_name, &
                                       model_clock,  &
                                       calendar,     &
                                       mpi )

    implicit none

    character(*),           intent(in)    :: program_name
    type(model_clock_type), intent(inout) :: model_clock
    class(calendar_type),   intent(in)    :: calendar
    class(mpi_type),        intent(inout) :: mpi

    type(inventory_by_mesh_type),  pointer :: chi_inventory => null()
    type(inventory_by_mesh_type),  pointer :: panel_id_inventory => null()
    procedure(filelist_populator), pointer :: files_init_ptr => null()

    character(len=*),   parameter   :: io_context_name = "shallow_water"
    character(str_def), allocatable :: base_mesh_names(:)
    logical(l_def)                  :: create_rdef_div_operators

    !-------------------------------------------------------------------------
    ! Initialise aspects of the infrastructure
    !-------------------------------------------------------------------------

    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    call set_derived_config( .true. )


    !-------------------------------------------------------------------------
    ! Work out which meshes are required
    !-------------------------------------------------------------------------

    ! Always have just one base mesh with shifted version
    allocate(base_mesh_names(1))
    base_mesh_names(1) = prime_mesh_name

    !-------------------------------------------------------------------------
    ! Initialise aspects of the grid
    !-------------------------------------------------------------------------

    ! TODO Stencil depth needs to be taken from configuration options
    ! Create the mesh
    call init_mesh( mpi%get_comm_rank(),  mpi%get_comm_size(), &
                    base_mesh_names,                                         &
                    required_stencil_depth = get_required_stencil_depth() )

    ! Create FEM specifics (function spaces and chi field)
    chi_inventory => get_chi_inventory()
    panel_id_inventory => get_panel_id_inventory()
    call init_fem(mesh_collection, chi_inventory, panel_id_inventory)

    !-------------------------------------------------------------------------
    ! Initialise aspects of output
    !-------------------------------------------------------------------------

    ! If using XIOS for diagnostic output or checkpointing, then set up XIOS
    ! domain and context

    files_init_ptr => init_shallow_water_files
    call init_io( io_context_name, mpi%get_comm(),   &
                  chi_inventory, panel_id_inventory, &
                  model_clock, calendar,             &
                  populate_filelist=files_init_ptr )

    !-------------------------------------------------------------------------
    ! Setup constants
    !-------------------------------------------------------------------------

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    create_rdef_div_operators = .true.
    call create_runtime_constants(mesh_collection, chi_inventory,  &
                                  panel_id_inventory, model_clock, &
                                  create_rdef_div_operators)

    deallocate(base_mesh_names)
    nullify(chi_inventory, panel_id_inventory)

  end subroutine initialise_infrastructure

  !=============================================================================
  !> @brief Initialises the shallow_water application.
  !> @param[in]     mesh       The primary mesh
  !> @param[in,out] model_data The working data set for the model run
  subroutine initialise_model( mesh,    &
                               model_data )

    use swe_timestep_alg_mod, only: swe_timestep_alg_init

    implicit none

    type(mesh_type),       pointer, intent(in)    :: mesh
    type(model_data_type), target,  intent(inout) :: model_data

    type(field_collection_type), pointer :: prognostic_fields => null()
    type(field_collection_type), pointer :: diagnostic_fields => null()

    ! Prognostic fields
    type(field_type), pointer :: wind => null()
    type(field_type), pointer :: geopot => null()
    type(field_type), pointer :: buoyancy => null()
    type(field_type), pointer :: q => null()

    call log_event( 'shallow_water: Initialising miniapp ...', LOG_LEVEL_INFO )

    prognostic_fields => model_data%prognostic_fields
    diagnostic_fields => model_data%diagnostic_fields

    call prognostic_fields%get_field("wind", wind)
    call prognostic_fields%get_field("buoyancy", buoyancy)
    call prognostic_fields%get_field("geopot", geopot)
    call prognostic_fields%get_field("q", q)

    ! Initialise transport and shallow water model
    call swe_timestep_alg_init( mesh, wind, geopot, buoyancy, q )

    call log_event( 'shallow_water: Miniapp initialised', LOG_LEVEL_INFO )

  end subroutine initialise_model


  !=============================================================================
  !> @brief Finalises infrastructure and constants used by the model.
  subroutine finalise_infrastructure()

    implicit none

    !-------------------------------------------------------------------------
    ! Finalise aspects of the grid
    !-------------------------------------------------------------------------
    call final_fem()

  end subroutine finalise_infrastructure

  !=============================================================================
  !> @brief Finalise the shallow_water application.
  !> @param[in,out] model_data   The working data set for the model run
  !> @param[in]     program_name An identifier given to the model begin run
  subroutine finalise_model( model_data, &
                             program_name )

    use swe_timestep_alg_mod, only: swe_timestep_alg_final

    implicit none

    type(model_data_type), target, intent(inout) :: model_data
    character(*),                  intent(in)    :: program_name

    type( field_collection_type ), pointer :: prognostic_fields => null()

    type(field_type), pointer :: wind => null()
    type(field_type), pointer :: geopot => null()
    type(field_type), pointer :: q => null()

    ! Pointer for setting I/O handlers on fields
    procedure(write_interface), pointer :: tmp_write_ptr => null()

    prognostic_fields => model_data%prognostic_fields

    call prognostic_fields%get_field('wind', wind)
    call prognostic_fields%get_field('geopot', geopot)
    call prognostic_fields%get_field('q', q)

    ! Checksums
    call checksum_alg( program_name, wind, 'wind', geopot, 'geopot', q, 'q' )

    ! Finalise transport
    call swe_timestep_alg_final()

  end subroutine finalise_model

end module shallow_water_model_mod
