!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the coupled miniapp.
!>
module coupled_driver_mod

  use add_mesh_map_mod,           only : assign_mesh_maps
  use base_mesh_config_mod,       only : GEOMETRY_SPHERICAL, &
                                         GEOMETRY_PLANAR
  use calendar_mod,               only : calendar_type
  use constants_mod,              only : i_def, str_def, str_longlong, &
                                         r_def, r_second
  use convert_to_upper_mod,       only : convert_to_upper
  use coupled_alg_mod,            only : coupled_alg
  use create_mesh_mod,            only : create_extrusion, create_mesh
  use driver_mesh_mod,            only : init_mesh
  use driver_modeldb_mod,         only : modeldb_type
  use driver_fem_mod,             only : init_fem, final_fem
  use extrusion_mod,              only : extrusion_type,         &
                                         uniform_extrusion_type, &
                                         PRIME_EXTRUSION, TWOD
  use field_collection_mod,       only : field_collection_type
  use field_mod,                  only : field_type
  use init_coupled_mod,           only : init_coupled
  use inventory_by_mesh_mod,      only : inventory_by_mesh_type
  use log_mod,                    only : log_event, log_scratch_space, &
                                         LOG_LEVEL_ALWAYS,             &
                                         LOG_LEVEL_ERROR,              &
                                         LOG_LEVEL_INFO
  use mesh_mod,                   only : mesh_type
  use mesh_collection_mod,        only : mesh_collection

  use sci_checksum_alg_mod,       only : checksum_alg

  implicit none

  private
  public initialise, step, finalise

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !> @param [in]     program_name Identifier given to the model being run
  !> @param [in,out] modeldb      The structure that holds model state
  !> @param [in]     calendar     The model calendar
  subroutine initialise( program_name, modeldb, calendar )

    implicit none

    character(*),            intent(in)    :: program_name
    type(modeldb_type),      intent(inout) :: modeldb
    class(calendar_type),    intent(in)    :: calendar

    ! Coordinate field
    type(field_type),             pointer :: chi(:)
    type(field_type),             pointer :: panel_id
    type(mesh_type),              pointer :: mesh
    type(mesh_type),              pointer :: mesh_twod
    type(inventory_by_mesh_type)          :: chi_inventory
    type(inventory_by_mesh_type)          :: panel_id_inventory

    character(str_def),    allocatable :: base_mesh_names(:)
    character(str_def),    allocatable :: twod_names(:)

    class(extrusion_type),        allocatable :: extrusion
    type(uniform_extrusion_type), allocatable :: extrusion_2d

    character(str_def) :: prime_mesh_name

    integer(i_def) :: stencil_depth
    integer(i_def) :: geometry
    integer(i_def) :: method
    integer(i_def) :: number_of_layers
    real(r_def)    :: domain_bottom
    real(r_def)    :: domain_height
    real(r_def)    :: scaled_radius
    logical        :: apply_partition_check

    integer(i_def) :: i
    integer(i_def), parameter :: one_layer = 1_i_def


    ! Extract namelist variables
    prime_mesh_name  = modeldb%config%base_mesh%prime_mesh_name()
    geometry         = modeldb%config%base_mesh%geometry()
    method           = modeldb%config%extrusion%method()
    domain_height    = modeldb%config%extrusion%domain_height()
    number_of_layers = modeldb%config%extrusion%number_of_layers()
    scaled_radius    = modeldb%config%planet%scaled_radius()

    ! Initialise mesh
    ! Determine the required meshes
    allocate(base_mesh_names(1))
    base_mesh_names(1) = prime_mesh_name

    ! Create the required extrusions
    select case (geometry)
    case (GEOMETRY_PLANAR)
      domain_bottom = 0.0_r_def
    case (GEOMETRY_SPHERICAL)
      domain_bottom = scaled_radius
    case default
      call log_event("Invalid geometry for mesh initialisation", &
                      LOG_LEVEL_ERROR)
    end select
    allocate( extrusion, source=create_extrusion( method,           &
                                                  domain_height,    &
                                                  domain_bottom,    &
                                                  number_of_layers, &
                                                  PRIME_EXTRUSION ) )

    extrusion_2d = uniform_extrusion_type( domain_bottom, &
                                           domain_bottom, &
                                           one_layer, TWOD )

    ! Create the required meshes
    stencil_depth = 1
    apply_partition_check = .false.
    call init_mesh( modeldb%configuration,       &
                    modeldb%mpi%get_comm_rank(), &
                    modeldb%mpi%get_comm_size(), &
                    base_mesh_names, extrusion,  &
                    stencil_depth, apply_partition_check )

    allocate( twod_names, source=base_mesh_names )
    do i=1, size(twod_names)
      twod_names(i) = trim(twod_names(i))//'_2d'
    end do
    call create_mesh( base_mesh_names, extrusion_2d, &
                      alt_name=twod_names )
    call assign_mesh_maps( twod_names )


    ! Build the FEM function spaces and coordinate fields
    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )

    ! Create and initialise prognostic fields
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    mesh_twod => mesh_collection%get_mesh(mesh, TWOD)
    call chi_inventory%get_field_array(mesh, chi)
    call panel_id_inventory%get_field(mesh, panel_id)
    call init_coupled( mesh_twod, chi, panel_id, modeldb )

    nullify(mesh, chi, panel_id)
    deallocate(base_mesh_names)

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs a time step.
  !> @param [in]     program_name An identifier given to the model being run
  !> @param [in,out] modeldb      The structure that holds model state
  subroutine step( program_name, modeldb )

    implicit none

    character(*),       intent(in)    :: program_name
    type(modeldb_type), intent(inout) :: modeldb

    ! Call an algorithm
    call coupled_alg(modeldb)

    ! Write out output file
    call log_event(program_name//": Writing diagnostic output", LOG_LEVEL_INFO)

  end subroutine step

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !> @param [in]     program_name An identifier given to the model being run
  !> @param [in,out] modeldb      The structure that holds model state
  subroutine finalise( program_name, modeldb )

    implicit none

    character(*),       intent(in)    :: program_name
    type(modeldb_type), intent(inout) :: modeldb

    type( field_collection_type ), pointer :: depository
    type( field_type ),            pointer :: field_2

    character(str_longlong), pointer       :: cpl_component_name

    depository => modeldb%fields%get_field_collection("depository")
    call depository%get_field("field_2", field_2)

    ! Write checksum of coupled field (on incoming component) to file 
    call modeldb%values%get_value("cpl_name", cpl_component_name)
    if (trim(cpl_component_name) == "lfric_i") then
      call checksum_alg(program_name, field_2, 'coupled_field_2')
    end if

    call log_event( program_name//': Miniapp completed', LOG_LEVEL_INFO )

    call final_fem()

  end subroutine finalise

end module coupled_driver_mod
