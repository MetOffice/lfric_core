!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the lbc_demo application.
!>
!>
module lbc_demo_driver_mod

  use add_mesh_map_mod,         only: assign_mesh_maps
  use sci_checksum_alg_mod,     only: checksum_alg
  use constants_mod,            only: i_def, str_def, &
                                      r_def
  use create_mesh_mod,          only: create_mesh, create_extrusion
  use driver_mesh_mod,          only: init_mesh
  use driver_modeldb_mod,       only: modeldb_type
  use driver_fem_mod,           only: init_fem, final_fem
  use driver_io_mod,            only: init_io, final_io
  use extrusion_mod,            only: extrusion_type,         &
                                      uniform_extrusion_type, &
                                      twod, prime_extrusion
  use field_collection_mod,     only: field_collection_type
  use field_mod,                only: field_type
  use apply_lbc_fields_alg_mod, only: apply_lbc_fields

  use inventory_by_mesh_mod,    only: inventory_by_mesh_type
  use log_mod,                  only: log_event,         &
                                      log_scratch_space, &
                                      log_level_info,    &
                                      log_level_error,   &
                                      log_level_trace
  use mesh_mod,                 only: mesh_type
  use mesh_collection_mod,      only: mesh_collection
  use namelist_mod,             only: namelist_type

  !------------------------------------
  ! Configuration modules
  !------------------------------------
  use base_mesh_config_mod, only: geometry_spherical, &
                                  geometry_planar,    &
                                  topology_non_periodic
  use lbc_demo_config_mod,  only: lbc_source_file, &
                                  field_type_real
  implicit none

  private

  type(inventory_by_mesh_type) :: chi_inventory
  type(inventory_by_mesh_type) :: panel_id_inventory

  public initialise, step, finalise

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Sets up required state in preparation for run.
!> @param [in]     program_name An identifier given to the model being run
!> @param [in,out] modeldb      The structure that holds model state
subroutine initialise( program_name, modeldb)

  use init_lam_fields_alg_mod, only: init_lam_fields
  use init_lbc_fields_alg_mod, only: init_lbc_fields

  implicit none

  character(*),         intent(in)    :: program_name
  type(modeldb_type),   intent(inout) :: modeldb

  ! Coordinate fields
  type(field_type), pointer :: chi(:)
  type(field_type), pointer :: panel_id
  type(mesh_type),  pointer :: mesh

  character(str_def), allocatable :: base_mesh_names(:)
  character(str_def), allocatable :: twod_names(:)

  class(extrusion_type),        allocatable :: extrusion
  type(uniform_extrusion_type), allocatable :: extrusion_2d

  type(namelist_type), pointer :: base_mesh_nml
  type(namelist_type), pointer :: planet_nml
  type(namelist_type), pointer :: extrusion_nml
  type(namelist_type), pointer :: io_nml
  type(namelist_type), pointer :: lbc_demo_nml

  character(str_def) :: prime_mesh_name
  character(str_def) :: lbc_mesh_name
  character(str_def) :: output_mesh_name
  character(str_def) :: context_name

  integer(i_def) :: stencil_depth
  integer(i_def) :: geometry
  integer(i_def) :: method
  integer(i_def) :: number_of_layers
  real(r_def)    :: domain_bottom
  real(r_def)    :: domain_height
  real(r_def)    :: scaled_radius
  logical        :: check_partitions

  logical :: write_diag

  logical :: enable_lbc
  logical :: apply_lbc
  logical :: write_lbc
  integer :: lbc_source
  integer :: topology

  integer :: i

  integer(i_def), parameter :: one_layer = 1_i_def

  !=======================================================================
  ! Extract configuration variables
  !=======================================================================
  base_mesh_nml => modeldb%configuration%get_namelist('base_mesh')
  planet_nml    => modeldb%configuration%get_namelist('planet')
  extrusion_nml => modeldb%configuration%get_namelist('extrusion')
  io_nml        => modeldb%configuration%get_namelist('io')
  lbc_demo_nml  => modeldb%configuration%get_namelist('lbc_demo')

  call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
  call base_mesh_nml%get_value( 'geometry', geometry )
  call base_mesh_nml%get_value( 'topology', topology )

  call extrusion_nml%get_value( 'method', method )
  call extrusion_nml%get_value( 'domain_height', domain_height )
  call extrusion_nml%get_value( 'number_of_layers', number_of_layers )
  call planet_nml%get_value( 'scaled_radius', scaled_radius )
  call io_nml%get_value( 'write_diag', write_diag)

  call lbc_demo_nml%get_value( 'enable_lbc', enable_lbc )
  call lbc_demo_nml%get_value( 'apply_lbc',  apply_lbc )
  call lbc_demo_nml%get_value( 'write_lbc',  write_lbc )
  call lbc_demo_nml%get_value( 'lbc_source', lbc_source )

  !=======================================================================
  ! Mesh setup
  !=======================================================================

  ! Two meshes required, the LAM and it's LBC counterpart. Construct
  ! the mesh names required.

  if (enable_lbc) then

    ! In order for LBC fields be enabled, the LBC mesh
    ! associated with the model LAM mesh must also
    ! be loaded from file.

    ! Perform check to see if the model mesh is suitable
    if (topology /= topology_non_periodic) then
      write(log_scratch_space,'(A)')                        &
          'Invalid configuration: Only non-periodic LAM '// &
          'domains currently suppported for LBC driving data.'
      call log_event(log_scratch_space, log_level_error)
    end if

    allocate(base_mesh_names(2))
    lbc_mesh_name = trim(prime_mesh_name)//'-lbc'
    base_mesh_names(2) = lbc_mesh_name

  else

    allocate(base_mesh_names(1))

  end if

  base_mesh_names(1) = prime_mesh_name

  !-------------------------------------------------------------------------
  ! Generate required extrusions
  ! Full depth and 2D extrusions of both are required test potential use cases.
  !-------------------------------------------------------------------------
  select case (geometry)

  case (geometry_planar)
    domain_bottom = 0.0_r_def

  case (geometry_spherical)
    domain_bottom = scaled_radius

  case default
    call log_event("Invalid geometry for mesh initialisation", &
                   log_level_error)
  end select

  allocate( extrusion, source=create_extrusion( method,           &
                                                domain_height,    &
                                                domain_bottom,    &
                                                number_of_layers, &
                                                prime_extrusion ) )

  extrusion_2d = uniform_extrusion_type( domain_bottom, &
                                         domain_bottom, &
                                         one_layer, twod )

  !-------------------------------------------------------------------------
  ! Initialise mesh objects and assign InterGrid maps
  !-------------------------------------------------------------------------
  stencil_depth = 1
  check_partitions = .false.

  call init_mesh( modeldb%configuration,       &
                  modeldb%mpi%get_comm_rank(), &
                  modeldb%mpi%get_comm_size(), &
                  base_mesh_names, extrusion,  &
                  stencil_depth, check_partitions )

  allocate( twod_names, source=base_mesh_names )

  do i=1, size(twod_names)
    twod_names(i) = trim(twod_names(i))//'_2d'
  end do

  call create_mesh( base_mesh_names, extrusion_2d, alt_name=twod_names )
  call assign_mesh_maps( twod_names )


  !=======================================================================
  ! Build the FEM function spaces and coordinate fields
  !=======================================================================
  call init_fem( mesh_collection, chi_inventory, panel_id_inventory )

  !=======================================================================
  ! Setup general I/O system.
  !=======================================================================
  ! Initialise I/O context
  if (write_diag) then

    if (enable_lbc .and. write_lbc) then
      ! Output LBCs
      context_name     = trim(program_name)//':lbcs'
      output_mesh_name = trim(prime_mesh_name)//'-lbc'
    else
      ! Diagnostics
      context_name     = trim(program_name)//':diags'
      output_mesh_name = trim(prime_mesh_name)
    end if

    write(log_scratch_space,'(A)')                     &
        'Setting up context ' // trim(context_name) // &
        ' for ' // trim(output_mesh_name)
    call log_event(log_scratch_space, log_level_info)

    call init_io( context_name, output_mesh_name, modeldb, &
                  chi_inventory, panel_id_inventory )
  end if

  !=======================================================================
  ! Create and initialise prognostic fields
  !=======================================================================
  mesh => mesh_collection%get_mesh(prime_mesh_name)
  call init_lam_fields( mesh, modeldb )

  !=======================================================================
  ! Create and initialise LBC fields
  !=======================================================================
  if (enable_lbc .and. apply_lbc) then
    call init_lbc_fields( mesh, chi_inventory, modeldb )
  end if

  nullify(mesh, chi, panel_id)
  deallocate(base_mesh_names)

end subroutine initialise

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Performs a time step.
!> @param [in]     program_name An identifier given to the model being run
!> @param [in,out] modeldb      The structure that holds model state
subroutine step( program_name, modeldb )

  use write_field_set_mod,    only: write_field_set
  use lfric_xios_context_mod, only: lfric_xios_context_type

  implicit none

  character(*),       intent(in)    :: program_name
  type(modeldb_type), intent(inout) :: modeldb

  type( field_collection_type ), pointer :: output_diags
  type(namelist_type), pointer :: io_nml
  type(namelist_type), pointer :: lbc_demo_nml

  logical :: apply_lbc, write_diag, write_lbc, enable_lbc
  character(str_def) :: suffix

  io_nml       => modeldb%configuration%get_namelist('io')
  lbc_demo_nml => modeldb%configuration%get_namelist('lbc_demo')

  call io_nml%get_value( 'write_diag',   write_diag)

  call lbc_demo_nml%get_value( 'apply_lbc',  apply_lbc )
  call lbc_demo_nml%get_value( 'enable_lbc', enable_lbc )
  call lbc_demo_nml%get_value( 'write_lbc',  write_lbc )

  if (apply_lbc) then
    ! Update prognostic with LBC fields
    call apply_lbc_fields(modeldb)
  end if

  if (write_diag) then
    ! Write out output file

    if (enable_lbc .and. write_lbc) then
      output_diags => modeldb%fields%get_field_collection('lbc')
      suffix=':lbc'
    else
      output_diags => modeldb%fields%get_field_collection('depository')
      suffix=''
    end if

    call write_field_set(output_diags, suffix)

  end if

end subroutine step


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Tidies up after a run.
!> @param [in]     program_name An identifier given to the model being run
!> @param [in,out] modeldb      The structure that holds model state
subroutine finalise( program_name, modeldb )

  implicit none

  character(*),       intent(in)    :: program_name
  type(modeldb_type), intent(inout) :: modeldb

  type(field_collection_type), pointer :: depository
  type(namelist_type),         pointer :: lbc_demo_nml

  integer(i_def) :: lbc_field_type

  lbc_demo_nml => modeldb%configuration%get_namelist('lbc_demo')
  call lbc_demo_nml%get_value( 'field_type',  lbc_field_type )

  !-------------------------------------------------------------------------
  ! Checksum output - Only for real fields
  !-------------------------------------------------------------------------
  if (lbc_field_type == field_type_real) then
    depository => modeldb%fields%get_field_collection('depository')
    call checksum_alg( program_name, field_collection=depository )
  end if

  call log_event( program_name//': model completed', log_level_trace )

  !-------------------------------------------------------------------------
  ! Driver layer finalise
  !-------------------------------------------------------------------------
  ! Finalise IO
  call final_io(modeldb)
  call final_fem()

end subroutine finalise

end module lbc_demo_driver_mod
