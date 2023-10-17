!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the simple_diffusion miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module simple_diffusion_driver_mod
  use base_mesh_config_mod,       only : prime_mesh_name
  use calendar_mod,               only : calendar_type
  use checksum_alg_mod,           only : checksum_alg
  use constants_mod,              only : i_def, str_def, &
                                         r_def, r_second
  use convert_to_upper_mod,       only : convert_to_upper
  use driver_mesh_mod,            only : init_mesh
  use driver_modeldb_mod,         only : modeldb_type
  use driver_fem_mod,             only : init_fem, final_fem
  use driver_io_mod,              only : init_io, final_io
  use field_collection_mod,       only : field_collection_type
  use field_mod,                  only : field_type
  use init_simple_diffusion_mod,          only : init_simple_diffusion
  use inventory_by_mesh_mod,      only : inventory_by_mesh_type
  use io_config_mod,              only : write_diag
  use log_mod,                    only : log_event, log_scratch_space, &
                                         LOG_LEVEL_ALWAYS, LOG_LEVEL_INFO, &
                                         LOG_LEVEL_TRACE
  use mesh_mod,                   only : mesh_type
  use mesh_collection_mod,        only : mesh_collection
  use model_clock_mod,            only : model_clock_type
  use mpi_mod,                    only : mpi_type
  use simple_diffusion_alg_mod,   only : simple_diffusion_alg

  implicit none

  private

  public initialise, step, finalise

contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !> @param [in]     program_name An identifier given to the model being run
  !> @param [in,out] modeldb      The structure that holds model state
  !> @param [in]     calendar     The model calendar
    subroutine initialise( program_name, modeldb, calendar )

    implicit none

    character(*),            intent(in)    :: program_name
    type(modeldb_type),      intent(inout) :: modeldb
    class(calendar_type),    intent(in)    :: calendar

    ! Coordinate field
    type(field_type),             pointer :: chi(:) => null()
    type(field_type),             pointer :: panel_id => null()
    type(mesh_type),              pointer :: mesh => null()
    type(inventory_by_mesh_type)          :: chi_inventory
    type(inventory_by_mesh_type)          :: panel_id_inventory
    character(str_def),       allocatable :: base_mesh_names(:)

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------
    ! Create the mesh
    allocate(base_mesh_names(1))
    base_mesh_names(1) = prime_mesh_name
    call init_mesh( modeldb%mpi%get_comm_rank(), &
                    modeldb%mpi%get_comm_size(), &
                    base_mesh_names )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )

    ! Initialise I/O context
    call init_io( program_name, modeldb%mpi%get_comm(), chi_inventory, &
                  panel_id_inventory, modeldb%clock, calendar )

    ! Create and initialise prognostic fields
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    call chi_inventory%get_field_array(mesh, chi)
    call panel_id_inventory%get_field(mesh, panel_id)
    call init_simple_diffusion( mesh, chi, panel_id, modeldb )

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
    type( field_collection_type ), pointer :: depository
    type( field_type ),            pointer :: diffusion_field

    depository => modeldb%fields%get_field_collection("depository")
    call depository%get_field("diffusion_field", diffusion_field)

    ! Call an algorithm
    call log_event(program_name//": Calculating diffusion", LOG_LEVEL_INFO)
    call simple_diffusion_alg(diffusion_field)

    if (write_diag ) then
        ! Write out output file
        call log_event(program_name//": Writing diagnostic output", LOG_LEVEL_INFO)
        call diffusion_field%write_field('diffusion_field')
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
    type( field_collection_type ), pointer :: depository
    type( field_type ),            pointer :: diffusion_field

    depository => modeldb%fields%get_field_collection("depository")
    call depository%get_field("diffusion_field", diffusion_field)

    !--------------------------------------------------------------------------
    ! Model finalise
    !--------------------------------------------------------------------------
    ! Write checksums to file
    call checksum_alg(program_name, diffusion_field, 'simple_diffusion_field_1')
    call log_event( program_name//': Miniapp completed', LOG_LEVEL_TRACE )

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------
    ! Finalise IO
    call final_io()
    call final_fem()

  end subroutine finalise

end module simple_diffusion_driver_mod
