!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Controls infrastructure related information used by the model

module gravity_wave_infrastructure_mod

  use base_mesh_config_mod,       only : prime_mesh_name
  use calendar_mod,               only : calendar_type
  use check_configuration_mod,    only : get_required_stencil_depth
  use constants_mod,              only : i_def,           &
                                         PRECISION_REAL,  &
                                         r_def, r_second, &
                                         l_def, str_def
  use convert_to_upper_mod,       only : convert_to_upper
  use derived_config_mod,         only : set_derived_config
  use geometric_constants_mod,    only : get_chi_inventory,  &
                                         get_panel_id_inventory
  use inventory_by_mesh_mod,      only : inventory_by_mesh_type
  use log_mod,                    only : log_event,          &
                                         log_scratch_space,  &
                                         LOG_LEVEL_ALWAYS,   &
                                         LOG_LEVEL_INFO
  use mesh_collection_mod,        only : mesh_collection
  use model_clock_mod,            only : model_clock_type
  use mpi_mod,                    only : mpi_type
  use field_mod,                  only : field_type
  use driver_fem_mod,             only : init_fem, init_function_space_chains
  use driver_io_mod,              only : init_io, final_io
  use driver_mesh_mod,            only : init_mesh
  use runtime_constants_mod,      only : create_runtime_constants
  use formulation_config_mod,     only : l_multigrid
  use multigrid_config_mod,       only : chain_mesh_tags

  implicit none

  private
  public initialise_infrastructure, finalise_infrastructure

  character(*), public, parameter   :: xios_ctx = 'gravity_wave'

contains

  !> @brief Initialises the infrastructure used by the model
  !> @param [out]    model_clock  Time within the model
  !> @param [in]     mpi          Communication object
  !> @param [in]     caledar      Interprets date strings.
  !>
  subroutine initialise_infrastructure( model_clock, mpi, calendar )

    use log_mod,                 only : log_level_error
    use step_calendar_mod,       only : step_calendar_type

    implicit none

    type(model_clock_type), intent(inout) :: model_clock
    class(mpi_type),        intent(inout) :: mpi
    class(calendar_type),   intent(in)    :: calendar

    type(inventory_by_mesh_type), pointer :: chi_inventory => null()
    type(inventory_by_mesh_type), pointer :: panel_id_inventory => null()

    logical(l_def)                  :: mesh_already_exists
    integer(i_def)                  :: i, j, mesh_ctr, max_num_meshes
    character(str_def), allocatable :: base_mesh_names(:)
    character(str_def), allocatable :: tmp_mesh_names(:)
    logical(l_def)                  :: create_rdef_div_operators

    write(log_scratch_space,'(A)')                        &
        'Application built with '//trim(PRECISION_REAL)// &
        '-bit real numbers'
    call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

    call set_derived_config( .false. )

    !-------------------------------------------------------------------------
    ! Work out which meshes are required
    !-------------------------------------------------------------------------

    ! Gather together names of meshes
    max_num_meshes = 1
    if ( l_multigrid ) max_num_meshes = max_num_meshes + SIZE(chain_mesh_tags)
    ! Don't know the full number of meshes, so allocate maximum for now
    allocate(tmp_mesh_names(max_num_meshes))

    ! Prime extrusions -------------------------------------------------------
    ! Always have the prime mesh
    mesh_ctr = 1
    tmp_mesh_names(1) = prime_mesh_name
    ! Multigrid meshes
    if ( l_multigrid ) then
      do i = 1, SIZE(chain_mesh_tags)
        ! Only add mesh if it has not already been added
        mesh_already_exists = .false.
        do j = 1, mesh_ctr
          if ( tmp_mesh_names(j) == chain_mesh_tags(i) ) then
            mesh_already_exists = .true.
            exit
          end if
        end do
        if ( .not. mesh_already_exists ) then
          mesh_ctr = mesh_ctr + 1
          tmp_mesh_names(mesh_ctr) = chain_mesh_tags(i)
        end if
      end do
    end if
    ! Transfer mesh names from temporary array to an array of appropriate size
    allocate(base_mesh_names(mesh_ctr))
    do i = 1, mesh_ctr
      base_mesh_names(i) = tmp_mesh_names(i)
    end do
    deallocate(tmp_mesh_names)

    !-------------------------------------------------------------------------
    ! Initialise aspects of the grid
    !-------------------------------------------------------------------------
    ! Create the mesh
    call init_mesh( mpi%get_comm_rank(), mpi%get_comm_size(),             &
                    base_mesh_names,                                      &
                    required_stencil_depth = get_required_stencil_depth() )

    ! Create FEM specifics (function spaces and chi field)
    chi_inventory => get_chi_inventory()
    panel_id_inventory => get_panel_id_inventory()
    call init_fem( mesh_collection, chi_inventory, panel_id_inventory )
    if ( l_multigrid ) then
      call init_function_space_chains( mesh_collection, chain_mesh_tags )
    end if

    !-------------------------------------------------------------------------
    ! Initialise aspects of output
    !-------------------------------------------------------------------------
    call init_io( xios_ctx, mpi%get_comm(), chi_inventory, panel_id_inventory, &
                  model_clock, calendar )

    !-------------------------------------------------------------------------
    ! Setup constants
    !-------------------------------------------------------------------------

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field and limited area masks.
    create_rdef_div_operators = .true.
    call create_runtime_constants( mesh_collection, chi_inventory,  &
                                   panel_id_inventory, model_clock, &
                                   create_rdef_div_operators )


    nullify(chi_inventory, panel_id_inventory)
    deallocate(base_mesh_names)

  end subroutine initialise_infrastructure


  !> @brief Finalises infrastructure used by the model
  subroutine finalise_infrastructure()

    implicit none

    ! Finalise I/O
    call final_io()

  end subroutine finalise_infrastructure

end module gravity_wave_infrastructure_mod
