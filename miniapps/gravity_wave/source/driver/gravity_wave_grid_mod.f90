!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Controls grid related information used by the model

module gravity_wave_grid_mod

  use base_mesh_config_mod,           only : prime_mesh_name
  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type
  use create_mesh_mod,                only : init_mesh
  use create_fem_mod,                 only : init_fem
  use fs_continuity_mod,              only : W3
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection_type, &
                                             function_space_collection
  use function_space_chain_mod,       only : function_space_chain_type
  use local_mesh_collection_mod,      only : local_mesh_collection, &
                                             local_mesh_collection_type
  use log_mod,                        only : log_event,          &
                                             log_scratch_space,  &
                                             LOG_LEVEL_INFO
  use formulation_config_mod,         only : l_multigrid
  use mesh_collection_mod,            only : mesh_collection, &
                                             mesh_collection_type
  use mpi_mod,                        only : get_comm_size, get_comm_rank

  implicit none

  private
  public initialise_grid

contains

  !> @brief Initialises grid related information used by the model
  !> @param [in,out] mesh_id The identifier of the primary mesh
  !> @param [in,out] twod_mesh_id The identifier of the primary 2d mesh
  !> @param [in,out] chi_xyz (X,Y,Z) coordinate fields
  !> @param [in,out] chi_sph Spherical coordinate fields
  !> @param [in,out] panel_id 2D field giving cubed sphere panel ids
  !> @param [in,out] multigrid_function_space_chain An ordered  list of function
  !>                                               spaces used by multigrid
  subroutine initialise_grid(mesh_id, twod_mesh_id,      &
                             chi_xyz, chi_sph, panel_id, &
                             multigrid_function_space_chain)

    implicit none

    integer(i_def),   intent(inout) :: mesh_id, twod_mesh_id
    type(field_type), intent(inout) :: chi_xyz(3)
    type(field_type), intent(inout) :: chi_sph(3)
    type(field_type), intent(inout) :: panel_id
    type(function_space_chain_type), intent(inout) :: &
                                       multigrid_function_space_chain

    type(function_space_type), pointer :: function_space

    integer(i_def)     :: total_ranks, local_rank
    integer(i_def)     :: i

    integer(i_def), allocatable :: multigrid_mesh_ids(:)

    allocate( local_mesh_collection, &
              source = local_mesh_collection_type() )
    allocate( mesh_collection, &
              source=mesh_collection_type() )

    ! Get the rank information from the virtual machine
    total_ranks = get_comm_size()
    local_rank  = get_comm_rank()

    ! Create the mesh
    call init_mesh( local_rank, total_ranks, mesh_id,     &
                    twod_mesh_id=twod_mesh_id,            &
                    multigrid_mesh_ids=multigrid_mesh_ids,&
                    use_multigrid=l_multigrid )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem(mesh_id, chi_xyz, chi_sph, panel_id)

    if (l_multigrid) then

      do i=1, size(multigrid_mesh_ids)

        ! Make sure this function_space is in the collection
        function_space => &
            function_space_collection%get_fs( multigrid_mesh_ids(i), 0, W3 )

        write( log_scratch_space,"(A,I0,A)")                       &
             'Adding function_space id ', function_space%get_id(), &
             ' to multigrid function_space chain'
        call log_event( log_scratch_space, LOG_LEVEL_INFO )

        call multigrid_function_space_chain%add( function_space )

      end do
    end if

    nullify( function_space )

  end subroutine initialise_grid

end module gravity_wave_grid_mod
