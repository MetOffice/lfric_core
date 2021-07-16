!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Initialises model, Controls runtime constants, FEM and mesh related
!>        information used by the model

module init_multires_coupling_model_mod

  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type
  use create_mesh_mod,                only : init_mesh, final_mesh
  use create_fem_mod,                 only : init_fem, final_fem
  use runtime_constants_mod,          only : create_runtime_constants, &
                                             final_runtime_constants
  use fs_continuity_mod,              only : W3
  use local_mesh_collection_mod,      only : local_mesh_collection, &
                                             local_mesh_collection_type
  use mesh_collection_mod,            only : mesh_collection, &
                                             mesh_collection_type
  use mpi_mod,                        only : get_comm_size, get_comm_rank
  use multires_coupling_config_mod,   only : multires_coupling_mesh_tags
  use formulation_config_mod,         only : l_multigrid, &
                                             use_multires_coupling

  implicit none

  private
  public initialise_multires_coupling_model, final_multires_coupling_model

contains

  !> @brief Initialises runtime constants, FEM and meshes used by the model
  !> @param [in,out] prime_mesh_id The identifier of the primary mesh
  !> @param [in,out] prime_2D_mesh_id The identifier of the primary 2d mesh
  !> @param [in,out] prime_shifted_mesh_id The identifier of the shifted primary mesh
  !> @param [in,out] prime_double_level_mesh_id The identifier of the primary double
  !>                 level mesh
  !> @param [in,out] multires_coupling_mesh_ids A list of identifiers for all meshes in the
  !>                 multires_coupling Miniapp
  !> @param [in,out] multires_coupling_2D_mesh_ids A list of identifiers for all 2d
  !>                 meshes in the multires_coupling Miniapp
  !> @param [in,out] multigrid_mesh_ids A list of identifiers for all meshes used for
  !>                 multigrid
  !> @param [in,out] multigrid_2D_mesh_ids A list of identifiers for all 2d
  !>                 meshes used for multigrid
  !> @param [in,out] prime_chi_xyz (X,Y,Z) coordinate fields
  !> @param [in,out] prime_chi_sph Spherical coordinate fields
  !> @param [in,out] prime_panel_id 2D field giving cubed sphere panel ids
  !> @param [in,out] prime_shifted_chi_xyz (X,Y,Z) coordinate fields on shifted mesh
  !> @param [in,out] prime_shifted_chi_sph Spherical coordinate fields on shifted mesh
  !> @param [in,out] prime_double_level_chi_xyz (X,Y,Z) coordinate fields on
  !>                 double level mesh
  !> @param [in,out] prime_double_level_chi_sph Spherical coordinate fields on
  !>                 double level mesh
  !> @param [in,out] chi_xyz_fields All (X,Y,Z) coordinate fields for the multires_coupling Miniapp
  !> @param [in,out] chi_sph_fields All spherical coordinate fields for the multires_coupling Miniapp
  !> @param [in,out] panel_id_fields 2D field giving cubed sphere panel id fields
  !>                 for multires_coupling Miniapp
  !> @param [in,out] chi_mg_sph Spherical coordinate fields for multigrid
  !> @param [in,out] panel_id_mg 2D fields giving cubed sphere panel ids for multigrid
  subroutine initialise_multires_coupling_model(                                     &
                                       prime_mesh_id, prime_2D_mesh_id,              &
                                       prime_shifted_mesh_id,                        &
                                       prime_double_level_mesh_id,                   &
                                       multires_coupling_mesh_ids,                   &
                                       multires_coupling_2D_mesh_ids,                &
                                       multigrid_mesh_ids, multigrid_2D_mesh_ids,    &
                                       prime_chi_xyz, prime_chi_sph, prime_panel_id, &
                                       prime_shifted_chi_xyz,                        &
                                       prime_shifted_chi_sph,                        &
                                       prime_double_level_chi_xyz,                   &
                                       prime_double_level_chi_sph,                   &
                                       chi_xyz_fields, chi_sph_fields,               &
                                       panel_id_fields, chi_mg_sph, panel_id_mg )

    implicit none

    integer(kind=i_def),              intent(inout) :: prime_mesh_id, &
                                                       prime_2D_mesh_id
    integer(kind=i_def),              intent(inout) :: prime_shifted_mesh_id, &
                                                       prime_double_level_mesh_id
    integer(kind=i_def), allocatable, intent(inout) :: multires_coupling_mesh_ids(:)
    integer(kind=i_def), allocatable, intent(inout) :: multires_coupling_2D_mesh_ids(:)
    integer(kind=i_def), allocatable, intent(inout) :: multigrid_mesh_ids(:)
    integer(kind=i_def), allocatable, intent(inout) :: multigrid_2D_mesh_ids(:)

    type(field_type),              intent(inout) :: prime_chi_xyz(3)
    type(field_type),              intent(inout) :: prime_chi_sph(3)
    type(field_type),              intent(inout) :: prime_shifted_chi_xyz(3)
    type(field_type),              intent(inout) :: prime_shifted_chi_sph(3)
    type(field_type),              intent(inout) :: prime_double_level_chi_xyz(3)
    type(field_type),              intent(inout) :: prime_double_level_chi_sph(3)
    type(field_type),              intent(inout) :: prime_panel_id
    type(field_type), allocatable, intent(inout) :: chi_xyz_fields(:,:)
    type(field_type), allocatable, intent(inout) :: chi_sph_fields(:,:)
    type(field_type), allocatable, intent(inout) :: panel_id_fields(:)
    type(field_type), allocatable, intent(inout) :: chi_mg_sph(:,:)
    type(field_type), allocatable, intent(inout) :: panel_id_mg(:)

    integer(kind=i_def) :: total_ranks, local_rank

    allocate( local_mesh_collection, &
              source = local_mesh_collection_type() )

    allocate( mesh_collection, &
              source=mesh_collection_type() )

    ! Get the rank information from the virtual machine
    total_ranks = get_comm_size()
    local_rank  = get_comm_rank()

    ! Create the meshes
    call init_mesh( local_rank, total_ranks,         &
                    prime_mesh_id, prime_2D_mesh_id, &
                    prime_shifted_mesh_id,           &
                    prime_double_level_mesh_id,      &
                    multigrid_mesh_ids,              &
                    multigrid_2D_mesh_ids,           &
                    l_multigrid,                     &
                    multires_coupling_mesh_ids,      &
                    multires_coupling_2D_mesh_ids,   &
                    multires_coupling_mesh_tags,     &
                    use_multires_coupling )

    ! Create FEM specifics (function spaces and chi field)
    call init_fem( prime_mesh_id, prime_chi_xyz, prime_chi_sph,            &
                   prime_panel_id, prime_shifted_mesh_id,                  &
                   prime_shifted_chi_xyz, prime_shifted_chi_sph,           &
                   prime_double_level_mesh_id, prime_double_level_chi_xyz, &
                   prime_double_level_chi_sph,                             &
                   multigrid_mesh_ids, multigrid_2D_mesh_ids,              &
                   chi_mg_sph, panel_id_mg,                                &
                   l_multigrid,                                            &
                   multires_coupling_mesh_ids,                             &
                   multires_coupling_2D_mesh_ids,                          &
                   chi_sph_fields, chi_xyz_fields, panel_id_fields,        &
                   use_multires_coupling )

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the fem algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_runtime_constants( prime_mesh_id, prime_2D_mesh_id,           &
                                   prime_chi_xyz, prime_chi_sph,              &
                                   prime_panel_id, prime_shifted_mesh_id,     &
                                   prime_shifted_chi_xyz,                     &
                                   prime_shifted_chi_sph,                     &
                                   prime_double_level_mesh_id,                &
                                   prime_double_level_chi_xyz,                &
                                   prime_double_level_chi_sph,                &
                                   multigrid_mesh_ids, multigrid_2D_mesh_ids, &
                                   chi_mg_sph, panel_id_mg,                   &
                                   multires_coupling_mesh_ids,                &
                                   multires_coupling_2D_mesh_ids,             &
                                   chi_sph_fields, panel_id_fields )

  end subroutine initialise_multires_coupling_model

  !> @brief Finalise FEM, runtime constants and meshes
  subroutine final_multires_coupling_model()

    implicit none

    call final_runtime_constants()
    call final_fem()
    call final_mesh()


  end subroutine final_multires_coupling_model


end module init_multires_coupling_model_mod
