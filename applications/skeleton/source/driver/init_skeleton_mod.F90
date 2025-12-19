!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Initialisation functionality for the skeleton miniapp

!> @details Handles init of prognostic fields and through the call to
!>          runtime_contants the coordinate fields and fem operators

module init_skeleton_mod

  use constants_mod,                 only: i_def
  use driver_modeldb_mod,            only: modeldb_type
  use field_collection_mod,          only: field_collection_type
  use field_mod,                     only: field_type
  use function_space_collection_mod, only: function_space_collection
  use function_space_mod,            only: function_space_type
  use fs_continuity_mod,             only: W3
  use log_mod,                       only: log_event, log_level_info
  use mesh_mod,                      only: mesh_type
  use skeleton_constants_mod,        only: create_skeleton_constants

  implicit none

  contains

  !> @details Initialises everything needed to run the skeleton miniapp
  !> @param[in,out] modeldb  The structure that holds model state
  !> @param[in]     mesh     Representation of the mesh the code will run on
  !> @param[in,out] chi      The co-ordinate field
  !> @param[in,out] panel_id 2d field giving the id for cubed sphere panels

  subroutine init_skeleton(modeldb, mesh, chi, panel_id)

    implicit none

    type(modeldb_type), target, intent(inout) :: modeldb
    type(mesh_type),   pointer, intent(in)    :: mesh

    ! Coordinate field
    type(field_type), intent(inout) :: chi(:)
    type(field_type), intent(inout) :: panel_id

    type(field_type) :: field_1

    type(field_collection_type),   pointer :: depository
    type(function_space_type),     pointer :: fs

    integer(i_def) :: order_h, order_v

    call log_event('skeleton: Initialising miniapp ...', log_level_info)

    depository => modeldb%fields%get_field_collection("depository")

    order_h = modeldb%config%finite_element%element_order_h()
    order_v = modeldb%config%finite_element%element_order_v()

    fs => function_space_collection%get_fs(mesh, order_h, order_v, W3)

    ! Create prognostic fields
    ! Creates a field in the W3 function space (fully discontinuous field)
    call field_1%initialise(fs, name="field_1")

    ! Add field to modeldb
    call depository%add_field(field_1)

    ! Create skeleton runtime constants. This creates various things
    ! needed by the fem algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_skeleton_constants(modeldb, mesh, chi, panel_id)

    call log_event('skeleton: Miniapp initialised', log_level_info)

  end subroutine init_skeleton

end module init_skeleton_mod
