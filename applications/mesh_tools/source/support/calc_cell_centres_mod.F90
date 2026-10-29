!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Module to hold routines for calculating cell centre coordinates
!           for meshes used with global/regional models.
!> @details There are separate routines for global/regional models.
!>
!> @todo    In future there should only be one routine to calculate cell
!>          centres. The two routines are maintianed here for legacy
!>          reasons. This issue should be addressed by ticket #3494.
!-----------------------------------------------------------------------------
module calc_cell_centres_mod

  use constants_mod,         only: i_def, r_def
  use coord_transform_mod,   only: ll2xyz, xyz2ll
  use reference_element_mod, only: NWB

  implicit none

  private
  public :: calc_cell_centres_global_model
  public :: calc_cell_centres_regional_model

  integer(i_def), parameter :: NW = NWB

contains

!-------------------------------------------------------------------------------
!> @brief   Calculates the mesh cell centres for meshes representing
!>          global models.
!> @details The face centres for the mesh are calculated based on the current
!>          node coordinates. The node_cordinates are assumed to be in [lon, lat].
!>          The resulting face centre coordinates are also in [lon, lat].
!>
!> @param[in]  nodes_on_cells  IDs of nodes connected to given the cell IDs.
!> @param[in]  node_coords     Lon,lat in radians.
!> @param[out] cell_coords     Lon,lat in radians.
!>
!-------------------------------------------------------------------------------
subroutine calc_cell_centres_global_model( nodes_on_cells, &
                                           node_coords,    &
                                           cell_coords )

  implicit none

  integer(i_def), intent(in)  :: nodes_on_cells(:,:)
  real(r_def),    intent(in)  :: node_coords(:,:)
  real(r_def),    intent(out) :: cell_coords(:,:)

  real(r_def) :: radius_ratio

  real(r_def), allocatable :: cell_node_coords_xyz(:,:)
  real(r_def), allocatable :: cell_node_coords_ll(:,:)
  real(r_def), allocatable :: cell_centre_xyz(:)

  integer(i_def), allocatable :: nodes_on_cell(:)

  integer(i_def) :: nnodes_per_cell
  integer(i_def) :: ncells

  ! Counters
  integer(i_def) :: cell, node

  nnodes_per_cell = size( nodes_on_cells, 1)
  ncells          = size( nodes_on_cells, 2)

  allocate( nodes_on_cell(nnodes_per_cell) )
  allocate( cell_node_coords_xyz(3,nnodes_per_cell) )
  allocate( cell_node_coords_ll(2,nnodes_per_cell) )
  allocate( cell_centre_xyz(3) )

  cell_coords(:,:) = 0.0_r_def

  do cell=1, ncells
    cell_centre_xyz(:) = 0.0_r_def

    ! Get the node ids on the cell
    nodes_on_cell(:) = nodes_on_cells(:,cell)

    do node=1, nnodes_per_cell
      ! Get the node coordinates (in radians)
      cell_node_coords_ll(:,node) = node_coords(:,nodes_on_cell(node))

      ! Get node coordinates as cartesian (x,y,z)
      call ll2xyz( cell_node_coords_ll(1,node),  &
                   cell_node_coords_ll(2,node),  &
                   cell_node_coords_xyz(1,node), &
                   cell_node_coords_xyz(2,node), &
                   cell_node_coords_xyz(3,node) )

      cell_centre_xyz(:) = cell_centre_xyz(:) + cell_node_coords_xyz(:,node)

    end do

    radius_ratio = 1.0_r_def/sqrt(  cell_centre_xyz(1)**2 &
                                  + cell_centre_xyz(2)**2 &
                                  + cell_centre_xyz(3)**2 )

    cell_centre_xyz(:) = cell_centre_xyz(:) * radius_ratio

    ! Convert cell centre back to lat long
    call xyz2ll( cell_centre_xyz(1),   & ! x
                 cell_centre_xyz(2),   & ! y
                 cell_centre_xyz(3),   & ! z
                 cell_coords(1,cell),  & ! longitude
                 cell_coords(2,cell) )   ! latititude
  end do

end subroutine calc_cell_centres_global_model


!-------------------------------------------------------------------------------
!> @brief   Calculates the cell centre coordinates for meshes representing
!>          regional models.
!> @details The face centres for the mesh are calculated based on the current
!>          node coordinates of the node on the NW corner of the cell.
!>
!>          The node_cordinates are assumed to be cartesian in the
!>          x,y plane. Resulting face centre coordinates are in equally
!>          disributed in [m] or [radians], depending on the units the node
!>          coordinates are given.
!>
!-------------------------------------------------------------------------------
subroutine calc_cell_centres_regional_model( dx, dy,        &
                                             nodes_on_cell, &
                                             node_coords,   &
                                             cell_coords )

  implicit none

  real(r_def), intent(in) :: dx
  real(r_def), intent(in) :: dy

  integer(i_def), intent(in)  :: nodes_on_cell(:,:)
  real(r_def),    intent(in)  :: node_coords(:,:)
  real(r_def),    intent(out) :: cell_coords(:,:)

  real(r_def)    :: offset_x
  real(r_def)    :: offset_y
  integer(i_def) :: ncells

  ! Counters
  integer(i_def) :: cell, base_node

  ncells = size(nodes_on_cell,2)

  ! 1.0 Initialise the face centres
  cell_coords(:,:) = 0.0_r_def

  ! 2.0 Open cells have `ghost` nodes/edges and are located
  !     on the Eastern/Southern edges of the panel. All cells,
  !     including the open cells have a unique NW node. Use this
  !     NW node and self%dx, self%dy to calculate the face centre,
  !     assuming the cells are parallel in the x and y directions.
  offset_x = dx*0.5_r_def
  offset_y = dy*0.5_r_def
  do cell=1, ncells
    base_node = nodes_on_cell(NW, cell)
    cell_coords(1,cell) = node_coords(1, base_node) + offset_x
    cell_coords(2,cell) = node_coords(2, base_node) - offset_y
  end do

end subroutine calc_cell_centres_regional_model

end module calc_cell_centres_mod
