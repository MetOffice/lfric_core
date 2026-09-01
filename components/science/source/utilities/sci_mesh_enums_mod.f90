!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------

!> @brief Module to access mesh enumerations.
module sci_mesh_enums_mod

  use constants_mod, only: i_def, imdi
  use log_mod,       only: log_event, log_level_error
  use mesh_mod,      only: mesh_type

  use base_mesh_config_mod, only:                                 &
          sci_geometry_spherical      => geometry_spherical,      &
          sci_geometry_planar         => geometry_planar,         &
          sci_topology_fully_periodic => topology_fully_periodic, &
          sci_topology_non_periodic   => topology_non_periodic

  implicit none

  private

  public :: geometry_spherical, geometry_planar
  public :: topology_fully_periodic, topology_non_periodic

  public :: get_mesh_enums

  ! These will get switched to something hardcoded for the science components
  ! at a later date to break dependence on base_mesh_config_mod
  integer(i_def), parameter :: geometry_spherical      = sci_geometry_spherical      ! 157
  integer(i_def), parameter :: geometry_planar         = sci_geometry_planar         ! 358
  integer(i_def), parameter :: topology_fully_periodic = sci_topology_fully_periodic ! 492
  integer(i_def), parameter :: topology_non_periodic   = sci_topology_non_periodic   ! 157

contains

!---------------------------------------------------------------------------
!> @brief      Returns mesh enumerations in line with science component values
!> @param[in]  mesh     Mesh object to query
!> @param[out] geometry [optional] Science component geometry enumeration
!> @param[out] topology [optional] Science component topology enumeration
!>
subroutine get_mesh_enums(mesh, geometry, topology)

  implicit none

  type(mesh_type), intent(in) :: mesh

  integer(i_def), optional, intent(out) :: geometry
  integer(i_def), optional, intent(out) :: topology

  if (present(geometry)) then
    geometry = imdi

    if (mesh%is_geometry_spherical()) then
      geometry = geometry_spherical

    else if (mesh%is_geometry_planar()) then
      geometry = geometry_planar

    else
      call log_event('Unsupported mesh geometry', log_level_error)

    end if
  end if

  if (present(topology)) then
    topology = imdi

    if (mesh%is_topology_periodic()) then
      topology = topology_fully_periodic

    else if (mesh%is_topology_non_periodic()) then
      topology = topology_non_periodic

    else
      call log_event('Unsupported mesh topology', log_level_error)

    end if
  end if

end subroutine get_mesh_enums

end module sci_mesh_enums_mod
