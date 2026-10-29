!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Module to define a coordinate transformation for a stretched
!!          regional mesh.
!> @details The coordinate transformation is defined using a polynomial
!!          function and applied to the unit mesh coordinates (u_coord)
!!          to give the transformed coordinates (t_coord), so t=f(u).
!!          The transformation is formed by the separate application to
!!          the positive and negative coordinates of each axis (N-S) and
!!          (E-W). i.e. working on single axis at a time, radiating from the
!!          centre to each of the 4 boundaries.
!>
module polynomial_stretching_mod

  use constants_mod,         only: r_def, i_def, l_def
  use stretch_transform_config_mod, &
                             only : cell_size_outer,      &
                                    cell_size_inner,      &
                                    n_cells_stretch_wsen, &
                                    n_cells_outer_wsen,   &
                                    poly_power
  implicit none

  public :: associated_axis_direction, &
            calculate_offset,          &
            polynomial_stretch,        &
            polynomial_parameters

  integer(i_def), public, parameter :: axis_ns = 2 ! Latitude  (North-South)
  integer(i_def), public, parameter :: axis_ew = 1 ! Longitude (East-West)
  integer(i_def), public, parameter :: boundary_w = 1 ! West
  integer(i_def), public, parameter :: boundary_s = 2 ! South
  integer(i_def), public, parameter :: boundary_e = 3 ! East
  integer(i_def), public, parameter :: boundary_n = 4 ! North

contains

!> @brief   Determine the axis associated with boundary
!> @details axis_ns is associated with boundary_n and boundary_s
!!          axis_ew is associated with boundary_e and boundary_w
function associated_axis_direction( boundary ) result(axis_direction)

  implicit none

  integer(i_def) :: boundary
  integer(i_def) :: axis_direction

  if (boundary == boundary_n .or. boundary == boundary_s) then
    ! North-South
    axis_direction = axis_ns
  else
    ! East-West
    axis_direction = axis_ew
 end if

end function associated_axis_direction

!> @brief Calculate the offset to apply.
!> @details If the number of cells in the outer/stretch region on one boundary
!!          e.g. the North, is not the same as the number of cells on the other
!!          boundary (the South) then calculate the offset so that the high
!!          resolution interior will be centred at (0,0).
!> @param axis_direction Axis (N-S) or (E-W)
function calculate_offset( axis_direction ) result(offset)

  integer(i_def), intent(in) :: axis_direction
  real(r_def) :: offset
  real(r_def) :: du

  du = cell_size_inner(axis_direction)
  if ( axis_direction == axis_ns ) then
    ! North-South
    offset = (n_cells_outer_wsen(boundary_N) + n_cells_stretch_wsen(boundary_N)) - &
             (n_cells_outer_wsen(boundary_S) + n_cells_stretch_wsen(boundary_S))
  else
    ! East-West
    offset = (n_cells_outer_wsen(boundary_E) + n_cells_stretch_wsen(boundary_E)) - &
             (n_cells_outer_wsen(boundary_W) + n_cells_stretch_wsen(boundary_W))
 end if
 offset = 0.5_r_def * du * offset

end function calculate_offset

!> @brief Calculate the polynomial stretching parameters
!> @details Stretching function t=f(u) applied to separate regions:
!!          Inner region:   t = b u
!!          Stretch region: t = a (u - ui) ^n + b u
!!          Outer region:   t = to + c (u - uo).
!!          This subroutine returns the parameters a, b and c
!!          used in these functions.
!> @param param_a   Parameter a, used in stretch region
!> @param param_b   Parameter b, used in inner region
!> @param param_c   Parameter c, used in outer region
!> @param u_domain  u coordinate of domain boundary (uniform-mesh)
!> @param u_inner   u coordinate of boundary between inner and stretch regions (uniform-mesh)
!> @param u_outer   u coordinate of boundary between stretch and outer regions (uniform-mesh)
!> @param du        Cell size (uniform-mesh)
!> @param boundary  Boundary enumeration which maps to a particular
!!                  direction [W|S|E|N]
subroutine polynomial_parameters( param_a, param_b, param_c, &
                                  u_domain, u_inner, u_outer, du, boundary )

  implicit none

  real(r_def), intent(inout) :: param_a, param_b, param_c, u_inner, u_outer
  real(r_def),    intent(in) :: u_domain, du
  integer(i_def), intent(in) :: boundary

  real(r_def) :: l_stretch
  integer(i_def) :: axis_direction

  ! Given the coordinates (u) with mesh size (du),
  ! define new coordinates (t) such that in the outer and inner regions,
  ! the spacing is cell_size_outer and cell_size_inner and in the
  ! stretch region (in between the inner and outer) the coordinates
  ! satisfy t = a ( u - ui) ^n + b u where ui is the boundary between the
  ! inner and stretch region.

  axis_direction = associated_axis_direction(boundary)

  ! We only consider the region [0,u_domain]
  ! | INNER    | STRETCH   |    OUTER   |
  !         u_inner     u_outer      u_domain

  ! Define the edges of the stretch region
  u_outer = u_domain - ( n_cells_outer_wsen(boundary) * du )
  u_inner = u_domain - ( ( n_cells_outer_wsen(boundary) + &
                           n_cells_stretch_wsen(boundary) ) * du )

  ! Define the total size or length of the stretch region
  l_stretch = ( u_outer - u_inner )

  ! In outer region t = c (u -uo)
  ! First derivative t' = c so c = target cell_size / du

  param_c = cell_size_outer(axis_direction) / du

  ! In inner region and at u = ui (between inner and stretch)
  ! First derivative t' = b so b = target cell_size /du

  param_b = cell_size_inner(axis_direction) / du

  ! In stretch region t = a (u - ui) ^n + bu
  ! Derivative t' = n a (u - ui) ^(n-1) + b
  ! At u = uo (between stretch and outer), where uo - ui = l
  ! Set n a (u - ui) ^(n-1) + b = c
  ! So a = (c - b) / ( n l ^(n-1) )

  param_a = ( param_c - param_b ) / &
            ( poly_power * l_stretch ** (poly_power - 1_i_def) )

end subroutine polynomial_parameters

!> @brief Apply a polynomial stretching transformation to a given coordinate
!> @details Stretching function t=f(u) applied to separate regions:
!!          Inner region:   t = b u
!!          Stretch region: t = a (u - ui) ^n + b u
!!          Outer region:   t = to + c (u - uo).
!!          This subroutine calculates the value t, given u.
!> @param u_coord   The input coordinate (unit-mesh)
!> @param param_a   Parameter a, used in stretch region
!> @param param_b   Parameter b, used in inner region
!> @param param_c   Parameter c, used in outer region
!> @param u_inner   u coordinate of boundary between inner and stretch regions (unit-mesh)
!> @param u_outer   u coordinate of boundary between stretch and outer regions (unit-mesh)
!> @param t_coord   The transformed output coordinate
function polynomial_stretch( u_coord, param_a, param_b, param_c, &
                             u_inner, u_outer ) &
                             result( t_coord )

  implicit none

  real(r_def), intent(in) :: u_coord
  real(r_def), intent(in) :: param_a, param_b, param_c, u_inner, u_outer

  real(r_def) :: t_coord, t_outer, l_stretch, new_u_coord

  logical(l_def) :: use_symmetry

  ! Define the total size or length of the stretch region
  l_stretch =  u_outer - u_inner

  ! Define a useful constant that describes the new coordinate at the
  ! point between the stretch and outer regions.
  t_outer = ( param_a * l_stretch ** poly_power ) + &
            ( param_b * u_outer )

  ! Use symmetry to define coords < 0
  if ( u_coord < 0.0_r_def ) then
    use_symmetry = .true.
    new_u_coord = -1.0_r_def * u_coord
  else
    use_symmetry= .false.
    new_u_coord = u_coord
  end if

  ! Assign new coordinates using transform y=f(x)
  if ( new_u_coord < u_inner ) then
    ! In inner t = b u
    t_coord = param_b * new_u_coord

  else if ( new_u_coord >= u_inner .and. new_u_coord < u_outer ) then
    ! In stretch t = a (u - ui) ^n + bu where a (u - ui) ^n >0
    t_coord = param_b * new_u_coord + &
              param_a * ( new_u_coord - u_inner ) ** poly_power

  else
    ! In outer t = c (u - uo) + to
    t_coord = param_c * ( new_u_coord - u_outer ) + t_outer
  end if

  ! To define coords <0
  if ( use_symmetry ) then
    t_coord = -1.0_r_def * t_coord
  end if

  return

end function polynomial_stretch

end module polynomial_stretching_mod
