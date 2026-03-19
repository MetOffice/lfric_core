!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Module to define a coordinate transformation for a stretched
!!          regional mesh.
!> @details The coordinate transformation is defined using a polynomial
!!          function.
!>
module polynomial_stretching_mod

  use constants_mod,         only: r_def, i_def, l_def
  use stretch_transform_config_mod, &
                             only : cell_size_outer,      &
                                    cell_size_inner,      &
                                    n_cells_stretch_nsew, &
                                    n_cells_outer_nsew,   &
                                    poly_power
  implicit none

  public :: associated_direction, &
            calculate_offset,     &
            polynomial_stretch,   &
            polynomial_parameters

contains

function associated_direction( boundary ) result(direction)

  implicit none

  integer(i_def) :: boundary
  integer(i_def) :: direction

  if (boundary == 1 .or. boundary == 2) then
    ! North-South
    direction = 2
  else
    ! East-West
    direction = 1
 end if

end function associated_direction

!> @brief Calculate the offset to apply before the polynomial stretching
!> @details If the number of cells in the outer/stretch region on one boundary
!!          e.g. the North, is not the same as the number of cells on the other
!!          boundary (the South) then calculate the offset so that the high
!!          resolution interior will be centred at (0,0).
!> @param direction Direction (N-S) or (E-W)
function calculate_offset( direction ) result(offset)

  integer(i_def), intent(in) :: direction
  real(r_def) :: offset
  real(r_def) :: dx

  dx = cell_size_inner(direction)
  if ( direction == 2 ) then
    ! North-South
    offset = (n_cells_outer_nsew(1) + n_cells_stretch_nsew(1)) - &
             (n_cells_outer_nsew(2) + n_cells_stretch_nsew(2))
    offset = 0.5_r_def * dx * offset
  else
    ! East-West
    offset = (n_cells_outer_nsew(3) + n_cells_stretch_nsew(3)) - &
             (n_cells_outer_nsew(4) + n_cells_stretch_nsew(4))
    offset = 0.5_r_def * dx * offset
 end if

end function calculate_offset

!> @brief Calculate the polynomial stretching parameters
!> @details In inner y = b x, in stretch y = a (x - xi) ^n + b x
!!          and in outer y = yo + c (x - xo).
!> @param param_a   Parameter a
!> @param param_b   Parameter b
!> @param param_c   Parameter c
!> @param x_domain  Uniform mesh coordinate at the end of the mesh
!> @param x_inner   Uniform mesh coordinate betwen inner and stretch
!> @param x_outer   Uniform mesh coordinate between stretch and outer
!> @param dx        Uniform mesh cell size
!> @param boundary  1 North, 2 South, 3 East or 4 West
subroutine polynomial_parameters( param_a, param_b, param_c, &
                                  x_domain, x_inner, x_outer, dx, boundary )

  implicit none

  real(r_def), intent(inout) :: param_a, param_b, param_c, x_inner, x_outer
  real(r_def),    intent(in) :: x_domain, dx
  integer(i_def), intent(in) :: boundary

  real(r_def) :: l_stretch
  integer(i_def) :: direction

  ! Given the coordinates x with mesh size dx,
  ! define new coordinates y such that in the outer and inner regions,
  ! the spacing is cell_size_outer and cell_size_inner and in the
  ! stretch region (in between the inner and outer) the coordinates
  ! satisfy y = a ( x - xi) ^n + b x where xi is the boundary between the
  ! inner and stretch region.

  direction = associated_direction(boundary)

  ! We only consider the region [0,x_domain]
  ! | INNER    | STRETCH   |    OUTER   |
  !         x_inner     x_outer      x_domain

  ! Define the edges of the stretch region
  x_outer = x_domain - ( n_cells_outer_nsew(boundary) * dx )
  x_inner = x_domain - ( ( n_cells_outer_nsew(boundary) + &
                           n_cells_stretch_nsew(boundary) ) * dx )

  ! Define the total size or length of the stretch region
  l_stretch = ( x_outer - x_inner )

  ! In outer region y = c (x -xo)
  ! y' = c so c = target cell_size / dx

  param_c = cell_size_outer(direction) / dx

  ! In inner region and at x = xi (between inner and stretch)
  ! y' = b so b = target cell_size /dx

  param_b = cell_size_inner(direction) / dx

  ! In stretch region y = a (x - xi) ^n + bx
  ! Derivative y' = n a (x - xi) ^(n-1) + b
  ! At x = xo (between stretch and outer), where xo - xi = l
  ! Set n a (x - xi) ^(n-1) + b = c
  ! So a = (c - b) / ( n l ^(n-1) )

  param_a = ( param_c - param_b ) / &
            ( poly_power * l_stretch ** (poly_power - 1_i_def) )

end subroutine polynomial_parameters

!> @brief Apply a polynomial stretching transformation to a given coordinate
!> @details In inner y = b x, in stretch y = a (x - xi) ^n + b x
!!          and in outer y = yo + c x
!> @param param_a   Parameter a
!> @param param_b   Parameter b
!> @param param_c   Parameter c
!> @param x_inner   Unit mesh coordinate betwen inner and stretch
!> @param x_outer   Unit mesh coordinate between stretch and outer
!> @param dx        Unit mesh cell size
!> @param direction North-south or East-west
function polynomial_stretch( x_coord, param_a, param_b, param_c, &
                             x_inner, x_outer ) &
                             result( y_coord )

  implicit none

  real(r_def), intent(in) :: x_coord
  real(r_def), intent(in) :: param_a, param_b, param_c, x_inner, x_outer

  real(r_def) :: y_coord, y_outer, l_stretch, new_x_coord

  logical(l_def) :: use_symmetry

  ! Define the total size or length of the stretch region
  l_stretch =  x_outer - x_inner

  ! Define a useful constant that describes the new coordinate at the
  ! point between the stretch and outer regions.
  y_outer = ( param_a * l_stretch ** poly_power ) + &
            ( param_b * x_outer )

  ! Use symmetry to define coords < 0
  if ( x_coord < 0.0_r_def ) then
    use_symmetry = .true.
    new_x_coord = -1.0_r_def * x_coord
  else
    use_symmetry= .false.
    new_x_coord = x_coord
  end if

  ! Assign new coordinates using transform y=f(x)
  if ( new_x_coord < x_inner ) then
    ! In inner y = b x
    y_coord = param_b * new_x_coord

  else if ( new_x_coord >= x_inner .and. new_x_coord < x_outer ) then
    ! In stretch y = a (x - xi) ^n + bx where a (x - xi) ^n >0
    y_coord = param_b * new_x_coord + &
              param_a * ( new_x_coord - x_inner ) ** poly_power

  else
    ! In outer y = c (x - xo) + yo
    y_coord = param_c * ( new_x_coord - x_outer ) + y_outer
  end if

  ! To define coords <0
  if ( use_symmetry ) then
    y_coord = -1.0_r_def * y_coord
  end if

  return

end function polynomial_stretch

end module polynomial_stretching_mod
