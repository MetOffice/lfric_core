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
                             only : cell_size_outer, &
                                    cell_size_inner, &
                                    n_cells_stretch, &
                                    n_cells_outer,   &
                                    poly_power
  implicit none

  public :: polynomial_stretch, &
            polynomial_parameters

contains

!> @brief Calculate the polynomial stretching parameters
!> @details In inner y = b x, in stretch y = a (x - xi) ^n + b x
!!          and in outer y = yo + c (x - xo).
!> @param param_a   Parameter a
!> @param param_b   Parameter b
!> @param param_c   Parameter c
!> @param x_inner   Unit mesh coordinate betwen inner and stretch
!> @param x_outer   Unit mesh coordinate between stretch and outer
!> @param dx        Unit mesh cell size
!> @param direction North-south or East-west
subroutine polynomial_parameters( param_a, param_b, param_c, &
                                  x_inner, x_outer, dx, direction )

  implicit none

  real(r_def), intent(inout) :: param_a, param_b, param_c, x_inner, x_outer
  real(r_def),    intent(in) :: dx
  integer(i_def), intent(in) :: direction

  real(r_def) :: l_stretch

  ! Given the coordinates x defined on [-1,1] with mesh size dx,
  ! define new coordinates y such that in the outer and inner regions,
  ! the spacing is cell_size_outer and cell_size_inner and in the
  ! stretch region (in between the inner and outer) the coordinates
  ! satisfy y = a ( x - xi) ^n + b x where xi is the boundary between the
  ! inner and stretch region.

  ! We assume that the mesh is symmetrical and centred on (0,0)
  ! | OUTER | STRETCH | INNER | STRETCH | OUTER |

  ! Considering the region [0,1], define the edges of the stretch region
  x_outer = 1.0_r_def - ( n_cells_outer(direction) * dx )
  x_inner = 1.0_r_def - ( ( n_cells_outer(direction) + &
                            n_cells_stretch(direction) ) * dx )

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
