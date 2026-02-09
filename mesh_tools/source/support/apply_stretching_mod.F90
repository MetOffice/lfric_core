module apply_stretching_mod

  use constants_mod  :: r_def, i_def 
  use gen_planar_mod :: gen_planar_type
  
  implicit none

  private :: cubic_stretch, &
             cubic_parameters

  public :: apply_uniform_resolution, &
            apply_cubic_stretch,      &
            apply_cosine_stretch

contains

!> @brief Apply the uniform resolution stretching to the unit mesh coordinates
subroutine apply_uniform_resolution(self)

  class(gen_planar_type), intent(inout) :: self

  implicit none

  real(r_def) :: dx, param_a
  integer(i_def) :: nverts, vert, direction

  do direction = 1, 2

    if ( direction == 1 ) then
      dx = 2.0_r_def / self%edge_cells_x
      param_a = self%dx / dx
    else
      dx = 2.0_r_def / self%edge_cells_y
      param_a = self%dy / dy
    end if
      
    nverts = size(self%vert_coords(direction, :))
      
    do vert = 1, nverts
 
      self%vert_coords(direction, vert) = param_a * &
                                          self%vert_coords(direction, vert)
         
    end do
      
  end do
   
  return

end subroutine apply_uniform_resolution

!> @brief Apply the cubic stretching to the unit mesh coordinates
subroutine apply_cubic_stretch(self)

  implicit none

  class(gen_planar_type), intent(inout)  :: self

  real(r_def) :: dx, param_a, param_b, param_c, x_inner, x_outer
  integer(i_def) :: nverts, vert, direction

  do direction = 1, 2

    ! Calculate the cell spacing and numeber of points of unit mesh
    if ( direction == 1 ) then
      dx = 2.0_r_def / self%edge_cells_x
    else
      dx = 2.0_r_def / self%edge_cells_y
    end if

    nverts = size(self%vert_coords(direction, :))

    ! Calculate the parameters required for the stretching transform
    call cubic_parameters( param_a, param_b, param_c, &
                           x_inner, x_outer, dx, direction )
     
    ! Apply the transformation to each coordinate
    do vert = 1, nverts
 
       self%vert_coords(direction, vert) = &
            cubic_stretch(self%vert_coords(direction, vert), &
                          param_a, param_b, param_c, x_inner, x_outer )
         
    end do
      
  end do
  
  return
  
end subroutine apply_cubic_stretch

!> @brief Apply the cosine stretching to the unit mesh coordinates
subroutine apply_cosine_stretch(self)
end subroutine apply_cosine_stretch

!> @brief Calculate the cubic stretching parameters
!> @details In inner y = b x, in stretch y = a x^3 + b x
!!          and in outer y = c x.
!> @param param_a   Parameter a
!> @param param_b   Parameter b
!> @param param_c   Parameter c
!> @param x_inner   Unit mesh coordinate betwen inner and stretch
!> @param x_outer   Unit mesh coordinate between stretch and outer
!> @param dx        Unit mesh cell size
!> @param direction North-south or East-west
subroutine cubic_parameters( param_a, param_b, param_c, &
                             x_inner, x_outer, dx, direction )

  use stretch_transform_config_mod, &
                                  only : cell_size_outer,           &
                                         cell_size_inner,           &
                                         n_cells_stretch,           &
                                         n_cells_outer
  implicit none

  real(r_def), intent(inout) :: param_a, param_b, param_c, x_inner, x_outer
  real(r_def), intent(in) :: dx
  integer(i_def), intent(in) :: direction
  
  real(r_def) :: l_stretch

  ! Given the coordinates x defined on [-1,1] with mesh size dx,
  ! define new coordinates y such that in the outer and inner regions,
  ! the spacing is cell_size_outer and cell_size_inner and in the
  ! stretch region (in between the inner and outer) the coordinates
  ! satisfy y = a x^3 + bx

  ! We assume that the mesh is symmetrical and centred on (0,0)
  ! | OUTER | STRETCH | INNER | STRETCH | OUTER |

  ! Considering the region [0,1], define the edges of the stretch region
  x_outer = 1.0_r_def - ( n_cells_outer(direction) * dx )
  x_inner = 1.0_r_def - ( ( n_cells_outer(direction) + &
                            n_cells_stretch(direction) ) * dx )

  ! Define the total size or length of the stretch region
  l_stretch = x_outer - x_inner

  ! y = a x^3 + bx
  ! Derivatives
  ! y'   = 3 a x^2 + b
  ! y''  = 6 a x
  ! y''' = 6 a
  
  ! In inner region and at x = 0 (between inner and stretch)
  ! dy/dx =b so b = target cell_size /dx
  
  param_b = cell_size_inner(direction) / dx
  
  ! At x = l_stretch (between stretch and outer)
  ! Use Taylor series expansion
  ! y(x + dx) = y(x) + dx y'(x) + dx^2 y''(x)/2 + dx^3 y'''(x)/6
  ! y(x + dx) - y(x) = dx (3 a x^2 + b) + dx^2 (6ax)/2 + dx^3(6a)/6
  ! dy_outer         = a dx ( 3x^2 + 3 dx x + dx^2 ) + b dx
  ! a = (dy - b dx) / dx( 3x^2 + 3 dx x + dx^2)
  
  param_a = ( cell_size_outer(direction) - param_b * dx ) / &
            ( dx * ( 3 * l_stretch ** 2  + 3 * l_stretch * dx + dx ** 2 ))

  ! In outer region y = cx
  
  param_c = cell_size_outer(direction) / dx

end subroutine cubic_parameters

!> @brief Apply cubic stretching transformation to a given coordinate
!> @details In inner y = b x, in stretch y = a x^3 + b x
!!          and in outer y = c x
!> @param param_a   Parameter a
!> @param param_b   Parameter b
!> @param param_c   Parameter c
!> @param x_inner   Unit mesh coordinate betwen inner and stretch
!> @param x_outer   Unit mesh coordinate between stretch and outer
!> @param dx        Unit mesh cell size
!> @param direction North-south or East-west
function cubic_stretch( x_coord, param_a, param_b, param_c, x_inner, x_outer ) &
                        result( y_coord )

  implicit none
  
  real(r_def), intent(in) :: x_coord
  real(r_def), intent(in) :: param_a, param_b, param_c, x_inner, x_outer
  
  real(r_def) :: y_coord, outer_constant, l_stretch, new_x_coord
  
  logical(l_def) :: use_symmetry

  ! Define the total size or length of the stretch region
  l_stretch = x_outer - x_inner

  ! Define a useful constant
  outer_constant = ( param_a * l_stretch **3 )  + ( param_b * x_outer )
      
  ! Use symmetry to define coords < 0
  if ( x_coord < 0.0_r_def ) then
    use_symmetry = .true.
    new_x_coord = -1.0_r_def * x_coord
  else
    use_symmetry= .false.
    new_x_coord = x_coord
  end if

  if ( new_x_coord < x_outer ) then
    ! In inner y = bx and in stretch y = ax^3 + bx
    y_coord = param_b * new_x_coord + &
         max ( param_a * ( new_x_coord - x_inner )**3, 0.0_r_def )
  else
    ! In outer y = cx
    y_coord = param_c * ( new_x_coord - x_outer ) + outer_constant  
  end if

  ! To define coords <0
  if ( use_symmetry ) then
     y_coord = -1.0_r_def * y_coord
  end if
 
  return
  
end function cubic_stretch

end module apply_stretching_mod
