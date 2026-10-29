!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Module to define functions for a stretched regional mesh.
!> @details The coordinates of a stretched regional mesh are defined
!>          from a prescribed function.
module stretch_transform_mod

  use constants_mod,              only : r_def, i_def, l_def
  use planar_mesh_config_mod,     only : domain_centre
  use stretch_transform_config_mod, &
                                  only : cell_size_outer,           &
                                         cell_size_inner,           &
                                         n_cells_stretch,           &
                                         n_cells_outer,             &
                                         stretching_on,             &
                                         stretching_on_cell_nodes,  &
                                         stretching_on_p_points,    &
                                         stretching_on_cell_centres
  use log_mod,                    only : log_event,                 &
                                         log_scratch_space,         &
                                         LOG_LEVEL_INFO,            &
                                         LOG_LEVEL_ERROR,           &
                                         LOG_LEVEL_DEBUG

  implicit none

  private :: calculate_inflation

  public :: calculate_settings, &
            stretch_transform

  ! Prefix for error messages.
  character(len=*), parameter :: PREFIX = "[Stretch Transform] "

contains

!> @brief   Calculate the inflation factor in the stretch region.
!> @details In the stretch region, the grid box size is gradually increased from
!>          delta_inner to delta_outer, with the grid size being multiplied by
!>          (1 + inflation) on every grid box.
!> @param[in] n_stretch    Number of cells in the stretch region.
!> @param[in] delta_inner  Grid spacing in the inner region.
!> @param[in] delta_outer  Grid spacing in the outer region.
!> @return    inflation    Inflation factor
function calculate_inflation( n_stretch, delta_inner, delta_outer ) &
                              result( inflation )

  implicit none

  real(r_def), intent(in) :: delta_inner
  real(r_def), intent(in) :: delta_outer
  integer(i_def), intent(in) :: n_stretch

  real(r_def) :: inflation

  ! We require:
  !
  ! delta_outer = delta_inner * (1 + inflation) ^ n_stretch
  !
  ! where ^ indicates to the power of.

  inflation  = ( delta_outer / delta_inner ) ** &
               ( 1.0_r_def / n_stretch ) - 1.0_r_def

  return

end function calculate_inflation

!> @brief Calculate the intermediate variables required to define the stretched
!!        grid mesh.
!> @details Define the number of cells and coordinates of each region.
!> @param[in] axis_direction 1 = x or longitude, 2 = y or latitude
!> @param[in] total_n        Total number of cells along mesh edge
!> @param[out] n_inner       Number of cells in the inner region.
!> @param[out] n_stretch     Number of cells in the stretch region.
!> @param[out] inflation_factor Inflation factor
!> @param[out] outer_ends_l  Coordinates at the left edge of the outer region
!!                           of the stretching function
!> @param[out] outer_ends_r  Coordinates at the right edge of the outer region
!!                           of the stretching function
!> @param[out] inner_ends_l  Coordinates at the left edge of the inner region
!!                           of the stretching function
!> @param[out] inner_ends_r  Coordinates at the right edge of the inner region
!!                           of the stretching function
subroutine calculate_settings( axis_direction, total_n,    &
                               n_inner, n_stretch,         &
                               inflation_factor,           &
                               outer_ends_l, outer_ends_r, &
                               inner_ends_l, inner_ends_r )

  integer(i_def), intent(in) :: axis_direction
  integer(i_def), intent(in) :: total_n
  integer(i_def), intent(inout) :: n_inner, n_stretch
  real(r_def), intent(inout) :: outer_ends_l, outer_ends_r
  real(r_def), intent(inout) :: inner_ends_l, inner_ends_r
  real(r_def), intent(inout) :: inflation_factor

  real(r_def) :: stretch_ends_l, stretch_ends_r, stretch_depth_real
  real(r_def) :: delta_inner, delta_outer
  integer(i_def) :: n_outer, i

  ! Extract the input data for the given direction from the config data
  delta_inner = cell_size_inner( axis_direction )
  delta_outer = cell_size_outer( axis_direction )

  ! Calculate the inner, outer and stretch region parameters for the
  ! stretching function, depending on what function space the stretching
  ! is applied. This ensures that the stretching function is created in such
  ! a way to ensure that when it is mapped to the node coordinates, the resulting
  ! number of cells in the stretch region is n_cells_stretch, and the number of cells
  ! in the outer region is n_cells_outer.
  select case (stretching_on)
    case( stretching_on_cell_centres )
      n_stretch = n_cells_stretch( axis_direction ) + 2
      n_outer = n_cells_outer( axis_direction ) - 2
      n_inner = total_n - ( 2 * n_stretch ) - ( 2 * n_outer )

    case( stretching_on_p_points )
      n_stretch = n_cells_stretch( axis_direction )
      n_outer = n_cells_outer( axis_direction ) - 1
      n_inner = total_n - ( 2 * n_stretch ) - ( 2 * n_outer )

    case( stretching_on_cell_nodes )
      n_stretch = n_cells_stretch( axis_direction ) + 1
      n_outer = n_cells_outer( axis_direction ) - 2
      n_inner = total_n - 1 - ( 2 * n_stretch ) - ( 2 * n_outer )

    case default
      call log_event( PREFIX//"Unrecognised value of stretching_on", &
                      LOG_LEVEL_ERROR )
  end select

  if ( n_stretch < 0 ) then

    write(log_scratch_space,'(A)') &
    PREFIX//'n_stretch is negative'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)

  end if

  if ( n_outer < 0 ) then

    write(log_scratch_space,'(A)') &
    PREFIX//'n_outer is negative'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)

  end if

  if ( n_inner < 0 ) then

    write(log_scratch_space,'(A)') &
    PREFIX//'n_inner is negative'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)

  end if

  ! Calculate the stretch factor, and number of cells in stretch region.
  inflation_factor = calculate_inflation( n_stretch, delta_inner, delta_outer )

  ! Subtract half the inner domain to find the left side of the inner,
  ! and similarly add half to find the right side.
  inner_ends_l = domain_centre( axis_direction ) - &
                 0.5_r_def * ( n_inner - 1 ) * delta_inner

  inner_ends_r = domain_centre( axis_direction ) + &
                 0.5_r_def * ( n_inner - 1 ) * delta_inner

  ! Calculate the coordinates of the boundaries between the inner, stretch
  ! and outer regions.

  !
  ! outer       stretch      inner         inner       stretch        outer
  ! ends_l      ends_l       ends_l        ends_r      ends_r         ends_r
  !  |            |             |            |            |             |
  !  |            |             |            |            |             |
  !  ---Outer----   --Stretch--    --Inner--   --Stretch-- ---Outer----
  !     Left          Left                        Right       Right
  !    n_outer      n_stretch       n_inner     n_stretch    n_outer
  !  delta_outer                  delta_inner              delta_outer
  !

  ! To calculate stretch_ends_l, start at inner_ends_l and keep subtracting
  ! cells until stretch_ends_l is reached. Similarly for stretch_ends_r.

  stretch_depth_real = 0.0_r_def
  do i=1, n_stretch
    stretch_depth_real = stretch_depth_real &
                       + ( 1.0_r_def + inflation_factor ) ** i * delta_inner
  end do
  stretch_ends_r = inner_ends_r + stretch_depth_real
  stretch_ends_l = inner_ends_l - stretch_depth_real

  ! outer_ends_l is stretch_ends_l minus n_outer cells
  outer_ends_l = stretch_ends_l - delta_outer * n_outer

  ! outer_ends_r is stretch_ends_r plus n_outer cells
  outer_ends_r = stretch_ends_r + delta_outer * n_outer

  !---------------------- Log derived values---------------------------

  if (axis_direction == 1) then
    write(log_scratch_space,'(A)') &
    'x-direction or Longitude'
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    write(log_scratch_space,'(A, 4(F16.11), A)') &
    'ends = [', outer_ends_l, inner_ends_l, inner_ends_r, outer_ends_r, ' ]'
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    write(log_scratch_space,'(A, 2(F16.11), A)') &
    'deltas = [', delta_outer, delta_inner, ']'
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    write(log_scratch_space,'(A, F16.11 )') &
    'infl_target = ', inflation_factor
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

  elseif (axis_direction == 2 ) then
    write(log_scratch_space,'(A)') &
    'y-direction or Latitude'
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    write(log_scratch_space,'(A, 4(F16.11), A)') &
    'ends = [', outer_ends_l, inner_ends_l, inner_ends_r, outer_ends_r, ']'
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    write(log_scratch_space,'(A, 2(F16.11), A)') &
    'deltas = [', delta_outer, delta_inner, ']'
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    write(log_scratch_space,'(A, F16.11 )') &
    'infl_target = ', inflation_factor
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

  else
    write(log_scratch_space,'(A)') &
    PREFIX//'axis_direction value unknown'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

end subroutine calculate_settings

!> @brief Stretch transform definition, as a function of the cell index.
!> @details This is used to create a variable-resolution regional mesh.
!>        For ease, variable names (l=left and r=right) are used to describe
!>        the function in the x-direction. However, the same function is used
!>        in y-direction.
!> @param[in]     cell     The number of cells from the mesh edge
!> @param[in,out] node_coord The coordinate of the cell node
!>                         (in - previous, out - current)
!> @param[in,out] prev_stretch_coord   The coordinate of the stretching function
!>                         (in - previous, out - current) when using
!>                         stretching on p_points
!> @param[in] axis_direction  1=x (Longitude), 2=y (Latitude)
!> @param[in] n_inner      Number of cells in the inner region.
!> @param[in] n_stretch    Number of cells in the stretch region.
!> @param[in] inflation_factor Inflation factor
!> @param[in] outer_ends_l Coordinates at the left edge of the outer region
!> @param[in] outer_ends_r Coordinates at the right edge of the outer region
!> @param[in] inner_ends_l Coordinates at the left edge of the inner region
!> @param[in] inner_ends_r Coordinates at the right edge of the inner region
subroutine stretch_transform( cell, node_coord,           &
                              prev_stretch_coord,         &
                              axis_direction,             &
                              n_inner, n_stretch,         &
                              inflation_factor,           &
                              outer_ends_l, outer_ends_r, &
                              inner_ends_l, inner_ends_r )

  implicit none

  integer(i_def), intent(in) :: cell
  integer(i_def), intent(in) :: axis_direction
  real(r_def), intent(inout) :: node_coord
  real(r_def), intent(inout) :: prev_stretch_coord

  integer(i_def) :: n_inner, n_outer, n_stretch
  real(r_def) :: delta_inner, delta_outer, inflation_factor
  real(r_def) :: outer_ends_l, outer_ends_r
  real(r_def) :: inner_ends_l, inner_ends_r

  integer(i_def) :: c, stretch_cell
  real(r_def) :: stretching_coord

  delta_inner = cell_size_inner( axis_direction )
  delta_outer = cell_size_outer( axis_direction )

  ! Adjustments to allow the stretching function
  ! to be applied to different function spaces.
  select case (stretching_on)

    case( stretching_on_cell_centres )
      stretch_cell = cell
      n_outer = n_cells_outer( axis_direction ) - 2

    case( stretching_on_p_points )
      stretch_cell = cell + 1
      n_outer = n_cells_outer( axis_direction ) - 1

    case( stretching_on_cell_nodes )
      stretch_cell = cell
      n_outer = n_cells_outer( axis_direction ) - 2

    case default
      call log_event( PREFIX//"Unrecognised value of stretching_on", &
                      LOG_LEVEL_ERROR )
  end select

  ! First cell
  if ( stretch_cell == 0 ) then

    stretching_coord = outer_ends_l - delta_outer

  ! Outer region - left
  else if ( stretch_cell > 0 .and. &
            stretch_cell <= n_outer + 1 ) then

    stretching_coord = outer_ends_l &
          + ( stretch_cell - 1 ) * delta_outer

  ! Stretch region - left
  else if ( stretch_cell > n_outer + 1 .and. &
            stretch_cell <= ( n_outer + n_stretch + 1 ) ) then

    stretching_coord = inner_ends_l
    do c = 1, ( n_outer + n_stretch + 1 ) - stretch_cell
      stretching_coord = stretching_coord - ( 1.0_r_def + inflation_factor )** (c) * delta_inner
    end do

  ! Inner region
  else if ( stretch_cell >  ( n_outer + n_stretch + 1 ) .and. &
            stretch_cell < ( n_outer + n_stretch + n_inner ) ) then

    stretching_coord = inner_ends_l &
          + (stretch_cell - ( n_outer + n_stretch + 1 ) ) * delta_inner

  ! Stretch region - right
  else if ( stretch_cell >= ( n_outer + n_stretch + n_inner ) .and. &
            stretch_cell < ( n_outer + 2*n_stretch + n_inner ) ) then

    stretching_coord = inner_ends_r
    do c = 1, stretch_cell - ( n_outer + n_stretch + n_inner )
      stretching_coord = stretching_coord + (1.0_r_def + inflation_factor)** (c) * delta_inner
    end do

  ! Outer region - right
  else if ( stretch_cell >= ( n_outer + 2*n_stretch + n_inner ) .and. &
            stretch_cell < ( 2*n_outer + 2*n_stretch + n_inner + 1 ) ) then

    stretching_coord = outer_ends_r - ( n_outer * delta_outer )  &
          + ( stretch_cell - ( n_outer + 2*n_stretch + n_inner ) ) * delta_outer

  ! Final cells
  else if ( stretch_cell == ( 2*n_outer + 2*n_stretch + n_inner + 1 ) ) then

    stretching_coord = outer_ends_r + delta_outer

  ! For spacing on cell-nodes
  else if ( stretch_cell == ( 2*n_outer + 2*n_stretch + n_inner + 2 ) ) then

    stretching_coord = outer_ends_r + 2.0_r_def * delta_outer

  else
    write(log_scratch_space,'(A, I0)') &
    PREFIX//"Invalid cell number ", cell
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  ! Use the stretching function coordinate to define the coordinates
  ! of the cell nodes, depending on which function space the stretching
  ! is to be defined.

  select case (stretching_on)

    case( stretching_on_cell_centres )
      ! Generate the cell-node coordinates such that the
      ! stretching_function_coord s_i is half-way between the cell node
      ! coords n_i.
      ! i.e. The stretching function is on the cell centres.
      !
      ! s_i = (n_i+1 + n_i)/2 which can be rearranged as
      ! n_i+1 = n_i + 2(s_i - n_i)
      !
      ! At the boundary to the inner and outer domains
      ! n_i is defined directly from delta_inner and delta_outer,
      ! to make sure that there is an even spacing in the
      ! regular parts of the domain.

      if ( cell == 0 ) then

        node_coord = stretching_coord + 0.5_r_def * delta_outer

      else if ( cell ==  ( n_outer + n_stretch + 1 ) ) then

        node_coord = inner_ends_l + 0.5_r_def * delta_inner

      else if ( cell == ( n_outer + 2*n_stretch + n_inner - 2  ) ) then

        node_coord = outer_ends_r - ( n_outer + 1.5_r_def ) * delta_outer

      else if ( cell == ( 2*n_outer + 2*n_stretch + n_inner  ) ) then

        node_coord = outer_ends_r + 0.5_r_def * delta_outer

      else
        node_coord = node_coord + 2.0_r_def * ( stretching_coord - node_coord )

      end if

    case( stretching_on_p_points )
      ! Generate the cell-node coordinates such that the cell-nodes
      ! are half-way between the stretching_function_coord. The stretching
      ! function is neither on the cell nodes or cell centres - it is on
      ! ficticious points such that the nodes are half-way between the p-points.

       if ( cell == 0 ) then

         node_coord = stretching_coord - 0.5_r_def * delta_outer

       else if ( cell == ( 2*n_outer + 2*n_stretch + n_inner ) ) then

         node_coord = stretching_coord - 0.5_r_def * delta_outer

       else

         node_coord = ( stretching_coord + prev_stretch_coord ) / 2.0_r_def

       end if

       prev_stretch_coord = stretching_coord

   case( stretching_on_cell_nodes )
      ! Generate the cell-node coordinates such that the cell-nodes
      ! are the same as the stretching_function_coord. i.e.
      ! the stretching function is on the cell nodes.

        node_coord = stretching_coord

   case default

  end select

  write(log_scratch_space,'(A, I0, 2(A,F16.10))') 'cell = ', cell, &
                                                  ' stretch-coord = ', stretching_coord, &
                                                  ' cell-node-coord = ', node_coord
  call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

end subroutine stretch_transform

end module stretch_transform_mod
