!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!>  @brief     Contains routines for calculation of mass through a cell edge for
!!             the split advection scheme.
!!
!!  @details   The routine return_mass returns the mass which passes through a
!!             cell wall in a given timestep.
!!             The inputs to the routine are the departure distance and the rho
!!             values for the neighbouring cells (in 1D), and also the subgrid
!!             coefficients for rho.
!!             The scheme is Courant number unlimited which means that the mass
!!             from a number of cells need to be added.
!-------------------------------------------------------------------------------
module cosmic_flux_mod

use constants_mod, only : r_def

implicit none

private

! Public subroutines
public :: calc_stencil_ordering
public :: stencil_ordering_and_orientation
public :: frac_and_int_part
public :: calc_integration_limits
public :: populate_array
public :: map_cell_index
public :: eval_integral
public :: return_part_mass
public :: w2_dof
public :: reorientate_w2field
public :: dof_to_update
public :: calc_local_vertical_index


!--------------------------------------------------------------------------------
! Contained functions / subroutines
!--------------------------------------------------------------------------------
contains

  !--------------------------------------------------------------------------------
  !>  @brief  Input x, which is typically the departure distance
  !!          Returns int_x and frac_x such that:
  !!          if x_value > 0.0 then x_value = (int_x-1) + frac_x
  !!          if x_value < 0.0 then x_value = -(int_x-1) - frac_x
  !!          with 0 <= frac_x < 1.
  !!
  !!  @param[in]   x_value  x value
  !!  @param[out]  int_x    Positive x_value rounds towards +inf,
  !!                        Negative x_value rounds towards -inf
  !!  @param[out]  frac_x   Fractional part of x
  !--------------------------------------------------------------------------------
  subroutine frac_and_int_part(x_value,int_x,frac_x)

    implicit none

    real(kind=r_def), intent(in)  :: x_value
    integer, intent(out)          :: int_x
    real(kind=r_def), intent(out) :: frac_x

    frac_x = abs(x_value - int(x_value))
    int_x = abs(int(x_value))+1

  end subroutine frac_and_int_part


  !--------------------------------------------------------------------------------
  !>  @brief  Populates an array which contains the indices of the cells to sum
  !!          for the Courant number unlimited split advection scheme.
  !!
  !!  @param[in]   ncells           Number of cells to sum
  !!  @param[out]  index_array      Indices of cells to sum
  !!  @param[in]   departure_dist   Departure distance
  !!  @param[in]   edge_option      Cell edge option for determining fluxes
  !--------------------------------------------------------------------------------
  subroutine populate_array(ncells,index_array,departure_dist,edge_option)

    implicit none

    integer, intent(in)           :: ncells
    integer, intent(out)          :: index_array(ncells)
    real(kind=r_def), intent(in)  :: departure_dist
    integer, intent(in)           :: edge_option

    real(kind=r_def) :: frac_distance
    integer          :: ii
    integer          :: int_distance
    integer          :: istep, istart

    call frac_and_int_part(departure_dist,int_distance,frac_distance)

    if (departure_dist < 0.0_r_def) then
      istart = 0
      istep = 1
    else
      istart = -1
      istep = -1
    end if

    do ii = 1, ncells
      index_array(ii) = istart + istep*(ii-1) + edge_option
    end do

  end subroutine populate_array


  !--------------------------------------------------------------------------------
  !>  @brief  Returns the limits of integration of the "part" cell which comes from
  !!          the fractional part of the departure distance.
  !!          The values of x_left_limit and x_right_limit satisfy
  !!          0 <= x_left_limit, x_right_limit <= 1.
  !!
  !!  @param[in]   departure_dist   Departure distance
  !!  @param[in]   frac_x           Fractional part of departure distance
  !!  @param[out]  x_left_limit     Left limit of integral part of mass sum
  !!  @param[out]  x_right_limit    Right limit of integral part of mass sum
  !--------------------------------------------------------------------------------
  subroutine calc_integration_limits(departure_dist,frac_x,x_left_limit,x_right_limit)

    implicit none

    real(kind=r_def), intent(in)  :: departure_dist
    real(kind=r_def), intent(in)  :: frac_x
    real(kind=r_def), intent(out) :: x_left_limit
    real(kind=r_def), intent(out) :: x_right_limit

    if (departure_dist > 0.0_r_def) then
      x_left_limit = 1.0_r_def-frac_x
      x_right_limit = 1.0_r_def
    else
      x_left_limit = 0.0_r_def
      x_right_limit = frac_x
    end if

  end subroutine calc_integration_limits


  !--------------------------------------------------------------------------------
  !>  @brief  This will return the index assuming that we have the indexing of cells
  !!          with the form -4 | -3 | -2 | -1 | 0 | 1 | 2 | 3 | 4
  !!          and will return the index for the ordering assuming the form
  !!          1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9
  !!          Note that the stencil_length must of odd order
  !!
  !!  @param[in]     ii              Index value, takes integer values in the
  !!                          range [-(stencil_length-1)/2,(stencil_length-1)/2]
  !!  @param[in]     stencil_length  Length of the stencil
  !!  @param[return] cell_index      Index value out, between 1 and stencil_length
  !--------------------------------------------------------------------------------
  function map_cell_index(ii,stencil_length) result(cell_index)

    implicit none

    integer, intent(in) :: ii
    integer, intent(in) :: stencil_length
    integer             :: cell_index

    cell_index = ii+int((stencil_length+1)/2)

  end function map_cell_index


  !--------------------------------------------------------------------------------
  !>  @brief  Function which returns the mass from the fractional cell.
  !!
  !!  @param[in]      n_coeffs        number of coefficients for subgrid approximation
  !!  @param[in]      subgrid_coeffs  coefficients which approximate subgrid rho
  !!  @param[in]      x_left_limit    left integration limit
  !!  @param[in]      x_right_limit   right integration limit
  !!  @param[return]  part_mass       integrated mass
  !--------------------------------------------------------------------------------
  function return_part_mass(n_coeffs,subgrid_coeffs,x_left_limit,x_right_limit) result(part_mass)

    implicit none

    integer, intent(in) :: n_coeffs
    real(kind=r_def), intent(in)    :: subgrid_coeffs(1:n_coeffs)
    real(kind=r_def), intent(in)    :: x_left_limit
    real(kind=r_def), intent(in)    :: x_right_limit
    real(kind=r_def) :: part_mass

    part_mass = eval_integral(n_coeffs,subgrid_coeffs,x_right_limit) - &
                            eval_integral(n_coeffs,subgrid_coeffs,x_left_limit)

  end function return_part_mass


  !--------------------------------------------------------------------------------
  !>  @brief  Returns the value of the function rho(x) = a0*x+(1/2)*a1*x^2+(1/3)*a2*x^3
  !!
  !!  @param[in]      n_coeffs        Length of the subgrid coeffs array
  !!  @param[in]      subgrid_coeffs  Real array of length n_coeffs containing (a0,a1,a2)
  !!  @param[in]      xx              Ralue at which to evaluate the function
  !!  @param[return]  func_at_xx      Function evaluated at xx
  !--------------------------------------------------------------------------------
  function eval_integral(n_coeffs,subgrid_coeffs,xx) result(func_at_xx)

    implicit none

    integer, intent(in) :: n_coeffs
    real(kind=r_def), intent(in)    :: subgrid_coeffs(1:n_coeffs)
    real(kind=r_def), intent(in)    :: xx
    real(kind=r_def)                :: func_at_xx

    func_at_xx = subgrid_coeffs(1)*xx + 0.5_r_def*subgrid_coeffs(2)*xx**2 + &
                                  (1.0_r_def/3.0_r_def)*subgrid_coeffs(3)*xx**3

  end function eval_integral


  !--------------------------------------------------------------------------------
  !>  @brief  The ordering of the 1D stencils as defined in stencil_dofmap_mod.F90
  !!          is of the form | 4 | 3 | 2 | 1 | 5 | 6 | 7 |
  !!          If the integer input to this routine is 7 then the integer array
  !!          which is returned is (/ 4, 3, 2, 1, 5, 6, 7 /)
  !!
  !!  @param[in]   stencil_length     The length of the stencil
  !!  @param[out]  stencil_order_out  An integer array
  !--------------------------------------------------------------------------------
  subroutine calc_stencil_ordering(stencil_length,stencil_order_out)

    implicit none

    integer, intent(in)   :: stencil_length
    integer, intent(out)  :: stencil_order_out(1:stencil_length)

    integer :: ii, n

    n = (stencil_length+1)/2
    stencil_order_out(n) = 1

    do ii = n+1, stencil_length
      stencil_order_out(ii) = ii
    end do

    do ii = n-1, 1, -1
      stencil_order_out(ii) = stencil_order_out(ii+1)+1
    end do

  end subroutine calc_stencil_ordering


  !--------------------------------------------------------------------------------
  !>  @brief  Returns the ordering of the W3 dofs in a stencil dependent on cell
  !!          orientation and direction
  !!
  !!  @param[in]   stencil_length     The length of the stencil
  !!  @param[in]   orientation        The orientation of cell 1
  !!  @param[in]   direction          The direction of cosmic update
  !!  @param[out]  stencil_order_out  An integer array
  !--------------------------------------------------------------------------------
  subroutine stencil_ordering_and_orientation(stencil_length,orientation,direction,stencil_order_out)

    use flux_direction_mod,     only : x_direction, y_direction
    use log_mod,                only : log_event, LOG_LEVEL_ERROR, log_scratch_space, LOG_LEVEL_INFO

    implicit none

    integer, intent(in)   :: stencil_length
    integer, intent(in)   :: orientation
    integer, intent(out)  :: stencil_order_out(1:stencil_length)
    integer, intent(in)   :: direction

    integer :: ii, n

    ! Further details of this routine can be found in Ticket #868

    n = (stencil_length+1)/2
    stencil_order_out(n) = 1

    if (direction == x_direction ) then
      if (orientation == 1 .or. orientation == 2) then

        do ii = n+1, stencil_length
          stencil_order_out(ii) = ii
        end do

        do ii = n-1, 1, -1
          stencil_order_out(ii) = stencil_order_out(ii+1)+1
        end do

      elseif  (orientation == 3 .or. orientation == 4) then

        do ii = n+1, stencil_length
          stencil_order_out(ii) = ii - (n-1)
        end do

        do ii = n-1, 1, -1
          stencil_order_out(ii) = stencil_order_out(stencil_length+1 - ii) + n-1
        end do

      end if
    elseif (direction == y_direction ) then
      if (orientation == 1 .or. orientation == 4) then

        do ii = n+1, stencil_length
          stencil_order_out(ii) = ii
        end do

        do ii = n-1, 1, -1
          stencil_order_out(ii) = stencil_order_out(ii+1)+1
        end do

      elseif  (orientation == 2 .or. orientation == 3) then

        do ii = n+1, stencil_length
          stencil_order_out(ii) = ii - (n-1)
        end do

        do ii = n-1, 1, -1
          stencil_order_out(ii) = stencil_order_out(stencil_length+1 - ii) + n-1
        end do

      end if
    end if

  end subroutine stencil_ordering_and_orientation


  !--------------------------------------------------------------------------------
  !>  @brief  Returns the W2 dof to update dependent on cell orientation. If the
  !!          cell has orientation 1 then w2_dof = local_dof but if
  !!          cell has orientation 2 then local_dof = 1 returns w2_dof = 2
  !!          due to the orientation.
  !!
  !!  @param[in]   orientation    Orientation of the cell
  !!  @param[in]   local_dof      W2 dof assuming the cell has orientation
  !!  @param[out]  w2_dof         Actual W2 dof given the cell's orientation
  !--------------------------------------------------------------------------------
  function w2_dof(orientation,local_dof)

    implicit none

    integer             :: w2_dof
    integer, intent(in) :: local_dof, orientation

    if (local_dof >= 5) then
      w2_dof = local_dof
    else
      w2_dof = mod(orientation+local_dof-2,4)+1
    end if

  end function



  !--------------------------------------------------------------------------------
  !>  @brief  Re-orientates a cell of a W2 field given the cells orientation such
  !!          that the new field has orientation 1
  !!
  !!  @param[in]   orientation            Orientation of the cell
  !!  @param[in]   w2field                Horizontal dofs at lowest order for a W2 cell
  !!  @param[out]  reorientate_w2field    Horizontal dofs assuming orientation 1
  !--------------------------------------------------------------------------------
  function reorientate_w2field(w2field,orientation)

    implicit none

    integer :: orientation
    real(kind=r_def) :: w2field(1:4)
    real(kind=r_def) :: reorientate_w2field(1:4)

    integer :: local_dof

    do local_dof = 1, 4
      reorientate_w2field(local_dof) = w2field(w2_dof(orientation,local_dof))
    end do

  end function


  !--------------------------------------------------------------------------------
  !>  @brief  Returns the local W2 dof for opposite cell faces in either the x or
  !.          y-direction.
  !>          This is to handle cells in the halo which may have different
  !>          orientation due to orientation of panels in the cubed-sphere.
  !!
  !!  @param[in]   orientation    Orientation of the cell
  !!  @param[in]   direction      Either x or y horizontal directions
  !!  @param[out]  w2_dofs        Local W2 dof for opposite faces
  !--------------------------------------------------------------------------------
  function dof_to_update(orientation,direction) result(w2_dofs)
    use flux_direction_mod, only : x_direction, y_direction
    use log_mod,            only : log_event, LOG_LEVEL_ERROR

    implicit none

    integer               :: w2_dofs(1:2)
    integer, intent(in)   :: orientation
    integer, intent(in)   :: direction

    select case(direction)
      case(x_direction)
        w2_dofs = (/ orientation, mod(orientation+1,4)+1 /)
      case(y_direction)
        w2_dofs = (/ mod(orientation,4)+1, mod(orientation+2,4)+1 /)
      case default
        call log_event( "Error: direction not defined in dof_to_update", LOG_LEVEL_ERROR )
    end select

  end function

  !--------------------------------------------------------------------------------
  !>  @brief  Returns the layers which to sum in the vertical for calculating the
  !>          vertical mass flux in Cosmic. The Cosmic transport scheme is
  !>          designed to work for Courant numbers greater than 1.
  !!
  !!  @param[out]  local_index        Indexes of cells to sum to calculate mass flux
  !!  @param[in]   departure_dist     Departure distance
  !!  @param[in]   n_cells_to_sum     Number of cells to sum
  !!  @param[in]   k                  Cell for which mass flux is calculated for
  !!  @param[in]   nlayers            Number of vertical levels
  !!  @param[in]   edge               Each cell has two edges in the vertical,
  !!                                  edge=1 denotes the top edge, edge=0 denotes
  !!                                  the bottom edge
  !--------------------------------------------------------------------------------
  subroutine calc_local_vertical_index(local_index,departure_dist,n_cells_to_sum,k,nlayers,edge)

    use log_mod, only : log_event, LOG_LEVEL_ERROR, log_scratch_space

    implicit none

    integer, intent(in)           :: n_cells_to_sum
    integer, intent(in)           :: k
    integer, intent(in)           :: nlayers
    integer, intent(in)           :: edge
    integer                       :: local_index(n_cells_to_sum)
    real(kind=r_def), intent(in)  :: departure_dist

    integer :: jj

    if (departure_dist>0.0) then
      do jj=0,n_cells_to_sum-1
        local_index(jj+1) = -jj-1+k+edge
      end do
    else
      do jj=0,n_cells_to_sum-1
        local_index(jj+1) = jj+k+edge
      end do
    end if

    ! Now check that the values are in the limits [0,nlayers-1]
    do jj=1,n_cells_to_sum
      if (local_index(jj)<0 .or. local_index(jj)>(nlayers-1)) then
        write( log_scratch_space, '( A, E12.4E3, A, I6, A, 10I6 )' ) &
           "ERROR: Departure distance", departure_dist,              &
           " is too big. Layer of interest is ",k ,                  &
           " Index values are ", local_index
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
    end do

  end subroutine calc_local_vertical_index

end module cosmic_flux_mod
