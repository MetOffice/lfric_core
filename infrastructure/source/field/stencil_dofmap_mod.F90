!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief A type which holds information about the dofmap.
!> @details Types for storing the general stencil dofmaps of a general size
!> Inherits from linked_list_data_type so it can be stored in a linked_list
!> Allowable stencil types (the size is variable):
!>
!> \verbatim
!> POINT  -->       |1|
!>
!> 1DX    -->   |3|2|1|4|5|
!>
!>                  |5|
!>                  |4|
!> 1DY    -->       |1|
!>                  |2|
!>                  |3|
!>
!>                  |9|
!>                  |8|
!> CROSS  -->   |3|2|1|6|7|
!>                  |4|
!>                  |5|
!>
!>                |9|8|7|
!> REGION -->     |2|1|6|
!>                |3|4|5|
!> \endverbatim

module stencil_dofmap_mod

use constants_mod,        only: i_def, l_def
use master_dofmap_mod,    only: master_dofmap_type
use linked_list_data_mod, only: linked_list_data_type


implicit none

private
public :: generate_stencil_dofmap_id

type, extends(linked_list_data_type), public :: stencil_dofmap_type
  private
  integer(i_def) :: dofmap_shape
  integer(i_def) :: dofmap_extent
  integer(i_def) :: dofmap_size
  integer(i_def), allocatable :: local_size(:)
  integer(i_def), allocatable :: dofmap(:,:,:)
contains
  procedure :: get_dofmap
  procedure :: get_whole_dofmap
  procedure :: get_size
  procedure :: get_local_size
  procedure :: get_stencil_sizes

  !> Forced clear of this object from memory.
  !> This routine should not need to be called manually except
  !> (possibly) in pfunit tests
  procedure, public :: clear

  !> Finalizer routine, should be called automatically
  !> by code when the object is out of scope
  final :: stencil_dofmap_destructor

end type stencil_dofmap_type

integer(i_def), public, parameter :: STENCIL_POINT  = 1100
integer(i_def), public, parameter :: STENCIL_1DX    = 1200
integer(i_def), public, parameter :: STENCIL_1DY    = 1300
integer(i_def), public, parameter :: STENCIL_CROSS  = 1400
integer(i_def), public, parameter :: STENCIL_REGION = 1500

interface stencil_dofmap_type
  module procedure stencil_dofmap_constructor
end interface

contains

!-----------------------------------------------------------------------------
! Construct the stencil dofmap
!-----------------------------------------------------------------------------
!> Function to construct a stencil dofmap
!> @param[in] st_shape The shape of the required stencil
!> @param[in] st_depth Parameter to control how far the stencil is computed to
!> @param[in] ndf Number of degrees of freedom per cell for this function space
!> @param[in] mesh Object to build the stencil on
!> @param[in] master_dofmap The cell dofmap to create the stencil from
!> @return The dofmap object
function stencil_dofmap_constructor(st_shape, st_depth, ndf, mesh, master_dofmap) result(self)

  use log_mod,               only: log_event,         &
                                   log_scratch_space, &
                                   LOG_LEVEL_ERROR
  use mesh_mod,              only: mesh_type
  use reference_element_mod, only: W, E, N, S, &
                                   reference_element_type

  implicit none

  integer(i_def),           intent(in) :: st_shape
  integer(i_def),           intent(in) :: st_depth
  integer(i_def),           intent(in) :: ndf
  type(mesh_type), pointer, intent(in) :: mesh
  type(master_dofmap_type), intent(in) :: master_dofmap
  type(stencil_dofmap_type)            :: self

  integer(i_def) :: cell
  integer(i_def) :: ncells
  integer(i_def) :: st_size
  integer(i_def), pointer :: map(:) => null()
  integer(i_def) :: cells_in_stencil
  integer(i_def) :: last_halo_index
  integer(i_def) :: stencil_dofmap_id
  integer(i_def) :: i
  integer(i_def), allocatable, dimension(:) :: direction_map
  integer(i_def), allocatable, dimension(:) :: stencil_cells
  integer(i_def) :: number_of_neighbours

  logical(l_def) :: region

  class(reference_element_type), pointer :: reference_element => null()

  reference_element => mesh%get_reference_element()
  number_of_neighbours = reference_element%get_number_2d_edges()

  ! Since this routine is only valid for quadrilateral elements throw an error
  ! if the number of (horizontal) neighbours is not four
  if ( number_of_neighbours /= 4_i_def ) &
  call log_event( 'Stencil dofmaps only valid for quad elements', LOG_LEVEL_ERROR )

  ! Get direction map for the requested stencil shape
  allocate( direction_map(number_of_neighbours) )
  call get_direction_map( st_shape, st_depth, direction_map )

  ! Compute the maximum number of cells in a stencil
  if ( st_shape == STENCIL_REGION ) then
    st_size = ((2*st_depth) + 1)**2
    region = .true.
  else
    st_size = sum(direction_map) + 1
    region = .false.
  end if

  ! Allocate stencil cells array to be of max possible stencil size
  allocate( stencil_cells(st_size) )

  self%dofmap_shape  = st_shape
  self%dofmap_extent = st_depth
  self%dofmap_size   = st_size

  stencil_dofmap_id = generate_stencil_dofmap_id(st_shape, st_depth)

  ! call base class set_id()
  call self%set_id(stencil_dofmap_id)

  ! How far can stencil be computed before stencil goes off edge of halo region
  last_halo_index = mesh%get_halo_depth() - st_depth
  if (last_halo_index < 0) then
    write( log_scratch_space, '( A, I4, A, I4, A, I4 )' ) &
       'Attempting to create stencil: ', st_shape,' of extent ',st_depth, &
       ' when halo is depth is too small:',mesh%get_halo_depth()
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  else if (last_halo_index == 0) then
    ! Stencil extent same as halo depth, so compute stencil for all owned cells
    ncells = mesh%get_last_edge_cell()
  else
    ! Stencil extent smaller than halo depth, so compute stencil into halo
    ncells = mesh%get_last_halo_cell(last_halo_index)
  end if

  ! Allocate the dofmap array for maximum stencil size
  allocate( self%dofmap( ndf, self%dofmap_size, ncells ) )
  ! Allocate array for storing individual stencil sizes
  allocate( self%local_size(ncells) )

  ! Loop over all cells in mesh and calculate stencil dofmap for each
  do cell = 1, ncells
    ! Set default value for stencil cells
    stencil_cells = -1_i_def
    ! Add current cell to stencil and update number of cells in stencil
    stencil_cells(1) = cell
    cells_in_stencil = 1
    ! Call subroutine to get the rest of the cells in the stencil
    call get_stencil_cells( mesh, cell, st_depth, number_of_neighbours, &
                            direction_map, region, cells_in_stencil, &
                            stencil_cells )
    ! Set actual stencil size for cell
    self%local_size(cell) = cells_in_stencil
    ! Create pointer to master dofmap of each cell in stencil and add to
    ! local stencil dofmap
    do i=1, cells_in_stencil
      map => master_dofmap%get_master_dofmap(stencil_cells(i))
      self%dofmap(:,i,cell) = map(:)
    end do
  end do

  ! Clean up arrays and pointers
  deallocate( direction_map, stencil_cells )
  nullify(map, reference_element)

  end function stencil_dofmap_constructor

!-----------------------------------------------------------------------------
! Get the stencil dofmap for a single cell
!-----------------------------------------------------------------------------
!> Function Returns a pointer to the dofmap for the cell
!! @param[in] self The calling function_space
!! @param[in] cell Which cell
!! @return The pointer which points to a slice of the dofmap
function get_dofmap(self,cell) result(map)
  implicit none
  class(stencil_dofmap_type), target, intent(in) :: self
  integer(i_def),                     intent(in) :: cell
  integer(i_def), pointer                        :: map(:,:)

  map => self%dofmap(:,:,cell)
  return
end function get_dofmap
!-----------------------------------------------------------------------------
! Get the stencil dofmap for the whole mesh
!-----------------------------------------------------------------------------
!> Function Returns a pointer to the dofmap for the mesh
!! @param[in] self The calling function_space
!! @return The pointer which points to the dofmap
function get_whole_dofmap(self) result(map)
  implicit none
  class(stencil_dofmap_type), target, intent(in) :: self
  integer(i_def), pointer                        :: map(:,:,:)

  map => self%dofmap(:,:,:)
  return
end function get_whole_dofmap

!> Returns the size of the stencil in cells
!! @param[in] self The calling function_space
!! @return The size of the stencil in cells
function get_size(self) result(stencil_size)
  implicit none
  class(stencil_dofmap_type), target, intent(in) :: self
  integer(i_def)                                 :: stencil_size

  stencil_size = self%dofmap_size
  return
end function get_size

!> Returns the size of the stencil in cells for this cell
!! @param[in] self The calling function_space
!! @param[in] cell The cell at the centre of the stencil
!! @return The size of the stencil in cells
function get_local_size(self, cell) result(stencil_size)
  implicit none
  class(stencil_dofmap_type), target, intent(in) :: self
  integer(i_def),                     intent(in) :: cell
  integer(i_def)                                 :: stencil_size

  if ( self%dofmap_shape == STENCIL_REGION ) then
    stencil_size = self%local_size(cell)
  else
    stencil_size = self%dofmap_size
  end if
  return
end function get_local_size

!-----------------------------------------------------------------------------
! Get the local stencil sizes
!-----------------------------------------------------------------------------
!> Returns a pointer to the stencil local_size array
!! @param[in] self The calling space function
!! @return Pointer to the array of local stencil sizes
! TODO - Ticket #2019 Design discussion on how stencil sizes are passed to
!        kernels.
function get_stencil_sizes(self) result(stencil_sizes)
  use log_mod,                only: log_event, &
                                    log_scratch_space, &
                                    LOG_LEVEL_ERROR
  implicit none
  class(stencil_dofmap_type), target, intent(in) :: self
  integer(i_def), pointer                        :: stencil_sizes(:)


    stencil_sizes => self%local_size(:)

  return
end function get_stencil_sizes

!> Returns a stencil dofmap id using stencil shape and extent
!> @param[in] stencil_shape   Shape code of stencil
!> @param[in] stencil_extent  Extent of stencil
!> @return    stencil_id
!==============================================================================
function generate_stencil_dofmap_id( stencil_shape,   &
                                     stencil_extent ) &
                                     result( stencil_id )

  implicit none

  integer(i_def), intent(in) :: stencil_shape
  integer(i_def), intent(in) :: stencil_extent

  integer(i_def) :: stencil_id

  stencil_id = stencil_shape*100 + stencil_extent

  return
end function generate_stencil_dofmap_id

!> @brief Returns an array of values for depth of stencil in the different
!> directions
!> @param[in] st_shape The shape of the required stencil
!> @param[in] st_depth The extent of the stencil
!> @param[inout] direction_map An array of depths for each direction (W,E,N,S)
subroutine get_direction_map(st_shape, st_depth, direction_map)

  use log_mod,               only: log_event,         &
                                   log_scratch_space, &
                                   LOG_LEVEL_ERROR
  use reference_element_mod, only: W, E, N, S

  implicit none

  integer(i_def),              intent(in) :: st_shape
  integer(i_def),              intent(in) :: st_depth
  integer(i_def), allocatable, intent(inout) :: direction_map(:)

  ! Set direction map to depth zero in all directions
  direction_map = 0

  ! Set direction map depth in different directions depending on stencil shape
  select case ( st_shape )
  case ( STENCIL_POINT)
  case ( STENCIL_1DX )
    direction_map(W) = 1
    direction_map(E) = 1
  case ( STENCIL_1DY )
    direction_map(S) = 1
    direction_map(N) = 1
  case ( STENCIL_CROSS, STENCIL_REGION )
    direction_map(W) = 1
    direction_map(E) = 1
    direction_map(S) = 1
    direction_map(N) = 1
  case default
    ! Throw error if stencil shape invalid
    write( log_scratch_space, '( A, I4 )' ) &
       "Invalid stencil type: ", st_shape
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

  ! Multiply direction map values by the requested stencil depth
  direction_map = direction_map * st_depth

end subroutine get_direction_map

!> @brief Get the cell IDs for all the cells included in the stencil
!> @param[in] mesh Mesh object to find cells from
!> @param[in] origin Starting cell for stencil (centre cell)
!> @param[in] st_depth Stencil depth
!> @param[in] number_of_neighbours Number of cell neighbours
!> @param[in] direction_map Array of direction sizes for stencil
!> @param[in] region Boolean option for creating a region stencil
!> @param[inout] cells_in_stencil Number of cells currently in stencil
!> @param[inout] stencil_cells Array of cell IDs for stencil
!==============================================================================
subroutine get_stencil_cells( mesh,                  &
                              origin,                &
                              st_depth,              &
                              number_of_neighbours,  &
                              direction_map,         &
                              region,                &
                              cells_in_stencil,      &
                              stencil_cells )

  use mesh_mod,              only: mesh_type
  use master_dofmap_mod,     only: master_dofmap_type

  implicit none

  type(mesh_type), pointer, intent(in) :: mesh
  integer(i_def), intent(in) :: origin
  integer(i_def), intent(in) :: st_depth
  integer(i_def), intent(in) :: number_of_neighbours
  integer(i_def), intent(in), dimension(number_of_neighbours) :: direction_map
  logical(l_def), intent(in) :: region
  integer(i_def), intent(inout) :: cells_in_stencil
  integer(i_def), allocatable, intent(inout) :: stencil_cells(:)

  integer(i_def), dimension(number_of_neighbours) :: direction_map_rotate
  integer(i_def) :: i, j, k
  integer(i_def) :: cell
  integer(i_def) :: new_cell
  integer(i_def) :: branch_cell
  integer(i_def) :: branch_new_cell
  integer(i_def) :: direction
  integer(i_def) :: direction_branch

  do i=1, number_of_neighbours
    ! Starting cell is origin
    cell = origin
    direction = i
    if (direction_map(direction) > 0) then
      do j=1, direction_map(direction)
        new_cell = mesh%get_cell_next(direction, cell)
        ! Only add the cell to the stencil_cells list if it is not already
        ! added and is greater than 0
        if (.not.(any(stencil_cells == new_cell)) .and. new_cell > 0) then
          cells_in_stencil = cells_in_stencil + 1
          stencil_cells(cells_in_stencil) = new_cell
        end if
        ! If new_cell ID is greater than 0 then find the correct direction
        ! This is to make sure the direction is correct when moving panels
        if (new_cell > 0) then
          direction = direction_rotation( mesh, number_of_neighbours, &
                                          new_cell, cell, 180_i_def )
          ! If this is a region stencil then find the cells adjacent to the
          ! new_cell using the stencil depth and a direction 90 degrees from
          ! the direction of the origin cell
          if(region) then
            branch_cell = new_cell
            ! Get direction 90 degrees from cell
            direction_branch = direction_rotation( mesh, number_of_neighbours, &
                                                   branch_cell, cell, -90_i_def )
            do k=1, st_depth
              branch_new_cell = mesh%get_cell_next(direction_branch, branch_cell)
              ! Only add the cell to the stencil_cells list if it isn't already
              ! added and is greater than 0
              if ( .not.(any(stencil_cells == branch_new_cell)) &
                   .and. branch_new_cell > 0 ) then
                cells_in_stencil = cells_in_stencil + 1
                stencil_cells(cells_in_stencil) = branch_new_cell
              end if
              ! If new_cell ID is greater than 0 then find the correct direction
              ! This is to make sure the direction is correct when moving panels
              if (branch_new_cell > 0) then
                direction_branch = direction_rotation( mesh, &
                                                       number_of_neighbours, &
                                                       branch_new_cell, &
                                                       branch_cell, 180_i_def )
                branch_cell = branch_new_cell
              end if
            end do
          end if
          ! new_cell is now origin cell for next step
          cell = new_cell
        end if
      end do
    end if
  end do

end subroutine get_stencil_cells


!> @brief Returns a rotated direction from the origin_cell
!> @param[in] mesh Pointer to mesh
!> @param[in] number_of_neighbours Number of neighbouring cells
!> @param[in] current_cell Current cell to find direction from
!> @param[in] origin_cell Origin cell for start direction
!> @param[in] rotation Value for rotation
!> @return rotated_direction The rotated direction from the origin_cell
!==============================================================================
function direction_rotation( mesh,   &
                             number_of_neighbours, &
                             current_cell, &
                             origin_cell, &
                             rotation ) &
                             result( rotated_direction )

  use log_mod,               only: log_event,         &
                                   log_scratch_space, &
                                   LOG_LEVEL_ERROR
  use mesh_mod,              only: mesh_type
  use reference_element_mod, only: W, E, N, S

  implicit none

  type(mesh_type), pointer, intent(in) :: mesh
  integer(i_def), intent(in) :: number_of_neighbours
  integer(i_def), intent(in) :: current_cell
  integer(i_def), intent(in) :: origin_cell
  integer(i_def), intent(in) :: rotation

  integer(i_def) :: rotated_direction
  integer(i_def) :: direction
  integer(i_def), dimension(4) :: direction_ref = (/ W, S, E, N /)
  integer(i_def), dimension(4) :: rotated

  ! Rotate the array using cshift to change the order of the direction map
  select case (rotation)
    case (90_i_def, -270_i_def)  ! North becomes East
      rotated = cshift(direction_ref, shift = 1)
    case (-90_i_def, 270_i_def)  ! North becomes West
      rotated = cshift(direction_ref, shift = -1)
    case (180_i_def, -180_i_def) ! North becomes south
      rotated = cshift(direction_ref, shift = 2)
    case default
      write( log_scratch_space, '( A, I3 )' ) &
         'Invalid rotation value', rotation
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

  ! Find direction to origin cell and return
  do direction=1, number_of_neighbours
    if (mesh%get_cell_next(direction, current_cell) == origin_cell) then
      rotated_direction = rotated(direction)
    end if
  end do

  return
end function direction_rotation

!==============================================================================
!> @brief Routine to destroy object.
!> @param[in] self, The calling object type
subroutine clear(self)

  implicit none

  class (stencil_dofmap_type), intent(inout) :: self

  if (allocated(self%dofmap)) deallocate( self%dofmap)

  return
end subroutine clear

!==============================================================================
!> @brief Finalizer routine which should automatically call clear
!>        when object is out of scope.
!> @param[in] self, The calling object_type
subroutine stencil_dofmap_destructor(self)

  implicit none

  type (stencil_dofmap_type), intent(inout) :: self

  call self%clear()

  return
end subroutine stencil_dofmap_destructor

end module stencil_dofmap_mod

