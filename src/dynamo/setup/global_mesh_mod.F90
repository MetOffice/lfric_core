!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief  Describes the cell ordering within the global mesh

!> @details This object holds the connectivities that fully
!>          describe the 2D topology of the global mesh

module global_mesh_mod

use constants_mod, only: r_def, i_def
implicit none

private

type, public :: global_mesh_type
  private
!> Horizontal coords of vertices in full domain
  real(kind=r_def), allocatable :: vert_coords(:,:)
!> Full domain cell to cell connectivities
  integer, allocatable :: cell_next_2d(:,:)
!> Full domain vertices on a cell
  integer, allocatable :: vert_on_cell_2d(:,:)
!> Full domain cells that surround a vertex
  integer, allocatable :: cell_on_vert_2d(:,:)
!> Full domain edges on a cell
  integer, allocatable :: edge_on_cell_2d(:,:)
!> Full domain cells either side of an edge
  integer, allocatable :: cell_on_edge_2d(:,:)
!> Total number of vertices in the full domain
  integer              :: nverts
!> Total number of edges in the full domain
  integer              :: nedges
!> total number of cells in full domain
  integer              :: ncells
!> number of vertices on each cell
  integer              :: nverts_per_cell
!> number of edges on each cell
  integer              :: nedges_per_cell
!> maximum number of cells around a vertex
  integer              :: max_cells_per_vertex

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains
  !> @brief  Returns id of a cell in a mesh
  !> @details Gets the cell id that is x_cells across (E is +ve) and
  !> y_cells up/down (N is +ve) from given cell_number.
  !> NOTE: For a cubed -sphere mesh, this will only return correct cell ids if
  !> the offset cell remains on same the cubed-sphere "face" as the start cell
  !> @param[in] cell_number Start position in the global cell id array
  !> @param[in] x_cells Offset in the E/W direction
  !> @param[in] y_cells Offset in the N/S direction
  !> @return cell_id The cell id of the cell at the given offset to the start cell
  procedure, public :: get_cell_id

  !> @brief  Gets the cells that are incident on a particular vertex
  !> @param[in] vertex_number The number of the vertex being queried
  !> @param[out] cells The cells around the given vertex
  procedure, public :: get_cell_on_vert 

  !> @brief  Gets the cells that are incident on a particular edge
  !> @param[in] edge_number The number of the edge being queried
  !> @param[out] cells The cells either side of the given edge
  procedure, public :: get_cell_on_edge

  !> @brief  Gets the total number of cells in the global domain
  !> @return ncells The total number of cells in the global domain
  procedure, public :: get_ncells

  !> @brief  Returns maximum number of cells around a vertex.
  !> @details The maximum number of cells that can be incident with a vertex (the
  !> actual number of cells at a particular vertex could be less (eg. in a cubed
  !> sphere mesh, there are generally four cells incident with a vertex except
  !> for the "corner" vertices where there are three.
  !> @return max_cells_per_vertex The maximum number of cells that can be incident with a vertex
  procedure, public :: get_max_cells_per_vertex

  !> @brief  Gets the edges that are incident with a particular cell
  !> @param[in]  cell_gid The global id of the cell being queried 
  !> @param[out] The edges around the given cell
  procedure, public :: get_edge_on_cell

  !> @brief  Gets the vertices that are incident with a particular cell
  !> @param[in]  cell_gid The global id of the cell being queried 
  !> @param[out] The vertices around the given cell
  procedure, public :: get_vert_on_cell

  !> @brief  Returns number of vertices per 2D-cell
  !> @return nverts_per_cell Number of vertices per 2D-cell
  procedure, public :: get_nverts_per_cell

  !> @brief  Returns number of edges on each cell
  !> @return nedges_per_cell Number of edges per 2D-cell
  procedure, public :: get_nedges_per_cell

  !> @brief  Returns cells adjacent to a cell
  !> @param [in]  cell_gid   The global id of the requested cell
  !> @return An array containing the global ids of cells adjacent to
  !>         the cell with global id [cell_gid]
  procedure, public :: get_cell_next

  !> @brief  Returns vertex coordinates
  !> @param [in]  vert_gid   The global id of the requested vertex
  !> @return A three-element array containing the coordinates of a 
  !>         single vertex on the global mesh. Currently, these are in
  !>         spherical coords [long,lat,radius] 
  procedure, public :: get_vert_coords

end type global_mesh_type

interface global_mesh_type
  module procedure global_mesh_constructor
  module procedure global_mesh_constructor_unit_test_data
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------------
! Constructs a global mesh object
!-------------------------------------------------------------------------------
!>  @brief   Constructs a global mesh object
!>  @details Constructs a global mesh object to hold the connectivities
!!           that fully describe the 2D topology of the mesh
!-------------------------------------------------------------------------------
function global_mesh_constructor( filename ) result(self)

use constants_mod,  only: str_def
use ugrid_2d_mod,   only: ugrid_2d_type
use ugrid_file_mod, only: ugrid_file_type
use ncdf_quad_mod,  only: ncdf_quad_type

implicit none

character(len = str_def), intent(in) :: filename

type(global_mesh_type) :: self

type(ugrid_2d_type) :: ugrid_2d

class(ugrid_file_type), allocatable :: file_handler

! dimensions from file
integer :: nvert_in
integer :: nface_in
integer :: nedge_in
integer :: num_nodes_per_face
integer :: num_nodes_per_edge
integer :: num_edges_per_face
integer :: max_num_faces_per_node

allocate( ncdf_quad_type :: file_handler )
call ugrid_2d%set_file_handler( file_handler )
call ugrid_2d%read_from_file( trim(filename) )

call ugrid_2d%get_dimensions( num_nodes              = nvert_in, &
                              num_edges              = nedge_in, &
                              num_faces              = nface_in, &
                              num_nodes_per_face     = num_nodes_per_face, &
                              num_edges_per_face     = num_edges_per_face, &
                              num_nodes_per_edge     = num_nodes_per_edge, &
                              max_num_faces_per_node = max_num_faces_per_node )

self%nverts = nvert_in
self%nedges = nedge_in
self%ncells = nface_in

self%nverts_per_cell      = num_nodes_per_face
self%nedges_per_cell      = num_edges_per_face
self%max_cells_per_vertex = max_num_faces_per_node

allocate( self%vert_coords(2, nvert_in) )
call ugrid_2d%get_node_coords(self%vert_coords)

allocate( self%cell_next_2d( num_edges_per_face, nface_in ) )
call ugrid_2d%get_face_face_connectivity( self%cell_next_2d )

allocate( self%vert_on_cell_2d( num_nodes_per_face, nface_in ) )
call ugrid_2d%get_face_node_connectivity( self%vert_on_cell_2d )

allocate( self%edge_on_cell_2d( num_edges_per_face, nedge_in ) )
call ugrid_2d%get_face_edge_connectivity( self%edge_on_cell_2d )

allocate( self%cell_on_vert_2d( self%max_cells_per_vertex, nvert_in ) )
call calc_cell_on_vertex( self%vert_on_cell_2d, &
                          num_nodes_per_face, &
                          nface_in, &
                          self%cell_on_vert_2d, &
                          self%max_cells_per_vertex, &
                          nvert_in)

! Populate cells either side of each edge
allocate( self%cell_on_edge_2d(2,nedge_in) )  ! There can only ever be 2 cells incident on an edge (whatever the topography!)
call calc_cell_on_edge( self%edge_on_cell_2d, &
                        num_edges_per_face, &
                        nface_in, &
                        self%cell_on_edge_2d, &
                        nedge_in )

end function global_mesh_constructor

!-------------------------------------------------------------------------------
! Returns the cells on vertices. PRIVATE subroutine.
!-------------------------------------------------------------------------------
! Details: Calculates the cells that are incident on a vertex by looping through
!          the vertices on all cells (which we store) and filling an array based 
!          on vertex number with the cells around it.
! Input:   vert_on_cell    Array with indices of vertices on cells
!          verts_per_cell  Number of vertices per cell
!          ncell           Number of cells
!          cells_per_vert  Number of cells per vertex
!          nvert           Number of vertices
! Output:  cell_on_vert    Array with indices of cells on vertices
!-------------------------------------------------------------------------------
subroutine calc_cell_on_vertex(vert_on_cell, &
                               verts_per_cell, &
                               ncell, &
                               cell_on_vert, &
                               cells_per_vert,  &
                               nvert)
implicit none

integer, intent(in)  :: verts_per_cell, ncell
integer, intent(in)  :: vert_on_cell(verts_per_cell, ncell)
integer, intent(in)  :: cells_per_vert, nvert
integer, intent(out) :: cell_on_vert(cells_per_vert, nvert)

integer :: cell
integer :: vertno
integer :: cellno
integer :: vert

cell_on_vert = 0

do cell = 1,ncell
  do vertno = 1,verts_per_cell

    vert = vert_on_cell(vertno,cell)

    do cellno = 1,cells_per_vert
      if(cell_on_vert(cellno,vert) == cell)exit
      if(cell_on_vert(cellno,vert) == 0)then
        cell_on_vert(cellno,vert) = cell
        exit
      end if
    end do
  end do
end do

end subroutine calc_cell_on_vertex

!-------------------------------------------------------------------------------
! Returns the cells on edges. PRIVATE subroutine.
!-------------------------------------------------------------------------------
! Details: Calculates the cells that are either side of an edge by looping
!          through the edges on all cells (which we store) and filling an array 
!          based on edge number with the cells around it
! Input:   edge_on_cell    Array with indices of edges on cells
!          edges_per_cell  Number of edges per cell
!          ncell           Number of cells
!          nedge           Number of edges
! Output:  cell_on_edge    Array with indices of cells on edges
!-------------------------------------------------------------------------------
subroutine calc_cell_on_edge( edge_on_cell, &
                              edges_per_cell, &
                              ncell, &
                              cell_on_edge, &
                              nedge)
implicit none

integer, intent(in)  :: edges_per_cell, ncell
integer, intent(in)  :: edge_on_cell(edges_per_cell, ncell)
integer, intent(in)  :: nedge
integer, intent(out) :: cell_on_edge(2, nedge)

integer :: cell
integer :: edgeno
integer :: cellno
integer :: edge

cell_on_edge=0

do cell=1,ncell
  do edgeno=1,edges_per_cell

    edge=edge_on_cell(edgeno,cell)

    do cellno=1, 2
      if(cell_on_edge(cellno,edge) == cell)exit
      if(cell_on_edge(cellno,edge) == 0)then
        cell_on_edge(cellno,edge)=cell
        exit
      end if
    end do
  end do
end do

end subroutine calc_cell_on_edge

!-------------------------------------------------------------------------------
! Returns id of a cell in a mesh
!-------------------------------------------------------------------------------
function get_cell_id( self, cell_number, x_cells, y_cells ) result ( cell_id )

use reference_element_mod, only : W, S, E, N

implicit none

class(global_mesh_type), intent(in) :: self

integer, intent(in) :: cell_number
integer, intent(in) :: x_cells, y_cells

integer :: cell_id

integer :: index, dist, i

cell_id=cell_number

if (x_cells > 0 )then
  index = E
  dist = x_cells
else if (x_cells < 0 )then
  index = W
  dist = abs(x_cells)
else
  index = W
  dist = 0
endif
do i = 1,dist
  cell_id = self%cell_next_2d(index,cell_id)
end do

if (y_cells > 0 )then
  index = N
  dist = y_cells
else if (y_cells < 0 )then
  index = S
  dist = abs(y_cells)
else
  index = S
  dist = 0
endif
do i = 1,dist
  cell_id = self%cell_next_2d(index,cell_id)
end do

end function get_cell_id

!-------------------------------------------------------------------------------
! Gets the cells that are incident on a particular vertex
!-------------------------------------------------------------------------------
subroutine get_cell_on_vert( self, vertex_number, cells)

implicit none

class(global_mesh_type), intent(in) :: self

integer, intent(in)  :: vertex_number
integer, intent(out) :: cells(:)

cells = self%cell_on_vert_2d(:,vertex_number)

end subroutine get_cell_on_vert

!-------------------------------------------------------------------------------
! Gets the cells that are incident on a particular edge
!-------------------------------------------------------------------------------
subroutine get_cell_on_edge( self, edge_number, cells)

implicit none

class(global_mesh_type), intent(in) :: self
integer, intent(in) :: edge_number
integer, intent(out) :: cells(:)

cells=self%cell_on_edge_2d(:,edge_number)

end subroutine get_cell_on_edge

!-------------------------------------------------------------------------------
! Gets the total number of cells in the global domain
!-------------------------------------------------------------------------------
function get_ncells( self ) result (ncells)

class(global_mesh_type), intent(in) :: self

integer :: ncells

ncells = self%ncells

end function get_ncells

!-------------------------------------------------------------------------------
! Gets maximum number of cells around a vertex
!-------------------------------------------------------------------------------
function get_max_cells_per_vertex( self ) result (max_cells_per_vertex)

class(global_mesh_type), intent(in) :: self

integer :: max_cells_per_vertex

max_cells_per_vertex = self%max_cells_per_vertex

end function get_max_cells_per_vertex

!-------------------------------------------------------------------------------
! Gets edges that are incident with a particular cell
!-------------------------------------------------------------------------------
subroutine get_edge_on_cell(self, cell_gid, edges)

  implicit none
  class (global_mesh_type), intent(in)  :: self
  integer (i_def),          intent(in)  :: cell_gid
  integer (i_def),          intent(out) :: edges(:)

  edges(:) = self%edge_on_cell_2d(:,cell_gid)

end subroutine get_edge_on_cell

!-------------------------------------------------------------------------------
! Gets the vertices that are incident with a particular cell
!-------------------------------------------------------------------------------
subroutine get_vert_on_cell(self, cell_gid, verts)

  implicit none
  class (global_mesh_type), intent(in)  :: self
  integer (i_def),          intent(in)  :: cell_gid
  integer (i_def),          intent(out) :: verts(:)

  verts(:) = self%vert_on_cell_2d(:,cell_gid)

end subroutine get_vert_on_cell

!-------------------------------------------------------------------------------
! Gets the number of vertices per 2D-cell
!-------------------------------------------------------------------------------
function get_nverts_per_cell( self ) result (nverts_per_cell)

  implicit none
  class(global_mesh_type), intent(in) :: self
  integer :: nverts_per_cell

  nverts_per_cell = self%nverts_per_cell

end function get_nverts_per_cell

!-------------------------------------------------------------------------------
! Gets the number of edges on each cell
!-------------------------------------------------------------------------------
function get_nedges_per_cell( self ) result (nedges_per_cell)

  implicit none
  class(global_mesh_type), intent(in) :: self
  integer                             :: nedges_per_cell

  nedges_per_cell = self%nedges_per_cell

end function get_nedges_per_cell

!-------------------------------------------------------------------------------
! Gets ids of cells adjacent to a cell with global id
!-------------------------------------------------------------------------------
subroutine get_cell_next (self, cell_gid, cell_next)

  implicit none
  class (global_mesh_type), intent(in)  :: self
  integer (i_def),          intent(in)  :: cell_gid
  integer(i_def),           intent(out) :: cell_next(:)

  cell_next(:) = self%cell_next_2d(:,cell_gid)

end subroutine get_cell_next

!-------------------------------------------------------------------------------
! Gets vertex coordinates
!-------------------------------------------------------------------------------
subroutine get_vert_coords (self, vert_gid, vert_coords)

  implicit none
  class (global_mesh_type), intent(in)  :: self
  integer,                  intent(in)  :: vert_gid
  real(r_def),              intent(out) :: vert_coords(:)

  vert_coords(1:2) = self%vert_coords(1:2,vert_gid)

end subroutine get_vert_coords


!==============================================================================
! The following routines are only available when setting data for unit testing.
!==============================================================================
!> @brief   Stucture-Constructor (for unit testing)
!> @return  A 2D global mesh object based on a 9-cell global mesh. 
!>          3x3 cell arrangement of quadralateral cells
!============================================================================
function global_mesh_constructor_unit_test_data() result (self)

  implicit none

  type(global_mesh_type) :: self

  integer(i_def) :: nverts = 16
  integer(i_def) :: nedges = 24

  ! Returns global_mesh_object of size 3x3 quad reference cell.
  ! As per reference cell, direction of numbering is anti-clockwise
  ! Starting point is
  ! Vertices: Bottom Left (south-west)
  ! Edges:    Left        (west)
  ! Faces:    Left        (west)
  self%ncells      = 9

  self%nverts_per_cell = 4
  self%nedges_per_cell = 4

  self%max_cells_per_vertex  = 4

  allocate( self%cell_next_2d    (self%nedges_per_cell, self%ncells) )
  allocate( self%vert_on_cell_2d (self%nverts_per_cell, self%ncells) )
  allocate( self%edge_on_cell_2d (self%nedges_per_cell, self%ncells) )
  allocate( self%cell_on_vert_2d (self%max_cells_per_vertex, nverts) )
  allocate( self%cell_on_edge_2d (2, nedges) )
  allocate( self%vert_coords     (3, nverts) )

  ! Note: These test coordinates are in [Long, Lat] in units of Radians
  self%vert_coords(1:2,1)  = [-2.0_r_def, -2.0_r_def]
  self%vert_coords(1:2,2)  = [-1.0_r_def, -2.0_r_def]
  self%vert_coords(1:2,3)  = [ 1.0_r_def, -2.0_r_def]
  self%vert_coords(1:2,4)  = [ 2.0_r_def, -2.0_r_def]
  self%vert_coords(1:2,5)  = [-2.0_r_def, -1.0_r_def]
  self%vert_coords(1:2,6)  = [-1.0_r_def, -1.0_r_def]
  self%vert_coords(1:2,7)  = [ 1.0_r_def, -1.0_r_def]
  self%vert_coords(1:2,8)  = [ 2.0_r_def, -1.0_r_def]
  self%vert_coords(1:2,9)  = [-2.0_r_def,  1.0_r_def]
  self%vert_coords(1:2,10) = [-1.0_r_def,  1.0_r_def]
  self%vert_coords(1:2,11) = [ 1.0_r_def,  1.0_r_def]
  self%vert_coords(1:2,12) = [ 2.0_r_def,  1.0_r_def]
  self%vert_coords(1:2,13) = [-2.0_r_def,  2.0_r_def]
  self%vert_coords(1:2,14) = [-1.0_r_def,  2.0_r_def]
  self%vert_coords(1:2,15) = [ 1.0_r_def,  2.0_r_def]
  self%vert_coords(1:2,16) = [ 2.0_r_def,  2.0_r_def]

  self%cell_next_2d(:,1)     = [0, 0, 2, 4]
  self%cell_next_2d(:,2)     = [1, 0, 3, 5]
  self%cell_next_2d(:,3)     = [2, 0, 0, 6]
  self%cell_next_2d(:,4)     = [0, 1, 5, 7]
  self%cell_next_2d(:,5)     = [4, 2, 6, 8]
  self%cell_next_2d(:,6)     = [5, 3, 0, 9]
  self%cell_next_2d(:,7)     = [0, 4, 8, 0]
  self%cell_next_2d(:,8)     = [7, 5, 9, 0]
  self%cell_next_2d(:,9)     = [8, 6, 0, 0]

  self%vert_on_cell_2d(:,1)  = [ 1,  2,  6,  5]
  self%vert_on_cell_2d(:,2)  = [ 2,  3,  7,  6]
  self%vert_on_cell_2d(:,3)  = [ 3,  4,  8,  7]
  self%vert_on_cell_2d(:,4)  = [ 5,  6, 10,  9]
  self%vert_on_cell_2d(:,5)  = [ 6,  7, 11, 10]
  self%vert_on_cell_2d(:,6)  = [ 7,  8, 12, 11]
  self%vert_on_cell_2d(:,7)  = [ 9, 10, 14, 13]
  self%vert_on_cell_2d(:,8)  = [10, 11, 15, 14]
  self%vert_on_cell_2d(:,9)  = [11, 12, 16, 15]

  self%edge_on_cell_2d(:,1)  = [4,   1,  5,  8]
  self%edge_on_cell_2d(:,2)  = [5,   2,  6,  9]
  self%edge_on_cell_2d(:,3)  = [6,   3,  7, 10]
  self%edge_on_cell_2d(:,4)  = [11,  8, 12, 15]
  self%edge_on_cell_2d(:,5)  = [12,  9, 13, 16]
  self%edge_on_cell_2d(:,6)  = [13, 10, 14, 17]
  self%edge_on_cell_2d(:,7)  = [18, 15, 19, 22]
  self%edge_on_cell_2d(:,8)  = [19, 16, 20, 23]
  self%edge_on_cell_2d(:,9)  = [20, 17, 21, 24]

  self%cell_on_vert_2d(:,1)  = [0, 0, 1, 0]
  self%cell_on_vert_2d(:,2)  = [0, 0, 2, 1]
  self%cell_on_vert_2d(:,3)  = [0, 0, 3, 2]
  self%cell_on_vert_2d(:,4)  = [0, 0, 0, 3]
  self%cell_on_vert_2d(:,5)  = [0, 1, 4, 0]
  self%cell_on_vert_2d(:,6)  = [1, 2, 5, 4]
  self%cell_on_vert_2d(:,7)  = [2, 3, 6, 5]
  self%cell_on_vert_2d(:,8)  = [3, 0, 0, 6]
  self%cell_on_vert_2d(:,9)  = [0, 4, 7, 0]
  self%cell_on_vert_2d(:,10) = [4, 5, 8, 7]
  self%cell_on_vert_2d(:,11) = [5, 6, 9, 8]
  self%cell_on_vert_2d(:,12) = [6, 0, 0, 9]
  self%cell_on_vert_2d(:,13) = [0, 7, 0, 0]
  self%cell_on_vert_2d(:,14) = [7, 8, 0, 0]
  self%cell_on_vert_2d(:,15) = [8, 9, 0, 0]
  self%cell_on_vert_2d(:,16) = [9, 0, 0, 0]

  self%cell_on_edge_2d(:,1)  = [0, 1]
  self%cell_on_edge_2d(:,2)  = [0, 2]
  self%cell_on_edge_2d(:,3)  = [0, 3]
  self%cell_on_edge_2d(:,4)  = [1, 0]
  self%cell_on_edge_2d(:,5)  = [2, 1]
  self%cell_on_edge_2d(:,6)  = [3, 2]
  self%cell_on_edge_2d(:,7)  = [0, 3]
  self%cell_on_edge_2d(:,8)  = [1, 4]
  self%cell_on_edge_2d(:,9)  = [2, 5]
  self%cell_on_edge_2d(:,10) = [3, 6]
  self%cell_on_edge_2d(:,11) = [4, 0]
  self%cell_on_edge_2d(:,12) = [5, 4]
  self%cell_on_edge_2d(:,13) = [6, 5]
  self%cell_on_edge_2d(:,14) = [0, 6]
  self%cell_on_edge_2d(:,15) = [4, 7]
  self%cell_on_edge_2d(:,16) = [5, 8]
  self%cell_on_edge_2d(:,17) = [6, 9]
  self%cell_on_edge_2d(:,18) = [7, 0]
  self%cell_on_edge_2d(:,19) = [8, 7]
  self%cell_on_edge_2d(:,20) = [9, 8]
  self%cell_on_edge_2d(:,21) = [0, 9]
  self%cell_on_edge_2d(:,22) = [7, 0]
  self%cell_on_edge_2d(:,23) = [8, 0]
  self%cell_on_edge_2d(:,24) = [9, 0]

end function global_mesh_constructor_unit_test_data

end module global_mesh_mod
