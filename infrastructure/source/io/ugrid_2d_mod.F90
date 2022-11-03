!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief   Store 2-dimensional ugrid mesh data.
!> @details Holds all information necessary to define ugrid vn0.9 compliant
!>          storage of 2-dimensional meshes. Pulling data out is currently
!>          done with accessor routines; this may change as dynamo matures.
!-------------------------------------------------------------------------------

module ugrid_2d_mod

use constants_mod,  only : i_def, r_def, str_def, str_longlong, l_def, &
                           imdi, rmdi, cmdi
use ugrid_file_mod, only : ugrid_file_type
use global_mesh_map_collection_mod, only: global_mesh_map_collection_type

implicit none

private

!-------------------------------------------------------------------------------
! Module parameters.
!-------------------------------------------------------------------------------

integer(i_def), parameter :: TOPOLOGY_DIMENSION  = 2

!-------------------------------------------------------------------------------
!> @brief Stores 2-dimensional grid information.
!-------------------------------------------------------------------------------
type, public :: ugrid_2d_type
  private

  character(str_def) :: mesh_name

  character(str_def) :: geometry
  character(str_def) :: topology
  character(str_def) :: coord_sys

  integer(i_def) :: npanels = imdi

  logical(l_def) :: periodic_x = .false.   !< Periodic in E-W direction.
  logical(l_def) :: periodic_y = .false.   !< Periodic in N-S direction.

  integer(i_def) :: max_stencil_depth = 0
  integer(i_def) :: void_cell = 0

  character(str_longlong) :: constructor_inputs !< Inputs used to generate mesh.

  character(str_def) :: coord_units_xy(2) = cmdi

  integer(i_def) :: edge_cells_x !< Number of cells on panel edge x-axis.
  integer(i_def) :: edge_cells_y !< Number of cells on panel edge y-axis.

  ! Numbers of different entities.
  integer(i_def) :: num_cells                !< Number of cells.
  integer(i_def) :: num_nodes                !< Number of nodes.
  integer(i_def) :: num_edges                !< Number of edges.
  integer(i_def) :: num_faces                !< Number of faces.

  integer(i_def) :: num_nodes_per_face       !< Number of nodes surrounding each face.
  integer(i_def) :: num_nodes_per_edge       !< Number of nodes defining each edge.
  integer(i_def) :: num_edges_per_face       !< Number of edges bordering each face.
  integer(i_def) :: max_num_faces_per_node   !< Maximum number of faces surrounding each node.

  ! Variables for LBC mesh only.
  integer(i_def) :: rim_depth = imdi

  ! Variables for Regional mesh only.
  ! Domain size along x/y-axes.
  real(r_def)    :: domain_size(2) = rmdi

  ! Variables for Local meshes only.
  character(str_def) :: partition_of = cmdi

  integer(i_def) :: inner_depth      = imdi
  integer(i_def) :: halo_depth       = imdi
  integer(i_def) :: num_edge         = imdi
  integer(i_def) :: last_edge_cell   = imdi
  integer(i_def) :: num_ghost        = imdi
  integer(i_def) :: last_ghost_cell  = imdi
  integer(i_def) :: num_global_cells = imdi
  integer(i_def) :: num_faces_global = imdi

  integer(i_def), allocatable :: node_cell_owner(:)
  integer(i_def), allocatable :: edge_cell_owner(:)

  integer(i_def), allocatable :: num_inner(:)
  integer(i_def), allocatable :: num_halo(:)
  integer(i_def), allocatable :: last_inner_cell(:)
  integer(i_def), allocatable :: last_halo_cell(:)

  integer(i_def), allocatable :: cell_gid(:)
  integer(i_def), allocatable :: node_on_cell_gid(:,:)
  integer(i_def), allocatable :: edge_on_cell_gid(:,:)

  ! Coordinates
  real(r_def), allocatable :: node_coordinates(:,:) !< Coordinates of nodes
  real(r_def), allocatable :: face_coordinates(:,:) !< Coordinates of faces

  ! Connectivity
  integer(i_def), allocatable :: face_node_connectivity(:,:) !< Nodes belonging to each face
  integer(i_def), allocatable :: face_edge_connectivity(:,:) !< Edges belonging to each face
  integer(i_def), allocatable :: face_face_connectivity(:,:) !< Neighbouring faces of each face
  integer(i_def), allocatable :: edge_node_connectivity(:,:) !< Nodes belonging to each edge

  ! Target mesh map variables
  integer(i_def) :: nmaps = 0 !< Number of mesh maps for this mesh (as source)

  character(str_def), allocatable :: target_mesh_names(:)   !< Target mesh names
  integer(i_def),     allocatable :: target_edge_cells_x(:) !< Target meshes panel edge cells in x-axis
  integer(i_def),     allocatable :: target_edge_cells_y(:) !< Target meshes panel edge cells in y-axis

  ! Global mesh maps
  type(global_mesh_map_collection_type), pointer :: target_global_mesh_maps => null()

  ! Information about the domain orientation
  real(r_def)    :: north_pole(2)   !< [Longitude, Latitude] of northt pole used
                                    !< for the domain orientation (degrees)
  real(r_def)    :: null_island(2)  !< [Longitude, Latitude] of null island
                                    !< used for the domain orientation (degrees)

  ! File handler
  class(ugrid_file_type), allocatable :: file_handler

contains
  procedure :: get_n_meshes
  procedure :: get_mesh_names
  procedure :: get_dimensions
  procedure :: set_by_generator
  procedure :: set_file_handler
  procedure :: set_from_file_read
  procedure :: write_to_file
  procedure :: append_to_file
  procedure :: get_metadata
  procedure :: get_coord_units
  procedure :: get_node_coords
  procedure :: get_face_coords
  procedure :: get_face_node_connectivity
  procedure :: get_face_edge_connectivity
  procedure :: get_face_face_connectivity
  procedure :: get_edge_node_connectivity
  procedure :: write_coordinates

  procedure :: set_global_mesh_maps
  generic   :: set_mesh_maps => set_global_mesh_maps

  procedure :: get_global_mesh_maps
  generic   :: get_mesh_maps => get_global_mesh_maps


  !> Routine to destroy object
  procedure :: clear

  !> Object finalizer
  final     :: ugrid_2d_destructor

end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines.
!-------------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------------
!>  @brief Gets number of nodes, edges, faces etc.
!>
!>  @param[out]     num_nodes               Number of nodes.
!>  @param[out]     num_edges               Number of edges.
!>  @param[out]     num_faces               Number of faces.
!>  @param[out]     num_nodes_per_face      Number of nodes around each face.
!>  @param[out]     num_edges_per_face      Number of edges around each face.
!>  @param[out]     num_nodes_per_edge      Number of nodes defining each edge.
!>  @param[out]     max_num_faces_per_node  Maximum number of faces surrounding each node.
!-------------------------------------------------------------------------------

subroutine get_dimensions( self, num_nodes, num_edges, num_faces,  &
                           num_nodes_per_face, num_edges_per_face, &
                           num_nodes_per_edge, max_num_faces_per_node )

  implicit none

  ! Arguments
  class(ugrid_2d_type), intent(in) :: self

  integer(i_def), intent(out) :: num_nodes
  integer(i_def), intent(out) :: num_edges
  integer(i_def), intent(out) :: num_faces
  integer(i_def), intent(out) :: num_nodes_per_face
  integer(i_def), intent(out) :: num_edges_per_face
  integer(i_def), intent(out) :: num_nodes_per_edge
  integer(i_def), intent(out) :: max_num_faces_per_node

  num_nodes = self%num_nodes
  num_edges = self%num_edges
  num_faces = self%num_faces

  num_nodes_per_face = self%num_nodes_per_face
  num_edges_per_face = self%num_edges_per_face
  num_nodes_per_edge = self%num_nodes_per_edge
  max_num_faces_per_node = self%max_num_faces_per_node

  return
end subroutine get_dimensions

!-------------------------------------------------------------------------------
!> @brief Gets a list of the mesh names in a ugrid file.
!>
!> @param[in]  filename   The name of the file to query.
!> @param[out] mesh_names Character[:], Name of mesh topologies in the file.
!-------------------------------------------------------------------------------
subroutine get_mesh_names(self, filename, mesh_names)

  implicit none

  class(ugrid_2d_type), intent(inout) :: self
  character(len=*),     intent(in)    :: filename
  character(len=*),     intent(out)   :: mesh_names(:)

  call self%file_handler%file_open(trim(filename))
  call self%file_handler%get_mesh_names(mesh_names)
  call self%file_handler%file_close()

  return
end subroutine get_mesh_names

!-------------------------------------------------------------------------------
!> @brief Gets the number of the mesh topologies in a ugrid file.
!>
!> @param[in]  filename   The name of the file to query.
!> @param[out] n_meshes   Integer, Number of mesh topologies in the file.
!-------------------------------------------------------------------------------
subroutine get_n_meshes(self, filename, n_meshes)

  implicit none

  class(ugrid_2d_type), intent(inout) :: self
  character(len=*),     intent(in)    :: filename
  integer(i_def),       intent(out)   :: n_meshes

  call self%file_handler%file_open(trim(filename))
  n_meshes = self%file_handler%get_n_meshes()
  call self%file_handler%file_close()

  return
end subroutine get_n_meshes

!-------------------------------------------------------------------------------
!>  @brief   Allocates ugrid_2d internal storage, populated by ugrid generator.
!>  @details Allocates component arrays according to sizes obtained from the
!>           passed generator strategy.
!>
!>  @param[in]  generator_strategy The generator strategy in use.
!-------------------------------------------------------------------------------

subroutine allocate_arrays(self, generator_strategy)
  use ugrid_generator_mod, only: ugrid_generator_type
  implicit none

  ! Arguments
  type(ugrid_2d_type),         intent(inout) :: self
  class(ugrid_generator_type), intent(in)    :: generator_strategy

  call generator_strategy%get_dimensions(                      &
         num_nodes              = self%num_nodes,              &
         num_edges              = self%num_edges,              &
         num_faces              = self%num_faces,              &
         num_nodes_per_face     = self%num_nodes_per_face,     &
         num_edges_per_face     = self%num_edges_per_face,     &
         num_nodes_per_edge     = self%num_nodes_per_edge,     &
         max_num_faces_per_node = self%max_num_faces_per_node )



  if ( allocated( self%node_coordinates ) )       &
                                     deallocate( self%node_coordinates )
  if ( allocated( self%face_coordinates ) )       &
                                     deallocate( self%face_coordinates )

  if ( allocated( self%edge_node_connectivity ) ) &
                                     deallocate( self%edge_node_connectivity )

  if ( allocated( self%face_node_connectivity ) ) &
                                     deallocate( self%face_node_connectivity )

  if ( allocated( self%face_edge_connectivity ) ) &
                                     deallocate( self%face_edge_connectivity )

  if ( allocated( self%face_face_connectivity ) ) &
                                     deallocate( self%face_face_connectivity )

  allocate(self%node_coordinates(2, self%num_nodes))
  allocate(self%face_coordinates(2, self%num_faces))
  allocate(self%edge_node_connectivity(self%num_nodes_per_edge, self%num_edges))
  allocate(self%face_node_connectivity(self%num_nodes_per_face, self%num_faces))
  allocate(self%face_edge_connectivity(self%num_edges_per_face, self%num_faces))
  allocate(self%face_face_connectivity(self%num_edges_per_face, self%num_faces))

  if (self%nmaps > 0) then
    if (allocated(self%target_mesh_names))   deallocate (self%target_mesh_names)
    if (allocated(self%target_edge_cells_x)) deallocate (self%target_edge_cells_x)
    if (allocated(self%target_edge_cells_y)) deallocate (self%target_edge_cells_y)

    allocate(self%target_mesh_names(self%nmaps))
    allocate(self%target_edge_cells_x(self%nmaps))
    allocate(self%target_edge_cells_y(self%nmaps))
  end if

  return
end subroutine allocate_arrays

!-------------------------------------------------------------------------------
!>  @brief   Allocates ugrid_2d internal storage, populated by ugrid file.
!>  @details Allocates component arrays according to sizes already stored in
!>           the ugrid_2d object. These arrays are populated elsewhere
!>           by a ugrid file strategy.
!>
!-------------------------------------------------------------------------------

subroutine allocate_arrays_for_file(self)
  implicit none

  ! Arguments
  type(ugrid_2d_type),    intent(inout) :: self

  if ( allocated( self%node_coordinates ) )       &
                                     deallocate( self%node_coordinates )
  if ( allocated( self%face_coordinates ) )       &
                                     deallocate( self%face_coordinates )

  if ( allocated( self%edge_node_connectivity ) ) &
                                     deallocate( self%edge_node_connectivity )

  if ( allocated( self%face_node_connectivity ) ) &
                                     deallocate( self%face_node_connectivity )

  if ( allocated( self%face_edge_connectivity ) ) &
                                     deallocate( self%face_edge_connectivity )

  if ( allocated( self%face_face_connectivity ) ) &
                                     deallocate( self%face_face_connectivity )

  if ( allocated( self%target_mesh_names ) ) &
                                     deallocate( self%target_mesh_names )

  allocate(self%node_coordinates(2, self%num_nodes))
  allocate(self%face_coordinates(2, self%num_faces))

  allocate(self%edge_node_connectivity(self%num_nodes_per_edge, self%num_edges))
  allocate(self%face_node_connectivity(self%num_nodes_per_face, self%num_faces))
  allocate(self%face_edge_connectivity(self%num_edges_per_face, self%num_faces))
  allocate(self%face_face_connectivity(self%num_edges_per_face, self%num_faces))

  if (self%nmaps > 0) then
    allocate( self%target_mesh_names(self%nmaps) )
  end if

  return
end subroutine allocate_arrays_for_file

!---------------------------------------------------------------------------------
!>  @brief   Populate arrays according to passed generator strategy.
!>  @details Calls back to the passed generator strategy in order to populate
!>           the coordinate and connectivity arrays.
!>
!>  @param[in,out] generator_strategy The generator with which to generate the mesh.
!---------------------------------------------------------------------------------
subroutine set_by_generator(self, generator_strategy)

  use ugrid_generator_mod, only: ugrid_generator_type

  implicit none

  class(ugrid_2d_type),        intent(inout) :: self
  class(ugrid_generator_type), intent(inout) :: generator_strategy

  call generator_strategy%get_metadata                &
      ( mesh_name          = self%mesh_name,          &
        geometry           = self%geometry,           &
        topology           = self%topology,           &
        coord_sys          = self%coord_sys,          &
        periodic_x         = self%periodic_x,         &
        periodic_y         = self%periodic_y,         &
        edge_cells_x       = self%edge_cells_x,       &
        edge_cells_y       = self%edge_cells_y,       &
        constructor_inputs = self%constructor_inputs, &
        north_pole         = self%north_pole,         &
        null_island        = self%null_island,        &
        rim_depth          = self%rim_depth,          &
        domain_size        = self%domain_size,        &
        nmaps              = self%nmaps,              &
        void_cell          = self%void_cell )


  self%npanels = generator_strategy%get_number_of_panels()

  call generator_strategy%generate()

  call allocate_arrays(self, generator_strategy)

  if (self%nmaps > 0) then
    call generator_strategy%get_metadata                &
        ( target_mesh_names = self%target_mesh_names,   &
          maps_edge_cells_x = self%target_edge_cells_x, &
          maps_edge_cells_y = self%target_edge_cells_y )

    self%target_global_mesh_maps => generator_strategy%get_global_mesh_maps()
  end if

  call generator_strategy%get_coordinates          &
      ( node_coordinates = self%node_coordinates,  &
        cell_coordinates = self%face_coordinates,  &
        coord_units_x    = self%coord_units_xy(1), &
        coord_units_y    = self%coord_units_xy(2) )

  call generator_strategy%get_connectivity                    &
      ( face_node_connectivity = self%face_node_connectivity, &
        edge_node_connectivity = self%edge_node_connectivity, &
        face_edge_connectivity = self%face_edge_connectivity, &
        face_face_connectivity = self%face_face_connectivity )

  return
end subroutine set_by_generator

!-------------------------------------------------------------------------------
!> @brief   Sets the file read/write strategy.
!> @details Receives a file-handler object and moves its allocation into
!>          the appropriate type component. On exit, the file_handler
!>          dummy argument will no longer be allocated.
!>
!> @param[in,out] file_handler The file handler to use for IO.
!-------------------------------------------------------------------------------

subroutine set_file_handler(self, file_handler)
  implicit none

  ! Arguments
  class(ugrid_2d_type),                intent(inout) :: self
  class(ugrid_file_type), allocatable, intent(inout) :: file_handler

  call move_alloc(file_handler, self%file_handler)

  return
end subroutine set_file_handler

!-------------------------------------------------------------------------------
!> @brief   Reads ugrid information and populates internal arrays.
!> @details Calls back to the file handler strategy (component) in order to
!>          read the ugrid mesh data and populate internal arrays with data
!>          from a file.
!>
!> @param[in]     mesh_name Name of the mesh topology.
!> @param[in]     filename  Name of the ugrid file.
!-------------------------------------------------------------------------------

subroutine set_from_file_read(self, mesh_name, filename)

  implicit none

  ! Arguments
  class(ugrid_2d_type), intent(inout) :: self
  character(len=*),     intent(in)    :: mesh_name
  character(len=*),     intent(in)    :: filename

  self%mesh_name = trim(mesh_name)

  call self%file_handler%file_open(trim(filename))

  call self%file_handler%get_dimensions(                    &
         mesh_name              = self%mesh_name,           &
         num_nodes              = self%num_nodes,           &
         num_edges              = self%num_edges,           &
         num_faces              = self%num_faces,           &
         num_nodes_per_face     = self%num_nodes_per_face,  &
         num_edges_per_face     = self%num_edges_per_face,  &
         num_nodes_per_edge     = self%num_nodes_per_edge,  &
         max_num_faces_per_node = self%max_num_faces_per_node )

  call allocate_arrays_for_file(self)

  call self%file_handler%read_mesh(                         &
      mesh_name              = self%mesh_name,              &
      geometry               = self%geometry,               &
      topology               = self%topology,               &
      coord_sys              = self%coord_sys,              &
      periodic_x             = self%periodic_x,             &
      periodic_y             = self%periodic_y,             &
      npanels                = self%npanels,                &
      max_stencil_depth      = self%max_stencil_depth,      &
      constructor_inputs     = self%constructor_inputs,     &
      node_coordinates       = self%node_coordinates,       &
      face_coordinates       = self%face_coordinates,       &
      coord_units_x          = self%coord_units_xy(1),      &
      coord_units_y          = self%coord_units_xy(2),      &
      void_cell              = self%void_cell,              &
      face_node_connectivity = self%face_node_connectivity, &
      edge_node_connectivity = self%edge_node_connectivity, &
      face_edge_connectivity = self%face_edge_connectivity, &
      face_face_connectivity = self%face_face_connectivity, &
      num_targets            = self%nmaps,                  &
      target_mesh_names      = self%target_mesh_names,      &
      north_pole             = self%north_pole,             &
      null_island            = self%null_island )

  call self%file_handler%file_close()

  return
end subroutine set_from_file_read

!-------------------------------------------------------------------------------
!> @brief   Writes stored ugrid information to data file.
!> @details Calls back to the file handler strategy (component) in order to
!>          read the ugrid mesh data and populate internal arrays with data
!>          from a file.
!>
!> @param[in] filename  Filenameto write contents to.
!-------------------------------------------------------------------------------

subroutine write_to_file(self, filename)

  use ugrid_generator_mod, only: ugrid_generator_type

  implicit none

  ! Arguments
  class(ugrid_2d_type), intent(inout) :: self
  character(len=*),     intent(in)    :: filename

  call self%file_handler%file_new(trim(filename))

  call self%file_handler%write_mesh(                         &
       mesh_name              = self%mesh_name,              &
       geometry               = self%geometry,               &
       topology               = self%topology,               &
       coord_sys              = self%coord_sys,              &
       periodic_x             = self%periodic_x,             &
       periodic_y             = self%periodic_y,             &
       npanels                = self%npanels,                &
       north_pole             = self%north_pole,             &
       null_island            = self%null_island,            &
       max_stencil_depth      = self%max_stencil_depth,      &
       constructor_inputs     = self%constructor_inputs,     &

       num_nodes              = self%num_nodes,              &
       num_edges              = self%num_edges,              &
       num_faces              = self%num_faces,              &
       node_coordinates       = self%node_coordinates,       &
       face_coordinates       = self%face_coordinates,       &
       coord_units_x          = self%coord_units_xy(1),      &
       coord_units_y          = self%coord_units_xy(2),      &
       void_cell              = self%void_cell,              &
       face_node_connectivity = self%face_node_connectivity, &
       edge_node_connectivity = self%edge_node_connectivity, &
       face_edge_connectivity = self%face_edge_connectivity, &
       face_face_connectivity = self%face_face_connectivity, &

       ! Intergrid maps
       num_targets             = self%nmaps,             &
       target_mesh_names       = self%target_mesh_names, &
       target_global_mesh_maps = self%target_global_mesh_maps )



  call self%file_handler%file_close()


  return
end subroutine write_to_file

!-------------------------------------------------------------------------------
!> @brief   Appends stored ugrid information to existing ugrid data file.
!> @details Calls back to the file handler strategy (component) in order to
!>          read the ugrid_2d_type mesh data and populate the file_handlers
!>          internal arrays which the file handler will append to the ugrid file.
!>
!> @param[in] filename  Output file to open to add data.
!-------------------------------------------------------------------------------

subroutine append_to_file(self, filename)

  use ugrid_generator_mod, only: ugrid_generator_type

  implicit none

  ! Arguments
  class(ugrid_2d_type), intent(inout) :: self
  character(len=*),     intent(in) :: filename

  call self%file_handler%file_open(trim(filename))

  call self%file_handler%append_mesh(                        &
       mesh_name              = self%mesh_name,              &
       geometry               = self%geometry,               &
       topology               = self%topology,               &
       coord_sys              = self%coord_sys,              &
       periodic_x             = self%periodic_x,             &
       periodic_y             = self%periodic_y,             &
       npanels                = self%npanels,                &
       north_pole             = self%north_pole,             &
       null_island            = self%null_island,            &
       max_stencil_depth      = self%max_stencil_depth,      &
       constructor_inputs     = self%constructor_inputs,     &
       num_nodes              = self%num_nodes,              &
       num_edges              = self%num_edges,              &
       num_faces              = self%num_faces,              &
       node_coordinates       = self%node_coordinates,       &
       face_coordinates       = self%face_coordinates,       &
       coord_units_x          = self%coord_units_xy(1),      &
       coord_units_y          = self%coord_units_xy(2),      &
       void_cell              = self%void_cell,              &
       face_node_connectivity = self%face_node_connectivity, &
       edge_node_connectivity = self%edge_node_connectivity, &
       face_edge_connectivity = self%face_edge_connectivity, &
       face_face_connectivity = self%face_face_connectivity, &

       ! InterGrid Maps
       num_targets             = self%nmaps,             &
       target_mesh_names       = self%target_mesh_names, &
       target_global_mesh_maps = self%target_global_mesh_maps )

  call self%file_handler%file_close()

  return
end subroutine append_to_file


!-------------------------------------------------------------------------------
!> @brief    Gets metadata of the current mesh in this object which was
!>           set by the ugrid_generator_type or read in from NetCDF (UGRID) file.
!> @details  The ugrid_2d object may hold a global or local mesh of differing
!>           topologies. As such, some of the metadata from this routine
!>           will not apply in all cases.
!>
!> @param[out] mesh_name          Name of the current mesh topology.
!> @param[out] geometry           Domain geometry enumeration key.
!> @param[out] topology           Domain topology enumeration key.
!> @param[out] coord_sys          Co-ordinate sys enumeration key.
!> @param[out] npanels            Number of panels used in mesh topology.
!> @param[out] periodic_x         Periodic in E-W direction.
!> @param[out] periodic_y         Periodic in N-S direction.
!> @param[out] max_stencil_depth  Maximum stencil depth supported (Local meshes).
!> @param[out] void_cell          Values to mark a cell as external
!>                                to domain.
!> @param[out] edge_cells_x       Number of panel edge cells (x-axis).
!> @param[out] edge_cells_y       Number of panel edge cells (y-axis).
!> @param[out] ncells_global      Number of cells in the global mesh.
!> @param[out] domain_size        Domain size ix x/y-axes.
!> @param[out] rim_depth          Rim depth (in cells) for LBC meshes.
!> @param[out] inner_depth        Number of inner halo layers.
!> @param[out] halo_depth         Number of outer halo layers.
!> @param[out] partition_of       Specifies the global mesh name, if this mesh a local mesh.
!> @param[out] num_edge           Number of cells in the edge layer of the partition.
!> @param[out] last_edge_cell     Array index of the last cell in the edge layer.
!> @param[out] num_ghost          Number of cells in the ghost layer of the partition.
!> @param[out] last_ghost_cell    Array index of the last cell in the ghost layer.
!> @param[out] constructor_inputs Input arguments use to create this mesh.
!> @param[out] nmaps              The number of intergrid maps from this mesh.
!> @param[out] target_mesh_names  Names of target mesh topologies in this file
!>                                which this mesh possesses cell-cell maps for.
!> @param[out] north_pole         Optional, [Longitude, Latitude] of north pole
!>                                used for domain orientation (degrees).
!> @param[out] null_island        Optional, [Longitude, Latitude] of null
!>                                island used for domain orientation (degrees).
!-------------------------------------------------------------------------------
subroutine get_metadata( self, mesh_name,                 &
                         geometry, topology, coord_sys,   &
                         npanels, periodic_x, periodic_y, &
                         max_stencil_depth,               &
                         void_cell,                       &
                         edge_cells_x, edge_cells_y,      &
                         ncells_global, domain_size,      &
                         rim_depth, inner_depth,          &
                         halo_depth, partition_of,        &
                         num_edge, last_edge_cell,        &
                         num_ghost, last_ghost_cell,      &
                         constructor_inputs, nmaps,       &
                         target_mesh_names,               &
                         north_pole, null_island )

  implicit none

  class(ugrid_2d_type), intent(in) :: self
  character(str_def),   optional, intent(out) :: mesh_name
  character(str_def),   optional, intent(out) :: geometry
  character(str_def),   optional, intent(out) :: topology
  character(str_def),   optional, intent(out) :: coord_sys
  integer(i_def),       optional, intent(out) :: npanels
  logical(l_def),       optional, intent(out) :: periodic_x
  logical(l_def),       optional, intent(out) :: periodic_y

  integer(i_def),       optional, intent(out) :: max_stencil_depth
  integer(i_def),       optional, intent(out) :: void_cell
  integer(i_def),       optional, intent(out) :: edge_cells_x
  integer(i_def),       optional, intent(out) :: edge_cells_y

  character(str_longlong),  optional, intent(out) :: constructor_inputs

  integer(i_def),       optional, intent(out) :: nmaps
  character(str_def),   optional, intent(out), &
                                  allocatable :: target_mesh_names(:)

  character(str_def),   optional, intent(out) :: partition_of

  real(r_def),    optional, intent(out) :: north_pole(2)
  real(r_def),    optional, intent(out) :: null_island(2)
  integer(i_def), optional, intent(out) :: inner_depth
  integer(i_def), optional, intent(out) :: halo_depth
  integer(i_def), optional, intent(out) :: num_edge
  integer(i_def), optional, intent(out) :: last_edge_cell
  integer(i_def), optional, intent(out) :: num_ghost
  integer(i_def), optional, intent(out) :: last_ghost_cell
  integer(i_def), optional, intent(out) :: ncells_global
  real(r_def),    optional, intent(out) :: domain_size(2)
  integer(i_def), optional, intent(out) :: rim_depth

  if (present(constructor_inputs)) constructor_inputs = self%constructor_inputs
  if (present(domain_size))        domain_size        = self%domain_size
  if (present(rim_depth))          rim_depth          = self%rim_depth
  if (present(partition_of))       partition_of       = self%partition_of
  if (present(inner_depth))        inner_depth        = self%inner_depth
  if (present(halo_depth))         halo_depth         = self%halo_depth
  if (present(num_edge))           num_edge           = self%num_edge
  if (present(last_edge_cell))     last_edge_cell     = self%last_edge_cell
  if (present(num_ghost))          num_ghost          = self%num_ghost
  if (present(last_ghost_cell))    last_ghost_cell    = self%last_ghost_cell
  if (present(north_pole))         north_pole         = self%north_pole
  if (present(null_island))        null_island        = self%null_island
  if (present(ncells_global))      ncells_global      = self%num_global_cells

  if (present(mesh_name))          mesh_name          = self%mesh_name
  if (present(geometry))           geometry           = self%geometry
  if (present(topology))           topology           = self%topology
  if (present(coord_sys))          coord_sys          = self%coord_sys
  if (present(npanels ))           npanels            = self%npanels
  if (present(periodic_x))         periodic_x         = self%periodic_x
  if (present(periodic_y))         periodic_y         = self%periodic_y
  if (present(max_stencil_depth))  max_stencil_depth  = self%max_stencil_depth
  if (present(void_cell))          void_cell          = self%void_cell
  if (present(constructor_inputs)) constructor_inputs = self%constructor_inputs
  if (present(edge_cells_x))       edge_cells_x       = self%edge_cells_x
  if (present(edge_cells_y))       edge_cells_y       = self%edge_cells_y
  if (present(nmaps))              nmaps              = self%nmaps
  if (self%nmaps > 0 .and. present(target_mesh_names)) then
    target_mesh_names = self%target_mesh_names
  end if
  if (present(north_pole))         north_pole(:)      = self%north_pole(:)
  if (present(null_island))        null_island(:)     = self%null_island(:)

end subroutine get_metadata

!-------------------------------------------------------------------------------
!> @brief  Gets node coordinate units in [x,y] or [longitude, latitude] directions.
!>
!> @param[out]  coord_units_xy  String for coordinate units along x/y-axes.
!-------------------------------------------------------------------------------
subroutine get_coord_units(self, coord_units_xy)
  implicit none

  class(ugrid_2d_type), intent(in)  :: self
  character(str_def),   intent(out) :: coord_units_xy(2)

  coord_units_xy = self%coord_units_xy

  return
end subroutine get_coord_units

!-------------------------------------------------------------------------------
!> @brief   Gets node coordinates with ugrid array index ordering.
!> @details Returns a rank-two array of node coordinates, with the
!>          coordinate dimension index innermost, and the node number
!>          outermost. Format: [long, lat, radius].
!>
!> @param[out]  node_coords  Node coordinate array.
!-------------------------------------------------------------------------------

subroutine get_node_coords(self, node_coords)
  implicit none

  class(ugrid_2d_type), intent(in)  :: self

  real(r_def),          intent(out) :: node_coords(:,:)

  integer(i_def) :: i

  do i = 1, self%num_nodes
    node_coords(1:2,i) = self%node_coordinates(1:2,i)
  end do

  return
end subroutine get_node_coords


!-------------------------------------------------------------------------------
!> @brief   Gets face coordinates with ugrid array index ordering.
!> @details Returns a rank-two array of node coordinates, with the
!>          coordinate dimension index innermost, and the face number
!>          outermost. Format: [long, lat].
!>
!> @return  face_coords Face coordinate array.
!-------------------------------------------------------------------------------

function get_face_coords(self) result(face_coords)

  implicit none

  class(ugrid_2d_type), intent(in)  :: self

  real(r_def), allocatable :: face_coords(:,:)

  face_coords = self%face_coordinates(:,:)

  return
end function get_face_coords

!-------------------------------------------------------------------------------
!> @brief   Gets an array of node indices surrounding each face.
!> @details Returns a rank-two array of nodes surrounding each face, with
!>          the nodes surrounding any single face being contiguous.
!>
!> @return  face_node_connectivity  Indices of nodes adjacent to faces.
!-------------------------------------------------------------------------------
function get_face_node_connectivity( self ) &
                             result( face_node_connectivity )
  implicit none

  class(ugrid_2d_type), intent(in) :: self
  integer(i_def), allocatable      :: face_node_connectivity(:,:)

  if ( allocated(face_node_connectivity) ) then
    deallocate(face_node_connectivity)
  end if

  allocate( face_node_connectivity, &
            source=self%face_node_connectivity )
  return
end function get_face_node_connectivity

!-------------------------------------------------------------------------------
!> @brief   Gets an array of edge indices surrounding each face.
!> @details Returns a rank-two array of edges surrounding each face, with
!>          the edges surrounding any single face being contiguous.
!>
!> @return  face_edge_connectivity  Indices of edges adjacent to faces.
!-------------------------------------------------------------------------------
function get_face_edge_connectivity( self ) &
                             result( face_edge_connectivity )
  implicit none

  class(ugrid_2d_type), intent(in) :: self
  integer(i_def), allocatable      :: face_edge_connectivity(:,:)

  if ( allocated(face_edge_connectivity) ) then
    deallocate(face_edge_connectivity)
  end if

  allocate( face_edge_connectivity, &
            source=self%face_edge_connectivity )
  return
end function get_face_edge_connectivity


!-------------------------------------------------------------------------------
!> @brief   Gets an array of face indices surrounding each face.
!> @details Returns a rank-two array of faces surrounding each face, with
!>          the faces surrounding any single face being contiguous.
!>
!> @return  face_face_connectivity  Indices of faces adjacent to faces.
!-------------------------------------------------------------------------------
function get_face_face_connectivity( self ) &
                             result( face_face_connectivity )
  implicit none

  class(ugrid_2d_type), intent(in) :: self
  integer(i_def), allocatable      :: face_face_connectivity(:,:)

  if ( allocated(face_face_connectivity) ) then
    deallocate(face_face_connectivity)
  end if

  allocate( face_face_connectivity, &
            source=self%face_face_connectivity )
  return
end function get_face_face_connectivity


!-------------------------------------------------------------------------------
!> @brief   Gets an array of node indices attached to each edge.
!> @details Returns a rank-two array of nodes attached to each edge.
!>
!> @return  edge_node_connectivity  Indices of nodes attached to edges.
!-------------------------------------------------------------------------------
function get_edge_node_connectivity( self ) &
                             result( edge_node_connectivity )
  implicit none

  class(ugrid_2d_type), intent(in) :: self
  integer(i_def), allocatable      :: edge_node_connectivity(:,:)

  if ( allocated(edge_node_connectivity) ) then
    deallocate(edge_node_connectivity)
  end if

  allocate( edge_node_connectivity, &
            source=self%edge_node_connectivity )
  return
end function get_edge_node_connectivity


!-------------------------------------------------------------------------------
!> @brief     Assigns a global_mesh_map_collection to the ugrid_2d object.
!> @param[in] global_maps  Global mesh map collection.
!-------------------------------------------------------------------------------
subroutine set_global_mesh_maps( self, global_maps  )

  implicit none

  class (ugrid_2d_type), intent(inout) :: self

  type(global_mesh_map_collection_type), &
                         intent(in), target :: global_maps

  self%target_global_mesh_maps => global_maps

end subroutine set_global_mesh_maps

!-------------------------------------------------------------------------------
!> @brief     Retrieves pointer to current global_mesh_map_collection
!>            assigned to the ugrid_2d object.
!> @param[out] global_maps  Pointer to global mesh map collection in ugrid_2d
!>                          object.
!-------------------------------------------------------------------------------
subroutine get_global_mesh_maps( self, global_maps  )

  implicit none

  class (ugrid_2d_type), intent(in), target :: self

  type(global_mesh_map_collection_type), &
                         intent(out), pointer :: global_maps

  global_maps => self%target_global_mesh_maps

end subroutine get_global_mesh_maps

!-------------------------------------------------------------------------------
!> @brief   Writes coordinates to a .dat file in the units they are held in
!>          within the UGRID file.
!> @details Produces a file with rough output of coordinates. Intended to be
!>          a temporary routine only: would be better implemented for the
!>          long-term as a plain-text ugrid file strategy. Hence the
!>          good-enough-for-now hardwired unit numbers and file names.
!>
!-------------------------------------------------------------------------------
subroutine write_coordinates(self)
  implicit none

  ! Arguments
  class(ugrid_2d_type), intent(in) :: self

  integer(i_def) :: inode

  open(56, file='nodes.dat')
  do inode = 1, self%num_nodes
    write(56,*) self%node_coordinates(:,inode)
  end do
  close(56)

  return
end subroutine write_coordinates

!-------------------------------------------------------------------------------
!> @brief Routine to destroy object.
!-------------------------------------------------------------------------------
subroutine clear(self)

  implicit none

  class (ugrid_2d_type), intent(inout) :: self

  if (allocated(self%node_coordinates))       deallocate( self%node_coordinates )
  if (allocated(self%face_coordinates))       deallocate( self%face_coordinates )

  if (allocated(self%face_node_connectivity)) deallocate( self%face_node_connectivity )
  if (allocated(self%edge_node_connectivity)) deallocate( self%edge_node_connectivity )
  if (allocated(self%face_edge_connectivity)) deallocate( self%face_edge_connectivity )
  if (allocated(self%face_face_connectivity)) deallocate( self%face_face_connectivity )

  if (allocated(self%target_mesh_names))      deallocate( self%target_mesh_names )
  if (allocated(self%target_edge_cells_x))    deallocate( self%target_edge_cells_x )
  if (allocated(self%target_edge_cells_y))    deallocate( self%target_edge_cells_y )

  if (allocated(self%file_handler))           deallocate( self%file_handler )

  self%mesh_name  = cmdi
  self%geometry   = cmdi
  self%topology   = cmdi
  self%coord_sys  = cmdi
  self%periodic_x = .false.
  self%periodic_y = .false.

  self%max_stencil_depth  = imdi
  self%npanels            = imdi
  self%constructor_inputs = cmdi

  self%coord_units_xy = cmdi
  self%edge_cells_x   = imdi
  self%edge_cells_y   = imdi

  self%num_cells = imdi
  self%num_nodes = imdi
  self%num_edges = imdi
  self%num_faces = imdi

  self%num_nodes_per_face = imdi
  self%num_nodes_per_edge = imdi
  self%num_edges_per_face = imdi
  self%max_num_faces_per_node = imdi

  self%nmaps            = 0
  self%rim_depth        = imdi
  self%domain_size(2)   = rmdi
  self%partition_of     = cmdi
  self%inner_depth      = imdi
  self%halo_depth       = imdi
  self%num_edge         = imdi
  self%last_edge_cell   = imdi
  self%num_ghost        = imdi
  self%last_ghost_cell  = imdi
  self%num_faces_global = imdi

  self%target_global_mesh_maps => null()

end subroutine clear

!-------------------------------------------------------------------------------
!> @brief Finalizer routine which should automatically call clear
!>        when object is out of scope.
!-------------------------------------------------------------------------------
subroutine ugrid_2d_destructor(self)

  implicit none

  type (ugrid_2d_type), intent(inout) :: self

  call self%clear()

end subroutine ugrid_2d_destructor

end module ugrid_2d_mod
