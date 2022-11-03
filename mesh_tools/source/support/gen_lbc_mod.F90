!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Module to define the ugrid generation strategy object (gen_lbc_type)
!>          for a LBC (rim) mesh.
!> @details This type is a subclass of the ugrid_generator_type. It describes
!>          a LBC mesh format suitable for storage as a ugrid file.
!>          All required connectivity is calculated via the generate method
!>          and made available to the ugrid writer.
!>
!-------------------------------------------------------------------------------
module gen_lbc_mod
!-------------------------------------------------------------------------------

  use constants_mod,                  only: r_def, i_def, l_def, str_def,   &
                                            str_longlong, imdi, rmdi, emdi, &
                                            i_native, radians_to_degrees
  use gen_planar_mod,                 only: gen_planar_type
  use global_mesh_map_mod,            only: global_mesh_map_type
  use global_mesh_map_collection_mod, only: global_mesh_map_collection_type
  use log_mod,                        only: log_event, log_scratch_space, &
                                            LOG_LEVEL_ERROR, LOG_LEVEL_WARNING
  use reference_element_mod,          only: W, S, E, N, &
                                            SWB, SEB, NWB, NEB
  use ugrid_generator_mod,            only: ugrid_generator_type

  use mesh_config_mod,                only: geometry_from_key,  &
                                            coord_sys_from_key, &
                                            key_from_geometry,  &
                                            key_from_topology,  &
                                            key_from_coord_sys, &
                                            topology_non_periodic


  implicit none

  private

  ! Module private functions
  private :: calc_panel_adjacency,         &
             calc_panel_node_connectivity, &
             calc_panel_edge_connectivity, &
             calc_panel_lbc_lam_map

  ! Module private subroutines
  private :: calc_adjacency,               &
             calc_cell_node_connectivity,  &
             calc_cell_edge_connectivity,  &
             calc_edge_node_connectivity,  &
             calc_lbc_lam_map,             &
             extract_coords


  ! Mesh Vertex directions: local aliases for reference_element_mod values
  integer(i_def), parameter :: NW = NWB
  integer(i_def), parameter :: NE = NEB
  integer(i_def), parameter :: SE = SEB
  integer(i_def), parameter :: SW = SWB

  ! Set to -9999 to mark null (void) cell connectivity.
  integer(i_def), parameter :: VOID_ID = -9999

  integer(i_def), parameter :: horz_npanels = 2
  integer(i_def), parameter :: vert_npanels = 2

  ! For a lbc meshes there is only one panel
  integer(i_def), parameter :: NPANELS = 1


!===================================================================
! A: METHODOLOGY
!===================================================================
! NOTE: The meshes referred to in this module as 2D reference to
!       "cell" are the equivalent of the a "face" and may in some
!       parts be used interchangeably.
!
! For the generation of a LBC mesh around the "rim" of a
! planar mesh domain, this strategy breaks down the "rim" into
! 4 generic retangular panels as shown below. The "base" panels
! are the same (horizontalx2, verticalx2) except for length and
! the IDs assigned to them.
!
! The LBC panels are numbered and arranged in relation to the
! parent LAM as shown below.
!
! LAM NW      LBC Panel #1      LAM NE
!        O---------------+----O
!        |      HP#1     |    |
!        +----+----------+    |
!        |    |          |    |
!  LBC   |    |          |VP#1| LBC
!  Panel |    |          |    | Panel
!   #4   |VP#2|          |    |  #2
!        |    |          |    |
!        |    +----------+----+
!        |    |      HP#2     |
!        O----+---------------O
! LAM SW      LBC Panel #2      LAM SE
!
! 1. Cells for each base panel are numbered in rows, when orientated
!    with the row as the longest edge. IDs on each base panel
!    continue from where the last panel ended.
! 2. LBC-LAM integrid cell map can be acquired as the
!    top-left corner of each base panel maps to a corner cell on
!    the parent LAM. This information along with the LAM "cell-next"
!    connectivity is used to create the integrid cell map.
! 3. Base panel connectivities are generated treating each
!    panel as being non-periodic at the base panel boundaries.
! 4. As shown, in the schematic, the "east" end of each base panel
!    connects to left most southern cells of the adjacent panel. This
!    property is used to "stitch" the base panels together as shown.
! 5, Each of the base panels connectivity information is cyclically
!    rotated so that all base panels are using the same frame of
!    reference.
! 6. Using the LBC-LAM intergrid map and connectivity information,
!    the node coordinates of the LBC mesh can be extracted from
!    the LAM parent.
!
! Note: The LBC mesh/integrid map produced does not require the cell
!       connectivity at the LAM domain boundaries as the LBC mesh
!       base panels are always considered to have no periodicity.
!       As such, the generation of a LBC mesh should work for a
!       planar mesh ugrid generator object regardless of periodicity
!       at the LAM domain boundaries


!===================================================================
! B: DERIVED TYPE DEFINITIONS
!===================================================================

!-------------------------------------------------------------------
! B1: Unrotated base panel for the construction of
!     LBC mesh ugrid generation strategy. (base_lbc_panel_type)
!     (PRIVATE OBJECT)
!-------------------------------------------------------------------
!
!     NW       Panel length
! Reference  O--------------+
!  Corner    |    Base      |   Panel width
!            |    Panel     | (LBC rim depth)
!            +--------------+
!
! Cell/Node/Edge id numbering ordered along rows beginning
! from NW corner of panel. Base panels treated in horizontal
! orientation until assembled into LBC mesh.
!-------------------------------------------------------------------
  type, private :: base_lbc_panel_type

    private

    integer(i_def) :: length  !> Length of base panel (cells)
    integer(i_def) :: depth   !> Depth of base panel (cells)
    integer(i_def) :: n_cells !> Number of cells in base panel
    integer(i_def) :: n_edges !> Number of edges in base panel
    integer(i_def) :: n_nodes !> Number of nodes in base panel

    integer(i_def) :: edges_per_face
    integer(i_def) :: nodes_per_face

    !> Cell-Cell connectivity array
    integer(i_def), allocatable :: cell_next(:,:)

    !> Local cell ids from NW-NE panel corners
    integer(i_def), allocatable :: north_cells(:)

    !> Local cell ids from SW-SE panel corners
    integer(i_def), allocatable :: south_cells(:)

    !> Local cell ids from NE-SE panel corners
    integer(i_def), allocatable :: east_cells(:)

    !> Local cell ids from NW-SW panel corners
    integer(i_def), allocatable :: west_cells(:)

  contains

    procedure :: base_lbc_panel_clear
    final     :: base_lbc_panel_final

  end type base_lbc_panel_type

!-------------------------------------------------------------------
! B2: Derived Type: LBC mesh ugrid generation strategy (gen_lbc_type)
!-------------------------------------------------------------------
! Ugrid mesh generator for LBC mesh topology / co-ordinates
! configured as a ring of cells from the boundaries of a
! parent LAM domain.
!
! LAM NW      LBC Panel #1      LAM NE
!        O---------------+----O
!        |      HP#1     |    |
!        +----+----------+    |
!        |    |          |    |
!  LBC   |    |          |VP#1| LBC
!  Panel |    |          |    | Panel
!   #4   |VP#2|          |    |  #2
!        |    |          |    |
!        |    +----------+----+
!        |    |      HP#2     |
!        O----+---------------O
! LAM SW      LBC Panel #2      LAM SE
!
!---------------------------------------------------------------
  type, extends(ugrid_generator_type), public :: gen_lbc_type

    private

    character(str_def) :: mesh_name

    integer(i_native)  :: geometry  = emdi
    integer(i_native)  :: topology  = emdi
    integer(i_native)  :: coord_sys = emdi
    character(str_def) :: coord_units_x
    character(str_def) :: coord_units_y

    integer(i_def) :: npanels = NPANELS
    integer(i_def) :: nmaps = 1      !> LBC meshes will only map
                                     !> to parent mesh
    real(r_def)    :: domain_size(2) !> Global domain size of LAM parent
    integer(i_def) :: outer_cells_x  !> Max number of cells in x-direction
    integer(i_def) :: outer_cells_y  !> Max number of cells in y-direction

    integer(i_def) :: rim_depth      !> Number of cells from LAM boundary to
                                     !> include in the LBC mesh

    character(str_longlong) :: constructor_inputs

    integer(i_def) :: n_nodes !> Number of nodes in this mesh
    integer(i_def) :: n_edges !> Number of edges in this mesh
    integer(i_def) :: n_faces !> Number of faces in this mesh

    !> Mesh names which this mesh has maps for
    character(str_def),                    allocatable :: target_mesh_names(:)
    type(global_mesh_map_collection_type), allocatable :: global_mesh_maps

    integer(i_def), allocatable :: cell_next(:,:)     ! (4, edge_cells_x*edge_cells_y)
    integer(i_def), allocatable :: nodes_on_cell(:,:) ! (4, edge_cells_x*edge_cells_y)
    integer(i_def), allocatable :: edges_on_cell(:,:) ! (4, edge_cells_x*edge_cells_y)
    integer(i_def), allocatable :: nodes_on_edge(:,:) ! (2, edge_cells_x*edge_cells_y)

    real(r_def),    allocatable :: node_coords(:,:)   ! (2, edge_cells_x*edge_cells_y)
    real(r_def),    allocatable :: cell_coords(:,:)   ! (2, edge_cells_x*edge_cells_y)

    ! Hold variables from the reference element.
    ! Done because the available Cray compiler has internal compiler errors
    ! when attempting to include the reference_element_type.
    integer(i_def) :: nodes_per_face
    integer(i_def) :: nodes_per_edge
    integer(i_def) :: edges_per_face

    integer(i_def) :: max_num_faces_per_node

    ! Used to generate LBC
    type(gen_planar_type), pointer :: lam_strategy => null()

    type(base_lbc_panel_type) :: base_h_panel
    type(base_lbc_panel_type) :: base_v_panel

    ! Information about the domain orientation
    real(r_def)    :: north_pole(2) = [rmdi,rmdi]  ! [Longitude, Latitude] of
                                                   ! north pole used
                                                   ! for domain orientation (degrees)
    real(r_def)    :: null_island(2) = [rmdi,rmdi] ! [Longitude, Latitude] of
                                                   ! null island used
                                                   ! for domain orientation (degrees)

  contains

    procedure :: generate
    procedure :: get_number_of_panels
    procedure :: get_metadata
    procedure :: get_dimensions
    procedure :: get_coordinates
    procedure :: get_connectivity
    procedure :: get_global_mesh_maps

    procedure :: clear
    final     :: gen_lbc_final

  end type gen_lbc_type

!-----------------------------------------------------------------------------

  interface gen_lbc_type
    module procedure gen_lbc_constructor
  end interface gen_lbc_type

  interface base_lbc_panel_type
    module procedure base_lbc_panel_constructor
  end interface base_lbc_panel_type

!-----------------------------------------------------------------------------



contains

!=============================================================================
! C: Module subroutines functions
!=============================================================================

!=============================================================================
! C1: Constructor/Finaliser for LBC mesh ugrid generator (gen_lbc_type)
!=============================================================================
! LAM NW      LBC Panel #1      LAM NE
!        O---------------+----O
!        |      HP#1     |    |
!        +----+----------+    |
!        |    |          |    |
!  LBC   |    |          |VP#1| LBC
!  Panel |    |          |    | Panel
!   #4   |VP#2|          |    |  #2
!        |    |          |    |
!        |    +----------+----+
!        |    |      HP#2     |
!        O----+---------------O
! LAM SW      LBC Panel #2      LAM SE
!-----------------------------------------------------------------------------
!> @brief   Constructor for LBC ugrid mesh generator (gen_lbc_type)
!> @details Constructs a LBC mesh based on the parent LAM mesh strategy and
!>          specified rim depth in cells.
!>
!> @param[in] lam_strategy  Planar mesh type object from which to generate a
!>                          LBC mesh ugrid generator.
!> @param[in] rim_depth     The depth of LBC rim (in cells).
!-----------------------------------------------------------------------------
function gen_lbc_constructor( lam_strategy, rim_depth ) result( self )

  implicit none

  type(gen_planar_type), intent(in), target :: lam_strategy
  integer(i_def),        intent(in) :: rim_depth

  type(gen_lbc_type) :: self

  character(str_def) :: mesh_parent
  character(str_def) :: geometry_key
  character(str_def) :: coord_sys_key

  integer(i_def) :: lam_n_nodes
  integer(i_def) :: lam_n_edges
  integer(i_def) :: lam_n_faces

  self%lam_strategy => lam_strategy
  self%rim_depth    =  rim_depth
  self%nmaps        =  1

  allocate( self%target_mesh_names(self%nmaps))
  allocate( self%global_mesh_maps, source=global_mesh_map_collection_type())

  ! 1.0 Interrogate LAM Generator strategy for LBC dimensions
  !----------
  call self%lam_strategy%get_dimensions                                   &
                              ( num_nodes = lam_n_nodes,                  &
                                num_edges = lam_n_edges,                  &
                                num_faces = lam_n_faces,                  &
                                num_nodes_per_face = self%nodes_per_face, &
                                num_edges_per_face = self%edges_per_face, &
                                num_nodes_per_edge = self%nodes_per_edge, &
                                max_num_faces_per_node =                  &
                                             self%max_num_faces_per_node )

  call lam_strategy%get_metadata                                    &
                        ( mesh_name    = self%target_mesh_names(1), &
                          edge_cells_x = self%outer_cells_x,        &
                          edge_cells_y = self%outer_cells_y,        &
                          north_pole   = self%north_pole,           &
                          null_island  = self%null_island,          &
                          geometry     = geometry_key,              &
                          coord_sys    = coord_sys_key,             &
                          domain_size  = self%domain_size )

  mesh_parent    = trim(self%target_mesh_names(1))
  self%mesh_name = trim(mesh_parent)//'-lbc'

  self%coord_sys   =  coord_sys_from_key(coord_sys_key)
  self%geometry    =  geometry_from_key(geometry_key)
  self%topology    =  topology_non_periodic



  ! 2.0 Check the rim depth is valid
  !----------
  if (self%rim_depth <= 0) then
    write(log_scratch_space,'(A)') &
        'LBC mesh requires rim depth > 0 cells.'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  else if ( (self%outer_cells_x - 2*self%rim_depth < 1) .or. &
            (self%outer_cells_y - 2*self%rim_depth < 1) ) then
    write(log_scratch_space,'(A,I0,A)')      &
        'LBC rim depth of ', self%rim_depth, &
        ' will erode entire LAM inner domain.'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  self%base_h_panel =                                           &
      base_lbc_panel_type( self%outer_cells_x - self%rim_depth, &
                           self%rim_depth )
  self%base_v_panel =                                           &
      base_lbc_panel_type( self%outer_cells_y - self%rim_depth, &
                           self%rim_depth )

  ! 3.0 Calculate the number of cells in the LBC mesh
  !----------
  self%n_faces = 2*(self%base_h_panel%n_cells + self%base_v_panel%n_cells)
  self%n_edges = 2*(self%base_h_panel%n_edges + self%base_v_panel%n_edges)
  self%n_nodes = 2*(self%base_h_panel%n_nodes + self%base_v_panel%n_nodes)

  write(self%constructor_inputs,'(A,I0)')       &
      'lam_strategy=<gen_planar_type,"'//       &
      trim(self%target_mesh_names(1))//'">;'//  &
      'rim_depth=', rim_depth

  return
end function gen_lbc_constructor


!-----------------------------------------------------------------------------
!> @brief Finaliser for the LBC ugrid mesh generator (gen_lbc_type)
!-----------------------------------------------------------------------------
subroutine gen_lbc_final(self)

  implicit none

  type(gen_lbc_type), intent(inout) :: self

  call self%clear()

  return
end subroutine gen_lbc_final


!=============================================================================
! C2: Constructor/Finaliser for LBC base panel (base_lbc_panel_type)
!=============================================================================
!> @brief Constructor for LBC mesh base panel object
!> @param[in] panel_length  Base panel length (in cells)
!> @param[in] panel_depth   Base panel depth (in cells).
!-----------------------------------------------------------------------------
function base_lbc_panel_constructor(panel_length, panel_depth) result(self)

  implicit none

  integer(i_def), intent(in) :: panel_length
  integer(i_def), intent(in) :: panel_depth

  type(base_lbc_panel_type) :: self

  integer(i_def) :: base_id
  integer(i_def) :: cell

  self%length = panel_length
  self%depth  = panel_depth
  self%edges_per_face = 4
  self%nodes_per_face = 4

  self%n_cells = self%length * self%depth
  self%n_edges = self%length * ( 2*self%depth + 1 )
  self%n_nodes = self%length * (   self%depth + 1 )

  allocate( self%north_cells ( self%length ) )
  allocate( self%south_cells ( self%length ) )
  allocate( self%east_cells  ( self%depth  ) )
  allocate( self%west_cells  ( self%depth  ) )
  allocate( self%cell_next   ( 4, self%n_cells) )

  ! 1.0 Capture cells on edge of panel
  !----------
  do cell=1, self%length
    self%north_cells(cell) = cell
    self%south_cells(cell) = cell + (self%depth-1)*self%length
  end do

  do cell=1, self%depth
    self%east_cells(cell) = cell*self%length
    self%west_cells(cell) = (cell-1)*self%length + 1
  end do

  ! 2.0 Assign default cell numbering
  !----------
  base_id = 0
  self%cell_next = calc_panel_adjacency( self, &
                                         base_id )

  return
end function base_lbc_panel_constructor

!-----------------------------------------------------------------------------
!> @brief Finaliser for LBC mesh base panel object
!-----------------------------------------------------------------------------
subroutine base_lbc_panel_final(self)
  implicit none

  type(base_lbc_panel_type), intent(inout) :: self
  call self%base_lbc_panel_clear()

  return
end subroutine base_lbc_panel_final



!=============================================================================
! D1: Methods for LBC mesh ugrid generator (gen_lbc_type)
!=============================================================================

!-----------------------------------------------------------------------------
!> @brief Populates the arguments with the dimensions defining
!>        the planar mesh.
!>
!> @param[out]  num_nodes              The number of nodes on the mesh.
!> @param[out]  num_edges              The number of edges on the mesh.
!> @param[out]  num_faces              The number of faces on the mesh.
!> @param[out]  num_nodes_per_face     The number of nodes around each face.
!> @param[out]  num_edges_per_face     The number of edges around each face.
!> @param[out]  num_nodes_per_face     The number of nodes around each edge.
!> @param[out]  max_num_faces_per_node The maximum number of faces surrounding
!>                                     each node.
!-----------------------------------------------------------------------------
subroutine get_dimensions( self,                                       &
                           num_nodes, num_edges, num_faces,            &
                           num_nodes_per_face, num_edges_per_face,     &
                           num_nodes_per_edge, max_num_faces_per_node )
  implicit none

  class(gen_lbc_type), intent(in) :: self

  integer(i_def), intent(out) :: num_nodes
  integer(i_def), intent(out) :: num_edges
  integer(i_def), intent(out) :: num_faces
  integer(i_def), intent(out) :: num_nodes_per_face
  integer(i_def), intent(out) :: num_edges_per_face
  integer(i_def), intent(out) :: num_nodes_per_edge
  integer(i_def), intent(out) :: max_num_faces_per_node

  num_nodes = self%n_nodes
  num_edges = self%n_edges
  num_faces = self%n_faces

  num_nodes_per_face = self%nodes_per_face
  num_edges_per_face = self%edges_per_face
  num_nodes_per_edge = self%nodes_per_edge

  max_num_faces_per_node = self%max_num_faces_per_node

  return
end subroutine get_dimensions

!-----------------------------------------------------------------------------
!> @brief   Populates the argument array with the coordinates of the
!>          mesh's nodes.
!> @details Exposes the instance's node_coords array to the caller.
!>
!> @param[out]  node_coordinates  Node coordinates.
!> @param[out]  cell_coordinates  Cell centre coordinates.
!> @param[out]  coord_units_x     Units for x-coordinates
!> @param[out]  coord_units_y     Units for y-coordinates
!-----------------------------------------------------------------------------
subroutine get_coordinates( self,             &
                            node_coordinates, &
                            cell_coordinates, &
                            coord_units_x,    &
                            coord_units_y )

  implicit none

  class(gen_lbc_type), intent(in)  :: self

  real(r_def),         intent(out) :: node_coordinates(:,:)
  real(r_def),         intent(out) :: cell_coordinates(:,:)
  character(str_def),  intent(out) :: coord_units_x
  character(str_def),  intent(out) :: coord_units_y

  node_coordinates = self%node_coords
  cell_coordinates = self%cell_coords
  coord_units_x    = self%coord_units_x
  coord_units_y    = self%coord_units_y

  return
end subroutine get_coordinates


!-----------------------------------------------------------------------------
!> @brief   Populates the argument arrays with the corresponding mesh
!>          connectivity information.
!> @details Implements the connectivity-providing interface required
!>          by the ugrid writer.
!>
!> @param[out]  face_node_connectivity  Face-node connectivity.
!> @param[out]  face_edge_connectivity  Face-edge connectivity.
!> @param[out]  face_face_connectivity  Face-face connectivity.
!> @param[out]  edge_node_connectivity  Edge-node connectivity.
!-----------------------------------------------------------------------------
subroutine get_connectivity( self,                   &
                             face_node_connectivity, &
                             face_edge_connectivity, &
                             face_face_connectivity, &
                             edge_node_connectivity )

  implicit none

  class(gen_lbc_type), intent(in) :: self

  integer(i_def), intent(out) :: face_node_connectivity(:,:)
  integer(i_def), intent(out) :: face_edge_connectivity(:,:)
  integer(i_def), intent(out) :: face_face_connectivity(:,:)
  integer(i_def), intent(out) :: edge_node_connectivity(:,:)

  face_node_connectivity = self%nodes_on_cell
  face_edge_connectivity = self%edges_on_cell
  face_face_connectivity = self%cell_next
  edge_node_connectivity = self%nodes_on_edge

  return
end subroutine get_connectivity


!-----------------------------------------------------------------------------
!> @brief   Generates the mesh and connectivity.
!> @details Calls each of the instance methods which calculate the
!>          specified mesh and populate its arrays.
!-----------------------------------------------------------------------------
subroutine generate(self)

  ! Treat each panel as a non-periodic panel
  ! To be reconnected in the following orientation
  !==============================================================
  ! LAM NW      LBC Panel #1      LAM NE
  !        O---------------+----O
  !        |      HP#1     |    |
  !        +----+----------+    |
  !        |    |          |    |
  !  LBC   |    |          |VP#1| LBC
  !  Panel |    |          |    | Panel
  !   #4   |VP#2|          |    |  #2
  !        |    |          |    |
  !        |    +----------+----+
  !        |    |      HP#2     |
  !        O----+---------------O
  ! LAM SW      LBC Panel #2      LAM SE

  implicit none

  class(gen_lbc_type), intent(inout) :: self

  integer(i_def) :: lam_n_nodes
  integer(i_def) :: lam_n_edges
  integer(i_def) :: lam_n_faces

  integer(i_def), allocatable :: lam_face_node(:,:)
  integer(i_def), allocatable :: lam_face_edge(:,:)
  integer(i_def), allocatable :: lam_edge_node(:,:)
  integer(i_def), allocatable :: lam_face_face(:,:)

  real(r_def), allocatable :: lam_node_coords(:,:)
  real(r_def), allocatable :: lam_cell_coords(:,:)

  ! 0.0 Has the LAM strategy been generated.
  !----------
  if ( .not. self%lam_strategy%is_generated() ) then
    call self%lam_strategy%generate()
  end if

  ! 1.0 Extract data from the lam_strategy to be
  !     passed to each generation routine
  !----------
  ! 1.1 LAM dimensions
  call self%lam_strategy%get_dimensions(                               &
                             num_nodes = lam_n_nodes,                  &
                             num_edges = lam_n_edges,                  &
                             num_faces = lam_n_faces,                  &
                             num_nodes_per_face = self%nodes_per_face, &
                             num_edges_per_face = self%edges_per_face, &
                             num_nodes_per_edge = self%nodes_per_edge, &
                             max_num_faces_per_node =                  &
                                          self%max_num_faces_per_node )

  ! 1.2 LAM coordinates
  allocate( lam_node_coords( 2, lam_n_nodes) )
  allocate( lam_cell_coords( 2, lam_n_faces) )
  call self%lam_strategy%get_coordinates(                      &
                             node_coordinates=lam_node_coords, &
                             cell_coordinates=lam_cell_coords, &
                             coord_units_x=self%coord_units_x, &
                             coord_units_y=self%coord_units_y )

  ! 1.3 LAM connectivity
  allocate( lam_face_node(self%nodes_per_face, lam_n_faces) )
  allocate( lam_face_face(self%edges_per_face, lam_n_faces) )
  allocate( lam_face_edge(self%edges_per_face, lam_n_faces) )
  allocate( lam_edge_node(self%nodes_per_edge, lam_n_edges) )

  call self%lam_strategy%get_connectivity(                           &
                             face_node_connectivity = lam_face_node, &
                             face_edge_connectivity = lam_face_edge, &
                             face_face_connectivity = lam_face_face, &
                             edge_node_connectivity = lam_edge_node )

  ! 2.0 Generate the LBC
  !----------
  call calc_adjacency( self )
  call calc_cell_edge_connectivity( self )
  call calc_cell_node_connectivity( self )
  call calc_edge_node_connectivity( self )

  ! 2.1 Generate LBC-LAM intergrid map
  call calc_lbc_lam_map(self, lam_face_face)

  ! 2.2 Extract coordinates for LBC mesh entities
  call extract_coords( self,            &
                       lam_node_coords, &
                       lam_cell_coords, &
                       lam_face_node )

  ! 3.0 Clear used objects once LBC has been generated
  !----------
  call self%base_h_panel%base_lbc_panel_clear()
  call self%base_v_panel%base_lbc_panel_clear()

  return
end subroutine generate

!-----------------------------------------------------------------------------
!> @brief Returns the number of panels in the mesh topology.
!> @description Panels are a subset of cells in the mesh domain which may
!>              exhibit common properties.
!> @return answer Number of panels resulting from this generation strategy.
!-----------------------------------------------------------------------------
function get_number_of_panels( self ) result( answer )

  implicit none

  class(gen_lbc_type), intent(in) :: self

  integer(i_def) :: answer

  answer = NPANELS

end function get_number_of_panels


!-----------------------------------------------------------------------------
!> @brief Returns mesh metadata information.
!> @details This subroutine is provided as a means to request specific metadata
!>          from the current mesh configuration.
!>
!> @param[out] mesh_name           Optional, Name of mesh instance to generate
!> @param[out] geometry            Optional, Mesh domain surface type.
!> @param[out] topology            Optional, Mesh boundary/connectivity type
!> @param[out] coord_sys           Optional, Coordinate system to position nodes.
!> @param[out] periodic_x          Optional, Periodic in E-W direction.
!> @param[out] periodic_y          Optional, Periodic in N-S direction.
!> @param[out] edge_cells_x        Optional, Number of panel edge cells (x-axis).
!> @param[out] edge_cells_y        Optional, Number of panel edge cells (y-axis).
!> @param[out] constructor_inputs  Optional, Inputs used to create this mesh from
!>                                           the mesh_generator
!> @param[out] nmaps               Optional, Number of maps to create with this mesh
!>                                           as source mesh
!> @param[out] rim_depth           Optional, Rim depth of LBC mesh (LAMs).
!> @param[out] domain_size         Optional, Size of global model domain.
!> @param[out] void_cell           Optional, Cell ID used for null connectivity.
!> @param[out] target_mesh_names   Optional, Mesh names of the target meshes that
!>                                           this mesh has maps for.
!> @param[out] maps_edge_cells_x   Optional, Number of panel edge cells (x-axis) of
!>                                           target mesh(es) to create map(s) for.
!> @param[out] maps_edge_cells_y   Optional, Number of panel edge cells (y-axis) of
!>                                           target mesh(es) to create map(s) for.
!> @param[out] north_pole          Optional, [Longitude, Latitude] of north pole
!>                                           used for domain orientation (degrees)
!> @param[out] null_island         Optional, [Longitude, Latitude] of null
!>                                           island used for domain orientation (degrees)
!-----------------------------------------------------------------------------
subroutine get_metadata( self,               &
                         mesh_name,          &
                         geometry,           &
                         topology,           &
                         coord_sys,          &
                         periodic_x,         &
                         periodic_y,         &
                         edge_cells_x,       &
                         edge_cells_y,       &
                         constructor_inputs, &
                         nmaps,              &
                         rim_depth,          &
                         domain_size,        &
                         void_cell,          &
                         target_mesh_names,  &
                         maps_edge_cells_x,  &
                         maps_edge_cells_y,  &
                         north_pole,         &
                         null_island )
  implicit none

  class(gen_lbc_type), intent(in) :: self

  character(str_def), optional, intent(out) :: mesh_name
  character(str_def), optional, intent(out) :: geometry
  character(str_def), optional, intent(out) :: topology
  character(str_def), optional, intent(out) :: coord_sys

  logical(l_def), optional, intent(out) :: periodic_x
  logical(l_def), optional, intent(out) :: periodic_y
  integer(i_def), optional, intent(out) :: edge_cells_x
  integer(i_def), optional, intent(out) :: edge_cells_y
  integer(i_def), optional, intent(out) :: nmaps
  integer(i_def), optional, intent(out) :: rim_depth
  integer(i_def), optional, intent(out) :: void_cell
  real(r_def),    optional, intent(out) :: domain_size(2)

  character(str_longlong), optional, intent(out) :: constructor_inputs

  character(str_def), optional, allocatable, intent(out) :: target_mesh_names(:)
  integer(i_def),     optional, allocatable, intent(out) :: maps_edge_cells_x(:)
  integer(i_def),     optional, allocatable, intent(out) :: maps_edge_cells_y(:)

  real(r_def), optional, intent(out) :: north_pole(2)
  real(r_def), optional, intent(out) :: null_island(2)


  if (present(mesh_name))    mesh_name    = self%mesh_name
  if (present(geometry))     geometry     = key_from_geometry(self%geometry)
  if (present(topology))     topology     = key_from_topology(self%topology)
  if (present(coord_sys))    coord_sys    = key_from_coord_sys(self%coord_sys)
  if (present(periodic_x))   periodic_x   = .false.
  if (present(periodic_y))   periodic_y   = .false.
  if (present(edge_cells_x)) edge_cells_x = self%outer_cells_x
  if (present(edge_cells_y)) edge_cells_y = self%outer_cells_y
  if (present(nmaps))        nmaps        = self%nmaps
  if (present(rim_depth))    rim_depth    = self%rim_depth
  if (present(void_cell))    void_cell    = VOID_ID

  if (present(constructor_inputs)) constructor_inputs = self%constructor_inputs

  if (present(edge_cells_x))    edge_cells_x = self%outer_cells_x
  if (present(edge_cells_y))    edge_cells_y = self%outer_cells_y
  if (present(nmaps))           nmaps        = self%nmaps

  if (present(void_cell)) void_cell = VOID_ID

  if (self%nmaps > 0) then
    if (present(target_mesh_names)) then
      allocate(target_mesh_names(1), source=self%target_mesh_names)
    end if

    if (present(maps_edge_cells_x)) then
      allocate(maps_edge_cells_x(1), source=self%outer_cells_x)
    end if

    if (present(maps_edge_cells_y)) then
      allocate(maps_edge_cells_y(1), source=self%outer_cells_y)
    end if
  end if

  ! These are inherited from LAM so should be in degrees for cf-compliance
  if (present(domain_size)) domain_size(:) = self%domain_size(:)
  if (present(north_pole))  north_pole(:)  = self%north_pole(:)
  if (present(null_island)) null_island(:) = self%null_island(:)

  return
end subroutine get_metadata


!-------------------------------------------------------------------------------
!> @brief  Gets the global mesh map collection which uses
!>         this mesh as the source mesh.
!>
!> @return global_mesh_maps global_mesh_map_collection_type
!-------------------------------------------------------------------------------
function get_global_mesh_maps(self) result(global_mesh_maps)

  implicit none

  class(gen_lbc_type), target, intent(in) :: self

  type(global_mesh_map_collection_type), pointer :: global_mesh_maps

  global_mesh_maps => self%global_mesh_maps

  return
end function get_global_mesh_maps


!-------------------------------------------------------------------------------
!> @brief Subroutine to manually deallocate any memory used by the object.
!-------------------------------------------------------------------------------
subroutine clear(self)

  implicit none

  class(gen_lbc_type), intent(inout) :: self

  if (allocated(self%cell_next))     deallocate( self%cell_next     )
  if (allocated(self%nodes_on_cell)) deallocate( self%nodes_on_cell )
  if (allocated(self%edges_on_cell)) deallocate( self%edges_on_cell )
  if (allocated(self%nodes_on_edge)) deallocate( self%nodes_on_edge )
  if (allocated(self%node_coords))   deallocate( self%node_coords   )
  if (allocated(self%cell_coords))   deallocate( self%cell_coords   )

  if (allocated(self%target_mesh_names)) deallocate( self%target_mesh_names )

  if (allocated(self%global_mesh_maps)) then
    call self%global_mesh_maps%clear()
    deallocate( self%global_mesh_maps )
  end if

  call self%base_h_panel%base_lbc_panel_clear()
  call self%base_v_panel%base_lbc_panel_clear()

  self%lam_strategy => null()

  return
end subroutine clear



!===============================================================================
! Module private subroutines LBC mesh ugrid generator (gen_lbc_type)
!===============================================================================
! These are private subroutines which are only called by the generate method
! when the object data is being generated.
!-------------------------------------------------------------------------------
!> @brief Generates the LBC mesh connectivity (adjacent cells)
!>        (PRIVATE SUBROUTINE)
!-------------------------------------------------------------------------------
subroutine calc_adjacency( self )

  implicit none

  type(gen_lbc_type) :: self

  integer(i_def), allocatable :: h_panel_cell_next(:,:,:)
  integer(i_def), allocatable :: v_panel_cell_next(:,:,:)
  integer(i_def), allocatable :: v_panel_cell_ids(:,:)
  integer(i_def), allocatable :: h_panel_cell_ids(:,:)

  integer(i_def) :: base_id, cell_id
  integer(i_def) :: cell_1, cell_2
  integer(i_def) :: i, j, cell, panel

  ! Treat each panel as a non-periodic panel
  ! To be reconnected in the following orientation
  !--------------------------------------------------------
  !
  !   +---------------+----+
  !   |      HP#1     |    |
  !   +----+----------+    |
  !   |    |          |    |
  !   |    |          |VP#1|
  !   |    |          |    |
  !   |VP#2|          |    |
  !   |    |          |    |
  !   |    +----------+----+
  !   |    |      HP#2     |
  !   +----+---------------+

  allocate( h_panel_cell_next( self%edges_per_face,       &
                               self%base_h_panel%n_cells, &
                               horz_npanels ) )
  allocate( v_panel_cell_next( self%edges_per_face,       &
                               self%base_v_panel%n_cells, &
                               vert_npanels ) )

  allocate( h_panel_cell_ids ( self%base_h_panel%n_cells, &
                               horz_npanels ) )
  allocate( v_panel_cell_ids ( self%base_v_panel%n_cells, &
                               vert_npanels ) )

  ! 1.0 Group cell ids on each panel to make it easier
  !     to connect panels later on.
  !----------
  base_id = 0
  do panel=1, 2
    ! 1.1 Horizontal panels
    do cell=1, self%base_h_panel%n_cells
      base_id = base_id + 1
      h_panel_cell_ids(cell, panel) = base_id
    end do

    ! 1.2 Vertical panels
    do cell=1, self%base_v_panel%n_cells
      base_id = base_id + 1
      v_panel_cell_ids(cell, panel) = base_id
    end do
  end do


  ! 2.0 Generate the adjacency for each panel
  !----------
  base_id = 0
  do i=1, 2
    ! 2.1 Horizontal panels
    h_panel_cell_next(:,:,i) =                   &
        calc_panel_adjacency( self%base_h_panel, &
                              base_id )

    base_id = base_id + self%base_h_panel%n_cells

    ! 2.2 Vertical panels
    v_panel_cell_next(:,:,i) =                   &
        calc_panel_adjacency( self%base_v_panel, &
                              base_id )

    base_id = base_id + self%base_v_panel%n_cells

  end do

  ! 3.0 Connect east edge of panels.
  !----------
  do j=1, self%rim_depth

    ! 3.1 Vertical joins in LBC
    cell_1 = self%base_h_panel%east_cells(j)
    cell_2 = self%base_v_panel%south_cells(j)

    ! Panels h1:v1
    h_panel_cell_next(E,cell_1,1) = v_panel_cell_ids(cell_2, 1)
    v_panel_cell_next(S,cell_2,1) = h_panel_cell_ids(cell_1, 1)

    ! Panels h2:v2
    h_panel_cell_next(E,cell_1,2) = v_panel_cell_ids(cell_2, 2)
    v_panel_cell_next(S,cell_2,2) = h_panel_cell_ids(cell_1, 2)

    ! 3.2 Horizontal joins in LBC
    cell_1 = self%base_v_panel%east_cells(j)
    cell_2 = self%base_h_panel%south_cells(j)

    ! Panels v1:h2
    v_panel_cell_next(E,cell_1,1) = h_panel_cell_ids(cell_2, 2)
    h_panel_cell_next(S,cell_2,2) = v_panel_cell_ids(cell_1, 1)

    ! Panels v2:h1
    v_panel_cell_next(E,cell_1,2) = h_panel_cell_ids(cell_2, 1)
    h_panel_cell_next(S,cell_2,1) = v_panel_cell_ids(cell_1, 2)

  end do


  ! 4.0 Rotate panels to use the same reference for North
  !----------
  ! Cyclic rotate panels 2:4, hp#2 and the vertcal panels
  do i=1, self%base_v_panel%n_cells
    v_panel_cell_next(:, i, 1) = cshift(v_panel_cell_next(:,i,1), 1, 1)
    v_panel_cell_next(:, i, 2) = cshift(v_panel_cell_next(:,i,2), 3, 1)
  end do

  do i=1, self%base_h_panel%n_cells
    h_panel_cell_next(:, i, 2) = cshift(h_panel_cell_next(:,i,2), 2, 1)
  end do


  ! 5.0 Reshape cell next arrays to the familar format [entity, n_cells]
  !----------
  if (allocated(self%cell_next)) deallocate(self%cell_next)
  allocate(self%cell_next(self%edges_per_face, self%n_faces))

  cell_id = 0
  do panel=1, 2
    do i=1, self%base_h_panel%n_cells
      cell_id = cell_id+1
      self%cell_next(:,cell_id) = h_panel_cell_next(:,i,panel)
    end do
    do i=1, self%base_v_panel%n_cells
      cell_id = cell_id+1
      self%cell_next(:,cell_id) = v_panel_cell_next(:,i,panel)
    end do
  end do

  deallocate( h_panel_cell_next, &
              v_panel_cell_next, &
              v_panel_cell_ids,  &
              h_panel_cell_ids )

  return
end subroutine calc_adjacency



!-------------------------------------------------------------------------------
!> @brief Generates the LBC mesh cell-node connectivity (nodes on a given cell)
!>        (PRIVATE SUBROUTINE)
!-------------------------------------------------------------------------------
subroutine calc_cell_node_connectivity(self)

  implicit none

  type(gen_lbc_type) :: self

  integer(i_def) :: panel, node_id, cell_id, i

  integer(i_def), allocatable :: hp_node_on_cells(:,:,:)
  integer(i_def), allocatable :: vp_node_on_cells(:,:,:)

  allocate( hp_node_on_cells( self%nodes_per_face,       &
                              self%base_h_panel%n_cells, &
                              horz_npanels) )

  allocate( vp_node_on_cells( self%nodes_per_face,       &
                              self%base_v_panel%n_cells, &
                              vert_npanels ) )

  hp_node_on_cells = imdi
  vp_node_on_cells = imdi

  node_id = 0
  do panel=1, 2
    hp_node_on_cells(:,:,panel) = &
        calc_panel_node_connectivity( self%base_h_panel, base_id=node_id )
    vp_node_on_cells(:,:,panel) = &
        calc_panel_node_connectivity( self%base_v_panel, base_id=node_id )
  end do


  ! 1.0 Connect the nodes on the east edge of the panels
  !----------
  ! 1.1 Horizontal panels
  do panel =1, horz_npanels
    do i=1, self%rim_depth
      hp_node_on_cells(NE, self%base_h_panel%east_cells(i), panel) = &
          vp_node_on_cells(SW, self%base_v_panel%south_cells(i), panel )
      hp_node_on_cells(SE, self%base_h_panel%east_cells(i), panel) = &
          vp_node_on_cells(SE, self%base_v_panel%south_cells(i), panel )
    end do
  end do

  ! 1.2 Vertical panels
  do i=1, self%rim_depth
    vp_node_on_cells(NE, self%base_v_panel%east_cells(i), 1) = &
        hp_node_on_cells(SW, self%base_h_panel%south_cells(i), 2 )
    vp_node_on_cells(SE, self%base_v_panel%east_cells(i), 1) = &
        hp_node_on_cells(SE, self%base_h_panel%south_cells(i), 2 )
    vp_node_on_cells(NE, self%base_v_panel%east_cells(i), 2) = &
        hp_node_on_cells(SW, self%base_h_panel%south_cells(i), 1 )
    vp_node_on_cells(SE, self%base_v_panel%east_cells(i), 2) = &
        hp_node_on_cells(SE, self%base_h_panel%south_cells(i), 1 )
  end do


  ! 2.0 Rotate panels to use the same reference for North
  !----------
  ! Cyclic rotate panels 2:4, hp#2, vertical panels
  do i=1, self%base_v_panel%n_cells
    vp_node_on_cells(:, i, 1)  = cshift(vp_node_on_cells(:,i,1), 1, 1)
    vp_node_on_cells(:, i, 2)  = cshift(vp_node_on_cells(:,i,2), 3, 1)
  end do

  do i=1, self%base_h_panel%n_cells
    hp_node_on_cells(:, i, 2)  = cshift(hp_node_on_cells(:,i,2), 2, 1)
  end do


  ! 3.0 Reshape cell next arrays to the familar format [entity, n_cells]
  !----------
  if (allocated(self%nodes_on_cell)) deallocate(self%nodes_on_cell)
  allocate( self%nodes_on_cell(self%nodes_per_face, self%n_faces) )
  cell_id=0
  do panel=1, 2
    do i=1, self%base_h_panel%n_cells
      cell_id = cell_id+1
      self%nodes_on_cell(:,cell_id) = hp_node_on_cells(:,i, panel)
    end do
    do i=1, self%base_v_panel%n_cells
      cell_id = cell_id+1
      self%nodes_on_cell(:,cell_id) = vp_node_on_cells(:,i, panel)
    end do
  end do

  return
end subroutine calc_cell_node_connectivity



!-------------------------------------------------------------------------------
!> @brief Generates the LBC mesh cell-edge connectivity, (edges on a given cell)
!>        (PRIVATE SUBROUTINE)
!-------------------------------------------------------------------------------
subroutine calc_cell_edge_connectivity(self)

  implicit none

  type(gen_lbc_type) :: self

  integer(i_def) :: panel, edge_id, cell_id,i

  integer(i_def), allocatable :: hp_edge_on_cells(:,:,:)
  integer(i_def), allocatable :: vp_edge_on_cells(:,:,:)


  allocate( hp_edge_on_cells( self%edges_per_face,       &
                              self%base_h_panel%n_cells, &
                              horz_npanels) )

  allocate( vp_edge_on_cells( self%edges_per_face,       &
                              self%base_v_panel%n_cells, &
                              vert_npanels ) )

  hp_edge_on_cells = imdi
  vp_edge_on_cells = imdi

  edge_id = 0
  do panel=1, 2
    hp_edge_on_cells(:,:,panel) = &
        calc_panel_edge_connectivity(self%base_h_panel, base_id=edge_id)
    vp_edge_on_cells(:,:,panel) = &
        calc_panel_edge_connectivity(self%base_v_panel, base_id=edge_id)
  end do


  ! 1.0 Connect the edges on the east edge of the panels
  !----------
  ! 1.1 Horizontal panels
  do panel=1, horz_npanels
    do i=1, self%rim_depth
      hp_edge_on_cells(E, self%base_h_panel%east_cells(i), panel) = &
          vp_edge_on_cells(S, self%base_v_panel%south_cells(i), panel )
    end do
  end do

  ! 1.2 Vertical panels
  do i=1, self%rim_depth
    vp_edge_on_cells(E, self%base_v_panel%east_cells(i), 1) = &
        hp_edge_on_cells(S, self%base_h_panel%south_cells(i), 2 )
    vp_edge_on_cells(E, self%base_v_panel%east_cells(i), 2) = &
        hp_edge_on_cells(S, self%base_h_panel%south_cells(i), 1 )
  end do


  ! 2.0 Rotate panels to use the same reference for North
  !----------
  ! Cyclic rotate panels 2:4, hp#2, vertical panels
  do i=1, self%base_v_panel%n_cells
    vp_edge_on_cells(:,i,1) = cshift( vp_edge_on_cells(:,i,1),1,1 )
    vp_edge_on_cells(:,i,2) = cshift( vp_edge_on_cells(:,i,2),3,1 )
  end do

  do i=1, self%base_h_panel%n_cells
    hp_edge_on_cells(:,i,2) = cshift( hp_edge_on_cells(:,i,2),2,1 )
  end do


  ! 3.0 Reshape cell next arrays to the familar format [entity, n_cells]
  !----------
  if (allocated(self%edges_on_cell)) deallocate(self%edges_on_cell)
  allocate(self%edges_on_cell(4,self%n_faces))

  cell_id=0
  do panel=1, 2
    do i=1, self%base_h_panel%n_cells
      cell_id = cell_id+1
      self%edges_on_cell(:,cell_id) = hp_edge_on_cells(:,i, panel)
    end do
    do i=1, self%base_v_panel%n_cells
      cell_id = cell_id+1
      self%edges_on_cell(:,cell_id) = vp_edge_on_cells(:,i, panel)
    end do
  end do

  return
end subroutine calc_cell_edge_connectivity



!-------------------------------------------------------------------------------
!> @brief Generates the LBC mesh edge-node connectivity
!>        (nodes on a given edge)  (PRIVATE SUBROUTINE)
!-------------------------------------------------------------------------------
subroutine calc_edge_node_connectivity(self)

  implicit none

  type(gen_lbc_type), intent(inout) :: self

  integer(i_def) :: cell_id, edge_id
  integer(i_def) :: i

  integer(i_def) :: edge_order(4), edge_index
  integer(i_def) :: edge_node_link(2,4)

  edge_order = [W,S,E,N]
  edge_node_link(:,W) = [NW,SW]
  edge_node_link(:,S) = [SW,SE]
  edge_node_link(:,E) = [NE,SE]
  edge_node_link(:,N) = [NW,NE]

  if (allocated(self%nodes_on_edge)) deallocate(self%nodes_on_edge)
  allocate(self%nodes_on_edge(2,self%n_edges))

  self%nodes_on_edge = imdi
  do cell_id=1, self%n_faces
    do i=1, self%edges_per_face

      edge_index = edge_order(i)
      edge_id    = self%edges_on_cell(edge_index, cell_id)
      if (any(self%nodes_on_edge(:,edge_id) == imdi)) then
        self%nodes_on_edge(1,edge_id) =                       &
            self%nodes_on_cell( edge_node_link(1,edge_index), &
                                cell_id )
        self%nodes_on_edge(2,edge_id) =                       &
            self%nodes_on_cell( edge_node_link(2,edge_index), &
                                cell_id )
      end if

    end do
  end do

  return
end subroutine calc_edge_node_connectivity



!-----------------------------------------------------------------------------
!> @brief Generates the LBC-LAM intergrid map. (PRIVATE SUBROUTINE)
!> @param[in] lam_face_face  The cell connectivity of the parent LAM
!-----------------------------------------------------------------------------
subroutine calc_lbc_lam_map(self, lam_face_face)

  implicit none

  type(gen_lbc_type), intent(inout) :: self
  integer(i_def),     intent(in)    :: lam_face_face(:,:)

  integer(i_def), allocatable :: cell_map(:,:,:)
  integer(i_def), allocatable :: h_lbc_lam_map(:,:)
  integer(i_def), allocatable :: v_lbc_lam_map(:,:)

  integer(i_def) :: panel, cell_id, lbc_cell_id

  integer(i_def) :: h_panel_lam_corners(2)
  integer(i_def) :: v_panel_lam_corners(2)
  integer(i_def) :: ref_corner, lam_ref_corner_id

  integer(i_def) :: source_id
  integer(i_def) :: target_id

  integer(i_def) :: ntarget_per_source_x = 1
  integer(i_def) :: ntarget_per_source_y = 1

  integer(i_def), allocatable :: lam_cell_next(:,:)

  allocate( lam_cell_next, source=lam_face_face )
  allocate( h_lbc_lam_map(self%base_h_panel%n_cells,2) )
  allocate( v_lbc_lam_map(self%base_v_panel%n_cells,2) )

  h_panel_lam_corners = [NW,SE]
  v_panel_lam_corners = [NE,SW]

  ! 1.0 Assign the LBC-LAM map to each LBC panel
  !     Starting from the panel reference corner and
  !     using tha LAM cell connectivity. The LAM connectivity
  !     is cyclically rotated to that it aligns with each base
  !     panel.
  !----------
  do panel=1, 2

    ! 1.1 Horizontal panel LBC-LAM cell map assignment
    if (panel /= 1) lam_cell_next(:,:) = cshift(lam_cell_next(:,:), -1, 1)
    ref_corner = h_panel_lam_corners(panel)
    lam_ref_corner_id = &
        self%lam_strategy%get_corner_gid( h_panel_lam_corners(panel) )

    h_lbc_lam_map(:,panel) =                       &
        calc_panel_lbc_lam_map( self%base_h_panel, &
                                lam_ref_corner_id, &
                                lam_cell_next )



    ! 1.2 Vertical panel LBC-LAM cell map assignment
    ref_corner = v_panel_lam_corners(panel)
    lam_cell_next(:,:) = cshift(lam_cell_next(:,:), -1, 1)
    lam_ref_corner_id = &
        self%lam_strategy%get_corner_gid( v_panel_lam_corners(panel) )
    v_lbc_lam_map(:,panel) =                       &
        calc_panel_lbc_lam_map( self%base_v_panel, &
                                lam_ref_corner_id, &
                                lam_cell_next )
  end do ! panel

  deallocate(lam_cell_next)

  ! 2.0 Construct the LBC-LAM cell map from the base panels
  !----------
  allocate(cell_map(ntarget_per_source_x, ntarget_per_source_y, self%n_faces))

  lbc_cell_id = 0
  do panel=1, 2
    do cell_id=1, self%base_h_panel%n_cells
      lbc_cell_id = lbc_cell_id + 1
      cell_map(1,1,lbc_cell_id) = h_lbc_lam_map(cell_id,panel)
    end do
    do cell_id=1, self%base_v_panel%n_cells
      lbc_cell_id = lbc_cell_id + 1
      cell_map(1,1,lbc_cell_id) = v_lbc_lam_map(cell_id,panel)
    end do
  end do


  ! 3.0 Add the map to the LBC mesh's global mesh-map collection
  !----------
  source_id = 1
  target_id = 2
  call self%global_mesh_maps%add_global_mesh_map( source_id, &
                                                  target_id, &
                                                  cell_map )
  deallocate(cell_map)

  return
end subroutine calc_lbc_lam_map

!-------------------------------------------------------------------------------
!> @brief   Extracts the LBC mesh node coordinates. (PRIVATE SUBROUTINE)
!> @details Node coordinates are extracted from the parent LAM and assigned to
!>          the corresponding LBC mesh entitites.
!> @param[in] lam_node_coords  Node coordinates of parent LAM
!> @param[in] lam_cell_coords  Face coordinates of parent LAM
!> @param[in] lam_face_node    Nodes on cace connectivity of parent LAM
!-------------------------------------------------------------------------------
subroutine extract_coords(self, lam_node_coords, lam_cell_coords, lam_face_node )

  implicit none

  type(gen_lbc_type), intent(inout) :: self

  real(r_def),    intent(in) :: lam_cell_coords(:,:)
  real(r_def),    intent(in) :: lam_node_coords(:,:)
  integer(i_def), intent(in) :: lam_face_node(:,:)

  integer(i_def) :: cell_id
  integer(i_def) :: i

  integer(i_def) :: node_order(4), node_index

  integer(i_def) :: n_cells, n_nodes
  integer(i_def) :: lam_cell_id, lbc_node_id, lam_node_id
  integer(i_def) :: source_id=1
  integer(i_def) :: target_id=2
  integer(i_def), allocatable :: cell_map(:,:,:)

  type(global_mesh_map_type), pointer :: lbc_lam_map => null()

  n_cells = self%n_faces
  n_nodes = self%n_nodes

  node_order = [SW,SE,NE,NW]

  if (allocated(self%node_coords)) deallocate(self%node_coords)
  if (allocated(self%cell_coords)) deallocate(self%cell_coords)

  allocate(self%node_coords(2,n_nodes))
  allocate(self%cell_coords(2,n_cells))

  allocate(cell_map(1,1,1))

  self%node_coords = imdi
  self%cell_coords = imdi

  lbc_lam_map => self%global_mesh_maps%get_global_mesh_map( source_id, &
                                                            target_id )

  do cell_id=1, n_cells
    call lbc_lam_map%get_cell_map( [cell_id], cell_map)

    lam_cell_id = cell_map(1, 1, 1)
    self%cell_coords(:,cell_id) = lam_cell_coords(:, lam_cell_id)

    do i=1, self%nodes_per_face

      node_index = node_order(i)

      lbc_node_id = self%nodes_on_cell(node_index, cell_id)
      lam_node_id = lam_face_node( node_index, lam_cell_id )

      if (any(self%node_coords(:,lbc_node_id) == imdi)) then
        self%node_coords(:,lbc_node_id) = lam_node_coords(:, lam_node_id)
      end if

    end do
  end do

  return
end subroutine extract_coords


!=============================================================================
! D2: Methods for LBC mesh base panel object (base_lbc_panel_type)
!=============================================================================
! These are private subroutine which perform functions on a base panel.
! Depending on the base id th entity numbering is offset.
!=============================================================================

!-----------------------------------------------------------------------------
!> @brief  Generate cell connectivity for a base panel, each panel is treated
!>         as having no periodicity. (PRIVATE SUBROUTINE)
!> @param[in] panel      A base panel type from the LBC mesh
!> @param[in] base_id    The number from which to begin numbering cells
!> @return    cell_next  Cell next array for the panel
!-----------------------------------------------------------------------------
function calc_panel_adjacency( panel, base_id ) &
                       result( cell_next )

  implicit none

  type(base_lbc_panel_type), intent(in) :: panel
  integer(i_def),            intent(in) :: base_id

  integer(i_def), allocatable :: cell_next(:,:)

  integer(i_def) :: base_cell
  integer(i_def) :: cell, i

  if (allocated(cell_next)) deallocate(cell_next)
  allocate( cell_next( panel%edges_per_face, &
                       panel%n_cells ) )

  base_cell = base_id

  ! 1.0 Assign panel cell numbering
  !----------
  do cell=1, panel%n_cells
    cell_next(W, cell) = base_cell + cell - 1
    cell_next(S, cell) = base_cell + cell + panel%length
    cell_next(E, cell) = base_cell + cell + 1
    cell_next(N, cell) = base_cell + cell - panel%length
  end do

  ! 2.0 Set all cells on panel borders to have no neighbours
  !     external to the LBC mesh domain
  !----------
  do i=1, panel%length
    cell = panel%north_cells(i)
    cell_next(N, cell) = VOID_ID

    cell = panel%south_cells(i)
    cell_next(S, cell) = VOID_ID
  end do

  do i=1, panel%depth
    cell = panel%east_cells(i)
    cell_next(E, cell) = VOID_ID

    cell = panel%west_cells(i)
    cell_next(W, cell) = VOID_ID
  end do

  return
end function calc_panel_adjacency


!-------------------------------------------------------------------------------
!> @brief Generate cell-node connectivity for a base panel. (PRIVATE SUBROUTINE)
!> @param[in]    panel           A base panel type from the LBC mesh
!> @param[inout] base_id         The number from which to begin numbering nodes
!> @return       nodes_on_cells  Cell-node connectivty array for the panel
!-------------------------------------------------------------------------------
function calc_panel_node_connectivity( panel,    &
                                       base_id ) &
                               result( nodes_on_cells )

  implicit none

  type(base_lbc_panel_type), intent(in)    :: panel
  integer(i_def),            intent(inout) :: base_id

  integer(i_def), allocatable :: nodes_on_cells(:,:)

  integer(i_def) :: cell, start_id
  integer(i_def) :: i, j

  if (allocated(nodes_on_cells)) deallocate(nodes_on_cells)
  allocate( nodes_on_cells( panel%nodes_per_face, &
                            panel%n_cells ) )

  nodes_on_cells = imdi

  ! 1.0 Number the panel nodes starting
  !     from the provided base id.
  !----------
  ! 1.1 Number NorthWest node on cells
  do i=1, panel%depth
    do j=1, panel%length
      cell = (i-1)*panel%length + j
      base_id = base_id + 1
      nodes_on_cells(NW, cell) = base_id
    end do
  end do

  ! 1.2 Number SouthWest node on cells at the
  !     bottom of panel
  start_id = panel%n_cells - panel%length + 1
  do i=start_id, panel%n_cells
    base_id = base_id + 1
    nodes_on_cells(SW, i) = base_id
  end do

  ! 2.0  All unique edges numbered, now complete
  !      the edge connectivities for each cell.
  !----------
  ! 2.1 SouthWest node for all cells except the bottom row
  do i=1, panel%n_cells - panel%length
    nodes_on_cells(SW, i) = &
        nodes_on_cells(NW, panel%cell_next(S,i))
  end do

  ! 2.2 Eastern Nodes for all except eastern cells in panel
  do i=1, panel%n_cells
    if ( any(panel%east_cells == i)) cycle
    nodes_on_cells(NE, i) = &
        nodes_on_cells(NW, panel%cell_next(E,i))
    nodes_on_cells(SE, i) = &
        nodes_on_cells(SW, panel%cell_next(E,i))
  end do

  return
end function calc_panel_node_connectivity

!-------------------------------------------------------------------------------
!> @brief  Generate cell-edge connectivity for a base panel (PRIVATE SUBROUTINE)
!> @param[in]    panel          A base panel type from the LBC mesh
!> @param[inout] base_id        The number from which to begin numbering edges
!> @return       edge_on_cells  Cell-edge connectivty array for the panel
!-------------------------------------------------------------------------------
function calc_panel_edge_connectivity( panel,    &
                                       base_id ) &
                               result( edges_on_cells )

  implicit none

  type(base_lbc_panel_type), intent(in)    :: panel
  integer(i_def),            intent(inout) :: base_id

  integer(i_def), allocatable :: edges_on_cells(:,:)

  integer(i_def) :: cell, start_id
  integer(i_def) :: i, j

  if (allocated(edges_on_cells)) deallocate(edges_on_cells)
  allocate( edges_on_cells( panel%edges_per_face, &
                            panel%n_cells ) )

  edges_on_cells = imdi

  ! 1.0 Number the panel edges starting
  !     from the provided base id.
  !----------
  do i=1, panel%depth

    ! 1.1 Number North edge on cells
    do j=1, panel%length
      cell = (i-1)*panel%length + j
      base_id = base_id + 1
      edges_on_cells( N, cell) = base_id
    end do

    ! 1.2 Number West edge on cells
    do j=1, panel%length
      cell = (i-1)*panel%length + j
      base_id = base_id + 1
      edges_on_cells( W, cell) = base_id
    end do

  end do

  ! 1.3 Number South edge on bottom of panel
  start_id = panel%n_cells - panel%length + 1
  do i=start_id, panel%n_cells
    base_id = base_id + 1
    edges_on_cells( S, i) = base_id
  end do

  ! 2.0 All unique edges numbered, now complete
  !     the edge connectivities for each cell.
  !----------
  ! 2.1 South edge for all cells except the bottom row
  do i=1, panel%n_cells - panel%length
    edges_on_cells( S, i) = &
        edges_on_cells( N, panel%cell_next(S,i))
  end do

  ! 2.2 East edge for all cells in panel
  do i=1, panel%n_cells
    if ( any(panel%east_cells == i)) cycle
    edges_on_cells( E, i ) = &
        edges_on_cells( W, panel%cell_next(E,i))
  end do

  return
end function calc_panel_edge_connectivity


!-------------------------------------------------------------------------------
!> @brief   Generate an integrid map for a LBC mesh base panel (PRIVATE SUBROUTINE)
!> @details A base panel is provided along with its ref corner, this locates it
!>          with respect to the parent LAM. This along with the parents LAM's cell
!>          connectivity is used to generate a LBC-LAM map for the panel.
!>
!> @return  lbc_lam_map An array which provides a one-to-one map of LBC mesh ids
!>                      in the given panel to the corresponding parent LAM mesh
!>                      ids
!-------------------------------------------------------------------------------
function calc_panel_lbc_lam_map( panel,             &
                                 lam_ref_corner_id, &
                                 lam_cell_next )    &
                         result( lbc_lam_map )

  implicit none

  type(base_lbc_panel_type), intent(in) :: panel
  integer(i_def),            intent(in) :: lam_ref_corner_id
  integer(i_def),            intent(in) :: lam_cell_next(:,:)

  integer(i_def), allocatable :: lbc_lam_map(:)

  integer(i_def) :: lbc_cell_id
  integer(i_def) :: start_id
  integer(i_def) :: lam_cell_id
  integer(i_def) :: cell_id

  if (allocated(lbc_lam_map)) deallocate(lbc_lam_map)
  allocate( lbc_lam_map(panel%n_cells) )

  ! Map is generated considering the base panel
  ! in a horizontal orientation.

  lbc_cell_id = 1
  lbc_lam_map(1) = lam_ref_corner_id

  ! 1.0 Top row of panel
  !----------
  do cell_id=2, panel%length
    lbc_cell_id = lbc_cell_id + 1
    lam_cell_id = lbc_lam_map(cell_id-1)
    lbc_lam_map(lbc_cell_id) = lam_cell_next(E, lam_cell_id)
  end do

  ! 2.0 remaining rows if required
  !----------
  if (panel%depth > 1) then
    start_id = panel%length + 1
    do cell_id=start_id, panel%n_cells
      lbc_cell_id = lbc_cell_id + 1
      lam_cell_id = lbc_lam_map( panel%cell_next(N,lbc_cell_id) )
      lbc_lam_map(lbc_cell_id) = lam_cell_next(S, lam_cell_id)
    end do
  end if

end function calc_panel_lbc_lam_map


!-----------------------------------------------------------------------------
!> @brief Allows for manual recovery of base_lbc_panel
!>        object memory. (PRIVATE SUBROUTINE)
!-----------------------------------------------------------------------------
subroutine base_lbc_panel_clear(self)

  implicit none

  class(base_lbc_panel_type), intent(inout) :: self

  if ( allocated(self%cell_next)   ) deallocate( self%cell_next   )
  if ( allocated(self%north_cells) ) deallocate( self%north_cells )
  if ( allocated(self%south_cells) ) deallocate( self%south_cells )
  if ( allocated(self%east_cells)  ) deallocate( self%east_cells  )
  if ( allocated(self%west_cells)  ) deallocate( self%west_cells  )
  return

end subroutine base_lbc_panel_clear



!-----------------------------------------------------------------------------
end module gen_lbc_mod
