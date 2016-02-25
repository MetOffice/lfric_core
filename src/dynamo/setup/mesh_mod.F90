!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
!
!> Module for mesh_type which defines local 3D mesh object.
!>
!> This module provides details for a mesh_type which is generated
!> using a global mesh object and a partition object, along some with
!> additional inputs.
!>
!> It is also contains a static mesh object for unit testing. This
!> is returned if a mesh_object is instatiated with a single integer
!> argument.

module mesh_mod

  use global_mesh_mod,   only: global_mesh_type
  use partition_mod,     only: partition_type
  use constants_mod,     only: i_def, r_def, l_def, IMDI, PI
  use configuration_mod, only: l_spherical, VGRID_UNIFORM, VGRID_QUADRATIC &
                             , VGRID_GEOMETRIC, VGRID_DCMIP
  use log_mod,           only: log_event, log_scratch_space, LOG_LEVEL_DEBUG   &
                             , LOG_LEVEL_ERROR

  implicit none

  private

  ! Connectivities in global ids, 2d-base level only.
  ! These are used in the creation of connectivity arrays from the
  ! global mesh, but are only temporary as the same connectivity
  ! information is present in the 3D-mesh once it has been extruded
  !
  ! Arrays where the connected entity ids are local ids

  !> @}
  !> @name Private variables used in object construction
  !>
  integer(i_def), allocatable :: vert_on_cell_2d (:,:) 
  !< Vertices connected to local 2d cell.
  !>
  integer(i_def), allocatable :: edge_on_cell_2d (:,:)
  !< Edges connected to local 2d cell.
  !>
  integer(i_def), allocatable :: cell_next_2d (:,:)
  !< Local 2d cell connectivity.
  !>
  real(r_def), allocatable :: vertex_coords_2d(:,:)
  !< Surface Coordinates in, [long, lat, radius]
  !< Units in Radians/metres
  !>
  integer(i_def) :: ncells_2d
  integer(i_def) :: nverts_2d
  integer(i_def) :: nedges_2d
  integer(i_def) :: nlayers
  integer(i_def) :: ncells
  integer(i_def) :: nverts_3d
  integer(i_def) :: nedges_3d
  integer(i_def) :: nfaces_3d
  !>
  !> @}

  !============================================================================
  ! Declare type definitions in this module
  !============================================================================

  ! types definitions for computing/storing the domain size
  type, private :: coordinate
    real(r_def) :: x,y,z
  end type coordinate

  type, public :: domain_limits
    type(coordinate) :: minimum, maximum
  end type domain_limits

  type, public :: mesh_type

    private

    !> The partition object that describes this local
    !> mesh's partition of the global mesh
    type(partition_type) :: partition

    !> The domain limits (x,y,z) for Cartesian domains
    !>                   (long, lat, radius) for spherical
    type (domain_limits) :: domain_size

    !> Mesh id number
    integer(i_def) :: mesh_id

    !> Number of 3d-cell layers in mesh object
    integer(i_def) :: nlayers

    !> Top of atmosphere above surface
    real(r_def) :: domain_top

    !> Non-dimensional vertical coordinate eta[0,1], eta(0:nlayers)
    real(r_def), allocatable :: eta(:)

    !> Depth of 3d-cell layer [m], dz(nlayers)
    real(r_def), allocatable :: dz(:)

    !> Vertex Coordinates
    !> The x-, y- and z-coordinates of vertices on the mesh [m]
    real(r_def), allocatable :: vertex_coords(:,:)

    !==========================================================================
    ! Mesh properties:
    !==========================================================================
    ! Local partition base level
    integer(i_def) :: nverts_2d            !< Number of verts in partition
    integer(i_def) :: nedges_2d            !< Number of edges in partition
    integer(i_def) :: ncells_2d            !< Number of cells in partition
    integer(i_def) :: ncells_2d_with_ghost !< Number of cells in partition
                                           !< @b including ghost cells

    ! Local partition 3d-mesh
    integer(i_def) :: nverts               !< Total number of verts in mesh
    integer(i_def) :: nedges               !< Total number of edges in mesh
    integer(i_def) :: nfaces               !< Total number of faces in mesh
    integer(i_def) :: ncells               !< Total number of cells in mesh
                                           !< @b excluding ghost cells
    integer(i_def) :: ncells_with_ghost    !< Total number of cells in mesh
                                           !< @b including ghost cells

    ! 3D-element properties
    integer(i_def) :: nverts_per_cell      !< Number of verts on 3d-cell
    integer(i_def) :: nedges_per_cell      !< Number of edges on 3d-cell
    integer(i_def) :: nfaces_per_cell      !< Number of faces on 3d-cell


    !==========================================================================
    ! Local/Global id Maps
    !==========================================================================
    ! Map is only for the base level, i.e. the surface level.
    ! The index of the array is taken to be the entities local id.
    integer(i_def), allocatable :: cell_lid_gid_map(:)

    !==========================================================================
    ! Connectivities
    !==========================================================================
    ! All vertex, edge, face and cell id numbers for connectivitivies are
    ! LOCAL (lid) unless explicitly stated in the variable name,
    ! i.e. gid (GLOBAL id numbers)

    ! All connectivity arrays are in the form (ids of entities, cell local id)
    ! i.e. cell_next(5,6) would hold the local id of the cell adjacent to
    !      face-5 of 3d-cell with local id = 6. Face-5 of a 3d-cell is the
    !      bottom so would be the 3d-cell directly below.

    ! For 3d-mesh
    !> Cell ids of adjacent cells
    integer(i_def), allocatable :: cell_next    (:,:)

    !> Vertex ids on cell
    integer(i_def), allocatable :: vert_on_cell (:,:)

    !> Face ids on cell
    integer(i_def), allocatable :: face_on_cell (:,:)

    !> Edge ids on cell
    integer(i_def), allocatable :: edge_on_cell (:,:)

    !> The cell that "owns" the vertex entities around each cell
    integer(i_def), allocatable :: vert_cell_owner(:,:)

    !> The cell that "owns" the edge entities around each cell
    integer(i_def), allocatable :: edge_cell_owner(:,:)

    !> The rank of the partition that "owns" the
    !> vertex entities around each cell
    integer(i_def), allocatable :: vertex_ownership(:,:)

    !> The rank of the partition that "owns" the
    !> edge entities around each cell
    integer(i_def), allocatable :: edge_ownership(:,:)

    !==========================================================================
    ! Colouring storage: these form the arguments to set_colours().
    !==========================================================================
    integer(i_def),              private :: ncolours
    integer(i_def), allocatable, private :: ncells_per_colour(:)
    integer(i_def), allocatable, private :: cells_in_colour(:,:)

  contains

    procedure, public :: get_id
    procedure, public :: get_nlayers
    procedure, public :: get_ncells_2d
    procedure, public :: get_ncells_2d_with_ghost
    procedure, public :: get_nedges_2d
    procedure, public :: get_nverts_2d
    procedure, public :: get_ncells
    procedure, public :: get_nverts
    procedure, public :: get_nedges
    procedure, public :: get_nfaces
    procedure, public :: get_vert_coords
    procedure, public :: get_cell_coords
    procedure, public :: get_column_coords
    procedure, public :: get_nverts_per_cell
    procedure, public :: get_nedges_per_cell
    procedure, public :: get_nfaces_per_cell
    procedure, public :: get_cell_gid
    procedure, public :: get_cell_lid
    procedure, public :: get_cell_next
    procedure, public :: get_face_on_cell
    procedure, public :: get_edge_on_cell
    procedure, public :: get_vert_on_cell
    procedure, public :: get_domain_size
    procedure, public :: get_domain_top
    procedure, public :: get_dz
    procedure, public :: get_eta
    procedure, public :: get_vertex_cell_owner
    procedure, public :: get_edge_cell_owner
    procedure, public :: is_vertex_owned
    procedure, public :: is_edge_owned
    procedure, public :: is_cell_owned
    procedure, public :: get_num_cells_core
    procedure, public :: get_num_cells_owned
    procedure, public :: get_halo_depth
    procedure, public :: get_num_cells_halo
    procedure, public :: get_num_cells_ghost
    procedure, public :: get_gid_from_lid

    ! Colouring now accessed through function_space_type
    procedure, public :: set_colours
    procedure, public :: get_ncolours
    procedure, public :: get_colours
    procedure, public :: is_coloured

    ! Overloaded assigment operator
    procedure         :: mesh_type_assign

    ! Destructor frees colouring storage
    final :: mesh_destructor

    ! Override default assignment for mesh_type pairs.
    generic           :: assignment(=) => mesh_type_assign

  end type mesh_type

  interface mesh_type
    module procedure mesh_constructor
    module procedure mesh_constructor_unit_test_data
  end interface

  !============================================================================
  ! Options for horizontal pFUnit test meshes
  !============================================================================
  !
  !> @}
  !> @name Horizontal Grid Types for pFunit tests
  integer(i_def), parameter, public :: PLANE             = 1
  integer(i_def), parameter, public :: PLANE_BI_PERIODIC = 2
  !> @}

  !> Counter variable to keep track next mesh id number
  integer(i_def) :: mesh_id_counter = 0

contains


  !============================================================================
  !> @brief Stucture-Constructor
  !> @param [in] partition     Partition object to base 3D-Mesh on
  !> @param [in] global_mesh   Global mesh object on which the partition is
  !>                           applied
  !> @param [in] nlayers_in    Number of 3D-cell layers in the 3D-Mesh object
  !> @param [in] domain_top    Top of atmosphere above surface
  !> @param [in] vgrid_option  Choice of vertical grid
  !> @return                   3D-Mesh object based on the list of partitioned
  !>                           cells on the given global mesh
  !============================================================================
  function mesh_constructor ( partition,     &
                              global_mesh,   &
                              nlayers_in,    &
                              domain_top,    &
                              vgrid_option ) &
                            result( self )

    ! User structure constructor function for 3d-mesh object
    ! for given 2d-partition of global skin.

    use reference_element_mod, only: nfaces, nverts

    implicit none

    type (partition_type),   intent(in) :: partition
    type (global_mesh_type), intent(in) :: global_mesh

    integer(i_def), intent(in) :: nlayers_in
    integer(i_def), intent(in) :: vgrid_option
    real(r_def),    intent(in) :: domain_top

    type(mesh_type) :: self

    integer(i_def) :: i, j        ! loop counters

    integer(i_def) :: nverts_per_2d_cell
    integer(i_def) :: nedges_per_2d_cell


    ! Arrays used in entity ownership calculation - see their names
    ! for a descriptioin of what they actually contain
    integer(i_def), allocatable :: verts( : )
    integer(i_def), allocatable :: edges( : )

    nverts_per_2d_cell = global_mesh%get_nverts_per_cell()
    nedges_per_2d_cell = global_mesh%get_nedges_per_cell()

    mesh_id_counter = mesh_id_counter+1

    self%mesh_id    = mesh_id_counter
    self%partition  = partition
    self%ncells_2d  = partition%get_num_cells_in_layer()
    self%nlayers    = nlayers_in
    self%ncells     = self%ncells_2d * self%nlayers
    self%domain_top = domain_top
    self%ncolours   = -1     ! Initialise ncolours to error status
    self%ncells_2d_with_ghost = self%ncells_2d                                 &
                              + partition%get_num_cells_ghost()
    self%ncells_with_ghost    = self%ncells_2d_with_ghost * self%nlayers

    ncells_2d       = self%ncells_2d_with_ghost
    ncells          = self%ncells_with_ghost
    nlayers         = self%nlayers

    allocate( self%eta ( 0:nlayers ) )
    allocate( self%dz  ( nlayers   ) )

    ! Calculate vertical coordinates eta[0,1] and dz in a separate subroutine
    call set_vertical_coordinate(self, vgrid_option)

    ! Calculate next-to cells and vertices on cells
    allocate ( self % cell_next    ( nfaces, ncells ) )
    allocate ( self % vert_on_cell ( nverts, ncells ) )

    ! Need connectivities of 2d mesh cells that make up this
    ! partition from the global mesh connectivities, in lid.
    call get_partition_stats (self, partition, global_mesh)

    self%nverts = nverts_3d
    self%nedges = nedges_3d
    self%nfaces = nfaces_3d

    call mesh_extruder     (self)
    call mesh_connectivity (self)
    call set_domain_size   (self)


    ! Assign ownership of cell vertices and cell edges
    allocate( &
      self%vertex_ownership (nverts_per_2d_cell, ncells_2d) )
    allocate( &
      self%edge_ownership   (nedges_per_2d_cell, ncells_2d) )
    allocate( &
      self%vert_cell_owner  (nverts_per_2d_cell, ncells_2d) )
    allocate( &
      self%edge_cell_owner  (nedges_per_2d_cell, ncells_2d) )

    allocate( verts (nverts_per_2d_cell) )
    allocate( edges (nedges_per_2d_cell) )

    do i=1, ncells_2d
      ! Vertex ownership
      call global_mesh%get_vert_on_cell(partition%get_gid_from_lid(i), verts)
      do j=1, nverts_per_2d_cell
        self%vert_cell_owner(j,i) = partition%get_lid_from_gid( &
                                     global_mesh%get_vert_cell_owner( verts(j) ) &
                                                               )
        if (self%vert_cell_owner(j,i) > 0) then
          self%vertex_ownership(j,i) = partition%get_cell_owner( &
                                                     self%vert_cell_owner(j,i) &
                                                               )
        else
          self%vertex_ownership(j,i) = partition%get_total_ranks() + 1
        end if
      end do

      ! Edge ownership
      call global_mesh%get_edge_on_cell(partition%get_gid_from_lid(i), edges)
      do j=1, nedges_per_2d_cell
        self%edge_cell_owner(j,i) = partition%get_lid_from_gid( &
                                     global_mesh%get_edge_cell_owner( edges(j) ) &
                                                              )
        if (self%edge_cell_owner(j,i) > 0) then
          self%edge_ownership(j,i) = partition%get_cell_owner( &
                                                     self%edge_cell_owner(j,i) &
                                                             )
        else
          self%edge_ownership(j,i) = partition%get_total_ranks() + 1
        end if
      end do
    end do

    deallocate( verts )
    deallocate( edges )

  end function mesh_constructor


  !============================================================================
  ! Mesh Type Methods
  !============================================================================
  !> @details This subrotuine returns 3-element array of vertex coords
  !>          in cartesian coords [x,y,z] with units in [m].
  !> @param[in]   vert_lid      The local id of the requested vertex
  !> @param[out]  vertex_coords A three-element array containing the
  !>                            cartesian coordinates of a single vertex
  !>                            in the mesh object
  !============================================================================
  subroutine get_vert_coords(self, vert_lid, vertex_coords)

    ! Returns 3-element array of vertex coords in
    ! cartesian coords [x,y,z] with units in [m].

    implicit none
    class(mesh_type), intent(in)  :: self
    integer(i_def),   intent(in)  :: vert_lid
    real(r_def),      intent(out) :: vertex_coords(:)

    vertex_coords(:) = self%vertex_coords(:,vert_lid)

  end subroutine get_vert_coords


  !> @details This subroutine returns 3-element array of vertex coords for
  !>          each vertex on the request local cell id. Coords are in
  !>          cartesian coords [x,y,z] in [m] and in same order as the
  !>          vertex on cell connectivity array
  !> @param[in]  cell_lid     The local id of the requested cell
  !> @param[out  cell_coords  A 2-dimensional array which contains
  !>                          the cartesian coordinates of each vertex
  !>                          on the requested cell_lid. The returned
  !>                          array will have dimensions of
  !>                          [3, nVertices on cell]
  !============================================================================
  subroutine get_cell_coords(self, cell_lid, cell_coords)

    ! Returns 3-element array of vertex coords for each
    ! vertex on the request local cell id. Coords are
    ! in cartesian coords [x,y,z] in [m] and in same order
    ! as the vertex on cell connectivity array

    implicit none
    class(mesh_type), intent(in)  :: self
    integer(i_def),   intent(in)  :: cell_lid
    real(r_def),      intent(out) :: cell_coords(:,:)

    integer(i_def) :: ivert, vert_lid

    do ivert=1, self%nverts_per_cell
      vert_lid = self%vert_on_cell(ivert, cell_lid)
      call self%get_vert_coords(vert_lid, cell_coords(:,ivert))
    end do

  end subroutine get_cell_coords


  !> @details This subroutine returns 3-element array of vertex coords
  !>          for a vertical cell column in the local mesh. Any local
  !>          cell id within the column can be provided.
  !> @param[in]  cell_lid      The local id of the requested cell
  !> @param[out] column_coords A 3-dimensional array which contains the
  !>                           cartesian coordinates of each vertex on the
  !>                           column of 3D-cells which include the requested
  !>                           cell_lid. The returned array will have
  !>                           dimensions of [3, nVertices on cell, nlayers]
  !============================================================================
  subroutine get_column_coords(self, cell_lid, column_coords)

    ! Returns 3-element array of vertex coords for a vertical
    ! cell column in the local mesh. Any local cell id within
    ! the column can be provided.

    implicit none
    class(mesh_type), intent(in)  :: self
    integer(i_def),   intent(in)  :: cell_lid
    real(r_def),      intent(out) :: column_coords(:,:,:)

    integer(i_def) :: base_id, k, icell

    if (cell_lid > self%ncells_2d) then
      base_id = modulo(cell_lid,self%ncells_2d)
      if (base_id == 0) base_id = self%ncells_2d
    else
      base_id = cell_lid
    end if

    do k=1, self%nlayers
      icell = base_id + (k-1)*self%ncells_2d
      call self%get_cell_coords(icell,column_coords(:,:,k))
    end do

  end subroutine get_column_coords


  !> @details This function returns the id number of the mesh, this number is
  !>          assigned when the object is first instatiated.
  !> @return  Mesh object id number
  !============================================================================
  function get_id(self) result (mesh_id)

    ! Function returns the id number of the mesh, this
    ! number is assigned when the object is first instatiated.

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def)               :: mesh_id

    mesh_id = self%mesh_id

  end function get_id


  !> @details This function returns the number of 3d-cell layers in the mesh
  !>          object
  !> @return  Number of 3d-cell vertical layers in the mesh object
  !============================================================================
  function get_nlayers(self) result (nlayers)

    ! This function returns the number of 3d-cell layers
    ! in the mesh object

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def)               :: nlayers

    nlayers = self%nlayers

  end function get_nlayers


  !> @details This function returns the number of
  !>          vertices per 3d-cell on this mesh
  !> @return  Number of vertices per 3d-cell on mesh
  !============================================================================
  function get_nverts_per_cell(self) result (nverts_per_cell)

    ! Returns number of vertices per 3d-cell on this mesh

    implicit none
    class (mesh_type), intent(in) :: self
    integer(i_def)                :: nverts_per_cell

    nverts_per_cell = self%nverts_per_cell

  end function get_nverts_per_cell


  !> @details This function returns the number of
  !>          edges per 3d-cell on this mesh
  !> @return  Number of edges per 3d-cell on mesh
  !============================================================================
  function get_nedges_per_cell(self) result (nedges_per_cell)

    ! Returns number of edges per 3d-cell on this mesh

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def)               :: nedges_per_cell

    nedges_per_cell = self%nedges_per_cell

  end function get_nedges_per_cell


  !> @details This function returns the number
  !>          of faces per 3d-cell on this et-mesh
  !> @return  Number of faces per 3d-cell on mesh
  !============================================================================
  function get_nfaces_per_cell(self) result (nfaces_per_cell)

    ! Returns number of faces per 3d-cell on this mesh

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def)               :: nfaces_per_cell

    nfaces_per_cell = self%nfaces_per_cell

  end function get_nfaces_per_cell


  !> @details Returns local cell id of adjacent cell on the specified face
  !>          (iface) of the given cell (icell)
  !> @param[in] iface         The index (on the known cell with local id, 
  !>                          @c cell_lid ) of the face common to both cells
  !> @param[in] cell_lid      The local id of the known cell
  !> @return                  The local id of cell adjacent to the
  !>                          known cell  ( @c cell_lid ) with the
  !>                          common face ( @c iface )
  !============================================================================
  function get_cell_next(self, iface, cell_lid) result (cell_next_lid)

    ! Returns local cell id of adjacent cell on the
    ! specified face (iface) of the given cell (cell_lid)

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def),   intent(in) :: iface      ! Index of face required
    integer(i_def),   intent(in) :: cell_lid   ! Local cell id
    integer(i_def)               :: cell_next_lid

    cell_next_lid = self%cell_next(iface, cell_lid)

  end function get_cell_next


  !> @details This function returns the local face id on local cell
  !> @param [in]  iface     The index face of interest
  !> @param [in]  icell     The local id of the cell which the
  !>                        face is a member of
  !> @return                The local id of face on index iface
  !>                        of cell with local id icell
  !============================================================================
  function get_face_on_cell(self, iface, icell) result (face_lid)

    ! Returns local face id on local cell

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def),   intent(in) :: iface   ! Index of face required
    integer(i_def),   intent(in) :: icell   ! Local cell id
    integer(i_def)               :: face_lid

    face_lid = self % face_on_cell(iface, icell)

  end function get_face_on_cell


  !> @details This function returns the local edge id on the local cell
  !> @param [in] iedge    The index of the edge whose local id is requested
  !> @param [in] icell    The local id of the cell on which the edge is located
  !> @return              Local id of edge on index iedge of cell with
  !>                      local id icell
  !============================================================================
  function get_edge_on_cell(self, iedge, icell) result (edge_lid)

    ! Returns local edge id on local cell

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def),   intent(in) :: iedge   ! Index of edge required
    integer(i_def),   intent(in) :: icell   ! Local cell id
    integer(i_def)               :: edge_lid

    edge_lid = self % edge_on_cell(iedge, icell)

  end function get_edge_on_cell


  !> @details This function returns the local vertex id on the local cell
  !> @param [in] ivert    The index of the vertex whose local id is requested
  !> @param [in] icell    The local id of the cell on which the vertex
  !>                      is located
  !> @return              Local id of vertex on index ivert of cell with
  !>                      local id icell
  !============================================================================
  function get_vert_on_cell(self, ivert, icell) result (vert_lid)

    ! Returns local vertex id on local cell

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def),   intent(in) :: ivert   ! Index of vertex required
    integer(i_def),   intent(in) :: icell   ! Local cell id
    integer(i_def)               :: vert_lid

    vert_lid = self % vert_on_cell(ivert, icell)

  end function get_vert_on_cell


  !> @details This function returns the number of cells in the horizontal
  !>          (i.e. the number of cells in the partition), @b excluding any
  !>          ghost cells.
  !> @return  Number of cells on a horizontal layer.
  !============================================================================
  function get_ncells_2d(self) result (ncells_2d)

    ! This function returns the number of cells in the horizontal
    ! (i.e. the number of cells in the partition), EXCLUDING any
    ! ghost cells.

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def) :: ncells_2d

    ncells_2d = self%ncells_2d

  end function get_ncells_2d


  !> @details This function returns the number of cells in the horizontal
  !>          (i.e. the number of cells in the partition), @b including any
  !>          ghost cells.
  !> @return  Number of cells on a horizontal layer
  !>         @b plus ghost cells.
  !============================================================================
  function get_ncells_2d_with_ghost(self) result (ncells_2d_with_ghost)

    ! This function returns the number of cells in the horizontal
    ! (i.e. the number of cells in the partition), INCLUDING any
    ! ghost cells.

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def) :: ncells_2d_with_ghost

    ncells_2d_with_ghost = self%ncells_2d_with_ghost

  end function get_ncells_2d_with_ghost


  !> @details This function returns the number of edges on one horizontal level
  !>          of the mesh. This is the same as the number of edges on the local
  !>          partition.
  !> @return  Number of edges on one horizontal level in the
  !>                   mesh object i.e. one edge thick
  !============================================================================
  function get_nedges_2d(self) result (nedges_2d)

    ! Returns total number of horizontal edges on each 2d-level
    ! of this mesh. This is the same as the number of edges on the
    ! local partition.

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def)               :: nedges_2d

    nedges_2d = self%nedges_2d

  end function get_nedges_2d


  !> @details This function returns the number of vertices on one horizontal
  !>          level of the mesh. This is the same as the number of vertices
  !>          on the local partition.
  !> @return  Number of vertices on horizontal layer in the
  !>                   mesh object i.e. one vertex in thick
  !============================================================================
  function get_nverts_2d(self) result (nverts_2d)

    ! Returns the number of vertices on one horizontal level
    ! of the mesh. This is the same as the number of vertices
    ! on the local partition.

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def)               :: nverts_2d

    nverts_2d = self%nverts_2d

  end function get_nverts_2d


  !> @details This function returns the number of 3d-cells in the local 3d-mesh.
  !> @return  Total number of 3d-cells in the mesh object
  !============================================================================
  function get_ncells(self) result (ncells)

    ! Returns total number of 3d-cells in this mesh

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def)               :: ncells

    ncells = self%ncells

  end function get_ncells


  !> @details This function returns the number of vertices in the local 3d-mesh.
  !> @return  Total number of vertices in the mesh object
  !============================================================================
  function get_nverts(self) result (nverts)

    ! Returns total number of vertices in this 3d-mesh

    class(mesh_type), intent(in) :: self
    integer(i_def) :: nverts

    nverts = self%nverts

  end function get_nverts


  !> @details This function returns the number of edges in the local 3d-mesh.
  !> @return  Total number of edges in the mesh object
  !============================================================================
  function get_nedges(self) result (nedges)

    ! Returns total number of edges in this 3d-mesh

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def)               :: nedges

    nedges = self%nedges

  end function get_nedges


  !> @details This function returns the number of faces in the local 3d-mesh.
  !> @return  Total number of faces in the mesh object
  !============================================================================
  function get_nfaces(self) result (nfaces)

    ! Returns total number of faces in this 3d-mesh

    class(mesh_type), intent(in) :: self
    integer(i_def)               :: nfaces

    nfaces = self%nfaces

  end function get_nfaces


  !> @details This function returns the global (across all processors) cell id
  !>          of a 2d-cell (i.e. on the local partition) with a local id given
  !>          by @c cell_lid .
  !> @param [in] cell_lid  The local id of the requested cell
  !> @return Global id of cell specified by @c cell_lid
  !============================================================================
  function get_cell_gid(self, cell_lid) result (cell_gid)

    ! Returns global cell id of a local 2d-cell on the partition
    ! that this mesh is based upon

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def),   intent(in) :: cell_lid
    integer(i_def)               :: cell_gid

    cell_gid = self%cell_lid_gid_map(cell_lid)

  end function get_cell_gid


  !> @details This function returns the local (to partition) cell id
  !>          of a 2d-cell with the global (across all processors) cell
  !>          id given by @c cell_gid .
  !> @param [in] cell_gid  Global id of the requested 2d-cell
  !> @return The local id of cell with global id @c cell_gid
  !============================================================================
  function get_cell_lid(self, cell_gid) result (cell_lid)

    ! Returns local cell id on this mesh of a given global
    ! cell id. If the global cell does not exist on this mesh,
    ! IMDI is returned

    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def),   intent(in) :: cell_gid
    integer(i_def)               :: cell_lid

    integer(i_def) :: i

    cell_lid = IMDI
    do i=1, self%ncells_2d
      if (self%cell_lid_gid_map(i) == cell_gid) then
        cell_lid = i
        exit
      end if
    end do

  end function get_cell_lid

  !> @details Returns the height above surface of the top of the mesh domain
  !> @return  Height of top of mesh
  !============================================================================
  function get_domain_top(self) result (domain_top)

    ! Returns top of model height above surface

    implicit none
    class (mesh_type), intent(in) :: self
    real  (r_def)                 :: domain_top

    domain_top = self%domain_top

  end function get_domain_top


  !> @details This functions returns an array of 3d-layer thicknesses in
  !>          metres
  !> @param[out]  dz  Vertical thickness of layers in [m], array of
  !>                  length @c nlayers . The number of layers in the mesh
  !>                  object can be obtained using the @c get_nlayers
  !>                  type-bound function.
  !============================================================================
  subroutine get_dz(self,dz)

    ! Returns array of vertical layer thicknesses (in metres) for this 3d-mesh

    implicit none
    class (mesh_type), intent(in) :: self
    real(r_def),      intent(out) :: dz(:)

    ! Get the thickness of layers
    dz(:) = self%dz(:)

  end subroutine get_dz


  !> @param[out] eta Array of dimensions (0:nlayers) of non-dimensional
  !>                 vertical coordinate normalised using the @c domain_top
  !============================================================================
  subroutine get_eta(self,eta)

    ! Returns non-dimensional vertical coordinate eta[0,1]

    implicit none
    class (mesh_type), intent(in) :: self
    real(r_def),      intent(out) :: eta(:)

    ! Get the thickness of layers
    eta(:) = self%eta(:)

  end subroutine get_eta


  !> @return  Domain size of mesh as a <domain_type>
  !============================================================================
  function get_domain_size(self) result (domain_size)

    ! Returns local mesh domain limits as a domain_limits type

    implicit none
    class(mesh_type), intent(in) :: self
    type(domain_limits)          :: domain_size

    domain_size = self%domain_size

  end function get_domain_size


  !> @param[in] vertex_index  The index of the vertex entity on the cell with
  !>                          local_id of cell_lid.
  !> @param[in] cell_lid      The local cell id on which the vertex entity is
  !>                          associated with.
  !> @return                  Gets the local id of the cell that "owns" a
  !>                          particular vertex entity
  !============================================================================
  function get_vertex_cell_owner( self, vertex_index, cell_lid )               &
    result (cell_owner)

    class(mesh_type), intent(in) :: self
    integer(i_def),   intent(in) :: vertex_index
    integer(i_def),   intent(in) :: cell_lid
    integer(i_def)               :: cell_owner

    cell_owner = self%vert_cell_owner( vertex_index, cell_lid )

  end function get_vertex_cell_owner


  !> @param[in] edge_index  The index of the edge entity on the cell with
  !>                        local_id of cell_lid.
  !> @param[in] cell_lid    The local cell id on which the edge entity is
  !>                        associated with.
  !> @return                Gets the local id of the cell that "owns" a
  !>                        particular edge entity
  !============================================================================
  function get_edge_cell_owner( self, edge_index, cell_lid )                   &
    result (cell_owner)

    class(mesh_type), intent(in) :: self
    integer(i_def),   intent(in) :: edge_index
    integer(i_def),   intent(in) :: cell_lid
    integer(i_def)               :: cell_owner

    cell_owner = self%edge_cell_owner( edge_index, cell_lid )

  end function get_edge_cell_owner


  !> @param[in] vertex_index  The index of the vertex entity on the cell with
  !>                          local_id of cell_lid.
  !> @param[in] cell_lid      The local cell id on which the edge entity is
  !>                          associated with.
  !> @return                  Whether the vertex is owned by the
  !>                          local parition
  !============================================================================
  function is_vertex_owned( self, vertex_index, cell_lid ) result (owned)

    class(mesh_type), intent(in) :: self
    integer(i_def),   intent(in) :: vertex_index
    integer(i_def),   intent(in) :: cell_lid
    logical(l_def)               :: owned

    owned = .false.
    if ( self%vertex_ownership( vertex_index, cell_lid ) ==                    &
         self%partition%get_local_rank() ) owned = .true.

  end function is_vertex_owned


  !> @param[in] edge_index  The index of the edge entity on the cell with
  !>                        local_id of cell_lid.
  !> @param[in] cell_lid    The local cell id which the edge entity is
  !>                        associated with.
  !> @return                Whether the edge is owned by the local parition
  !============================================================================
  function is_edge_owned( self, edge_index, cell_lid ) result (owned)

    class(mesh_type), intent(in) :: self
    integer(i_def),   intent(in) :: edge_index
    integer(i_def),   intent(in) :: cell_lid
    logical(l_def)               :: owned

    owned = .false.
    if ( self%edge_ownership( edge_index, cell_lid ) ==                        &
         self%partition%get_local_rank() ) owned = .true.

  end function is_edge_owned


  !> @param[in] cell_lid  The local cell id on the local partition.
  !> @return              Whether the cell is owned by the local parition
  !============================================================================
  function is_cell_owned( self, cell_lid ) result (owned)

    class(mesh_type), intent(in) :: self
    integer(i_def),   intent(in) :: cell_lid
    logical(l_def)               :: owned

    owned = .false.
    if ( self%partition%get_cell_owner(cell_lid) ==                            &
         self%partition%get_local_rank())    owned = .true.

  end function is_cell_owned


  !> @return The total number of core cells on the local partition
  !>         from the partition object associated with this
  !>         mesh object
  !============================================================================
  function get_num_cells_core( self ) result ( core_cells )

    implicit none

    class(mesh_type), intent(in) :: self
    integer(i_def)               :: core_cells

    core_cells = self%partition%get_num_cells_core()

  end function get_num_cells_core


  !> @details Get the number of owned cells from the partition object
  !> @return  The total number of core cells on the local partition
  !============================================================================
  function get_num_cells_owned( self ) result ( owned_cells )

    implicit none

    class(mesh_type), intent(in) :: self
    integer(i_def)               :: owned_cells

    owned_cells = self%partition%get_num_cells_owned()

  end function get_num_cells_owned


  !> @details Returns the maximum depth of the halo from the partition object
  !> @return  The maximum depth of halo cells
  !============================================================================
  function get_halo_depth( self ) result ( halo_depth )
    implicit none

    class(mesh_type), intent(in) :: self
    integer(i_def)               :: halo_depth

    halo_depth = self%partition%get_halo_depth()

  end function get_halo_depth


  !> @details Returns the total number of halo cells in a particular depth
  !>          of halo in a 2d slice from the partition object
  !> @param[in] depth       The depth of the halo being queried
  !> @return                The total number of halo cells of the particular
  !>                        depth on the local partition
  !============================================================================
  function get_num_cells_halo( self, depth ) result ( halo_cells )

    implicit none

    class(mesh_type), intent(in) :: self
    integer(i_def),   intent(in) :: depth
    integer(i_def)               :: halo_cells

    if (depth > self%get_halo_depth()) then
      halo_cells = 0
    else
      halo_cells = self%partition%get_num_cells_halo(depth)
    end if

  end function get_num_cells_halo


  !> @details Get the total number of ghost cells in a slice around
  !>          the local partition
  !> @return  The total number of ghost cells around the local partition
  !============================================================================
  function get_num_cells_ghost( self ) result ( ghost_cells )

    implicit none

    class(mesh_type), intent(in) :: self
    integer(i_def)               :: ghost_cells

    ghost_cells = self%partition%get_num_cells_ghost()

  end function get_num_cells_ghost


  !> @details Returns the global index of the cell that corresponds
  !>          to the given local index on the local partition
  !> @param[in] cell_lid  The id of a cell in local index space
  !> @return              The id of a cell in global index space
  !============================================================================
  function get_gid_from_lid( self, cell_lid ) result ( cell_gid )

    implicit none

    class(mesh_type), intent(in) :: self
    integer(i_def),   intent(in) :: cell_lid  ! local index
    integer(i_def)               :: cell_gid  ! global index

    cell_gid = self%partition%get_gid_from_lid(cell_lid)

  end function get_gid_from_lid


  !> @details Returns count of colours used in colouring mesh.
  !> @return          Number of colours used to colour this mesh.
  !============================================================================
  function get_ncolours(self) result(ncolours)
    implicit none
    class(mesh_type), intent(in) :: self
    integer(i_def)               :: ncolours

    ncolours = self%ncolours

  end function get_ncolours

  !============================================================================
  !> @brief Populates args with colouring info.
  !> @param[in] self  The mesh_type instance.
  !> @param[out] ncolours  Number of colours used to colour this mesh. 
  !> @param[out] ncells_per_colour  Count of cells in each colour.
  !> @param[out] colour_map         Indices of cells in each colour.
  !============================================================================
  subroutine get_colours(self, ncolours, ncells_per_colour, colour_map)
    implicit none
    class(mesh_type), intent(in), target      :: self
    integer(i_def), intent(out)               :: ncolours
    integer(i_def), pointer, intent(out)  :: ncells_per_colour(:)
    integer(i_def), pointer, intent(out)  :: colour_map(:,:)
  

    ncolours = self%ncolours
    ncells_per_colour => self%ncells_per_colour
    colour_map => self%cells_in_colour

  end subroutine get_colours

  !============================================================================
  !> @brief  Returns state of colouring: has colouring yet been applied to
  !>         this mesh?
  !>
  !> @return  Logical true is mesh coloured, false if not. 
  !============================================================================
  function is_coloured(self) result(cstat)
    implicit none
    class(mesh_type), intent(in) :: self
    logical(l_def)               :: cstat

    if(self%ncolours <= 0) then
      cstat = .false.
    else
      cstat = .true.
    end if

  end function is_coloured 

  !============================================================================
  !> @brief  Invoke calculation of colouring for this mesh.
  !>
  !> @param[in] self  The mesh_type instance.
  !============================================================================
  subroutine set_colours(self)
    use mesh_colouring_mod, only : colour_mod_set_colours => set_colours
    implicit none
    class(mesh_type), intent(inout) :: self


    call colour_mod_set_colours(self%get_ncells_2d(), &
                                self%cell_next, &
                                self%ncolours, &
                                self%ncells_per_colour, &
                                self%cells_in_colour)

  end subroutine set_colours


!==============================================================================
! PRIVATE ROUTINES/FUNCTIONS - Internal to object
!==============================================================================

  !----------------------------------------------------------------------------
  !> @details Called during <code>mesh_object<code>_construction, this routine
  !>          calculates partition statistics using the list of cells in the
  !>          partition from the <code>partition_object<code> and the
  !>          connectivity information of the same cells in the
  !>          <code>global_mesh_object<code>. Some information is added to the
  !>          mesh_object, while other info is held in module and used by the
  !>          subroutine mesh_extruder
  !>
  !> @param [in]  partition   A partition_object specifying which 2d-cells in
  !>                          the global_mesh are to be included in this
  !>                          partition
  !> @param [in]  global_mesh A global mesh object that describes the layout
  !>                          of the global mesh
  !> @param [out] self        The mesh_object with updated the partition
  !>                          statistics
  subroutine get_partition_stats(self, partition, global_mesh)

    implicit none

    class (mesh_type)                    :: self
    type  (partition_type),   intent(in) :: partition
    type  (global_mesh_type), intent(in) :: global_mesh

    integer (i_def) :: nverts_per_2d_cell
    integer (i_def) :: nedges_per_2d_cell

    integer (i_def) :: max_num_vertices_2d
    integer (i_def) :: max_num_edges_2d
    integer (i_def) :: counter, i, j


    integer (i_def) :: cell_gid
    integer (i_def) :: edge_gid
    integer (i_def) :: vert_gid
    integer (i_def) :: n_uniq_verts
    integer (i_def) :: n_uniq_edges


    integer (i_def), allocatable :: tmp_list(:)
    integer (i_def) :: tmp_int
    logical (l_def) :: edge_gid_present
    logical (l_def) :: vert_gid_present

    ! Arrays where the connected entity ids are global ids
    integer(i_def), allocatable :: &
      vert_on_cell_2d_gid (:,:)    &! Vertices connected to local 2d cell.
    , edge_on_cell_2d_gid (:,:)    &! Edges connected to local 2d cell.
    , cell_next_2d_gid    (:,:)     ! Local 2d cell connectivity.

    integer(i_def), allocatable :: &
      vert_lid_gid_map(:)           ! Vertices local-to-global id map


    ! Get global mesh statistics to size connectivity arrays
    nverts_per_2d_cell   = global_mesh%get_nverts_per_cell()
    nedges_per_2d_cell   = global_mesh%get_nedges_per_cell()

    self%nverts_per_cell = 2*nedges_per_2d_cell
    self%nedges_per_cell = 2*nedges_per_2d_cell + nverts_per_2d_cell
    self%nfaces_per_cell = nedges_per_2d_cell + 2


    ! Get partition statistics
    max_num_vertices_2d  = ncells_2d*nverts_per_2d_cell
    max_num_edges_2d     = ncells_2d*nedges_per_2d_cell


    ! Allocate arrays to hold partition connectivities, as global ids
    ! These will be deallocated in mesh_extruder
    allocate( vert_on_cell_2d_gid (nverts_per_2d_cell, ncells_2d) )
    allocate( edge_on_cell_2d_gid (nedges_per_2d_cell, ncells_2d) )
    allocate( cell_next_2d_gid    (nedges_per_2d_cell, ncells_2d) )


    ! Get 2d cell lid/gid map
    allocate( self % cell_lid_gid_map(ncells_2d) )

    do i=1, ncells_2d

      ! Note: For a single partition, the global ids should be the same
      !       as the local cell ids
      cell_gid = partition % get_gid_from_lid(i)
      self % cell_lid_gid_map(i) = cell_gid

      call global_mesh % get_vert_on_cell (cell_gid, vert_on_cell_2d_gid(:,i))
      call global_mesh % get_edge_on_cell (cell_gid, edge_on_cell_2d_gid(:,i))
      call global_mesh % get_cell_next    (cell_gid, cell_next_2d_gid(:,i))

    end do


    !--------------------------------------------------------------------------
    ! Now get a list of all unique entities in partition (cells/vertices/edges)
    !--------------------------------------------------------------------------
    ! A. Generate cell_next for partition, in local ids
    allocate( cell_next_2d (nedges_per_2d_cell, ncells_2d) )

    ! Set default to 0, For missing cells, i.e. at the edge of a partition
    ! there is no cell_next. So default to zero if a cell_next is not found
    cell_next_2d(:,:) = 0

    do i=1, ncells_2d
      do j=1, nedges_per_2d_cell
        ! For a 2D cell in the parition, get the each global id of its
        ! adjacent cell.
        cell_gid = cell_next_2d_gid(j,i)

        do counter=1, ncells_2d
          ! Loop over all cells in the partition, getting the cell
          ! global id.
          tmp_int = partition % get_gid_from_lid(counter)

          ! Set the adjacent cell id to the local id if
          ! it's global appears in the partition
          if (tmp_int == cell_gid) then
            cell_next_2d(j,i) = counter
            exit
          end if
        end do

      end do
    end do

    deallocate( cell_next_2d_gid )


    !--------------------------------------------------------------------------
    ! B. Get global ids of vertices in partition
    ! Get global ids of all vertices on partition, by looping over the
    ! cell-vertex connectivity on the partition.
    allocate( tmp_list(max_num_vertices_2d) )
    allocate( vert_on_cell_2d (nverts_per_2d_cell, ncells_2d) )

    n_uniq_verts = 0

    do i=1, ncells_2d
      do j=1, nverts_per_2d_cell

        vert_gid = vert_on_cell_2d_gid(j,i)
        vert_gid_present = .false.

        do counter=1, n_uniq_verts
          if (tmp_list(counter) == vert_gid) then
            vert_gid_present = .true.
            vert_on_cell_2d(j,i) = counter
            exit
          end if
        end do

        if (vert_gid_present .eqv. .false.) then
          n_uniq_verts = n_uniq_verts + 1
          tmp_list(n_uniq_verts) = vert_gid
          vert_on_cell_2d(j,i) = n_uniq_verts
        end if

      end do
    end do

    deallocate (vert_on_cell_2d_gid)

    allocate(vert_lid_gid_map(n_uniq_verts))
    vert_lid_gid_map(:) = tmp_list(1:n_uniq_verts)

    nverts_2d        = n_uniq_verts
    self % nverts_2d = n_uniq_verts
    deallocate(tmp_list)


    !-------------------------------------------------------------------------
    ! C. Get global ids of edges in partition
    ! Get global ids of all edges on partition, by looping over the cell-edge
    ! connectivity on the partition.
    allocate( tmp_list(max_num_edges_2d) )
    allocate( edge_on_cell_2d (nedges_per_2d_cell, ncells_2d) )

    n_uniq_edges = 0

    do i=1, ncells_2d            ! cells in local order in partition
      do j=1, nedges_per_2d_cell ! edges from that cell

        ! Get the global id of the edge
        edge_gid = edge_on_cell_2d_gid(j,i)

        ! Initialise the flag, we assume this gedge global id has not been
        ! encountered yet
        edge_gid_present = .false.

        ! Loop over the number of existing uniq edges
        do counter=1, n_uniq_edges
          ! Test to see if this global edge id is in the list
          if (tmp_list(counter) == edge_gid) then
            ! Set flag if it's present
            edge_gid_present = .true.

            ! Mark this edge as the local id, which is it's index in the
            ! array of unique ids
            edge_on_cell_2d(j,i) = counter

            ! Found this edge already so stop looping
            exit
          end if
        end do

        ! If the edge gid was not found
        if (edge_gid_present .eqv. .false.) then
          ! This is a new edge gid. So increment the number of unique
          ! edges by 1
          n_uniq_edges = n_uniq_edges + 1

          ! Add this edge gid to the next element in the array
          tmp_list(n_uniq_edges) = edge_gid

          ! and map it's local id
          edge_on_cell_2d(j,i) = n_uniq_edges
        end if

      end do
    end do

    deallocate(edge_on_cell_2d_gid)

    nedges_2d        = n_uniq_edges
    self % nedges_2d = n_uniq_edges

    deallocate(tmp_list)

    !--------------------------------------------------------------------------
    ! Get partition vertices lat-lon-z coords, Note z=surface height
    !--------------------------------------------------------------------------
    allocate(vertex_coords_2d(3,n_uniq_verts))
    do i=1, n_uniq_verts
      ! Get coords of vertices
      vert_gid = vert_lid_gid_map(i)
      call global_mesh % get_vert_coords(vert_gid,vertex_coords_2d(:,i))
    end do

    deallocate( vert_lid_gid_map )

    nverts_3d = nverts_2d * (nlayers+1)
    nedges_3d = nedges_2d * (nlayers+1) + nverts_2d*nlayers
    nfaces_3d = ncells_2d * (nlayers+1) + nedges_2d*nlayers

  end subroutine get_partition_stats


  !----------------------------------------------------------------------------
  !> @details Called during <code>mesh_object<code>_construction, this routine
  !>          extrudes the a surface 2D-mesh into the 3D-mesh object.
  !>
  !> @param self The mesh_object being constructed on the partition
  subroutine mesh_extruder(self)

    use coord_transform_mod,   only: llr2xyz
    use reference_element_mod, only: W, S, E, N, B, T,                        &
                                     nverts_h, nedges_h, nfaces_h,            &
                                     nverts,   nedges,   nfaces,              &
                                     SWB, SEB, NEB, NWB, SWT, SET, NET, NWT
    use configuration_mod,     only: l_spherical, earth_radius
    implicit none

    class(mesh_type) :: self

    ! Loop indices
    integer(i_def) :: i, j, k, id, jd, iedge, ivert, base_id

    ! From reference element
    ! nverts_h = number of vertices on 2d horizontal cell
    ! nedges_h = number of edges    on 2d horizontal cell
    !
    ! nverts   = number of vertices on 3d cell
    ! nedges   = number of edges    on 3d cell
    ! nfaces   = number of faces    on 3d cell

    integer(i_def) :: nverts_per_3d_cell
    integer(i_def) :: nedges_per_3d_cell
    integer(i_def) :: nfaces_per_3d_cell

    integer(i_def) :: nverts_per_2d_cell
    integer(i_def) :: nedges_per_2d_cell

    ! lat/long coordinates
    real(r_def) :: long, lat, r

    ! The height of the lowest z-level
    real(r_def) :: base_z

    ! Reference element stats
    nverts_per_2d_cell = nverts_h
    nedges_per_2d_cell = nedges_h

    nverts_per_3d_cell = nverts
    nedges_per_3d_cell = nedges
    nfaces_per_3d_cell = nfaces

    ! Allocate mesh object arrays
    allocate ( self % vertex_coords ( 3, nverts_3d ) )

    ! Apply default cell_next values
    self%cell_next(:,:) = 0

    do i=1, ncells_2d
      do ivert=1, nverts_per_2d_cell
        self % vert_on_cell(ivert,i) = vert_on_cell_2d(ivert,i)
      end do

      do iedge=1, nedges_per_2d_cell
        self % cell_next(iedge,i) = cell_next_2d(iedge,i)
      end do
    end do

    ! Add connectivity for up/down
    ! index nfaces_h+1 is the bottom of the 3d cell
    !                  (set to zero as it's the surface)
    ! index nfaces_h+2 is the top of the 3d cell
    do j=1,ncells_2d
      self%cell_next(nfaces_h+2,j) = j + ncells_2d
    end do

    ! Perform vertical extrusion for connectivity
    do k=1, nlayers-1
      do i=1, ncells_2d
        id = i  + k*ncells_2d
        jd = id - ncells_2d

        do j=1, nfaces_h ! only over vertical faces
          if (self % cell_next(j,jd) /= 0)                                    &
                    self % cell_next(j,id) = self % cell_next(j,jd) + ncells_2d
        end do

        self%cell_next(nfaces_h+1,id) = id - ncells_2d
        self%cell_next(nfaces_h+2,id) = id + ncells_2d

        if (k==nlayers-1) self%cell_next(nfaces_h+2,id) = 0

      end do
    end do

    ! NOTE: dz and domain top will depend on
    !       the orography, vertical resolution file, top of model
    !       and how vertical and horizontal smoothing is applied
    !       after the application of orography. At present, the
    !       vertical depth is hard-coded uniformly to 1.0

    ! The assumption is that the global mesh coords are provided in
    ! [longitude, latitude, radius] (long/lat in rads)

    ! perform vertical extrusion for vertices
    if( l_spherical )then
      !> @todo We shouldn't be using earth_radius here - it should be
      !!       some form of scaled planet radius - but that is a much
      !!       bigger change for a different ticket.
      base_z = earth_radius
    else
      base_z = 0.0_r_def
    end if
    do j=1, nverts_2d
     ! k = 0
     self%vertex_coords(1,j) = vertex_coords_2d(1,j)
     self%vertex_coords(2,j) = vertex_coords_2d(2,j)
     self%vertex_coords(3,j) = base_z
     ! k = 1, nlayers
      do k=1, nlayers
        self%vertex_coords(1,j+k*nverts_2d) = vertex_coords_2d(1,j)
        self%vertex_coords(2,j+k*nverts_2d) = vertex_coords_2d(2,j)
        self%vertex_coords(3,j+k*nverts_2d) = base_z + real(k)*self%dz(k)
      end do
    end do

    deallocate(vertex_coords_2d)

    if( l_spherical )then
      ! Convert (long,lat,r) -> (x,y,z)
      do j=1, nverts_2d
        do k=0, nlayers
          long = self % vertex_coords(1,j+k*nverts_2d)
          lat  = self % vertex_coords(2,j+k*nverts_2d)
          r    = self % vertex_coords(3,j+k*nverts_2d)
          call llr2xyz(long,lat,r,self % vertex_coords(1,j+k*nverts_2d), &
                                  self % vertex_coords(2,j+k*nverts_2d), &
                                  self % vertex_coords(3,j+k*nverts_2d))
        end do
      end do
    end if

    ! assign vertices to cells
    ! Loop over lowest layer of cells first, to set the cell ids above
    ! the lowest layer
    do i=1, ncells_2d
      do ivert=1, nverts_per_2d_cell
        self % vert_on_cell(nverts_per_2d_cell + ivert, i) =                  &
          self % vert_on_cell(ivert,i) + nverts_2d
      end do
    end do

    ! Do vertical extrusion
    do base_id=1, ncells_2d
      do k=1, nlayers-1

        i = base_id + k*ncells_2d
        do ivert=1, nverts_per_3d_cell
          self % vert_on_cell(ivert,i) = self % vert_on_cell(ivert,base_id)   &
                                       + k*nverts_2d
        end do

      end do
    end do

    deallocate (vert_on_cell_2d)
    deallocate (edge_on_cell_2d)
    deallocate (cell_next_2d)


    ! Diagnostic information
    call log_event('grid connectivity', LOG_LEVEL_DEBUG)
    do i=1, ncells
      write(log_scratch_space,'(7i6)') i,                                     &
        self % cell_next(S,i), self % cell_next(E,i),                         &
        self % cell_next(N,i), self % cell_next(W,i),                         &
        self % cell_next(B,i), self % cell_next(T,i)
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end do

    call log_event('verts on cells', LOG_LEVEL_DEBUG)
    do i=1, ncells
      write(log_scratch_space, '(9i6)') i,                                    &
        self % vert_on_cell(SWB,i), self % vert_on_cell(SEB,i),               &
        self % vert_on_cell(NEB,i), self % vert_on_cell(NWB,i),               &
        self % vert_on_cell(SWT,i), self % vert_on_cell(SET,i),               &
        self % vert_on_cell(NET,i), self % vert_on_cell(NWT,i)
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end do

    call log_event('vert coords', LOG_LEVEL_DEBUG)
    do i=1, nverts_3d
      write(log_scratch_space, '(i6,3ES20.10E3)') i,                          &
        self%vertex_coords(1,i),                                              &
        self%vertex_coords(2,i),                                              &
        self%vertex_coords(3,i)
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end do

  end subroutine mesh_extruder


  !----------------------------------------------------------------------------
  !> @details Compute mesh connectivity.
  !>
  !> <pre>
  !> cells->faces (3,2)
  !> cells->edges (3,1)
  !> </pre>
  !>
  !> @param ncells Total number of cells in horizontal layer.
  !
  subroutine mesh_connectivity( self )

  use reference_element_mod, only: W,   S,  E,  N,  B,  T,                    &
                                   EB, ET, WB, WT, NB, NT, SB, ST,            &
                                   SW, SE, NW, NE,                            &
                                   nedges_h, nverts_h, nfaces_h,              &
                                   SWB, SEB, NEB, NWB, SWT, SET, NET, NWT

  implicit none

  class (mesh_type) :: self

  integer (i_def) :: &
    cell             &! cell loop index
  , face             &! face loop index
  , edge             &! edge loop index
  , vert             &! vert loop index
  , face_id          &! unique face index
  , edge_id          &! unique edge index
  , edge_upper       &! index of edge on top face
  , cell_nbr         &! cell index for a neightbour cell
  , face_nbr         &! face index for a face on a neighbour cell
  , edge_nbr         &! edge index for a edge on a neighbour cell
  , vert_nbr         &! vert index for a edge on a neighbour cell
  , cell_nbr_nbr     &! cell index for a neighbour cell of a neighbour cell
  , vert_nbr_nbr      ! vert index for a neighbour cell of a neighbour cell

  allocate( self % face_on_cell(self%nfaces_per_cell, ncells_2d) )
  allocate( self % edge_on_cell(self%nedges_per_cell, ncells_2d) )

  self%face_on_cell(:,:) = 0
  self%edge_on_cell(:,:) = 0


  !==================================================
  ! Compute index of faces on the cell
  !==================================================
  face_id = 1
  do cell=1, ncells_2d

    ! First do faces on E,S,W,N sides of cell
    do face=1, nfaces_h

      if ( self%face_on_cell(face,cell) == 0 ) then
        ! Assign face_id to this face as it is not already assigned

        self%face_on_cell(face,cell) = face_id

        ! Find matching face on the neighbouring face.
        ! Note: The orientations of cells may not be the same,
        !       e.g. Cubedsphere, for a so search all faces of
        !       cell_next on the neighbouring cell
        cell_nbr = self%cell_next(face,cell)
        if (cell_nbr /= 0) then
          do face_nbr=1, nfaces_h
            if ( self%cell_next(face_nbr,cell_nbr) == cell ) then
              ! Found matching face, assign in and increment count
              self%face_on_cell(face_nbr,cell_nbr) = face_id
            end if
          end do
        end if

        face_id = face_id + 1

      end if ! test for assigned face

    end do  ! n_faces_h

    ! Now do Faces on B and T of cell
    do face=nfaces_h+1, self%nfaces_per_cell
      self%face_on_cell(face, cell) = face_id
      face_id = face_id + 1
    end do

  end do ! ncells_2d

  ! Compute the index of edges on the cell
  edge_id = 1
  do cell=1, ncells_2d
    ! horizontal edges ( edges on Bottom and Top faces)
    ! This uses the fact that edges of the bottom face correspond to the
    ! vertical faces, i.e edge i is the bottom edge of face i and hence
    ! we can use cell_next array (which is addressed through faces) for
    ! edge indexes
    do edge=1, nedges_h
      if ( self%edge_on_cell(edge,cell) == 0 ) then

        ! Index of edge on bottom face
        self%edge_on_cell(edge,cell) = edge_id

        ! Corresponding index of edge on top face
        edge_upper = edge + nedges_h + nverts_h
        self%edge_on_cell(edge_upper,cell) = edge_id + 1

        ! find matching edge in neighbouring cell
        cell_nbr = self%cell_next(edge,cell)
        if (cell_nbr /= 0) then
          do edge_nbr = 1,nedges_h
            if ( self%cell_next(edge_nbr,cell_nbr) == cell ) then
              ! Found cell which shares edge
              self%edge_on_cell(edge_nbr,cell_nbr) = edge_id
              edge_upper = edge_nbr+nedges_h+nverts_h
              self%edge_on_cell(edge_upper,cell_nbr) = edge_id + 1
            end if
          end do
        end if

        edge_id = edge_id + 2

      end if
    end do

    ! vertical edges (edges not on Bottom and Top faces)
    ! This uses the fact that vertical edges correspond to vertices of
    ! the bottom face i.e edge i has the same index as vertex i and hence
    ! we can use vert_on_cell array (which is addressed through verts) for
    ! edge indexes
    do edge = nedges_h+1,nedges_h+nverts_h
      if ( self%edge_on_cell(edge,cell) == 0 ) then
        self%edge_on_cell(edge,cell) = edge_id

        ! Find matching edge on two neighbouring cells
        ! this edge is an extrusion of the corresponding vertex
        vert = self%vert_on_cell(edge-nedges_h,cell)
        do face = 1,nfaces_h
          cell_nbr = self%cell_next(face,cell)
          if (cell_nbr /= 0) then
            do vert_nbr = 1,nverts_h
              if ( self%vert_on_cell(vert_nbr,cell_nbr) == vert ) then
                ! Found matching vert in neighbour cell
                self%edge_on_cell(vert_nbr+nedges_h,cell_nbr) = edge_id
                do face_nbr = 1,nfaces_h
                  ! Now find matching vert in neighbour of neighbour
                  cell_nbr_nbr = self%cell_next(face_nbr,cell_nbr)
                  if (cell_nbr_nbr /= 0) then
                    do vert_nbr_nbr = 1,nverts_h
                      if ( self%vert_on_cell(vert_nbr_nbr,cell_nbr_nbr) ==    &
                           vert ) then
                        self%edge_on_cell(vert_nbr_nbr+nedges_h,              &
                                          cell_nbr_nbr) = edge_id
                      end if
                    end do
                  end if ! cell_nbr_nbr
                end do
              end if
            end do
          end if
        end do
        edge_id = edge_id + 1
      end if
    end do
  end do

  call log_event( 'faces on cells', LOG_LEVEL_DEBUG )
  do cell = 1, ncells_2d
    write( log_scratch_space, '(7i6)' ) cell,                                 &
           self%face_on_cell(S,cell), self%face_on_cell(E,cell),              &
           self%face_on_cell(N,cell), self%face_on_cell(W,cell),              &
           self%face_on_cell(B,cell), self%face_on_cell(T,cell)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do

  call log_event( 'edges on cells', LOG_LEVEL_DEBUG )
  do cell = 1, ncells_2d
    write( log_scratch_space, '(13i6)' ) cell,                                &
           self%edge_on_cell(SB,cell), self%edge_on_cell(EB,cell),            &
           self%edge_on_cell(NB,cell), self%edge_on_cell(WB,cell),            &
           self%edge_on_cell(SW,cell), self%edge_on_cell(SE,cell),            &
           self%edge_on_cell(NE,cell), self%edge_on_cell(NW,cell),            &
           self%edge_on_cell(ST,cell), self%edge_on_cell(ET,cell),            &
           self%edge_on_cell(NT,cell), self%edge_on_cell(WT,cell)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do

  end subroutine mesh_connectivity



  !----------------------------------------------------------------------------
  !> @details Compute the domain limits (x,y,z) for
  !           Cartesian domains and (lambda,phi,r) for spherical domains
  subroutine set_domain_size(self)

    implicit none

    class(mesh_type) :: self

    if ( l_spherical ) then
      self % domain_size%minimum%x =  0.0_r_def
      self % domain_size%maximum%x =  2.0_r_def*PI
      self % domain_size%minimum%y = -0.5_r_def*PI
      self % domain_size%maximum%y =  0.5_r_def*PI
      self % domain_size%minimum%z =  0.0_r_def
      self % domain_size%maximum%z =  self % domain_top
    else
      self % domain_size%minimum%x =  minval( self % vertex_coords(1,:))
      self % domain_size%maximum%x =  maxval( self % vertex_coords(1,:))
      self % domain_size%minimum%y =  minval( self % vertex_coords(2,:))
      self % domain_size%maximum%y =  maxval( self % vertex_coords(2,:))
      self % domain_size%minimum%z =  minval( self % vertex_coords(3,:))
      self % domain_size%maximum%z =  maxval( self % vertex_coords(3,:))
    end if

    !> @todo Need to do a global reduction of maxs and mins when the
    !> code is parallel

  end subroutine set_domain_size


  !============================================================================
  !> @details Set up and calculate vertical coordinate. Called during
  !>          @c mesh_type object construction, this routine
  !>          sets up the vertical grid layers and type of vertical grid
  !> @param [inout] self          The mesh_object with set up vertical
  !>                              coordinate eta and layer spacing dz
  !> @param [in]    vgrid_option  Choice of vertical grid
  subroutine set_vertical_coordinate(self, vgrid_option)

    implicit none

    type(mesh_type), intent(inout) :: self
    integer(i_def),  intent(in)    :: vgrid_option

    integer(i_def) :: k, nlayers
    real   (r_def) :: stretching_factor, phi_flatten, delta_eta, eta_uni

    ! Get the number of layers
    nlayers = self%nlayers

    ! Calculate eta depending on uniform/stretching option
    select case (vgrid_option)

      ! UNIFORM GRID (constant delta_eta)
      case (VGRID_UNIFORM)
        do k = 0, nlayers
          self%eta(k) = real(k,r_def)/real(nlayers,r_def)
        end do

      ! QUADRATIC GRID: eta(k) = (k/numlayers)^2
      case (VGRID_QUADRATIC)
        do k = 0, nlayers
          self%eta(k) = ( real(k,r_def)/real(nlayers,r_def) )**2_i_def
        end do

      ! GEOMETRIC GRID
      ! Source: John Thuburn's ENDGame code for staggered grid.
      !         deta = (stretch - 1.0d0)/(stretch**(2*nz) - 1.0d0)
      ! Here:   The grid is non-staggered grid so it must be
      !         deta = (stretch - 1.0d0)/(stretch**(nz) - 1.0d0)
      case (VGRID_GEOMETRIC)
        stretching_factor = 1.03_r_def
        self%eta(0) = 0.0_r_def
        delta_eta = ( stretching_factor - 1.0_r_def ) / &
                    ( stretching_factor**(nlayers) - 1.0_r_def )
        do k = 1, nlayers
          self%eta(k) = self%eta(k-1) + delta_eta
          delta_eta = delta_eta*stretching_factor
        end do

      ! DCMIP GRID
      ! Source: DCMIP-TestCaseDocument_v1.7.pdf, Appendix F.2. - Eq. 229)
      ! phi_flatten is a flattening parameter (usually phi_flatten = 15)
      case (VGRID_DCMIP)
        phi_flatten = 15.0_r_def
        do k = 0, nlayers
          eta_uni = real(k,r_def)/real(nlayers,r_def)
          self%eta(k) = ( sqrt(phi_flatten*(eta_uni**2_i_def) + 1.0_r_def) &
                          - 1.0_r_def ) / &
                        ( sqrt(phi_flatten + 1.0_r_def) - 1.0_r_def )
        end do

      ! Default case - automatically make the uniform grid
      case default
        do k = 0, nlayers
          self%eta(k) = real(k,r_def)/real(nlayers,r_def)
        end do
      end select

    ! Calculate dz
    do k = 1, nlayers
       self%dz(k) = ( self%eta(k) - self%eta(k-1) )*self%domain_top
    end do

  end subroutine set_vertical_coordinate

!-----------------------------------------------------------------------------
! Mesh destructor
!-----------------------------------------------------------------------------
  subroutine mesh_destructor(self)

    implicit none

    type (mesh_type) :: self


    if (allocated(self%cell_lid_gid_map)) deallocate(self%cell_lid_gid_map)
    if (allocated(self%cell_next))        deallocate(self%cell_next)
    if (allocated(self%vert_on_cell))     deallocate(self%vert_on_cell)
    if (allocated(self%face_on_cell))     deallocate(self%face_on_cell)
    if (allocated(self%edge_on_cell))     deallocate(self%edge_on_cell)
    if (allocated(self%vertex_coords))    deallocate(self%vertex_coords)
    if (allocated(self%vert_cell_owner))  deallocate(self%vert_cell_owner)
    if (allocated(self%edge_cell_owner))  deallocate(self%edge_cell_owner)
    if (allocated(self%edge_ownership))   deallocate(self%edge_ownership)
    if (allocated(self%vertex_ownership)) deallocate(self%vertex_ownership)
    if (allocated(self%eta))              deallocate(self%eta)
    if (allocated(self%dz))               deallocate(self%dz)

    ! NB Allocation tests necessary as dtor can be called on
    ! temporary unconstructed object; e.g. during assignment

    if(allocated(self%ncells_per_colour)) then
      deallocate(self%ncells_per_colour)
    end if

    if(allocated(self%cells_in_colour)) then
      deallocate(self%cells_in_colour)
    end if


  end subroutine mesh_destructor

  !> Assignment operator between mesh_type pairs.
  !>
  !> @param[out] dest   field_type lhs
  !> @param[in]  source field_type rhs

  subroutine mesh_type_assign ( dest, source)

    use reference_element_mod, only: nfaces, nverts


    implicit none

    class(mesh_type), intent(out)     :: dest
    class(mesh_type), intent(in)      :: source


    dest%mesh_id    = source%mesh_id
    dest%partition  = source%partition
    dest%ncells_2d  = source%ncells_2d
    dest%nlayers    = source%nlayers
    dest%ncells     = source%ncells
    dest%domain_top = source%domain_top
    dest%ncolours   = source%ncolours


    ! Copy vertical coordinates eta[0,1] and dz
    allocate( dest%eta ( 0:source%nlayers ) )
    allocate( dest%dz  ( source%nlayers   ) )

    dest%eta(:) = source%eta(:)
    dest%dz(:) = source%dz(:)

    ! Copy vertex coords
    allocate ( dest % vertex_coords ( 3, source%nverts ) )
    dest%vertex_coords(:,:) = source%vertex_coords(:,:)

    dest%nverts_2d = source%nverts_2d
    dest%nedges_2d = source%nedges_2d
    dest%ncells_2d = source%ncells_2d
    dest%ncells_2d_with_ghost = source%ncells_2d_with_ghost

    dest%nverts = source%nverts
    dest%nedges = source%nedges
    dest%nfaces = source%nfaces
    dest%ncells = source%ncells
    dest%ncells_with_ghost = source%ncells_with_ghost

    dest%nverts_per_cell = source%nverts_per_cell
    dest%nedges_per_cell = source%nedges_per_cell
    dest%nfaces_per_cell = source%nfaces_per_cell

    ! Copy 2d cell lid/gid map
    allocate( dest % cell_lid_gid_map(source%ncells_2d_with_ghost) )
    dest%cell_lid_gid_map(:) = source%cell_lid_gid_map(:)

    ! Calculate next-to cells and vertices on cells
    allocate ( dest % cell_next    ( nfaces, source%ncells_with_ghost ) )
    allocate ( dest % vert_on_cell ( nverts, source%ncells_with_ghost ) )

    dest%cell_next(:,:) = source%cell_next(:,:)
    dest%vert_on_cell(:,:) = source%vert_on_cell(:,:)


    allocate( dest % face_on_cell(source%nfaces_per_cell, source%ncells_2d_with_ghost) )
    allocate( dest % edge_on_cell(source%nedges_per_cell, source%ncells_2d_with_ghost) )

    dest%face_on_cell(:,:) = source%face_on_cell(:,:)
    dest%edge_on_cell(:,:) = source%edge_on_cell(:,:)


    ! copy domain size

    dest % domain_size%minimum%x =  source % domain_size%minimum%x
    dest % domain_size%maximum%x =  source % domain_size%maximum%x
    dest % domain_size%minimum%y =  source % domain_size%minimum%y
    dest % domain_size%maximum%y =  source % domain_size%maximum%y
    dest % domain_size%minimum%z =  source % domain_size%minimum%z
    dest % domain_size%maximum%z =  source % domain_size%maximum%z
    

    ! Copy ownership of cell vertices and cell edges
    allocate( &
      dest%vertex_ownership (lbound(source%vertex_ownership,1):ubound(source%vertex_ownership,1) &
                             , lbound(source%vertex_ownership,2):ubound(source%vertex_ownership,2) ) )
    allocate( &
      dest%edge_ownership (lbound(source%edge_ownership,1):ubound(source%edge_ownership,1) &
                           , lbound(source%edge_ownership,2):ubound(source%edge_ownership,2) ) )
    allocate( &
      dest%vert_cell_owner (lbound(source%vert_cell_owner,1):ubound(source%vert_cell_owner,1) &
                            , lbound(source%vert_cell_owner,2):ubound(source%vert_cell_owner,2) ) )
    allocate( &
      dest%edge_cell_owner (lbound(source%edge_cell_owner,1):ubound(source%edge_cell_owner,1) &
                           , lbound(source%edge_cell_owner,2):ubound(source%edge_cell_owner,2) ) )

    dest%vertex_ownership(:,:) = source%vertex_ownership(:,:)
    dest%edge_ownership(:,:) = source%edge_ownership(:,:)
    dest%vert_cell_owner(:,:) = source%vert_cell_owner(:,:)
    dest%edge_cell_owner(:,:) = source%edge_cell_owner(:,:)


  end subroutine mesh_type_assign


  !============================================================================
  ! This routine is only available when setting data for unit testing.
  !============================================================================
  !> @brief     Stucture-Constructor (for unit testing)
  !> @param[in] mesh_cfg Sets the type of test mesh returned.
  !>                     [PLANE|PLANE_BI_PERIODIC].
  !>                     PLANE - returns a 5-layer non-biperiodic mesh
  !>                     PLANE_BI_PERIODIC - returns a 3-layer bi-periodic mesh
  !> @return             A 3D-Mesh object based on a 3x3-cell global mesh
  !>                     with one partition.
  !============================================================================
  function mesh_constructor_unit_test_data(mesh_cfg) result (self)

    ! Mesh returned is based on a 3x3 partition with the following
    ! cell numbering.
    !
    !    +-----------+
    !    | 7 | 8 | 9 |
    !    |---+---+---|
    !    | 4 | 5 | 6 |
    !    |---+---+---|
    !    | 1 | 2 | 3 |
    !    +---+---+---+
    !

    implicit none

    integer(i_def), intent(in) :: mesh_cfg
    type(mesh_type), target    :: self

    integer(i_def), parameter :: nverts_per_cell    = 8
    integer(i_def), parameter :: nedges_per_cell    = 12
    integer(i_def), parameter :: nverts_per_2d_cell = 4
    integer(i_def), parameter :: nedges_per_2d_cell = 4


    mesh_id_counter = mesh_id_counter+1_i_def

    self%mesh_id         = mesh_id_counter
    self%partition       = partition_type()
    self%nverts_per_cell = 8
    self%nedges_per_cell = 12
    self%nfaces_per_cell = 6
    self%ncolours        = -1  ! Initialise ncolours to error status


    if (mesh_cfg == PLANE) then
      self%domain_top = 10000.0_r_def

      self%ncells_2d = 9
      self%nverts_2d = 16
      self%nedges_2d = 24

      self%nlayers = 5
      self%ncells  = 45
      self%nverts  = 96
      self%nfaces  = 156
      self%nedges  = 224

    else if (mesh_cfg == PLANE_BI_PERIODIC) then
      self%domain_top = 6000.0_r_def

      ! 3x3x3 mesh bi-periodic
      self%ncells_2d = 9
      self%nverts_2d = 9
      self%nedges_2d = 18

      self%nlayers = 3
      self%ncells  = 27
      self%nverts  = 36
      self%nfaces  = 90
      self%nedges  = 99
    end if


    self%ncells_2d_with_ghost = self%ncells_2d                          &
                              + self%partition%get_num_cells_ghost()
    self%ncells_with_ghost    = self%ncells_2d_with_ghost * self%nlayers


    allocate( self%cell_next         ( self%nfaces_per_cell, self%ncells) )
    allocate( self%cell_lid_gid_map  ( self%ncells_2d) )
    allocate( self%vert_on_cell      ( self%nverts_per_cell, self%ncells) )
    allocate( self%face_on_cell      ( self%nfaces_per_cell, self%ncells_2d) )
    allocate( self%edge_on_cell      ( self%nedges_per_cell, self%ncells_2d) )
    allocate( self%vertex_coords     ( 3, self%nverts) )

    allocate( self%vert_cell_owner   ( nverts_per_2d_cell, self%ncells_2d ) )
    allocate( self%edge_cell_owner   ( nedges_per_2d_cell, self%ncells_2d ) )
    allocate( self%edge_ownership    ( nedges_per_2d_cell, self%ncells_2d ) )
    allocate( self%vertex_ownership  ( nverts_per_2d_cell, self%ncells_2d ) )

    allocate( self%eta               ( 0:self%nlayers ) )
    allocate( self%dz                ( self%nlayers   ) )

    ! Calculate vertical coordinates eta[0,1] and dz in a separate subroutine
    ! for the unit tests.
    ! Uniform vertical grid is used for pFunit tests (hard-wired).
    call set_vertical_coordinate(self, VGRID_UNIFORM)

    self%vert_cell_owner (:,:) = reshape( [ &
         9, 8, 5, 6, &  ! Cell 1
         8, 9, 6, 5, &  ! Cell 2
         9, 9, 6, 6, &  ! Cell 3
         6, 5, 8, 9, &  ! Cell 4
         5, 6, 9, 8, &  ! Cell 5
         6, 6, 9, 9, &  ! Cell 6
         9, 8, 8, 9, &  ! Cell 7
         8, 9, 9, 8, &  ! Cell 8
         9, 9, 9, 9  &  ! Cell 9
         ], shape(self%vert_cell_owner) )

    self%edge_cell_owner (:,:) = reshape( [ &
         3, 7, 2, 4, &  ! Cell 1
         2, 8, 3, 5, &  ! Cell 2
         3, 9, 3, 6, &  ! Cell 3
         6, 4, 5, 7, &  ! Cell 4
         5, 5, 6, 8, &  ! Cell 5
         6, 6, 6, 9, &  ! Cell 6
         9, 7, 8, 7, &  ! Cell 7
         8, 8, 9, 8, &  ! Cell 8
         9, 9, 9, 9  &  ! Cell 9
         ], shape(self%edge_cell_owner) )

    self%edge_ownership   (:,:) = 0
    self%vertex_ownership (:,:) = 0


    ! Global ids of local cells (given by index) in the
    ! partition, in this case there is only 1 partition
    ! so it is a one-to-one mapping, i.e. array element 1(lid)
    ! contains tha data 1 (gid).
    self%cell_lid_gid_map = [1,2,3,4,5,6,7,8,9]



    if (mesh_cfg == PLANE) then
      !=========================================================
      ! Assign 3D cell local ids on adjacent to given cell
      !
      ! Index ordering follows:
      ! 1) West
      ! 2) South
      ! 3) East
      ! 4) North
      ! 5) Bottom
      ! 6) Top
      !=========================================================
      ! Layer 1
      self%cell_next(:, 1) = [ 0,  0,  2,  4,  0, 10]
      self%cell_next(:, 2) = [ 1,  0,  3,  5,  0, 11]
      self%cell_next(:, 3) = [ 2,  0,  0,  6,  0, 12]
      self%cell_next(:, 4) = [ 0,  1,  5,  7,  0, 13]
      self%cell_next(:, 5) = [ 4,  2,  6,  8,  0, 14]
      self%cell_next(:, 6) = [ 5,  3,  0,  9,  0, 15]
      self%cell_next(:, 7) = [ 0,  4,  8,  0,  0, 16]
      self%cell_next(:, 8) = [ 7,  5,  9,  0,  0, 17]
      self%cell_next(:, 9) = [ 8,  6,  0,  0,  0, 18]

      ! Layer 2
      self%cell_next(:,10) = [ 0,  0, 11, 13,  1, 19]
      self%cell_next(:,11) = [10,  0, 12, 14,  2, 20]
      self%cell_next(:,12) = [11,  0,  0, 15,  3, 21]
      self%cell_next(:,13) = [ 0, 10, 14, 16,  4, 22]
      self%cell_next(:,14) = [13, 11, 15, 17,  5, 23]
      self%cell_next(:,15) = [14, 12,  0, 18,  6, 24]
      self%cell_next(:,16) = [ 0, 13, 17,  0,  7, 25]
      self%cell_next(:,17) = [16, 14, 18,  0,  8, 26]
      self%cell_next(:,18) = [17, 15,  0,  0,  9, 27]

      ! Layer 3
      self%cell_next(:,19) = [ 0,  0, 20, 22, 10, 28]
      self%cell_next(:,20) = [19,  0, 21, 23, 11, 29]
      self%cell_next(:,21) = [20,  0,  0, 24, 12, 30]
      self%cell_next(:,22) = [ 0, 19, 23, 25, 13, 31]
      self%cell_next(:,23) = [22, 20, 24, 26, 14, 32]
      self%cell_next(:,24) = [23, 21,  0, 27, 15, 33]
      self%cell_next(:,25) = [ 0, 22, 26,  0, 16, 34]
      self%cell_next(:,26) = [25, 23, 27,  0, 17, 35]
      self%cell_next(:,27) = [26, 24,  0,  0, 18, 36]

      ! Layer 4
      self%cell_next(:,28) = [ 0,  0, 29, 31, 19, 37]
      self%cell_next(:,29) = [28,  0, 30, 32, 20, 38]
      self%cell_next(:,30) = [29,  0,  0, 33, 21, 39]
      self%cell_next(:,31) = [ 0, 28, 32, 34, 22, 40]
      self%cell_next(:,32) = [31, 29, 33, 35, 23, 41]
      self%cell_next(:,33) = [32, 30,  0, 36, 24, 42]
      self%cell_next(:,34) = [ 0, 31, 35,  0, 25, 43]
      self%cell_next(:,35) = [34, 32, 36,  0, 26, 44]
      self%cell_next(:,36) = [35, 33,  0,  0, 27, 45]

      ! Layer 5
      self%cell_next(:,37) = [ 0,  0, 38, 40, 18,  0]
      self%cell_next(:,38) = [37,  0, 39, 41, 19,  0]
      self%cell_next(:,39) = [38,  0,  0, 42, 20,  0]
      self%cell_next(:,40) = [ 0, 37, 41, 43, 40,  0]
      self%cell_next(:,41) = [40, 38, 42, 44, 32,  0]
      self%cell_next(:,42) = [41, 39,  0, 45, 33,  0]
      self%cell_next(:,43) = [ 0, 40, 44,  0, 34,  0]
      self%cell_next(:,44) = [43, 41, 45,  0, 35,  0]
      self%cell_next(:,45) = [44, 42,  0,  0, 36,  0]

      !=========================================================
      ! Assign vertex local ids on cell corners
      !
      ! Index ordering follows:
      ! 1) South-West Bottom
      ! 2) South-East Bottom
      ! 3) North-East Bottom
      ! 4) North-West Bottom
      ! 5) South-West Top
      ! 6) South-East Top
      ! 7) North-East Top
      ! 8) North-West Top
      !=========================================================
      ! Layer 1
      self%vert_on_cell(:, 1) = [ 1,  2,  3,  4, 17, 18, 19, 20]
      self%vert_on_cell(:, 2) = [ 2,  5,  6,  3, 18, 21, 22, 19]
      self%vert_on_cell(:, 3) = [ 5,  7,  8,  6, 21, 23, 24, 22]
      self%vert_on_cell(:, 4) = [ 4,  3,  9, 10, 20, 19, 25, 26]
      self%vert_on_cell(:, 5) = [ 3,  6, 11,  9, 19, 22, 27, 25]
      self%vert_on_cell(:, 6) = [ 6,  8, 12, 11, 22, 24, 28, 27]
      self%vert_on_cell(:, 7) = [10,  9, 13, 14, 26, 25, 29, 30]
      self%vert_on_cell(:, 8) = [ 9, 11, 15, 13, 25, 27, 31, 29]
      self%vert_on_cell(:, 9) = [11, 12, 16, 15, 27, 28, 32, 31]

      ! Layer 2
      self%vert_on_cell(:,10) = [17, 18, 19, 20, 33, 34, 35, 36]
      self%vert_on_cell(:,11) = [18, 21, 22, 19, 34, 37, 38, 35]
      self%vert_on_cell(:,12) = [21, 23, 24, 22, 37, 39, 40, 38]
      self%vert_on_cell(:,13) = [20, 19, 25, 26, 36, 35, 41, 42]
      self%vert_on_cell(:,14) = [19, 22, 27, 25, 35, 38, 43, 41]
      self%vert_on_cell(:,15) = [22, 24, 28, 27, 38, 40, 44, 43]
      self%vert_on_cell(:,16) = [26, 25, 29, 30, 42, 41, 45, 46]
      self%vert_on_cell(:,17) = [25, 27, 31, 29, 41, 43, 47, 45]
      self%vert_on_cell(:,18) = [27, 28, 32, 31, 43, 44, 48, 47]

      ! Layer 3
      self%vert_on_cell(:,19) = [33, 34, 35, 36, 49, 50, 51, 52]
      self%vert_on_cell(:,20) = [34, 37, 38, 35, 50, 53, 54, 51]
      self%vert_on_cell(:,21) = [37, 39, 40, 38, 53, 55, 56, 54]
      self%vert_on_cell(:,22) = [36, 35, 41, 42, 52, 51, 57, 58]
      self%vert_on_cell(:,23) = [35, 38, 43, 41, 51, 54, 59, 57]
      self%vert_on_cell(:,24) = [38, 40, 44, 43, 54, 56, 60, 59]
      self%vert_on_cell(:,25) = [42, 41, 45, 46, 58, 57, 61, 62]
      self%vert_on_cell(:,26) = [41, 43, 47, 45, 57, 59, 63, 61]
      self%vert_on_cell(:,27) = [43, 44, 48, 47, 59, 60, 64, 63]

      ! Layer 4
      self%vert_on_cell(:,28) = [49, 50, 51, 52, 65, 66, 67, 68]
      self%vert_on_cell(:,29) = [50, 53, 54, 51, 66, 69, 70, 67]
      self%vert_on_cell(:,30) = [53, 55, 56, 54, 69, 71, 72, 70]
      self%vert_on_cell(:,31) = [52, 51, 57, 58, 68, 67, 73, 74]
      self%vert_on_cell(:,32) = [51, 54, 59, 57, 67, 70, 75, 73]
      self%vert_on_cell(:,33) = [54, 56, 60, 59, 70, 71, 76, 75]
      self%vert_on_cell(:,34) = [58, 57, 61, 62, 74, 73, 77, 78]
      self%vert_on_cell(:,35) = [57, 59, 63, 61, 73, 75, 79, 77]
      self%vert_on_cell(:,36) = [59, 60, 64, 63, 75, 76, 80, 79]

      ! Layer 5
      self%vert_on_cell(:,37) = [65, 66, 67, 68, 81, 82, 83, 84]
      self%vert_on_cell(:,38) = [66, 69, 70, 67, 82, 85, 86, 83]
      self%vert_on_cell(:,39) = [69, 71, 72, 70, 85, 87, 88, 86]
      self%vert_on_cell(:,40) = [68, 67, 73, 74, 84, 83, 89, 90]
      self%vert_on_cell(:,41) = [67, 70, 75, 73, 83, 86, 91, 89]
      self%vert_on_cell(:,42) = [70, 71, 76, 75, 86, 87, 82, 91]
      self%vert_on_cell(:,43) = [74, 73, 77, 78, 90, 89, 93, 94]
      self%vert_on_cell(:,44) = [73, 75, 79, 77, 89, 91, 95, 93]
      self%vert_on_cell(:,45) = [75, 76, 80, 79, 91, 92, 96, 95]

      !=========================================================
      ! Assign edge local ids on cell edges
      !
      ! Index ordering follows:
      ! 1)  West  Bottom
      ! 2)  South Bottom
      ! 3)  East  Bottom
      ! 4)  North Bottom
      ! 5)  South-West Middle
      ! 6)  South-East Middle
      ! 7)  North-East Middle
      ! 8)  North-West Middle
      ! 9)  West  Top
      ! 10) South Top
      ! 11) East  Top
      ! 12) North Top
      !=========================================================
      ! Layer 1
      self%edge_on_cell(:,1) = [ 1,  3,  5,  7,  9, 10, 11, 12,  2,  4,  6,  8]
      self%edge_on_cell(:,2) = [ 5, 13, 15, 17, 10, 19, 20, 11,  6, 14, 16, 18]
      self%edge_on_cell(:,3) = [15, 21, 23, 25, 19, 27, 28, 20, 16, 22, 24, 26]
      self%edge_on_cell(:,4) = [29,  7, 31, 33, 12, 11, 35, 36, 30,  8, 32, 34]
      self%edge_on_cell(:,5) = [31, 17, 37, 39, 11, 20, 41, 35, 32, 18, 38, 40]
      self%edge_on_cell(:,6) = [37, 25, 42, 44, 20, 28, 46, 41, 38, 26, 43, 45]
      self%edge_on_cell(:,7) = [47, 33, 49, 51, 36, 35, 53, 54, 48, 34, 50, 52]
      self%edge_on_cell(:,8) = [49, 39, 55, 57, 35, 41, 59, 53, 50, 40, 56, 58]
      self%edge_on_cell(:,9) = [55, 44, 60, 62, 41, 46, 64, 59, 56, 45, 61, 63]

      !=========================================================
      ! Assign face local ids on cell sides
      !
      ! Index ordering follows:
      ! 1)  West
      ! 2)  South
      ! 3)  East
      ! 4)  North
      ! 5)  Bottom
      ! 6)  Top
      !=========================================================
      ! Layer 1
      self%face_on_cell(:,1) = [ 1,  2,  3,  4,  5,  6]
      self%face_on_cell(:,2) = [ 3,  7,  8,  9, 10, 11]
      self%face_on_cell(:,3) = [ 8, 12, 13, 14, 15, 16]
      self%face_on_cell(:,4) = [17,  4, 18, 19, 20, 21]
      self%face_on_cell(:,5) = [18,  9, 22, 23, 24, 25]
      self%face_on_cell(:,6) = [22, 14, 26, 27, 28, 29]
      self%face_on_cell(:,7) = [30, 19, 31, 32, 33, 34]
      self%face_on_cell(:,8) = [31, 23, 35, 36, 37, 38]
      self%face_on_cell(:,9) = [35, 27, 39, 40, 41, 42]

      !=========================================================
      ! Assign [x,y,z] vertex coords in (m), with [0,0,0] at centre
      ! of planet (radius=30000.0).
      ! Level 0
      self%vertex_coords (:, 1) = [  5195.345687_r_def,  11352.037430_r_def, &
                                    -27278.922805_r_def ]
      self%vertex_coords (:, 2) = [ -6745.352861_r_def,  10505.264651_r_def, &
                                    -27278.922805_r_def ]
      self%vertex_coords (:, 3) = [  8757.797452_r_def, -13639.461402_r_def, &
                                    -25244.129544_r_def ]
      self%vertex_coords (:, 4) = [ -6745.352861_r_def, -14738.864893_r_def, &
                                    -25244.129544_r_def ]
      self%vertex_coords (:, 5) = [ -6745.352861_r_def, -10505.264651_r_def, &
                                    -27278.922805_r_def ]
      self%vertex_coords (:, 6) = [  8757.797452_r_def,  13639.461402_r_def, &
                                    -25244.129544_r_def ]
      self%vertex_coords (:, 7) = [  5195.345687_r_def, -11352.037430_r_def, &
                                    -27278.922805_r_def ]
      self%vertex_coords (:, 8) = [ -6745.352861_r_def,  14738.864893_r_def, &
                                    -25244.129544_r_def ]
      self%vertex_coords (:, 9) = [  8757.797452_r_def, -13639.461402_r_def, &
                                     25244.129544_r_def ]
      self%vertex_coords (:,10) = [ -6745.352861_r_def, -14738.864893_r_def, &
                                     25244.129544_r_def ]
      self%vertex_coords (:,11) = [  8757.797452_r_def,  13639.461402_r_def, &
                                     25244.129544_r_def ]
      self%vertex_coords (:,12) = [ -6745.352861_r_def,  14738.864893_r_def, &
                                     25244.129544_r_def ]
      self%vertex_coords (:,13) = [ -6745.352861_r_def,  10505.264651_r_def, &
                                     27278.922805_r_def ]
      self%vertex_coords (:,14) = [  5195.345687_r_def,  11352.037430_r_def, &
                                     27278.922805_r_def ]
      self%vertex_coords (:,15) = [ -6745.352861_r_def, -10505.264651_r_def, &
                                     27278.922805_r_def ]
      self%vertex_coords (:,16) = [  5195.345687_r_def, -11352.037430_r_def, &
                                     27278.922805_r_def ]

      ! Level 1
      self%vertex_coords (:,17) = [  5541.702066_r_def,  12108.839925_r_def, &
                                    -29097.517658_r_def ]
      self%vertex_coords (:,18) = [ -7195.043052_r_def,  11205.615628_r_def, &
                                    -29097.517658_r_def ]
      self%vertex_coords (:,19) = [  9341.650615_r_def, -14548.758829_r_def, &
                                    -26927.071514_r_def ]
      self%vertex_coords (:,20) = [ -7195.043052_r_def, -15721.455886_r_def, &
                                    -26927.071514_r_def ]
      self%vertex_coords (:,21) = [ -7195.043052_r_def, -11205.615628_r_def, &
                                    -29097.517658_r_def ]
      self%vertex_coords (:,22) = [  9341.650615_r_def,  14548.758829_r_def, &
                                    -26927.071514_r_def ]
      self%vertex_coords (:,23) = [  5541.702066_r_def, -12108.839925_r_def, &
                                    -29097.517658_r_def ]
      self%vertex_coords (:,24) = [ -7195.043052_r_def,  15721.455886_r_def, &
                                    -26927.071514_r_def ]
      self%vertex_coords (:,25) = [  9341.650615_r_def, -14548.758829_r_def, &
                                     26927.071514_r_def ]
      self%vertex_coords (:,26) = [ -7195.043052_r_def, -15721.455886_r_def, &
                                     26927.071514_r_def ]
      self%vertex_coords (:,27) = [  9341.650615_r_def,  14548.758829_r_def, &
                                     26927.071514_r_def ]
      self%vertex_coords (:,28) = [ -7195.043052_r_def,  15721.455886_r_def, &
                                     26927.071514_r_def ]
      self%vertex_coords (:,29) = [ -7195.043052_r_def,  11205.615628_r_def, &
                                     29097.517658_r_def ]
      self%vertex_coords (:,30) = [  5541.702066_r_def,  12108.839925_r_def, &
                                     29097.517658_r_def ]
      self%vertex_coords (:,31) = [ -7195.043052_r_def, -11205.615628_r_def, &
                                     29097.517658_r_def ]
      self%vertex_coords (:,32) = [  5541.702066_r_def, -12108.839925_r_def, &
                                     29097.517658_r_def ]

      ! Level 2
      self%vertex_coords (:,33) = [  5888.058445_r_def,  12865.642420_r_def, &
                                    -30916.112512_r_def ]
      self%vertex_coords (:,34) = [ -7644.733242_r_def,  11905.966605_r_def, &
                                    -30916.112512_r_def ]
      self%vertex_coords (:,35) = [  9925.503779_r_def, -15458.056256_r_def, &
                                    -28610.013483_r_def ]
      self%vertex_coords (:,36) = [ -7644.733242_r_def, -16704.046879_r_def, &
                                    -28610.013483_r_def ]
      self%vertex_coords (:,37) = [ -7644.733242_r_def, -11905.966605_r_def, &
                                    -30916.112512_r_def ]
      self%vertex_coords (:,38) = [  9925.503779_r_def,  15458.056256_r_def, &
                                    -28610.013483_r_def ]
      self%vertex_coords (:,39) = [  5888.058445_r_def, -12865.642420_r_def, &
                                    -30916.112512_r_def ]
      self%vertex_coords (:,40) = [ -7644.733242_r_def,  16704.046879_r_def, &
                                    -28610.013483_r_def ]
      self%vertex_coords (:,41) = [  9925.503779_r_def, -15458.056256_r_def, &
                                     28610.013483_r_def ]
      self%vertex_coords (:,42) = [ -7644.733242_r_def, -16704.046879_r_def, &
                                     28610.013483_r_def ]
      self%vertex_coords (:,43) = [  9925.503779_r_def,  15458.056256_r_def, &
                                     28610.013483_r_def ]
      self%vertex_coords (:,44) = [ -7644.733242_r_def,  16704.046879_r_def, &
                                     28610.013483_r_def ]
      self%vertex_coords (:,45) = [ -7644.733242_r_def,  11905.966605_r_def, &
                                     30916.112512_r_def ]
      self%vertex_coords (:,46) = [  5888.058445_r_def,  12865.642420_r_def, &
                                     30916.112512_r_def ]
      self%vertex_coords (:,47) = [ -7644.733242_r_def, -11905.966605_r_def, &
                                     30916.112512_r_def ]
      self%vertex_coords (:,48) = [  5888.058445_r_def, -12865.642420_r_def, &
                                     30916.112512_r_def ]

      ! Level 3
      self%vertex_coords (:,49) = [  6234.414824_r_def,  13622.444916_r_def, &
                                    -32734.707366_r_def ]
      self%vertex_coords (:,50) = [ -8094.423433_r_def,  12606.317581_r_def, &
                                    -32734.707366_r_def ]
      self%vertex_coords (:,51) = [ 10509.356942_r_def, -16367.353683_r_def, &
                                    -30292.955453_r_def ]
      self%vertex_coords (:,52) = [ -8094.423433_r_def, -17686.637872_r_def, &
                                    -30292.955453_r_def ]
      self%vertex_coords (:,53) = [ -8094.423433_r_def, -12606.317581_r_def, &
                                    -32734.707366_r_def ]
      self%vertex_coords (:,54) = [ 10509.356942_r_def,  16367.353683_r_def, &
                                    -30292.955453_r_def ]
      self%vertex_coords (:,55) = [  6234.414824_r_def, -13622.444916_r_def, &
                                    -32734.707366_r_def ]
      self%vertex_coords (:,56) = [ -8094.423433_r_def,  17686.637872_r_def, &
                                    -30292.955453_r_def ]
      self%vertex_coords (:,57) = [ 10509.356942_r_def, -16367.353683_r_def, &
                                     30292.955453_r_def ]
      self%vertex_coords (:,58) = [ -8094.423433_r_def, -17686.637872_r_def, &
                                     30292.955453_r_def ]
      self%vertex_coords (:,59) = [ 10509.356942_r_def,  16367.353683_r_def, &
                                     30292.955453_r_def ]
      self%vertex_coords (:,60) = [ -8094.423433_r_def,  17686.637872_r_def, &
                                     30292.955453_r_def ]
      self%vertex_coords (:,61) = [ -8094.423433_r_def,  12606.317581_r_def, &
                                     32734.707366_r_def ]
      self%vertex_coords (:,62) = [  6234.414824_r_def,  13622.444916_r_def, &
                                     32734.707366_r_def ]
      self%vertex_coords (:,63) = [ -8094.423433_r_def, -12606.317581_r_def, &
                                     32734.707366_r_def ]
      self%vertex_coords (:,64) = [  6234.414824_r_def, -13622.444916_r_def, &
                                     32734.707366_r_def ]

      ! Level 4
      self%vertex_coords (:,65) = [  6580.771204_r_def,  14379.247411_r_def, &
                                    -34553.302219_r_def ]
      self%vertex_coords (:,66) = [ -8544.113624_r_def,  13306.668558_r_def, &
                                    -34553.302219_r_def ]
      self%vertex_coords (:,67) = [ 11093.210106_r_def, -17276.651110_r_def, &
                                    -31975.897423_r_def ]
      self%vertex_coords (:,68) = [ -8544.113624_r_def, -18669.228864_r_def, &
                                    -31975.897423_r_def ]
      self%vertex_coords (:,69) = [ -8544.113624_r_def, -13306.668558_r_def, &
                                    -34553.302219_r_def ]
      self%vertex_coords (:,70) = [ 11093.210106_r_def,  17276.651110_r_def, &
                                    -31975.897423_r_def ]
      self%vertex_coords (:,71) = [  6580.771204_r_def, -14379.247411_r_def, &
                                    -34553.302219_r_def ]
      self%vertex_coords (:,72) = [ -8544.113624_r_def,  18669.228864_r_def, &
                                    -31975.897423_r_def ]
      self%vertex_coords (:,73) = [ 11093.210106_r_def, -17276.651110_r_def, &
                                     31975.897423_r_def ]
      self%vertex_coords (:,74) = [ -8544.113624_r_def, -18669.228864_r_def, &
                                     31975.897423_r_def ]
      self%vertex_coords (:,75) = [ 11093.210106_r_def,  17276.651110_r_def, &
                                     31975.897423_r_def ]
      self%vertex_coords (:,76) = [ -8544.113624_r_def,  18669.228864_r_def, &
                                     31975.897423_r_def ]
      self%vertex_coords (:,77) = [ -8544.113624_r_def,  13306.668558_r_def, &
                                     34553.302219_r_def ]
      self%vertex_coords (:,78) = [  6580.771204_r_def,  14379.247411_r_def, &
                                     34553.302219_r_def ]
      self%vertex_coords (:,79) = [ -8544.113624_r_def, -13306.668558_r_def, &
                                     34553.302219_r_def ]
      self%vertex_coords (:,80) = [  6580.771204_r_def, -14379.247411_r_def, &
                                     34553.302219_r_def ]

      ! Level 5
      self%vertex_coords (:,81) = [  6927.127583_r_def,  15136.049906_r_def, &
                                    -36371.897073_r_def ]
      self%vertex_coords (:,82) = [ -8993.803815_r_def,  14007.019535_r_def, &
                                    -36371.897073_r_def ]
      self%vertex_coords (:,83) = [ 11677.063269_r_def, -18185.948537_r_def, &
                                    -33658.839392_r_def ]
      self%vertex_coords (:,84) = [ -8993.803815_r_def, -19651.819857_r_def, &
                                    -33658.839392_r_def ]
      self%vertex_coords (:,85) = [ -8993.803815_r_def, -14007.019535_r_def, &
                                    -36371.897073_r_def ]
      self%vertex_coords (:,86) = [ 11677.063269_r_def,  18185.948537_r_def, &
                                    -33658.839392_r_def ]
      self%vertex_coords (:,87) = [  6927.127583_r_def, -15136.049906_r_def, &
                                    -36371.897073_r_def ]
      self%vertex_coords (:,88) = [ -8993.803815_r_def,  19651.819857_r_def, &
                                    -33658.839392_r_def ]
      self%vertex_coords (:,89) = [ 11677.063269_r_def, -18185.948536_r_def, &
                                     33658.839392_r_def ]
      self%vertex_coords (:,90) = [ -8993.803815_r_def, -19651.819857_r_def, &
                                     33658.839392_r_def ]
      self%vertex_coords (:,91) = [ 11677.063269_r_def,  18185.948536_r_def, &
                                     33658.839392_r_def ]
      self%vertex_coords (:,92) = [ -8993.803815_r_def,  19651.819857_r_def, &
                                     33658.839392_r_def ]
      self%vertex_coords (:,93) = [ -8993.803815_r_def,  14007.019535_r_def, &
                                     36371.897073_r_def ]
      self%vertex_coords (:,94) = [  6927.127583_r_def,  15136.049906_r_def, &
                                     36371.897073_r_def ]
      self%vertex_coords (:,95) = [ -8993.803815_r_def, -14007.019535_r_def, &
                                     36371.897073_r_def ]
      self%vertex_coords (:,96) = [  6927.127582_r_def, -15136.049906_r_def, &
                                     36371.897073_r_def ]


      ! Domain limits
      self % domain_size%minimum%x =  0.0_r_def
      self % domain_size%maximum%x =  2.0_r_def*PI
      self % domain_size%minimum%y = -0.5_r_def*PI
      self % domain_size%maximum%y =  0.5_r_def*PI
      self % domain_size%minimum%z =  0.0_r_def
      self % domain_size%maximum%z =  self % domain_top

    else if (mesh_cfg == PLANE_BI_PERIODIC) then
      !=========================================================
      ! Assign 3D cell local ids on adjacent to given cell
      !
      ! Index ordering follows:
      ! 1) West
      ! 2) South
      ! 3) East
      ! 4) North
      ! 5) Bottom
      ! 6) Top
      !=========================================================
      ! Layer 1
      self%cell_next(:, 1) = [ 3,  7,  2,  4,  0, 10]
      self%cell_next(:, 2) = [ 1,  8,  3,  5,  0, 11]
      self%cell_next(:, 3) = [ 2,  9,  1,  6,  0, 12]
      self%cell_next(:, 4) = [ 6,  1,  5,  7,  0, 13]
      self%cell_next(:, 5) = [ 4,  2,  6,  8,  0, 14]
      self%cell_next(:, 6) = [ 5,  3,  4,  9,  0, 15]
      self%cell_next(:, 7) = [ 9,  4,  8,  1,  0, 16]
      self%cell_next(:, 8) = [ 7,  5,  9,  2,  0, 17]
      self%cell_next(:, 9) = [ 8,  6,  7,  3,  0, 18]

      ! Layer 2
      self%cell_next(:,10) = [12, 16, 11, 13,  1, 19]
      self%cell_next(:,11) = [10, 17, 12, 14,  2, 20]
      self%cell_next(:,12) = [11, 18, 10, 15,  3, 21]
      self%cell_next(:,13) = [15, 10, 14, 16,  4, 22]
      self%cell_next(:,14) = [13, 11, 15, 17,  5, 23]
      self%cell_next(:,15) = [14, 12, 13, 18,  6, 24]
      self%cell_next(:,16) = [18, 13, 17, 10,  7, 25]
      self%cell_next(:,17) = [16, 14, 18, 11,  8, 26]
      self%cell_next(:,18) = [17, 15, 16, 12,  9, 27]

      ! Layer 3
      self%cell_next(:,19) = [21, 25, 20, 22, 10, 28]
      self%cell_next(:,20) = [19, 26, 21, 23, 11, 29]
      self%cell_next(:,21) = [20, 27, 19, 24, 12, 30]
      self%cell_next(:,22) = [24, 19, 23, 25, 13, 31]
      self%cell_next(:,23) = [22, 20, 24, 26, 14, 32]
      self%cell_next(:,24) = [23, 21, 22, 27, 15, 33]
      self%cell_next(:,25) = [27, 22, 26, 19, 16, 34]
      self%cell_next(:,26) = [25, 23, 27, 20, 17, 35]
      self%cell_next(:,27) = [26, 24, 25, 21, 18, 36]



      !=========================================================
      ! Assign vertex local ids on cell corners
      !
      ! Index ordering follows:
      ! 1) South-West Bottom
      ! 2) South-East Bottom
      ! 3) North-East Bottom
      ! 4) North-West Bottom
      ! 5) South-West Top
      ! 6) South-East Top
      ! 7) North-East Top
      ! 8) North-West Top
      !=========================================================
      ! Layer 1
      self%vert_on_cell(:, 1) = [ 1,  2,  3,  4, 10, 11, 12, 13]
      self%vert_on_cell(:, 2) = [ 2,  5,  6,  3, 11, 14, 15, 12]
      self%vert_on_cell(:, 3) = [ 5,  1,  4,  6, 14, 10, 13, 15]
      self%vert_on_cell(:, 4) = [ 4,  3,  7,  8, 13, 12, 16, 17]
      self%vert_on_cell(:, 5) = [ 3,  6,  9,  7, 12, 15, 18, 16]
      self%vert_on_cell(:, 6) = [ 6,  4,  8,  9, 15, 13, 17, 18]
      self%vert_on_cell(:, 7) = [ 8,  7,  2,  1, 17, 16, 11, 10]
      self%vert_on_cell(:, 8) = [ 7,  9,  5,  2, 16, 18, 14, 11]
      self%vert_on_cell(:, 9) = [ 9,  8,  1,  5, 18, 17, 10, 14]


      ! Layer 2
      self%vert_on_cell(:,10) = [10, 11, 12, 13, 19, 20, 21, 22]
      self%vert_on_cell(:,11) = [11, 14, 15, 12, 20, 23, 24, 21]
      self%vert_on_cell(:,12) = [14, 10, 13, 15, 23, 19, 22, 24]
      self%vert_on_cell(:,13) = [13, 12, 16, 17, 22, 21, 25, 26]
      self%vert_on_cell(:,14) = [12, 15, 18, 16, 21, 24, 27, 25]
      self%vert_on_cell(:,15) = [15, 13, 17, 18, 24, 22, 26, 27]
      self%vert_on_cell(:,16) = [17, 16, 11, 10, 26, 25, 20, 19]
      self%vert_on_cell(:,17) = [16, 18, 14, 11, 25, 27, 23, 20]
      self%vert_on_cell(:,18) = [18, 17, 10, 14, 27, 26, 19, 23]

      ! Layer 3
      self%vert_on_cell(:,19) = [19, 20, 21, 22, 28, 29, 30, 31]
      self%vert_on_cell(:,20) = [20, 23, 24, 21, 29, 32, 33, 30]
      self%vert_on_cell(:,21) = [23, 19, 22, 24, 32, 28, 31, 33]
      self%vert_on_cell(:,22) = [22, 21, 25, 26, 31, 30, 34, 35]
      self%vert_on_cell(:,23) = [21, 24, 27, 25, 30, 33, 36, 34]
      self%vert_on_cell(:,24) = [24, 22, 26, 27, 33, 31, 35, 36]
      self%vert_on_cell(:,25) = [26, 25, 20, 19, 35, 34, 29, 28]
      self%vert_on_cell(:,26) = [25, 27, 23, 20, 34, 36, 32, 29]
      self%vert_on_cell(:,27) = [27, 26, 19, 23, 36, 35, 28, 32]

      !=========================================================
      ! Assign edge local ids on cell edges
      !
      ! Index ordering follows:
      ! 1)  West  Bottom
      ! 2)  South Bottom
      ! 3)  East  Bottom
      ! 4)  North Bottom
      ! 5)  South-West Middle
      ! 6)  South-East Middle
      ! 7)  North-East Middle
      ! 8)  North-West Middle
      ! 9)  West  Top
      ! 10) South Top
      ! 11) East  Top
      ! 12) North Top
      !=========================================================
      ! Layer 1
      self%edge_on_cell(:,1) = [ 1,  3,  5,  7,  9, 10, 11, 12,  2,  4,  6,  8]
      self%edge_on_cell(:,2) = [ 5, 13, 15, 17, 10, 19, 20, 11,  6, 14, 16, 18]
      self%edge_on_cell(:,3) = [15, 21,  1, 23, 19,  9, 12, 20, 16, 22,  2, 24]
      self%edge_on_cell(:,4) = [25,  7, 27, 29, 12, 11, 31, 32, 26,  8, 28, 30]
      self%edge_on_cell(:,5) = [27, 17, 33, 35, 11, 20, 37, 31, 28, 18, 34, 36]
      self%edge_on_cell(:,6) = [33, 23, 25, 38, 20, 12, 32, 37, 34, 24, 26, 39]
      self%edge_on_cell(:,7) = [40, 29, 42,  3, 32, 31, 10,  9, 41, 30, 43,  4]
      self%edge_on_cell(:,8) = [42, 35, 44, 13, 31, 37, 19, 10, 43, 36, 45, 14]
      self%edge_on_cell(:,9) = [44, 38, 40, 21, 37, 32,  9, 19, 45, 39, 41, 22]


      !=========================================================
      ! Assign face local ids on cell sides
      !
      ! Index ordering follows:
      ! 1)  West
      ! 2)  South
      ! 3)  East
      ! 4)  North
      ! 5)  Bottom
      ! 6)  Top
      !=========================================================
      self%face_on_cell(:,1) = [ 1,  2,  3,  4,  5,  6]
      self%face_on_cell(:,2) = [ 3,  7,  8,  9, 10, 11]
      self%face_on_cell(:,3) = [ 8, 12,  1, 13, 14, 15]
      self%face_on_cell(:,4) = [16,  4, 17, 18, 19, 20]
      self%face_on_cell(:,5) = [17,  9, 21, 22, 23, 24]
      self%face_on_cell(:,6) = [21, 13, 16, 25, 26, 27]
      self%face_on_cell(:,7) = [28, 18, 29,  2, 30, 31]
      self%face_on_cell(:,8) = [29, 22, 32,  7, 33, 34]
      self%face_on_cell(:,9) = [32, 25, 28, 12, 35, 36]

    end if

  end function mesh_constructor_unit_test_data


end module mesh_mod
