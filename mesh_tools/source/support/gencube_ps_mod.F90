!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief   Module to define the gencube_ps_type, a subclass of the
!>          ugrid_generator_type which generates a cubed-sphere mesh in a
!>          format suitable for storage as a ugrid file.
!> @details Type implements the ugrid_generator_type interface to
!>          construct a cubed-sphere mesh.  All required connectivity is
!>          calculated and made available to the ugrid writer.
!>
!>      +---+
!>      | 5 |
!>  +---+---+---+---+
!>  | 1 | 2 | 3 | 4 |
!>  +---+---+---+---+
!>      | 6 |
!>      +---+
!>
!>  The equator passes through the centres of panels 1-4, while the zero degrees
!>  longitude passes through the centre of cell 1 and the axis of the north pole
!>  through panels 5 and 6.
!-------------------------------------------------------------------------------
module gencube_ps_mod
!-------------------------------------------------------------------------------

  use calc_global_cell_map_mod,       only: calc_global_cell_map
  use constants_mod,                  only: r_def, i_def, str_def, l_def,     &
                                            str_long, str_longlong,           &
                                            PI, radians_to_degrees,           &
                                            degrees_to_radians, rmdi, imdi,   &
                                            emdi
  use coord_transform_mod,            only: ll2xyz, xyz2ll
  use global_mesh_map_collection_mod, only: global_mesh_map_collection_type
  use log_mod,                        only: log_event, log_scratch_space, &
                                            LOG_LEVEL_ERROR, LOG_LEVEL_INFO
  use reference_element_mod,          only: W, S, E, N, SWB, SEB, NWB, NEB
  use ugrid_generator_mod,            only: ugrid_generator_type


  use mesh_config_mod,                only: coord_sys_ll, coord_sys_xyz, &
                                            key_from_geometry,           &
                                            key_from_topology,           &
                                            key_from_coord_sys,          &
                                            geometry_spherical,          &
                                            topology_periodic

  use rotation_mod,                   only: rotate_mesh_coords, &
                                            TRUE_NORTH_POLE_LL, &
                                            TRUE_NULL_ISLAND_LL

  use panel_decomposition_mod, only: panel_decomposition_type
  implicit none

  private

  public :: set_partition_parameters

!-----------------------------------------------------------------------------
  ! Mesh Vertex directions: local aliases for reference_element_mod values
  integer(i_def), parameter :: NW = NWB
  integer(i_def), parameter :: NE = NEB
  integer(i_def), parameter :: SE = SEB
  integer(i_def), parameter :: SW = SWB

  ! For a cubesphere these panels are the 6 faces of the domain cube
  integer(i_def), parameter :: NPANELS = 6
  integer(i_def), parameter :: PANEL_ROTATIONS(NPANELS) = (/ 0, 0, 1, 1, -1, 0 /)

  ! Prefix for error messages
  character(*),       parameter :: PREFIX = "[Cubed-Sphere Mesh] "

  ! flag to print out mesh data for debugging purposes
  logical(l_def),     parameter :: DEBUG = .false.

  ! Set to -9999 to be used for fill value so meshes
  ! are more Cf-compliant.
  integer(i_def), parameter :: VOID_ID = -9999
!-------------------------------------------------------------------------------

  type, extends(ugrid_generator_type), public :: gencube_ps_type

    private

    character(str_def) :: mesh_name
    integer(i_def)     :: geometry  = geometry_spherical
    integer(i_def)     :: topology  = topology_periodic
    integer(i_def)     :: coord_sys = emdi
    character(str_def) :: coord_units_x
    character(str_def) :: coord_units_y
    integer(i_def)     :: edge_cells
    real(r_def)        :: domain_extents(2,4) =rmdi
    integer(i_def)     :: npanels   = NPANELS
    real(r_def)        :: north_pole(2)
    real(r_def)        :: null_island(2)

    character(str_longlong) :: constructor_inputs

    integer(i_def)     :: nsmooth
    real(r_def)        :: equatorial_latitude = 0.0_r_def
    logical(l_def)     :: rotate_mesh    = .false.

    ! Connectivity and coordinates
    integer(i_def), allocatable :: cell_next(:,:)
    integer(i_def), allocatable :: verts_on_cell(:,:)
    integer(i_def), allocatable :: edges_on_cell(:,:)
    integer(i_def), allocatable :: verts_on_edge(:,:)
    real(r_def),    allocatable :: vert_coords(:,:)
    real(r_def),    allocatable :: cell_coords(:,:)

    ! Intergrid maps
    integer(i_def) :: nmaps
    character(str_def), allocatable :: target_mesh_names(:)
    integer(i_def),     allocatable :: target_edge_cells(:)
    type(global_mesh_map_collection_type), allocatable :: global_mesh_maps

    integer(i_def) :: max_num_faces_per_node

  contains
    procedure :: generate
    procedure :: get_number_of_panels
    procedure :: get_metadata
    procedure :: get_dimensions
    procedure :: get_coordinates
    procedure :: get_connectivity
    procedure :: get_global_mesh_maps
    procedure :: write_mesh

    procedure :: clear
    final     :: gen_cubedsphere_final

  end type gencube_ps_type

!-------------------------------------------------------------------------------
  interface gencube_ps_type
    module procedure gencube_ps_constructor
  end interface gencube_ps_type
!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> @brief   Constructor for gencube_ps_type.
!> @details Accepts mesh dimension for initialisation and validation.
!>
!> @param[in] mesh_name          Name of this mesh topology
!> @param[in] edge_cells         Number of cells per panel edge of the cubed-sphere.
!>                               Each panel will contain edge_cells*edge_cells faces.
!> @param[in] nsmooth            Number of smoothing passes to be performed on mesh nodes.
!>                               Each panel will contain edge_cells*edge_cells faces.
!> @param[in] coord_sys          Coordinate system to position nodes.
!> @param[in, optional] rotate_mesh
!>                               Logical to indicate rotation of the resulting mesh
!> @param[in, optional] target_north_pole
!>                               If rotating the cubed-sphere, then these are the target
!>                               coordinates [longitude, latitude] (degrees)
!>                               to move the reference north pole.
!> @param[in, optional] target_null_island
!>                               If rotating the cubed-sphere, then these are the target
!>                               coordinates [longitude, latitude] (degrees)
!>                               to move the reference null island.
!> @param[in, optional] target_mesh_names
!>                               Names of meshes to map to.
!> @param[in, optional] target_edge_cells
!>                               Number of cells per panel edge of the meshes to map to.
!> @param[in, optional]  equatorial_latitude
!>                               Latitude in degrees of the equator of the mesh,
!>                               implying a stretching towards one of the poles.
!>                               This is currently implemented through the
!>                               Schmidt stretching transform.
!>
!> @return    self               Instance of gencube_ps_type
!-------------------------------------------------------------------------------
 function gencube_ps_constructor( mesh_name, edge_cells, nsmooth, coord_sys,   &
                                  rotate_mesh, target_north_pole,              &
                                  target_null_island,                          &
                                  target_mesh_names, target_edge_cells,        &
                                  equatorial_latitude )                        &
                                  result( self )

  implicit none

  character(str_def), intent(in) :: mesh_name
  integer(i_def),     intent(in) :: edge_cells
  integer(i_def),     intent(in) :: nsmooth
  integer(i_def),     intent(in) :: coord_sys

  logical,            optional, intent(in) :: rotate_mesh
  real(r_def),        optional, intent(in) :: target_north_pole(2)
  real(r_def),        optional, intent(in) :: target_null_island(2)
  real(r_def),        optional, intent(in) :: equatorial_latitude

  character(str_def), optional, intent(in) :: target_mesh_names(:)
  integer(i_def),     optional, intent(in) :: target_edge_cells(:)

  type( gencube_ps_type ) :: self

  integer(i_def) :: i
  integer(i_def) :: remainder = 0

  character(str_long) :: target_mesh_names_str
  character(str_long) :: target_edge_cells_str

  self%mesh_name  = trim(mesh_name)
  self%edge_cells = edge_cells
  self%nsmooth    = nsmooth
  self%nmaps      = 0
  self%coord_sys  = coord_sys

  self%domain_extents(:,1) = [ -180.0_r_def, -90.0_r_def ]
  self%domain_extents(:,2) = [  180.0_r_def, -90.0_r_def ]
  self%domain_extents(:,3) = [  180.0_r_def,  90.0_r_def ]
  self%domain_extents(:,4) = [ -180.0_r_def,  90.0_r_def ]

  ! There are a maximum of 4 faces around a node in this type of mesh
  self%max_num_faces_per_node = 4

  ! Construct string to report construction arguments
  ! This string should be updated if new arguments are
  ! applied.
  write(self%constructor_inputs,'(2(A,I0))')     &
      'edge_cells=',    self%edge_cells,  ';' // &
      'smooth_passes=', self%nsmooth

  ! Only attempt to stretch if equatorial lattiude given
  if ( present(equatorial_latitude) ) then

    if ( (equatorial_latitude <= -90.0_r_def &
          .or. equatorial_latitude >= 90.0_r_def) ) then
      call log_event( 'Invalid equatorial latitude. Must be between -90 and 90', &
                      LOG_LEVEL_ERROR )
    end if

    self%equatorial_latitude = degrees_to_radians * equatorial_latitude

  else
    ! Default value is zero
    self%equatorial_latitude = 0.0_r_def
  end if

  self%rotate_mesh = .false.
  if ( present(rotate_mesh) ) self%rotate_mesh = rotate_mesh

  if ( self%rotate_mesh )then
    ! Inputs are in degrees, so convert
    ! and store them as radians.
    self%north_pole  = degrees_to_radians * target_north_pole
    self%null_island = degrees_to_radians * target_null_island
  else
    ! Default value is also given in degrees so
    ! convert to radians
    self%north_pole  = degrees_to_radians * TRUE_NORTH_POLE_LL
    self%null_island = degrees_to_radians * TRUE_NULL_ISLAND_LL
  end if

  ! Constructor inputs for any target mesh maps
  if ( present(target_edge_cells) .and. &
       present(target_mesh_names) ) then

    if ( size(target_edge_cells) == size(target_mesh_names) ) then

      self%nmaps = size(target_mesh_names)
      allocate(self%target_edge_cells(self%nmaps))
      allocate(self%target_mesh_names(self%nmaps))

      self%target_mesh_names(:) = target_mesh_names(:)
      self%target_edge_cells(:) = target_edge_cells(:)

      do i=1, self%nmaps

        if (target_edge_cells(i) < self%edge_cells) then
          remainder = mod(self%edge_cells, target_edge_cells(i))
          write(log_scratch_space,'(2(A,I0),A)')           &
               'Target edge_cells[',target_edge_cells(i),  &
               '] must be a factor of source edge_cells[', &
               self%edge_cells, ']'
        else if (target_edge_cells(i) > self%edge_cells) then
          remainder = mod(target_edge_cells(i), self%edge_cells)
          write(log_scratch_space,'(2(A,I0),A)')           &
               'Source edge_cells[',target_edge_cells(i),  &
               '] must be a factor of target edge_cells[', &
               self%edge_cells, ']'
        else
          ! Don't map to oneself
        end if

        if (remainder == 0) then
          self%target_edge_cells(i) = target_edge_cells(i)
        else
          call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
        end if

      end do

      write(target_mesh_names_str, '(A)') "'"//trim(adjustl(target_mesh_names(1)))//"'"
      write(target_edge_cells_str,'(I0)') target_edge_cells(1)

      if (size(target_mesh_names) > 1) then
        do i=2, self%nmaps
          write(target_mesh_names_str,'(A)')                   &
               trim(adjustl(target_mesh_names_str)) // ",'" // &
               trim(adjustl(target_mesh_names(i)))  // "'"
          write(target_edge_cells_str,'(A,I0)')                &
               trim(adjustl(target_edge_cells_str)) // ',',    &
               target_edge_cells(i)

        end do
      end if

      write(target_mesh_names_str,'(A)')  &
          'target_mesh_names=['//trim(adjustl(target_mesh_names_str))//']'
      write(target_edge_cells_str,'(A)')  &
          'target_edge_cells=['//trim(adjustl(target_edge_cells_str))//']'

      write(self%constructor_inputs,'(A)')          &
          trim(self%constructor_inputs) // ';' //   &
          trim(target_mesh_names_str)   // ';' //   &
          trim(target_edge_cells_str)

    end if

  end if

  return
end function gencube_ps_constructor

!-------------------------------------------------------------------------------
!> @brief   For each cell, calculates the set of cells to which it is adjacent.
!>          (Private Routine)
!> @details Allocates and populates the instance's cell_next(:,:) array
!>          with the id of each cell to which the index cell is adjacent.
!>
!> @param[in]   gen_cube   Generator strategy for a cubed-sphere
!> @param[out]  cell_next  A rank 2 (4,ncells)-sized array containing the
!>                         adjacency map.
!-------------------------------------------------------------------------------
subroutine calc_adjacency(gen_cube, cell_next)

  implicit none

  class(gencube_ps_type),      intent(in)  :: gen_cube
  integer(i_def), allocatable, intent(out) :: cell_next(:,:)

  integer(i_def) :: edge_cells, ncells, cpp, i
  integer(i_def) :: cell, astat, panel_number

  integer(i_def), allocatable :: panel_edge_cells_west(:)
  integer(i_def), allocatable :: panel_edge_cells_south(:)
  integer(i_def), allocatable :: panel_edge_cells_east(:)
  integer(i_def), allocatable :: panel_edge_cells_north(:)

  integer(i_def), allocatable :: panel_next(:,:)
  integer(i_def), allocatable :: panel_edge_cells(:,:,:) ! (cell id, panel edge, panel number )


  edge_cells = gen_cube%edge_cells
  cpp        = edge_cells*edge_cells
  ncells     = cpp*NPANELS


  allocate(cell_next(4, ncells), stat=astat)
  if (astat /= 0)                                               &
      call log_event( PREFIX//"Failure to allocate cell_next.", &
                      LOG_LEVEL_ERROR )

  allocate(panel_edge_cells_west(edge_cells))
  allocate(panel_edge_cells_south(edge_cells))
  allocate(panel_edge_cells_east(edge_cells))
  allocate(panel_edge_cells_north(edge_cells))

  cell_next = 0

  ! Panels are arranged and numbered as indicated
  !==============================================
  !      +---+
  !      | 5 |
  !  +---+---+---+---+
  !  | 1 | 2 | 3 | 4 |
  !  +---+---+---+---+
  !      | 6 |
  !      +---+
  !
  ! Cells in each panel are numbered from NW panel corner
  ! i.e. 3x3 panel number #1 is:
  !
  ! +-----+
  ! |1|2|3|
  ! |-+-+-|
  ! |4|5|6|
  ! |-+-+-|
  ! |7|8|9|
  ! +-----+

  allocate(panel_next(4, NPANELS))
  allocate(panel_edge_cells (edge_cells,4,NPANELS))

  ! Ordering : W,S,E,N
  panel_next(:,1) = [4,6,2,5]
  panel_next(:,2) = [1,6,3,5]
  panel_next(:,3) = [2,6,4,5]
  panel_next(:,4) = [3,6,1,5]
  panel_next(:,5) = [1,2,3,4]
  panel_next(:,6) = [1,4,3,2]

  call get_panel_edge_cell_ids( edge_cells, panel_edge_cells )

  ! Default settings
  do cell=1, ncells
    ! Default: W, S, E, N
    cell_next(:, cell) = (/ cell-1, cell+edge_cells, cell+1, cell-edge_cells /)
  end do



  ! Panel I
  !===========================================================================
  panel_number = 1
  do i=1, edge_cells
    ! Top edge
    cell = panel_edge_cells(i,N, panel_number)
    cell_next(N, cell) = panel_edge_cells(i, W, panel_next(N,panel_number) )

    ! Right edge
    cell = panel_edge_cells(i,E, panel_number)
    cell_next(E, cell) = panel_edge_cells(i, W, panel_next(E,panel_number) )

    ! Bottom edge
    cell = panel_edge_cells(i,S, panel_number)
    cell_next(S, cell) = panel_edge_cells(edge_cells+1-i, W, panel_next(S,panel_number) )

    ! Left edge
    cell = panel_edge_cells(i,W, panel_number)
    cell_next(W, cell) = panel_edge_cells(i, E, panel_next(W,panel_number) )
  end do

  ! Panel II
  !===========================================================================
  panel_number = 2
  do i=1, edge_cells
    ! Top edge
    cell = panel_edge_cells(i,N, panel_number)
    cell_next(N, cell) = panel_edge_cells(i, S, panel_next(N,panel_number) )

    ! Right edge
    cell = panel_edge_cells(i,E, panel_number)
    cell_next(E, cell) = panel_edge_cells(i, W, panel_next(E,panel_number) )

    ! Bottom edge
    cell = panel_edge_cells(i,S, panel_number)
    cell_next(S, cell) = panel_edge_cells(i, N, panel_next(S,panel_number) )

    ! Left edge
    cell = panel_edge_cells(i,W, panel_number)
    cell_next(W, cell) = panel_edge_cells(i, E, panel_next(W,panel_number) )
  end do

  ! Panel III
  !===========================================================================
  panel_number = 3
  do i=1, edge_cells

    ! Top edge
    cell = panel_edge_cells(i,N, panel_number)
    cell_next(N, cell) = panel_edge_cells(edge_cells+1-i, E, panel_next(N,panel_number) )

    ! Right edge
    cell = panel_edge_cells(i,E, panel_number)
    cell_next(E, cell) = panel_edge_cells(i, W, panel_next(E,panel_number) )

    ! Bottom edge
    cell = panel_edge_cells(i,S, panel_number)
    cell_next(S, cell) = panel_edge_cells(i, E, panel_next(S,panel_number) )

    ! Left edge
    cell = panel_edge_cells(i,W, panel_number)
    cell_next(W, cell) = panel_edge_cells(i, E, panel_next(W,panel_number) )

  end do


  ! Panel IV
  !===========================================================================
  panel_number = 4
  do i=1, edge_cells
    ! Top edge
    cell = panel_edge_cells(i,N, panel_number)
    cell_next(N, cell) = panel_edge_cells(edge_cells+1-i, N, panel_next(N,panel_number) )

    ! Right edge
    cell = panel_edge_cells(i,E, panel_number)
    cell_next(E, cell) = panel_edge_cells(i, W, panel_next(E,panel_number) )

    ! Bottom edge
    cell = panel_edge_cells(i,S, panel_number)
    cell_next(S, cell) = panel_edge_cells(edge_cells+1-i, S, panel_next(S,panel_number) )

    ! Left edge
    cell = panel_edge_cells(i,W, panel_number)
    cell_next(W, cell) = panel_edge_cells(i, E, panel_next(W,panel_number) )
  end do

  ! Panel V
  !===========================================================================
  panel_number = 5
  do i=1, edge_cells
    ! Top edge
    cell = panel_edge_cells(i,N, panel_number)
    cell_next(N, cell) = panel_edge_cells(edge_cells+1-i, N, panel_next(N,panel_number) )

    ! Right edge
    cell = panel_edge_cells(i,E, panel_number)
    cell_next(E, cell) = panel_edge_cells(edge_cells+1-i, N, panel_next(E,panel_number) )

    ! Bottom edge
    cell = panel_edge_cells(i,S, panel_number)
    cell_next(S, cell) = panel_edge_cells(i, N, panel_next(S,panel_number) )

    ! Left edge
    cell = panel_edge_cells(i,W, panel_number)
    cell_next(W, cell) = panel_edge_cells(i, N, panel_next(W,panel_number) )
  end do



  ! Panel VI
  !===========================================================================
  panel_number = 6
  do i=1, edge_cells
    ! Top edge
    cell = panel_edge_cells(i,N, panel_number)
    cell_next(N, cell) = panel_edge_cells(i, S, panel_next(N,panel_number) )

    ! Right edge
    cell = panel_edge_cells(i,E, panel_number)
    cell_next(E, cell) = panel_edge_cells(i, S, panel_next(E,panel_number) )

    ! Bottom edge
    cell = panel_edge_cells(i,S, panel_number)
    cell_next(S, cell) = panel_edge_cells(edge_cells+1-i, S, panel_next(S,panel_number) )

    ! Left edge
    cell = panel_edge_cells(i,W, panel_number)
    cell_next(W, cell) = panel_edge_cells(edge_cells+1-i, S, panel_next(W,panel_number) )
  end do

  return
end subroutine calc_adjacency

!-------------------------------------------------------------------------------
!> @brief   For each cell, calculates the four vertices which comprise it.
!>          (Private Routine)
!> @details Allocates and populates the instance's verts_on_cell(:,:) array with
!>          the vertices which form each cell.
!>
!> @param[in]   gen_cube       Generator strategy for a cubed-sphere
!> @param[out]  verts_on_cell  A rank 2 (4,ncells)-sized integer array of vertices
!>                             which constitute each cell.
!-------------------------------------------------------------------------------
subroutine calc_face_to_vert(gen_cube, verts_on_cell)

  implicit none

  class(gencube_ps_type),      intent(in)  :: gen_cube
  integer(i_def), allocatable, intent(out) :: verts_on_cell(:,:)

  integer(i_def) :: edge_cells, ncells, cpp
  integer(i_def) :: cell, idx, panel, nxf, astat

  edge_cells = gen_cube%edge_cells
  cpp        = edge_cells*edge_cells
  ncells     = cpp*NPANELS

  allocate(verts_on_cell(4, 6*cpp), stat=astat)

  if (astat /= 0)                                                   &
      call log_event( PREFIX//"Failure to allocate verts_on_cell.", &
                      LOG_LEVEL_ERROR )

  verts_on_cell = 0
  cell = 1
  nxf = 1

  ! NW vert of every cell in panels 1:4
  do cell = 1, 4*cpp
    verts_on_cell(NW, cell) = cell
  end do

  nxf = 4*cpp + 1

  ! Copy NE from E neighbour for every cell panels 1:4
  do cell = 1, 4*cpp
    verts_on_cell(NE, cell) = verts_on_cell(NW, gen_cube%cell_next(E, cell))
  end do

  ! Copy from S for non-S-edge rows of panels 1:4
  do panel = 1, 4
    do cell = (panel-1)*cpp+1, panel*cpp-edge_cells
      verts_on_cell(SW, cell) = verts_on_cell(NW, gen_cube%cell_next(S, cell))
      verts_on_cell(SE, cell) = verts_on_cell(NE, gen_cube%cell_next(S, cell))
    end do
  end do

  ! SE internals panel 5
  do idx = 0, edge_cells-2
    do cell = 4*cpp+(idx*edge_cells)+1, 4*cpp+(idx*edge_cells)+edge_cells-1
      verts_on_cell(SE, cell) = nxf
      nxf = nxf+1
      verts_on_cell(SW, gen_cube%cell_next(E, cell)) = verts_on_cell(SE, cell)
      verts_on_cell(NE, gen_cube%cell_next(S, cell)) = verts_on_cell(SE, cell)
      ! Transitive is valid here
      verts_on_cell(NW, gen_cube%cell_next(E, gen_cube%cell_next(S, cell))) &
                                                 = verts_on_cell(SE, cell)
    end do
  end do

  ! NW vert of every cell in panel 6
  do cell = 5*cpp+1, 6*cpp
    verts_on_cell(NW, cell) = nxf
    nxf = nxf+1
  end do

  ! Copy SW from S neighbour for non-S-edge rows of panel 6
  do cell = 5*cpp+1, 6*cpp-edge_cells
    verts_on_cell(SW, cell) = verts_on_cell(NW, gen_cube%cell_next(S, cell))
  end do

  ! SW verts of bottom row, panel 6
  do cell = 6*cpp-edge_cells+1, 6*cpp
    verts_on_cell(SW, cell) = nxf
    nxf = nxf+1
  end do

  ! SW verts of bottom row, panel 3
  do cell = 3*cpp-edge_cells+1, 3*cpp
    verts_on_cell(SW, cell) = nxf
    verts_on_cell(SE, gen_cube%cell_next(W, cell)) = verts_on_cell(SW, cell)
    nxf = nxf+1
  end do

  ! Vert at SE of panel 3
  cell = 3*cpp
  verts_on_cell(SE, cell) = nxf
  verts_on_cell(SW, gen_cube%cell_next(E, cell)) = verts_on_cell(SE, cell)

  ! Panel boundary joins...

  ! S=>W join, I=>VI
  do cell = cpp-edge_cells+1, cpp
    verts_on_cell(SW, cell) = verts_on_cell(SW, gen_cube%cell_next(S, cell))
    verts_on_cell(SE, cell) = verts_on_cell(NW, gen_cube%cell_next(S, cell))
  end do

  ! N=>W join, I=>V
  do cell = 1, edge_cells
    verts_on_cell(NW, gen_cube%cell_next(N, cell)) = verts_on_cell(NW, cell)
    verts_on_cell(SW, gen_cube%cell_next(N, cell)) = verts_on_cell(NE, cell)
  end do

  ! N=>E join, III=>V
  do cell = 2*cpp+1, 2*cpp+edge_cells
    verts_on_cell(SE, gen_cube%cell_next(N, cell)) = verts_on_cell(NW, cell)
    verts_on_cell(NE, gen_cube%cell_next(N, cell)) = verts_on_cell(NE, cell)
  end do

  ! E=>S join, III=>VI
  do cell = 3*cpp-edge_cells+1, 3*cpp
     verts_on_cell(NE, gen_cube%cell_next(S, cell)) = verts_on_cell(SW, cell)
     verts_on_cell(SE, gen_cube%cell_next(S, cell)) = verts_on_cell(SE, cell)
  end do

  ! N=>N join, IV=>V
  do cell = 3*cpp+1, 3*cpp+edge_cells
    verts_on_cell(NE, gen_cube%cell_next(N, cell)) = verts_on_cell(NW, cell)
    verts_on_cell(NW, gen_cube%cell_next(N, cell)) = verts_on_cell(NE, cell)
  end do

  ! Copy NE,SE from E neighbour for non-E-edge cells of panel 6
  do idx = 0, edge_cells-1
    do cell = 5*cpp+1+(idx*edge_cells), 5*cpp+1+(idx*edge_cells)+edge_cells-2
      verts_on_cell(NE, cell) = verts_on_cell(NW, gen_cube%cell_next(E, cell))
      verts_on_cell(SE, cell) = verts_on_cell(SW, gen_cube%cell_next(E, cell))
    end do
  end do

  ! S=>S join, VI=>IV
  do cell = 6*cpp-edge_cells+1, 6*cpp
    verts_on_cell(SW, gen_cube%cell_next(S, cell)) = verts_on_cell(SE, cell)
    verts_on_cell(SE, gen_cube%cell_next(S, cell)) = verts_on_cell(SW, cell)
  end do

  ! S=>N join, II=>VI
  do cell = 2*cpp-edge_cells+1, 2*cpp
    verts_on_cell(SW, cell) = verts_on_cell(NW, gen_cube%cell_next(S, cell))
    verts_on_cell(SE, cell) = verts_on_cell(NE, gen_cube%cell_next(S, cell))
  end do

  ! N=>S join, II=>V
  do cell = cpp+1, cpp+edge_cells
    verts_on_cell(SW, gen_cube%cell_next(N, cell)) = verts_on_cell(NW, cell)
    verts_on_cell(SE, gen_cube%cell_next(N, cell)) = verts_on_cell(NE, cell)
  end do

  return
end subroutine calc_face_to_vert

!-------------------------------------------------------------------------------
!> @brief   Calculates the edges which are found on each cell and the
!>          pair of vertices which are found on each edge.
!>          (Private Routine)
!> @details Allocates and populates both the edges_on_cell and
!>          verts_on_edge arrays for the instance.
!>
!> @param[in]   gen_cube       Generator strategy for a cubed-sphere
!> @param[out]  edges_on_cell  A rank-2 (4,ncells)-sized integer array of
!>                             the edges found on each cell.
!> @param[out]  verts_on_edge  A rank-2 (2,2*ncells)-sized integer array
!>                             of the vertices found on each edge.
!-------------------------------------------------------------------------------
subroutine calc_edges(gen_cube, edges_on_cell, verts_on_edge)

  implicit none

  class(gencube_ps_type),      intent(in)  :: gen_cube
  integer(i_def), allocatable, intent(out) :: edges_on_cell(:,:)
  integer(i_def), allocatable, intent(out) :: verts_on_edge(:,:)

  integer(i_def) :: edge_cells, ncells, cpp
  integer(i_def) :: cell, panel, idx, nxf, astat

  edge_cells = gen_cube%edge_cells
  cpp        = edge_cells*edge_cells
  ncells     = cpp*NPANELS

  allocate(edges_on_cell(4, ncells), stat=astat)

  if (astat /= 0)                                                   &
      call log_event( PREFIX//"Failure to allocate edges_on_cell.", &
                      LOG_LEVEL_ERROR )

  allocate(verts_on_edge(2, 2*ncells), stat=astat)

  if (astat /= 0)                                                   &
      call log_event( PREFIX//"Failure to allocate verts_on_edge.", &
                      LOG_LEVEL_ERROR )

  edges_on_cell = 0
  verts_on_edge = 0
  cell = 1
  nxf = 1

  do panel = 1, 4
    ! Top row of panel
    do cell = (panel-1)*cpp + 1, (panel-1)*cpp + edge_cells
      edges_on_cell(N, cell) = nxf
      edges_on_cell(W, cell) = nxf+1
      edges_on_cell(S, cell) = nxf+2
      verts_on_edge(1, nxf) = gen_cube%verts_on_cell(NW, cell)
      verts_on_edge(2, nxf) = gen_cube%verts_on_cell(NE, cell)
      verts_on_edge(1, nxf+1) = gen_cube%verts_on_cell(SW, cell)
      verts_on_edge(2, nxf+1) = gen_cube%verts_on_cell(NW, cell)
      verts_on_edge(1, nxf+2) = gen_cube%verts_on_cell(SE, cell)
      verts_on_edge(2, nxf+2) = gen_cube%verts_on_cell(SW, cell)
      nxf = nxf + 3
      ! Push W edge to W neighbour
      edges_on_cell(E, gen_cube%cell_next(W, cell)) = edges_on_cell(W, cell)
    end do
    ! Remainder of panel
    do cell = (panel-1)*cpp+1+edge_cells, panel*cpp
      edges_on_cell(W, cell) = nxf
      edges_on_cell(S, cell) = nxf+1
      verts_on_edge(1, nxf) = gen_cube%verts_on_cell(SW, cell)
      verts_on_edge(2, nxf) = gen_cube%verts_on_cell(NW, cell)
      verts_on_edge(1, nxf+1) = gen_cube%verts_on_cell(SE, cell)
      verts_on_edge(2, nxf+1) = gen_cube%verts_on_cell(SW, cell)
      nxf = nxf+2
      ! Copy N edge from N cell
      edges_on_cell(N, cell) = edges_on_cell(S, gen_cube%cell_next(N, cell))
      ! Push W edge to W neighbour
      edges_on_cell(E, gen_cube%cell_next(W, cell)) = edges_on_cell(W, cell)
    end do
  end do

  ! Panel V non-S-edge rows
  do cell = 4*cpp+1, 5*cpp-edge_cells
    edges_on_cell(S, cell) = nxf
    verts_on_edge(1, nxf) = gen_cube%verts_on_cell(SE, cell)
    verts_on_edge(2, nxf) = gen_cube%verts_on_cell(SW, cell)
    nxf = nxf + 1
    ! Push S edge to S neighbour
    edges_on_cell(N, gen_cube%cell_next(S, cell)) = edges_on_cell(S, cell)
  end do

  ! Panel V non-E-edge columns
  do idx = 0, edge_cells-1
    do cell = 4*cpp+1+idx*edge_cells, 4*cpp+(idx+1)*edge_cells-1
      edges_on_cell(E, cell) = nxf
      verts_on_edge(1, nxf) = gen_cube%verts_on_cell(NE, cell)
      verts_on_edge(2, nxf) = gen_cube%verts_on_cell(SE, cell)
      nxf = nxf + 1
      ! Push E edge to E neighbour
      edges_on_cell(W, gen_cube%cell_next(E, cell)) = edges_on_cell(E, cell)
    end do
  end do

  ! Panel VI non-S-edge rows
  do cell = 5*cpp+1, 6*cpp-edge_cells
    edges_on_cell(S, cell) = nxf
    verts_on_edge(1, nxf) = gen_cube%verts_on_cell(SE, cell)
    verts_on_edge(2, nxf) = gen_cube%verts_on_cell(SW, cell)
    nxf = nxf + 1
    ! Copy from N neighbour
    edges_on_cell(N, cell) = edges_on_cell(S, gen_cube%cell_next(N, cell))
  end do

  ! Panel VI non-E-edge columns
  do idx = 0, edge_cells-1
    do cell = 5*cpp+1+idx*edge_cells, 5*cpp+(idx+1)*edge_cells-1
      edges_on_cell(E, cell) = nxf
      verts_on_edge(1, nxf) = gen_cube%verts_on_cell(NE, cell)
      verts_on_edge(2, nxf) = gen_cube%verts_on_cell(SE, cell)
      nxf = nxf + 1
      ! Push E edge to E neighbour
      edges_on_cell(W, gen_cube%cell_next(E, cell)) = edges_on_cell(E, cell)
    end do
  end do

  ! Panel VI S-edge row copy in N
  do cell = 6*cpp-edge_cells+1, 6*cpp
    edges_on_cell(N, cell) = edges_on_cell(S, gen_cube%cell_next(N, cell))
  end do

  ! Join edges on panel boundaries...
  ! N=>W join, I=>V
  do cell = 1, edge_cells
    edges_on_cell(W, gen_cube%cell_next(N, cell)) = edges_on_cell(N, cell)
  end do

  ! S=>W join, I=>VI
  do cell = cpp-edge_cells+1, cpp
    edges_on_cell(W, gen_cube%cell_next(S, cell)) = edges_on_cell(S, cell)
  end do

  ! N=>E join, III=>V
  do cell = 2*cpp+1, 2*cpp+edge_cells
    edges_on_cell(E, gen_cube%cell_next(N, cell)) = edges_on_cell(N, cell)
  end do

  ! E=>S join, III=>VI
  do cell = 3*cpp-edge_cells+1, 3*cpp
    edges_on_cell(E, gen_cube%cell_next(S, cell)) = edges_on_cell(S, cell)
  end do

  ! N=>N join, IV=>V
  do cell = 3*cpp+1, 3*cpp+edge_cells
    edges_on_cell(N, gen_cube%cell_next(N, cell)) = edges_on_cell(N, cell)
  end do

  ! S=>N join, II=>VI
  do cell = 2*cpp-edge_cells+1, 2*cpp
    edges_on_cell(N, gen_cube%cell_next(S, cell)) = edges_on_cell(S, cell)
  end do

  ! N=>S join, II=>V
  do cell = cpp+1, cpp+edge_cells
    edges_on_cell(S, gen_cube%cell_next(N, cell)) = edges_on_cell(N, cell)
  end do

  ! S=>S join, IV=>VI
  do cell = 4*cpp-edge_cells+1, 4*cpp
    edges_on_cell(S, gen_cube%cell_next(S, cell)) = edges_on_cell(S, cell)
  end do

  return
end subroutine calc_edges

!-------------------------------------------------------------------------------
!> @brief   Calculates the coordinates of vertices in the mesh.
!> @details Assigns an (x,y) lat-long coordinate to each mesh
!>          vertex according to its Cartesian position in the mesh.
!>
!> @param[in]   gen_cube       Generator strategy for a cubed-sphere
!> @param[out]  vert_coords    A rank 2 (2,ncells)-sized real array of long and
!>                             lat coordinates (degrees) respectively for
!>                             each vertex.
!> @param[out]  coord_units_x  Units of x-coordinate.
!> @param[out]  coord_units_y  Units of y-coordinate.
!-------------------------------------------------------------------------------
subroutine calc_coords(gen_cube, vert_coords, coord_units_x, coord_units_y)

  implicit none

  class(gencube_ps_type),   intent(in)  :: gen_cube
  real(r_def), allocatable, intent(out) :: vert_coords(:,:)

  character(str_def), intent(out) :: coord_units_x
  character(str_def), intent(out) :: coord_units_y

  integer(i_def) :: ncells, edge_cells, nverts
  integer(i_def) :: cell, x, y, astat, cpp, vert, vert0
  real(r_def)    :: lat, long
  real(r_def)    :: x0, y0, z0
  real(r_def)    :: xs, ys, zs

  real(r_def)    :: dlambda
  real(r_def)    :: lambda1, lambda2
  real(r_def)    :: t1, t2

  real(r_def), parameter :: pio4 = PI/4.0_r_def

  edge_cells = gen_cube%edge_cells
  cpp        = edge_cells*edge_cells
  ncells     = cpp*NPANELS
  nverts     = ncells+2

  allocate(vert_coords(2, nverts), stat=astat)

  if (astat /= 0)                                                 &
      call log_event( PREFIX//"Failure to allocate vert_coords.", &
                      LOG_LEVEL_ERROR )

  vert_coords = 0.0_r_def
  dlambda = 0.5_r_def*PI/edge_cells  ! dlamba in radians
  vert = 1

! Panels I-IV
  do y=1,edge_cells
    lambda2 = (y-1)*dlambda - pio4
    t2 = tan(lambda2)
    do x=1,edge_cells
      lambda1 = (x-1)*dlambda - pio4
      t1 = tan(lambda1)

      ! Panel I
      xs = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2)
      ys = xs*t1
      zs = xs*t2

      call xyz2ll(xs, ys, zs, long, lat)

      vert_coords(1, vert) = long
      vert_coords(2, vert) = -lat

      ! Panel II
      vert0 = vert + cpp
      x0 = -ys
      y0 =  xs
      z0 =  zs

      call xyz2ll(x0, y0, z0, long, lat)
      vert_coords(1, vert0) = long
      vert_coords(2, vert0) = -lat

      ! Panel III
      vert0 = vert + 2*cpp
      x0 = -xs
      y0 = -ys
      z0 =  zs

      call xyz2ll(x0, y0, z0, long, lat)
      vert_coords(1, vert0) = long
      vert_coords(2, vert0) = -lat

      ! Panel IV
      vert0 = vert + 3*cpp
      x0 =  ys
      y0 = -xs
      z0 =  z0

      call xyz2ll(x0, y0, z0, long, lat)
      vert_coords(1, vert0) = long
      vert_coords(2, vert0) = -lat

      vert = vert + 1
    end do
  end do

! Panel V
! NB Change to cell-based vert lookup from here
  cell = 4*cpp

  do y=1, edge_cells-1
    lambda2 = y*dlambda - pio4 ! NB not y-1
    t2 = tan(lambda2)
    do x=1, edge_cells-1
      lambda1 = x*dlambda - pio4 ! NB not x-1
      t1 = tan(lambda1)

      xs = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2)
      ys = xs*t1
      zs = xs*t2
      ! Lookup vert with x offset
      vert0 = gen_cube%verts_on_cell(SE, cell+x)

      x0 = -ys
      y0 =  zs
      z0 = -xs

      call xyz2ll(x0, y0, z0, long, lat)
      vert_coords(1, vert0) = long
      vert_coords(2, vert0) = -lat

    end do
    cell = cell + edge_cells
  end do

! Panel VI
  cell = 5*cpp + 1

  do y=1, edge_cells
    lambda2 = (y-1)*dlambda - pio4
    t2 = tan(lambda2)
    do x=1, edge_cells
      lambda1 = (x-1)*dlambda - pio4
      t1 = tan(lambda1)

      xs = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2)
      ys = xs*t1
      zs = xs*t2
      ! Lookup vert
      vert0 = gen_cube%verts_on_cell(NW, cell)

      x0 = -ys
      y0 = -zs
      z0 =  xs

      call xyz2ll(x0, y0, z0, long, lat)
      vert_coords(1, vert0) = long
      vert_coords(2, vert0) = -lat

      cell = cell + 1
    end do
  end do

! Panel VI: Bottom edge
  cell = 6*cpp-edge_cells+1

  ! y constant
  lambda2 = edge_cells*dlambda - pio4
  t2 = tan(lambda2)
  do x=1, edge_cells
    lambda1 = (x-1)*dlambda - pio4
    t1 = tan(lambda1)

    xs = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2)
    ys = xs*t1
    zs = xs*t2

    vert0 = gen_cube%verts_on_cell(SW, cell)

    x0 = -ys
    y0 = -zs
    z0 =  xs

    call xyz2ll(x0, y0, z0, long, lat)
    vert_coords(1, vert0) = long
    vert_coords(2, vert0) = -lat

    cell = cell + 1
  end do

! Panel VI: Right edge
  cell = 5*cpp+edge_cells

  ! x constant
  lambda1 = edge_cells*dlambda - pio4
  t1 = tan(lambda1)
  do y=1, edge_cells
    lambda2 = (y-1)*dlambda - pio4
    t2 = tan(lambda2)

    xs = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2)
    ys = xs*t1
    zs = xs*t2

    vert0 = gen_cube%verts_on_cell(NE, cell)

    x0 = -ys
    y0 = -zs
    z0 =  xs

    call xyz2ll(x0, y0, z0, long, lat)
    vert_coords(1, vert0) = long
    vert_coords(2, vert0) = -lat

    cell = cell + edge_cells  ! NB Step size
  end do

! 6*edge_cells*edge_cells+2
  lambda1 = edge_cells*dlambda - pio4
  t1 = tan(lambda1)
  lambda2 = edge_cells*dlambda - pio4
  t2 = tan(lambda2)

  xs = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2)
  ys = xs*t1
  zs = xs*t2

  x0 = -ys
  y0 = -zs
  z0 =  xs


  call xyz2ll(x0, y0, z0, long, lat)

  vert0 = 6*cpp+2
  vert_coords(1, vert0) = long
  vert_coords(2, vert0) = -lat

  ! Output units from xyz2ll are in radians
  coord_units_x = 'radians'
  coord_units_y = 'radians'

  return
end subroutine calc_coords

!-------------------------------------------------------------------------------
!> @brief   Apply the Schmidt transform to the generation of the cubedsphere.
!> @details Attracts points to the north or south pole according
!>          to the value in self%equatorial_latitude. This gives a variable
!>          resolution mesh.
!>          (PRIVATE ROUTINE)
!>
!> @param[in,out] gen_cube  Generator strategy for a cubed-sphere
!-------------------------------------------------------------------------------
subroutine stretch_mesh(gen_cube)

  implicit none

  class(gencube_ps_type), intent(inout)  :: gen_cube

  real(r_def) :: stretching ! Holding variable for stretch function
  real(r_def) :: lat
  real(r_def) :: stretch_factor

  integer(i_def) :: nverts, vert

  nverts = size(gen_cube%vert_coords, dim=2)

  ! Find stretch factor from equatorial latitude
  stretch_factor = sqrt( (1.0_r_def - sin(gen_cube%equatorial_latitude)) &
                         / (1.0_r_def + sin(gen_cube%equatorial_latitude)) )

  ! Apply Schmidt stretching transformation
  if ( stretch_factor > 0.0_r_def ) then

    stretching = (1.0_r_def - stretch_factor**2) &
                /(1.0_r_def + stretch_factor**2)

    do vert = 1,nverts

      lat = gen_cube%vert_coords(2, vert)
      gen_cube%vert_coords(2, vert) = asin( (stretching + sin(lat)) &
                                           /(1.0_r_def + stretching*sin(lat)) )

    end do

  end if

end subroutine stretch_mesh

!-------------------------------------------------------------------------------
!> @brief Populates the arguments with the dimensions defining
!>        the mesh.
!>
!> @param[in]   self                   The gencube_ps_type instance reference.
!> @param[out]  num_nodes              The number of nodes on the mesh.
!> @param[out]  num_edges              The number of edges on the mesh.
!> @param[out]  num_faces              The number of faces on the mesh.
!> @param[out]  num_nodes_per_face     The number of nodes around each face.
!> @param[out]  num_edges_per_face     The number of edges around each face.
!> @param[out]  num_nodes_per_edge     The number of nodes around each edge.
!> @param[out]  max_num_faces_per_node The maximum number of faces surrounding
!>                                     each node.
!-------------------------------------------------------------------------------
subroutine get_dimensions(self, num_nodes, num_edges, num_faces,           &
                                num_nodes_per_face, num_edges_per_face,    &
                                num_nodes_per_edge, max_num_faces_per_node )

  implicit none

  class(gencube_ps_type), intent(in) :: self

  integer(i_def), intent(out) :: num_nodes
  integer(i_def), intent(out) :: num_edges
  integer(i_def), intent(out) :: num_faces
  integer(i_def), intent(out) :: num_nodes_per_face
  integer(i_def), intent(out) :: num_edges_per_face
  integer(i_def), intent(out) :: num_nodes_per_edge
  integer(i_def), intent(out) :: max_num_faces_per_node

  integer(i_def) :: edge_cells, cpp, ncells

  edge_cells = self%edge_cells
  cpp        = edge_cells*edge_cells
  ncells     = cpp*NPANELS

  num_faces = ncells
  num_nodes = ncells + 2
  num_edges = ncells * 2

  num_nodes_per_face = 4
  num_edges_per_face = 4
  num_nodes_per_edge = 2

  max_num_faces_per_node = self%max_num_faces_per_node

  return
end subroutine get_dimensions

!-------------------------------------------------------------------------------
!> @brief   Populates the argument array with the coordinates of
!>          the mesh's vertices.
!> @details Exposes the instance's vert_coords array to the caller.
!>
!> @param[out]  node_coordinates  The argument to receive the vert_coords data.
!> @param[out]  cell_coordinates  Cell centre coordinates
!> @param[out]  domain_extents    Principal coordiantes describing domain.
!> @param[out]  coord_units_x     Units of x-coordinate.
!> @param[out]  coord_units_y     Units of y-coordinate.
!-------------------------------------------------------------------------------
subroutine get_coordinates(self, node_coordinates, &
                                 cell_coordinates, &
                                 domain_extents,   &
                                 coord_units_x,    &
                                 coord_units_y)

  implicit none

  class(gencube_ps_type), intent(in)  :: self
  real(r_def),            intent(out) :: node_coordinates(:,:)
  real(r_def),            intent(out) :: cell_coordinates(:,:)
  real(r_def),            intent(out) :: domain_extents(:,:)
  character(str_def),     intent(out) :: coord_units_x
  character(str_def),     intent(out) :: coord_units_y

  node_coordinates = self%vert_coords
  cell_coordinates = self%cell_coords
  domain_extents   = self%domain_extents
  coord_units_x    = self%coord_units_x
  coord_units_y    = self%coord_units_y

  return
end subroutine get_coordinates

!-------------------------------------------------------------------------------
!> @brief   Populates the argument arrays with the corresponding mesh
!>          connectivity information.
!> @details Implements the connectivity-providing interface required
!>          by the ugrid writer.
!>
!>  @param[out]  face_node_connectivity  Face-node connectivity.
!>  @param[out]  face_edge_connectivity  Face-edge connectivity.
!>  @param[out]  face_face_connectivity  Face-face connectivity.
!>  @param[out]  edge_node_connectivity  Edge-node connectivity.
!-------------------------------------------------------------------------------
subroutine get_connectivity(self, face_node_connectivity, &
                                  face_edge_connectivity, &
                                  face_face_connectivity, &
                                  edge_node_connectivity )
  implicit none

  class(gencube_ps_type), intent(in) :: self

  integer(i_def), intent(out) :: face_node_connectivity(:,:)
  integer(i_def), intent(out) :: face_edge_connectivity(:,:)
  integer(i_def), intent(out) :: face_face_connectivity(:,:)
  integer(i_def), intent(out) :: edge_node_connectivity(:,:)

  face_node_connectivity = self%verts_on_cell
  face_edge_connectivity = self%edges_on_cell
  face_face_connectivity = self%cell_next
  edge_node_connectivity = self%verts_on_edge

  return
end subroutine get_connectivity

!-------------------------------------------------------------------------------
!> @brief  Gets the global mesh map collection which uses
!>         this mesh as the source mesh
!>
!> @return global_mesh_maps global_mesh_map_collection_type
!-------------------------------------------------------------------------------
function get_global_mesh_maps(self) result (global_mesh_maps)

  implicit none

  class(gencube_ps_type), target, intent(in) :: self

  type(global_mesh_map_collection_type), pointer  :: global_mesh_maps

  global_mesh_maps => self%global_mesh_maps

  return
end function get_global_mesh_maps

!-------------------------------------------------------------------------------
!> @brief   Generates the mesh and connectivity.
!> @details Calls each of the instance methods which calculate the
!>          specified mesh and populate the arrays.
!>
!> @param[in,out]  self  The gencube_ps_type instance reference.
!-------------------------------------------------------------------------------
subroutine generate(self)

  implicit none

  class(gencube_ps_type), intent(inout) :: self

  call calc_adjacency(self, self%cell_next)
  call calc_face_to_vert(self, self%verts_on_cell)
  call calc_edges(self, self%edges_on_cell, self%verts_on_edge)

  if (self%nmaps > 0_i_def) call calc_global_mesh_maps(self, PANEL_ROTATIONS)

  ! Co-ord output from calc_coords in radians
  call calc_coords(self, self%vert_coords,   &
                         self%coord_units_x, &
                         self%coord_units_y)

  call orient_lfric(self, PANEL_ROTATIONS)

  ! Stretching and smoothing of the mesh should be done before
  ! any rotations are done.
  if (self%nsmooth > 0_i_def) call smooth(self)

  if (self%equatorial_latitude /= 0.0_r_def) call stretch_mesh(self)

  if (self%rotate_mesh) call rotate_mesh_coords(self%vert_coords, &
                                                self%north_pole)

  call calc_cell_centres(self)

  ! Convert coordinate units to degrees to be CF compliant
  if (trim(self%coord_units_x) == 'radians') then
    self%vert_coords(1,:) = self%vert_coords(1,:) * radians_to_degrees
    self%cell_coords(1,:) = self%cell_coords(1,:) * radians_to_degrees
    self%coord_units_x = 'degrees_east'
  end if

  if (trim(self%coord_units_y) == 'radians') then
    self%vert_coords(2,:) = self%vert_coords(2,:) * radians_to_degrees
    self%cell_coords(2,:) = self%cell_coords(2,:) * radians_to_degrees
    self%coord_units_y = 'degrees_north'
  end if

  if (DEBUG) call write_mesh(self)

  return
end subroutine generate

!-------------------------------------------------------------------------------
!> @brief   PRIVATE subroutine to generate the reqeusted global mesh maps.
!> @details A map is generated for each requested target mesh based on this
!>          mesh objects mesh details, and those of the requested target
!>          meshes using calc_global_cell_map.
!>
!> @param[in,out]  self                  The gencube_ps_type instance reference.
!> @param[in]      panel_rotation_array  Cyclic rotations to apply to mesh
!>                                       connectivities
!-------------------------------------------------------------------------------
subroutine calc_global_mesh_maps(self, panel_rotation_array)

  implicit none

  class(gencube_ps_type), intent(inout) :: self
  integer(i_def),         intent(in)    :: panel_rotation_array(:)

  integer(i_def) :: source_id, source_cpp, source_ncells,         &
                    target_edge_cells, target_cpp, target_ncells, &
                    ntarget_per_source_x, ntarget_per_source_y, i
  integer(i_def), allocatable :: cell_map(:,:,:)



  allocate( self%global_mesh_maps, source=global_mesh_map_collection_type())

  source_id  = 1
  source_cpp = self%edge_cells*self%edge_cells
  source_ncells = source_cpp*NPANELS

  do i=1, size(self%target_mesh_names)

    target_edge_cells    = self%target_edge_cells(i)
    target_cpp           = target_edge_cells*target_edge_cells
    target_ncells        = target_cpp*NPANELS
    ntarget_per_source_x = max(1,target_edge_cells/self%edge_cells)
    ntarget_per_source_y = max(1,target_edge_cells/self%edge_cells)
    allocate(cell_map(ntarget_per_source_x,ntarget_per_source_y,source_ncells))
    call calc_global_cell_map(self, target_edge_cells, target_edge_cells, &
                              cell_map, panel_rotation_array )
    call self%global_mesh_maps%add_global_mesh_map( source_id, i+1, cell_map )

    deallocate(cell_map)

  end do

  return
end subroutine calc_global_mesh_maps


!-------------------------------------------------------------------------------
!> @brief Writes out the mesh and connectivity for debugging purposes.
!>
!> @param[in]  self  The gencube_ps_type instance reference.
!-------------------------------------------------------------------------------
subroutine write_mesh(self)

  use global_mesh_map_collection_mod, only: global_mesh_map_collection_type
  use global_mesh_map_mod, only: global_mesh_map_type

  use, intrinsic :: iso_fortran_env, only: stdout => output_unit

  implicit none

  class(gencube_ps_type), intent(in) :: self

  integer(i_def) :: i, j, ncells
  integer(i_def) :: cell, edge, vert

  type(global_mesh_map_type), pointer :: global_mesh_map => null()

  integer(i_def), allocatable :: cell_map (:,:)
  integer(i_def), allocatable :: tmp_map (:,:,:)
  integer(i_def) :: nsource
  integer(i_def) :: ntarget_per_cell
  integer(i_def) :: ntarget_cells_per_source_x
  integer(i_def) :: ntarget_cells_per_source_y

  character(str_long) :: tmp_str

  ncells = NPANELS*self%edge_cells*self%edge_cells

  write(stdout,'(A)')    "====DEBUG INFO===="
  write(stdout,'(A)')    "Mesh name: "// trim(self%mesh_name)
  write(stdout,'(A)')    "Geometry:  "// trim(key_from_geometry(self%geometry))
  write(stdout,'(A)')    "Topology:  "// trim(key_from_topology(self%topology))
  write(stdout,'(A,I0)') "Panels:    ", NPANELS
  write(stdout,'(A,I0)') "Panel edge cells: ", self%edge_cells
  write(stdout,'(A)')    "Coord_sys:  "// trim(key_from_coord_sys(self%coord_sys))
  write(stdout,'(A)')    "Co-ord (x) units: "// trim(self%coord_units_x)
  write(stdout,'(A)')    "Co-ord (y) units: "// trim(self%coord_units_y)
  write(stdout,'(A,I0)') "Smoothing iterations:     ", self%nsmooth
  write(stdout,'(A,I0)') "Mappings to other meshes: ", self%nmaps
  do i=1, self%nmaps
    write(stdout,'(T4,A,I0,A)') trim(self%target_mesh_names(i))// &
                                '(', self%target_edge_cells(i),')'
  end do


  write(stdout,'(A)') ''
  write(stdout,'(A)') '=========================='
  write(stdout,'(A)') ' Cell adjacency (W,S,E,N)'
  write(stdout,'(A)') '=========================='
  do cell=1, ncells
    tmp_str=''
    write(tmp_str,'(I07,A,4(I07,"  "))') cell,' => ', self%cell_next(:,cell)
    write(stdout,('(A)')) trim(tmp_str)
  end do


  write(stdout,'(A)') ''
  write(stdout,'(A)') '================================='
  write(stdout,'(A)') ' Vertices on cells (SW,SE,NW,NE)'
  write(stdout,'(A)') '================================='
  do cell=1, ncells
    tmp_str=''
    write(tmp_str,'(I07,A,4(I07,"  "))') cell,' => ', self%verts_on_cell(:,cell)
    write(stdout,('(A)')) trim(tmp_str)
  end do


  write(stdout,'(A)')  ''
  write(stdout,'(A)') '================================='
  write(stdout,'(A)') ' Edges on cells (W,S,E,N)'
  write(stdout,'(A)') '================================='
  do cell=1, ncells
    tmp_str=''
    write(tmp_str,'(I07,A,4(I07,"  "))') cell, ' => ', self%edges_on_cell(:,cell)
    write(stdout,('(A)')) trim(tmp_str)
  end do


  write(stdout,'(A)') ''
  write(stdout,'(A)') '================================='
  write(stdout,'(A)') ' Vertices on edges'
  write(stdout,'(A)') '================================='
  do edge=1, size(self%verts_on_edge, 2)
    tmp_str=''
    write(tmp_str,'(I07,A,I07,A,I07)')           &
        edge,' => ', self%verts_on_edge(1,edge), &
        ' -- ', self%verts_on_edge(2,edge)
    write(stdout,('(A)')) trim(tmp_str)
  end do


  write(stdout,'(A)')  ''
  write(stdout,'(A)') '================================='
  write(stdout,'(A)') ' Node coordinates (lon,lat)'
  write(stdout,'(A)') '================================='
  do vert=1, ncells+2
    tmp_str=''
    write(tmp_str,'(I07,A,2(F10.4,A))')          &
        vert,' => ( ', self%vert_coords(1,vert), &
        ',  ', self%vert_coords(2,vert), ' )'
    write(stdout,('(A)')) trim(tmp_str)
  end do

  write(stdout,'(A)')  ''
  write(stdout,'(A)') '================================='
  write(stdout,'(A)') ' Cell centre coordinates (lon,lat)'
  write(stdout,'(A)') '================================='
  do cell=1, ncells
    tmp_str=''
    write(tmp_str,'(I07,A,2(F10.4,A))')          &
        cell,' => ( ', self%cell_coords(1,cell), &
        ',  ', self%cell_coords(2,cell), ' )'
    write(stdout,('(A)')) trim(tmp_str)
  end do


  write(stdout,'(A)') ''
  write(stdout,'(A)') '================================='
  write(stdout,'(A)') ' Mappings to other meshes'
  write(stdout,'(A)') '================================='

  do i=1, self%nmaps
    global_mesh_map  => self%global_mesh_maps%get_global_mesh_map(1,i+1)
    nsource          = global_mesh_map%get_nsource_cells()
    ntarget_per_cell = global_mesh_map%get_ntarget_cells_per_source_cell()
    ntarget_cells_per_source_x = global_mesh_map%get_ntarget_cells_per_source_x()
    ntarget_cells_per_source_y = global_mesh_map%get_ntarget_cells_per_source_y()
    if (allocated(cell_map)) deallocate(cell_map)
    if (allocated(tmp_map)) deallocate(tmp_map)
    allocate(cell_map(ntarget_per_cell, nsource))
    allocate(tmp_map(ntarget_cells_per_source_x, ntarget_cells_per_source_y, 1))
    do j=1, nsource
      call global_mesh_map%get_cell_map([j], tmp_map)
      cell_map(:, j) = reshape(tmp_map(:,:,1), (/ ntarget_per_cell/) )
    end do
    write(stdout,'(2(A,I0),A)')                               &
        trim(self%mesh_name)//'(', self%edge_cells, ') => '// &
        trim(self%target_mesh_names(i))//'(', self%target_edge_cells(i), '):'
    do j=1, nsource
      write(stdout,'(I7,A,10(I0," "))') j,' => ' , cell_map(:, j)
    end do
  end do

  write(stdout,'(A)')    "====END DEBUG INFO===="
  return
end subroutine write_mesh

!-------------------------------------------------------------------------------
!> @brief   Reorients the cubed-sphere to be compatible with
!>          the orientation assumed by the LFRic infrastructure.
!>          (Private Routine)
!> @details Performs circular shifts on appropriate panels
!>
!> @param[in,out] gen_cube             Generator strategy for a cubed-sphere
!> @param[in]     panel_rotation_array Specifies the amount of rotation for
!>                                     each panel
!-------------------------------------------------------------------------------
subroutine orient_lfric(gen_cube, panel_rotation_array)

  implicit none

  class(gencube_ps_type), intent(inout) :: gen_cube
  integer(i_def),         intent(in)    :: panel_rotation_array(:)

  integer(i_def) :: cpp, p0, p1, i

  cpp = gen_cube%edge_cells*gen_cube%edge_cells

  ! Loop through panels
  do i = 1, SIZE(panel_rotation_array)

    ! Get indices of first and last cells of panel
    p0 = (i-1)*cpp + 1
    p1 = i*cpp

    ! Rotate left if panel rotation is 1
    if (panel_rotation_array(i) == 1_i_def) then
      ! verts
      gen_cube%verts_on_cell(:, p0:p1) = cshift(gen_cube%verts_on_cell(:, p0:p1), 1, 1)
      ! adj
      gen_cube%cell_next(:, p0:p1) = cshift(gen_cube%cell_next(:, p0:p1), 1, 1)
      ! edges
      gen_cube%edges_on_cell(:, p0:p1) = cshift(gen_cube%edges_on_cell(:, p0:p1), 1, 1)

      ! Rotate right if panel rotation is -1
    else if (panel_rotation_array(i) == -1_i_def) then
      ! verts
      gen_cube%verts_on_cell(:, p0:p1) = cshift(gen_cube%verts_on_cell(:, p0:p1), -1, 1)
      ! adj
      gen_cube%cell_next(:, p0:p1) = cshift(gen_cube%cell_next(:, p0:p1), -1, 1)
      ! edges
      gen_cube%edges_on_cell(:, p0:p1) = cshift(gen_cube%edges_on_cell(:, p0:p1), -1, 1)

    end if

  end do

  return
end subroutine orient_lfric

!-------------------------------------------------------------------------------
!> @brief   Smooth the cube grid (Private Routine)
!> @details Smooth the grid by iteratively computing the face centres as
!>          barycentres of the surrounding vertices and then the vertices
!>          as barycentres of the surrounding faces
!>
!> @param[in,out] gen_cube  Generator strategy for a cubed-sphere
!-------------------------------------------------------------------------------
subroutine smooth( gen_cube )

  implicit none

  class(gencube_ps_type), intent(inout) :: gen_cube

  integer(i_def) :: ncells
  integer(i_def) :: nverts
  integer(i_def) :: vert_id

  real(r_def)    :: x0(3)
  real(r_def)    :: xc(3)
  real(r_def)    :: radius_ratio
  real(r_def)    :: ll(2)

  integer(i_def), allocatable, dimension(:,:) :: cell_on_vert
  integer(i_def), allocatable, dimension(:)   :: ncell_on_vert
  real(r_def),    allocatable, dimension(:,:) :: cell_coords

  ! Counters
  integer(i_def) :: i, j, smooth_pass, cell, vert

  ncells = NPANELS*gen_cube%edge_cells*gen_cube%edge_cells
  nverts = ncells + 2

  allocate( cell_on_vert(4,nverts) )
  allocate( ncell_on_vert(nverts) )
  allocate( cell_coords(3,ncells) )

  ! Preliminary - Compute cell on vertices look up
  cell_on_vert(:,:) = -1_i_def
  ncell_on_vert(:)  =  0_i_def

  do cell=1, ncells
    do i=1, 4
      vert_id = gen_cube%verts_on_cell(i,cell)
      do j=1, 4
        if (cell_on_vert(j,vert_id) == -1_i_def ) then
          cell_on_vert(j,vert_id) = cell
          ncell_on_vert(vert_id)  = ncell_on_vert(vert_id) + 1_i_def
          exit
        end if
      end do
    end do
  end do

  ! Preliminary - Compute cell centre coordinates
  call calc_cell_centres( gen_cube )

  do cell=1, ncells
    call ll2xyz( gen_cube%cell_coords(1,cell), &
                 gen_cube%cell_coords(2,cell), &
                 cell_coords(1,cell),          &
                 cell_coords(2,cell),          &
                 cell_coords(3,cell) )

  end do

  do smooth_pass=1, gen_cube%nsmooth

    ! Compute vertices of barycentres of surrounding faces
    do vert=1, nverts
      xc(:) = 0.0_r_def
      do cell=1, ncell_on_vert(vert)
        xc(:) = xc(:) + cell_coords(:,cell_on_vert(cell, vert))
      end do
      radius_ratio = 1.0_r_def/sqrt( xc(1)**2 + xc(2)**2 + xc(3)**2)
      x0(:) =  xc(:)*radius_ratio
      call xyz2ll( x0(1), x0(2), x0(3),          &
                   gen_cube%vert_coords(1,vert), &
                   gen_cube%vert_coords(2,vert) )
    end do

    ! Compute faces as barycentres of surrounding vertices
    do cell=1, ncells
      xc(:) = 0.0_r_def
      do vert=1, 4
        ll = gen_cube%vert_coords(:,gen_cube%verts_on_cell(vert,cell))
        call ll2xyz(ll(1),ll(2),x0(1),x0(2),x0(3))
        xc(:) = xc(:) + x0(:)
      end do
      radius_ratio = 1.0_r_def/sqrt( xc(1)**2 + xc(2)**2 + xc(3)**2)
      cell_coords(:,cell) = xc(:)*radius_ratio
    end do

  end do

  return
end subroutine smooth

!-------------------------------------------------------------------------------
!> @brief   Calculates the mesh cell centres. (Private Routine)
!> @details The face centres for the mesh are calculated based on the current
!>          node coordinates. The node_cordinates are assumed to be in [lon, lat].
!>          Resulting face centre coordinates are in [lon, lat].
!>
!> @param[in,out] gen_cube  Generator strategy for a cubed-sphere
!-------------------------------------------------------------------------------
subroutine calc_cell_centres( gen_cube )

  implicit none

  class(gencube_ps_type), intent(inout) :: gen_cube

  integer(i_def) :: ncells

  real(r_def)    :: radius_ratio

  integer(i_def), allocatable :: verts_on_cell(:)

  real(r_def),    allocatable :: cell_vert_coords_xyz(:,:)
  real(r_def),    allocatable :: cell_vert_coords_ll(:,:)

  real(r_def), allocatable :: cell_centre_xyz(:)

  integer(i_def) :: nverts_per_cell = 4

  ! Counters
  integer(i_def) :: cell, vert

  ncells = NPANELS*gen_cube%edge_cells*gen_cube%edge_cells

  allocate( verts_on_cell(nverts_per_cell) )
  allocate( cell_vert_coords_xyz(3,nverts_per_cell) )
  allocate( cell_vert_coords_ll(2,nverts_per_cell) )
  allocate( cell_centre_xyz(3) )

  if (.not. allocated(gen_cube%cell_coords)) allocate( gen_cube%cell_coords(2,ncells) )

  gen_cube%cell_coords(:,:) = 0.0_r_def

  do cell=1, ncells
    cell_centre_xyz(:) = 0.0_r_def

    ! Get the vertex ids on the cell
    verts_on_cell(:) = gen_cube%verts_on_cell(:,cell)

    do vert=1, nverts_per_cell
      ! Get the vertex coords (in radians)
      cell_vert_coords_ll(:,vert) = gen_cube%vert_coords(:,verts_on_cell(vert))

      ! Get vertex coords as cartesian (x,y,z)
      call ll2xyz( cell_vert_coords_ll(1,vert),  &
                   cell_vert_coords_ll(2,vert),  &
                   cell_vert_coords_xyz(1,vert), &
                   cell_vert_coords_xyz(2,vert), &
                   cell_vert_coords_xyz(3,vert) )

      cell_centre_xyz(:) = cell_centre_xyz(:) + cell_vert_coords_xyz(:,vert)

    end do

    radius_ratio = 1.0_r_def/sqrt(  cell_centre_xyz(1)**2 &
                                  + cell_centre_xyz(2)**2 &
                                  + cell_centre_xyz(3)**2 )

    cell_centre_xyz(:) = cell_centre_xyz(:) * radius_ratio

    ! Convert cell centre back to lat long
    call xyz2ll( cell_centre_xyz(1),           & ! x
                 cell_centre_xyz(2),           & ! y
                 cell_centre_xyz(3),           & ! z
                 gen_cube%cell_coords(1,cell), & ! longitude
                 gen_cube%cell_coords(2,cell) )  ! latititude
  end do

end subroutine calc_cell_centres

!-----------------------------------------------------------------------------
!> @brief Returns the number of panels in the mesh topology.
!> @description Panels are a subset of cells in the mesh domain which may
!>              exhibit common properties.
!> @return answer Number of panels resulting from this generation strategy.
!-----------------------------------------------------------------------------
function get_number_of_panels( self ) result( answer )

  implicit none

  class(gencube_ps_type), intent(in) :: self

  integer(i_def) :: answer

  answer = self%npanels

end function get_number_of_panels


!-----------------------------------------------------------------------------
!> @brief Returns mesh metadata information.
!>
!> @param[out]  mesh_name          Optional, Name of mesh instance to generate
!> @param[out]  geometry           Optional, Mesh domain surface type.
!> @param[out]  topology           Optional, Mesh boundary/connectivity type
!> @param[out]  coord_sys          Optional, Coordinate system to position nodes.
!> @param[out]  periodic_x         Optional, Domain periodicity in x-axis
!> @param[out]  periodic_y         Optional, Domain periodicity in y-axis
!> @param[out]  edge_cells_x       Optional, Number of panel edge cells (x-axis).
!> @param[out]  edge_cells_y       Optional, Number of panel edge cells (y-axis).
!> @param[out]  constructor_inputs Optional, Inputs used to create this mesh from
!>                                           the this ugrid_generator_type
!> @param[out]  nmaps              Optional, Number of maps to create with this mesh
!>                                           as source mesh.
!> @param[out]  rim_depth          Optional, Depth of LBC mesh rim (in cells)
!> @param[out]  void_cell          Optional, Cell ID for null connectivity.
!> @param[out]  target_mesh_names  Optional, Mesh names of the target meshes that
!>                                           this mesh has maps for.
!> @param[out]  maps_edge_cells_x  Optional, Number of panel edge cells (x-axis) of
!>                                           target mesh(es) to create map(s) for.
!> @param[out]  maps_edge_cells_y  Optional, Number of panel edge cells (y-axis) of
!>                                           target mesh(es) to create map(s) for.
!> @param[out]  north_pole         Optional, [Longitude, Latitude] of north pole
!>                                           used for domain orientation (degrees)
!> @param[out]  null_island        Optional, [Longitude, Latitude] of null
!>                                           island used for domain orientation (degrees)
!> @param[out]  equatorial_latitude Optional, latitude of equator following
!>                                  stretching torwards pole (degrees)
!-----------------------------------------------------------------------------
subroutine get_metadata( self,               &
                         mesh_name,          &
                         geometry,           &
                         topology,           &
                         coord_sys,          &
                         periodic_xy,        &
                         edge_cells_x,       &
                         edge_cells_y,       &
                         constructor_inputs, &
                         nmaps,              &
                         rim_depth,          &
                         void_cell,          &
                         target_mesh_names,  &
                         maps_edge_cells_x,  &
                         maps_edge_cells_y,  &
                         north_pole,         &
                         null_island,        &
                         equatorial_latitude )

  implicit none

  class(gencube_ps_type), intent(in)  :: self

  character(str_def), optional, intent(out) :: mesh_name
  character(str_def), optional, intent(out) :: geometry
  character(str_def), optional, intent(out) :: topology
  character(str_def), optional, intent(out) :: coord_sys
  logical(l_def),     optional, intent(out) :: periodic_xy(2)
  integer(i_def),     optional, intent(out) :: edge_cells_x
  integer(i_def),     optional, intent(out) :: edge_cells_y
  integer(i_def),     optional, intent(out) :: nmaps
  integer(i_def),     optional, intent(out) :: rim_depth
  integer(i_def),     optional, intent(out) :: void_cell

  character(str_longlong), optional, intent(out) :: constructor_inputs

  character(str_def), allocatable, optional, intent(out) :: target_mesh_names(:)
  integer(i_def),     allocatable, optional, intent(out) :: maps_edge_cells_x(:)
  integer(i_def),     allocatable, optional, intent(out) :: maps_edge_cells_y(:)

  real(r_def),    optional, intent(out) :: north_pole(2)
  real(r_def),    optional, intent(out) :: null_island(2)
  real(r_def),    optional, intent(out) :: equatorial_latitude

  if (present(mesh_name))    mesh_name      = trim(self%mesh_name)
  if (present(geometry))     geometry       = key_from_geometry(self%geometry)
  if (present(topology))     topology       = key_from_topology(self%topology)
  if (present(coord_sys))    coord_sys      = key_from_coord_sys(self%coord_sys)
  if (present(periodic_xy))  periodic_xy    = .true.
  if (present(edge_cells_x)) edge_cells_x   = self%edge_cells
  if (present(edge_cells_y)) edge_cells_y   = self%edge_cells
  if (present(nmaps))        nmaps          = self%nmaps
  if (present(rim_depth))    rim_depth      = imdi
  if (present(void_cell))    void_cell      = VOID_ID

  if (present(north_pole))     north_pole(:)  = radians_to_degrees * self%north_pole(:)
  if (present(null_island))    null_island(:) = radians_to_degrees * self%null_island(:)
  if (present(equatorial_latitude)) &
                          equatorial_latitude = radians_to_degrees * self%equatorial_latitude

  if (present(constructor_inputs)) constructor_inputs = trim(self%constructor_inputs)

  if (self%nmaps > 0) then
    if (present(target_mesh_names)) target_mesh_names = self%target_mesh_names
    if (present(maps_edge_cells_x)) maps_edge_cells_x = self%target_edge_cells
    if (present(maps_edge_cells_x)) maps_edge_cells_y = self%target_edge_cells
  end if

  return
end subroutine get_metadata

!-------------------------------------------------------------------------------
!> @brief Subroutine to manually deallocate any memory used by the object.
!-------------------------------------------------------------------------------
subroutine clear(self)

  implicit none

  class (gencube_ps_type), intent(inout) :: self

  if (allocated(self%target_mesh_names)) deallocate( self%target_mesh_names )
  if (allocated(self%target_edge_cells)) deallocate( self%target_edge_cells )

  if (allocated(self%cell_next))     deallocate( self%cell_next     )
  if (allocated(self%verts_on_cell)) deallocate( self%verts_on_cell )
  if (allocated(self%edges_on_cell)) deallocate( self%edges_on_cell )
  if (allocated(self%verts_on_edge)) deallocate( self%verts_on_edge )
  if (allocated(self%vert_coords))   deallocate( self%vert_coords   )
  if (allocated(self%cell_coords))   deallocate( self%cell_coords   )

  if (allocated(self%global_mesh_maps)) then
    call self%global_mesh_maps%clear()
    deallocate( self%global_mesh_maps )
  end if


  return
end subroutine clear

!-----------------------------------------------------------------------------
!> @brief Finaliser for the cubedsphere ugrid mesh generator (gencube_ps_type)
!-----------------------------------------------------------------------------
subroutine gen_cubedsphere_final(self)

  implicit none

  type(gencube_ps_type), intent(inout) :: self

  call self%clear()

  return
end subroutine gen_cubedsphere_final


subroutine get_panel_edge_cell_ids( edge_cells, panel_edge_cells )

  implicit none

!
!          North
!      o---->>>----o         Cell ids on panel edges are
!      |           |         listed in direction shown:
! West Y   Panel   Y East
!      |           |         panel_edge_cells( cell_ids,
!      o---->>>----o                           side of the panel,
!          South                               panel number )
!

  integer(i_def), intent(in)  :: edge_cells
  integer(i_def), intent(out) :: panel_edge_cells(edge_cells,4,NPANELS)

  integer(i_def) :: i, cpp, panel

  cpp = edge_cells*edge_cells

  do panel=1, NPANELS
    ! Panel edge ordering W,S,E,N
    do i=1, edge_cells
      panel_edge_cells(i,W,panel) = (panel-1)*cpp + edge_cells*(i-1) + 1
      panel_edge_cells(i,S,panel) = (panel-1)*cpp + (cpp-edge_cells) + i
      panel_edge_cells(i,E,panel) = (panel-1)*cpp + edge_cells*(i)
      panel_edge_cells(i,N,panel) = (panel-1)*cpp + i
    end do
  end do
end subroutine get_panel_edge_cell_ids



!>==============================================================================
!> @brief Sets common partition parameters to be applied to global meshes.
!>
!> @param[out]  decomposition     Object containing decomposition parameters and
!>                                method
!> @param[out]  partitioner_ptr   Mesh partitioning strategy
!>==============================================================================
subroutine set_partition_parameters( decomposition, partitioner_ptr )

  use panel_decomposition_mod, only: panel_decomposition_type,           &
                                     custom_decomposition_type,          &
                                     auto_decomposition_type,            &
                                     row_decomposition_type,             &
                                     column_decomposition_type,          &
                                     auto_nonuniform_decomposition_type, &
                                     guided_nonuniform_decomposition_type

  use partition_mod,     only: partitioner_interface,   &
                               partitioner_cubedsphere, &
                               partitioner_cubedsphere_serial

  ! Configuration modules
  use partitions_config_mod, only: n_partitions,                        &
                                   panel_xproc, panel_yproc,            &
                                   panel_decomposition,                 &
                                   panel_decomposition_auto,            &
                                   panel_decomposition_row,             &
                                   panel_decomposition_column,          &
                                   panel_decomposition_custom,          &
                                   panel_decomposition_auto_nonuniform, &
                                   panel_decomposition_guided_nonuniform

  implicit none

  class(panel_decomposition_type),  intent(out), allocatable :: decomposition
  procedure(partitioner_interface), intent(out), pointer :: partitioner_ptr

  partitioner_ptr => null()

  if ( n_partitions == 1 ) then
    ! Use the serial cubed-sphere partitioner
    partitioner_ptr => partitioner_cubedsphere_serial
    call log_event( "Using serial cubed sphere partitioner", &
                    LOG_LEVEL_INFO )
  else if( mod(n_partitions, NPANELS) == 0 )then
    ! Use the parallel cubed-sphere partitioner
    partitioner_ptr => partitioner_cubedsphere

    select case(panel_decomposition)
      case( panel_decomposition_auto )
        decomposition = auto_decomposition_type()

      case( panel_decomposition_row )
        decomposition = row_decomposition_type()

      case( panel_decomposition_column )
        decomposition = column_decomposition_type()

      case( panel_decomposition_custom )
        decomposition = custom_decomposition_type( panel_xproc, panel_yproc )

      case ( panel_decomposition_auto_nonuniform )
        decomposition = auto_nonuniform_decomposition_type()

      case ( panel_decomposition_guided_nonuniform )
        decomposition = guided_nonuniform_decomposition_type( panel_xproc )

      case default
        call log_event( "Missing entry for panel decomposition, "// &
                        "specify 'auto' if unsure.", LOG_LEVEL_ERROR )

    end select

    call log_event( "Using parallel cubed sphere partitioner", &
                     LOG_LEVEL_INFO )

  else
    call log_event( "Number of partitions must be 1 "//      &
                    "or a multiple of the number of panels", &
                     LOG_LEVEL_ERROR )
  end if

end subroutine set_partition_parameters

end module gencube_ps_mod
