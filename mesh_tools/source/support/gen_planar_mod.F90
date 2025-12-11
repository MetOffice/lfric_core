!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief   Module to define the gen_planar_type
!> @details This type is a subclass of the ugrid_generator_type. It describes
!>          a planar mesh format suitable for storage as a ugrid file.
!>          All required connectivity is calculated via the generate method
!>          and made available to the ugrid writer.
!>
!-------------------------------------------------------------------------------
module gen_planar_mod
!-------------------------------------------------------------------------------

  use calc_global_cell_map_mod,       only: calc_global_cell_map
  use constants_mod,                  only: r_def, i_def, l_def, str_def, &
                                            str_long, imdi, rmdi, emdi,   &
                                            str_longlong,                 &
                                            radians_to_degrees,           &
                                            degrees_to_radians, PI
  use global_mesh_map_collection_mod, only: global_mesh_map_collection_type
  use global_mesh_map_mod,            only: generate_global_mesh_map_id
  use log_mod,                        only: log_event, log_scratch_space, &
                                            LOG_LEVEL_ERROR, LOG_LEVEL_INFO
  use mesh_config_mod,                only: key_from_coord_sys,          &
                                            key_from_geometry,           &
                                            key_from_topology,           &
                                            coord_sys_ll, coord_sys_xyz, &
                                            topology_non_periodic,       &
                                            geometry_spherical
  use planar_mesh_config_mod,         only: apply_stretch_transform
  use reference_element_mod,          only: reference_element_type, &
                                            reference_cube_type,    &
                                            W, S, E, N,             &
                                            SWB, SEB, NWB, NEB
  use ugrid_generator_mod,            only: ugrid_generator_type

  use rotation_mod,                   only: rotate_mesh_coords, &
                                            TRUE_NORTH_POLE_LL, &
                                            TRUE_NULL_ISLAND_LL
  use stretch_transform_mod,          only: stretch_transform,  &
                                            calculate_settings

  implicit none

  private

  public :: set_partition_parameters

  ! Mesh Vertex directions: local aliases for reference_element_mod
  ! values.
  integer(i_def), parameter :: NW = NWB
  integer(i_def), parameter :: NE = NEB
  integer(i_def), parameter :: SE = SEB
  integer(i_def), parameter :: SW = SWB

  integer(i_def), parameter :: OPEN     =  100
  integer(i_def), parameter :: CLOSED   = -100
  integer(i_def), parameter :: PERIODIC =  101

  ! Set to -9999 to be used for fill value so meshes
  ! are more Cf-compliant.
  integer(i_def), parameter :: VOID_ID = -9999

  ! For a planar meshes there is only one panel.
  integer(i_def), parameter :: NPANELS  = 1

  integer(i_def), parameter :: NBORDERS = 4

  ! Prefix for error messages.
  character(len=*), parameter :: PREFIX = "[Planar Mesh] "

  ! Flag to print out mesh data for debugging purposes.
  logical(l_def),   parameter :: DEBUG = .false.

  type, extends(ugrid_generator_type), public :: gen_planar_type

    private

    logical(l_def)     :: generated = .false.

    character(str_def) :: mesh_name
    integer(i_def)     :: geometry
    integer(i_def)     :: topology
    integer(i_def)     :: coord_sys = emdi
    character(str_def) :: coord_units_x
    character(str_def) :: coord_units_y
    integer(i_def)     :: edge_cells_x
    integer(i_def)     :: edge_cells_y
    integer(i_def)     :: fine_mesh_edge_cells_x
    integer(i_def)     :: fine_mesh_edge_cells_y
    integer(i_def)     :: npanels = NPANELS
    real(r_def)        :: domain_size(2)
    real(r_def)        :: domain_centre(2) = [0.0_r_def,0.0_r_def]
    real(r_def)        :: domain_extents(2,4)
    real(r_def)        :: north_pole(2)  = TRUE_NORTH_POLE_LL
    real(r_def)        :: null_island(2) = TRUE_NULL_ISLAND_LL
    real(r_def)        :: equatorial_latitude = TRUE_NULL_ISLAND_LL(2)

    character(str_longlong) :: constructor_inputs

    real(r_def)    :: dx, dy
    logical(l_def) :: periodic_xy(2)
    logical(l_def) :: rotate_mesh   = .false.

    ! Unique element types.
    integer(i_def) :: n_nodes
    integer(i_def) :: n_edges
    integer(i_def) :: n_faces

    ! Location index sets of cells adjacent to
    ! domain boundaries.
    integer(i_def), allocatable :: north_cells(:)
    integer(i_def), allocatable :: east_cells(:)
    integer(i_def), allocatable :: south_cells(:)
    integer(i_def), allocatable :: west_cells(:)

    ! Connectivity and coordinates.
    integer(i_def), allocatable :: cell_next(:,:)     ! (4, edge_cells_x*edge_cells_y)
    integer(i_def), allocatable :: verts_on_cell(:,:) ! (4, edge_cells_x*edge_cells_y)
    integer(i_def), allocatable :: edges_on_cell(:,:) ! (4, edge_cells_x*edge_cells_y)
    integer(i_def), allocatable :: verts_on_edge(:,:) ! (2, edge_cells_x*edge_cells_y)
    real(r_def),    allocatable :: vert_coords(:,:)   ! (2, edge_cells_x*edge_cells_y)
    real(r_def),    allocatable :: cell_coords(:,:)   ! (2, edge_cells_x*edge_cells_y)

    ! Intergrid maps.
    integer(i_def) :: nmaps
    character(str_def), allocatable :: target_mesh_names(:)
    integer(i_def),     allocatable :: target_edge_cells_x(:)
    integer(i_def),     allocatable :: target_edge_cells_y(:)
    type(global_mesh_map_collection_type), allocatable :: global_mesh_maps

    ! Hold variables from the reference element.
    ! Done because the available Cray compiler has internal compiler errors
    ! when attempting to include the reference_element_type.
    integer(i_def) :: nodes_per_face
    integer(i_def) :: edges_per_face
    integer(i_def) :: nodes_per_edge
    integer(i_def) :: max_num_faces_per_node

  contains

    procedure :: is_generated
    procedure :: get_corner_gid

    procedure :: generate
    procedure :: get_number_of_panels
    procedure :: get_metadata
    procedure :: get_dimensions
    procedure :: get_coordinates
    procedure :: get_connectivity
    procedure :: get_global_mesh_maps

    procedure :: write_mesh
    procedure :: clear
    final     :: gen_planar_final

  end type gen_planar_type

!-------------------------------------------------------------------------------
  interface gen_planar_type
    module procedure gen_planar_constructor
  end interface gen_planar_type
!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------
!> @brief   Constructor for gen_planar_type
!> @details Accepts inputs to configure a planar mesh. Intergrid mesh maps
!>          may be are optionally specified with varying numbers of grid cells.
!>
!>          Note: All target meshes specified for mapping to must be an integer
!>                factor of the main mesh, e.g. edge_cells_x=3;
!>                target_edge_cells_x=[6,9,12,15]
!>
!> @param[in] mesh_name       Name of this mesh topology.
!> @param[in] edge_cells_x    Number of cells in planar mesh x-axis.
!> @param[in] edge_cells_y    Number of cells in planar mesh y-axis.
!> @param[in] fine_mesh_edge_cells_x Number of cells in planar mesh x-axis.
!> @param[in] fine_mesh_edge_cells_y Number of cells in planar mesh y-axis.
!> @param[in] periodic_x      Logical for specifying periodicity in x-axis.
!> @param[in] periodic_y      Logical for specifying periodicity in y-axis.
!> @param[in] domain_size     Size of domain in x/y axes
!> @param[in] domain_centre   Co-ordinates of domain centre
!> @param[in] coord_sys       Coordinate system used to position nodes.
!> @param[in, optional] target_mesh_names
!>                            Names of mesh(es) to map to.
!> @param[in, optional] target_edge_cells_x
!>                            Number of cells in x axis of
!>                            target mesh(es) to map to.
!> @param[in, optional] target_edge_cells_y
!>                            Number of cells in y axis of
!>                            target mesh(es) to map to.
!> @param[in, optional] rotate_mesh
!>                            Logical to indicate whether to rotate the pole.
!> @param[in, optional] target_north_pole
!>                            The [longitude,latitude] co-ords for the new
!>                            north pole location.
!> @param[in, optional] target_null_island
!>                            The [longitude,latitude] co-ords for the new
!>                            null island location.
!-------------------------------------------------------------------------------
function gen_planar_constructor( reference_element,          &
                                 mesh_name,                  &
                                 geometry,                   &
                                 topology,                   &
                                 coord_sys,                  &
                                 edge_cells_x, edge_cells_y, &
                                 fine_mesh_edge_cells_x,     &
                                 fine_mesh_edge_cells_y,     &
                                 periodic_x, periodic_y,     &
                                 domain_size, domain_centre, &
                                 target_mesh_names,          &
                                 target_edge_cells_x,        &
                                 target_edge_cells_y,        &
                                 rotate_mesh,                &
                                 target_north_pole,          &
                                 target_null_island )        &
                                 result( self )

  implicit none

  class(reference_element_type), intent(in) :: reference_element

  character(str_def), intent(in) :: mesh_name
  integer(i_def),     intent(in) :: geometry
  integer(i_def),     intent(in) :: topology
  integer(i_def),     intent(in) :: coord_sys
  integer(i_def),     intent(in) :: edge_cells_x, edge_cells_y
  integer(i_def),     intent(in) :: fine_mesh_edge_cells_x
  integer(i_def),     intent(in) :: fine_mesh_edge_cells_y
  logical(l_def),     intent(in) :: periodic_x, periodic_y
  real(r_def),        intent(in) :: domain_size(2)
  real(r_def),        intent(in) :: domain_centre(2)

  logical,            optional, intent(in) :: rotate_mesh
  real(r_def),        optional, intent(in) :: target_north_pole(2)
  real(r_def),        optional, intent(in) :: target_null_island(2)

  character(str_def), optional, intent(in) :: target_mesh_names(:)
  integer(i_def),     optional, intent(in) :: target_edge_cells_x(:)
  integer(i_def),     optional, intent(in) :: target_edge_cells_y(:)

  type( gen_planar_type ) :: self

  character(str_long) :: target_mesh_names_str
  character(str_long) :: target_edge_cells_x_str
  character(str_long) :: target_edge_cells_y_str
  character(str_def)  :: rchar_x
  character(str_def)  :: rchar_y
  character(str_def)  :: lchar_periodic_x
  character(str_def)  :: lchar_periodic_y
  character(str_def)  :: lchar_coord_sys
  character(str_def)  :: lchar_geometry
  character(str_def)  :: lchar_topology

  character(str_def)  :: lon_str
  character(str_def)  :: lat_str
  character(str_def)  :: logic_str
  character(str_def)  :: temp_str

  character(str_def)  :: domain_size_str
  character(str_def)  :: domain_centre_str

  integer(i_def)      :: nodes_x
  integer(i_def)      :: nodes_y
  integer(i_def)      :: i
  integer(i_def)      :: n_edges
  integer(i_def)      :: remainder = 0_i_def

  integer(i_def) :: min_cells_x
  integer(i_def) :: min_cells_y

  ! At present this mesh generator strategy only supports 2d quad elements
  ! i.e. cube elements.
  select type(reference_element)
    type is (reference_cube_type)
      ! Carry on.
    class default
      call log_event( PREFIX//'Un-supported reference element type. ' // &
                      'Use reference_cube_type.', LOG_LEVEL_ERROR )
  end select

  self%nodes_per_face = reference_element%get_number_2d_vertices()
  self%edges_per_face = reference_element%get_number_2d_edges()
  self%nodes_per_edge = reference_element%get_number_verts_per_edge()

  ! There are a maximum of 4 faces around a node in this type of mesh.
  self%max_num_faces_per_node = 4

  min_cells_x = 1
  min_cells_y = 1
  if (periodic_x) min_cells_x = 2
  if (periodic_y) min_cells_y = 2

  if ( edge_cells_x < min_cells_x .or. &
       edge_cells_y < min_cells_y ) then
    call log_event( PREFIX//"For chosen periodic bounds:.", LOG_LEVEL_INFO )
    write(log_scratch_space,'(A,I0)') PREFIX//'minimum edge_cells_x = ', min_cells_x
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write(log_scratch_space,'(A,I0)') PREFIX//'minimum edge_cells_y = ', min_cells_y
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    call log_event( PREFIX//"Invalid dimension choices.", LOG_LEVEL_ERROR )
  end if

  self%mesh_name    = trim(mesh_name)
  self%geometry     = geometry
  self%topology     = topology
  self%edge_cells_x = edge_cells_x
  self%edge_cells_y = edge_cells_y
  self%fine_mesh_edge_cells_x = fine_mesh_edge_cells_x
  self%fine_mesh_edge_cells_y = fine_mesh_edge_cells_y

  self%nmaps          = 0_i_def
  self%periodic_xy(1) = periodic_x
  self%periodic_xy(2) = periodic_y
  self%coord_sys      = coord_sys
  self%domain_size    = domain_size
  self%domain_centre  = domain_centre
  self%domain_extents(:,:) = rmdi

  if ( ANY(self%domain_size(:) <= 0.0_r_def) ) then
    call log_event( PREFIX//" domain size values must be > 0.0", &
                    LOG_LEVEL_ERROR )
  end if

  if (present(rotate_mesh)) then
    self%rotate_mesh = rotate_mesh
  else
    self%rotate_mesh = .false.
  end if

  if (self%rotate_mesh) then
    if ( .not. (geometry  == geometry_spherical .and. &
                coord_sys == coord_sys_ll ) ) then
      write(log_scratch_space,'(A)')                             &
         'Rotated meshes are only supported for meshes with ' // &
         'spherical geometry and coordinates.'
      call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
    end if
  end if

  select case (self%coord_sys)
  case(coord_sys_xyz)
    self%dx = self%domain_size(1) / self%edge_cells_x
    self%dy = self%domain_size(2) / self%edge_cells_y
    self%coord_units_x = 'm'
    self%coord_units_y = 'm'

  case(coord_sys_ll)
    ! The namelist inputs were in degrees, so convert
    ! and store them as radians.
    self%domain_size   = degrees_to_radians * self%domain_size
    self%domain_centre = degrees_to_radians * self%domain_centre

    self%dx = self%domain_size(1) / self%edge_cells_x
    self%dy = self%domain_size(2) / self%edge_cells_y

    self%rotate_mesh   = rotate_mesh
    self%coord_units_x = 'radians'
    self%coord_units_y = 'radians'

    if ( self%rotate_mesh ) then
      self%north_pole  = degrees_to_radians * target_north_pole
      self%null_island = degrees_to_radians * target_null_island
    else
      ! Default value is also given in degrees so
      ! convert to radians.
      self%north_pole  = degrees_to_radians * TRUE_NORTH_POLE_LL
      self%null_island = degrees_to_radians * TRUE_NULL_ISLAND_LL
    end if

  case default
    write(log_scratch_space,'(A,I0)') &
        'Unset coordinate system enumeration: ',self%coord_sys
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )

  end select

  ! Construct input string.
  write(lchar_periodic_x,'(L8)')    periodic_x
  write(lchar_periodic_y,'(L8)')    periodic_y
  write(lchar_coord_sys, '(A)')     trim(key_from_coord_sys(self%coord_sys))
  write(lchar_geometry,  '(A)')     trim(key_from_geometry(self%geometry))
  write(lchar_topology,  '(A)')     trim(key_from_topology(self%topology))


  write(rchar_x, '(F10.2)') domain_size(1)
  write(rchar_y, '(F10.2)') domain_size(2)
  write(domain_size_str, '(A)') '['//trim(adjustl(rchar_x))// &
                                ','//trim(adjustl(rchar_y))//']'
  write(rchar_x, '(F10.2)') domain_centre(1)
  write(rchar_y, '(F10.2)') domain_centre(2)
  write(domain_centre_str, '(A)') '['//trim(adjustl(rchar_x))// &
                                  ','//trim(adjustl(rchar_y))//']'


  write(self%constructor_inputs,'(A,2(A,I0),(A))')         &
    'geometry='  // trim(adjustl(lchar_geometry))//  ';'// &
    'topology='  // trim(adjustl(lchar_topology))//  ';'// &
    'coord_sys=' // trim(adjustl(lchar_coord_sys))// ';',  &
    'edge_cells_x=', self%edge_cells_x,              ';'// &
    'edge_cells_y=', self%edge_cells_y,              ';'// &
    'periodic_x='//trim(adjustl(lchar_periodic_x))// ';'// &
    'periodic_y='//trim(adjustl(lchar_periodic_y))// ';'// &
    'domain_size='//trim(domain_size_str)//          ';'// &
    'domain_centre='//trim(domain_centre_str)



  if (self%coord_sys == coord_sys_ll) then

    ! Append rotate_mesh.
    write(logic_str,'(L8)') self%rotate_mesh
    write(temp_str,'(A)') 'rotate_mesh='//trim(adjustl(logic_str))
    write(self%constructor_inputs,'(A)') &
        trim(self%constructor_inputs) // ';' // trim(temp_str)

    if (self%rotate_mesh) then

      ! Append target pole coordinates.
      write(lon_str,'(F10.2)') target_north_pole(1)
      write(lat_str,'(F10.2)') target_north_pole(2)
      write(temp_str,'(A)')                &
          'north_pole=[' //               &
          trim(adjustl(lon_str)) // ',' // &
          trim(adjustl(lat_str)) // ']'

      write(self%constructor_inputs,'(A)') &
          trim(self%constructor_inputs) // ';' // trim(temp_str)

      ! Append null island coordinates.
      write(lon_str,'(F10.2)') target_null_island(1)
      write(lat_str,'(F10.2)') target_null_island(2)
      write(temp_str,'(A)')                &
          'null_island=[' //               &
          trim(adjustl(lon_str)) // ',' // &
          trim(adjustl(lat_str)) // ']'
      write(self%constructor_inputs,'(A)') &
          trim(self%constructor_inputs) // ';' // trim(temp_str)

    end if ! rotate_mesh

  end if ! coords_sys_ll


  ! Initialise as a plane with cyclic boundaries.
  nodes_x = self%edge_cells_x
  nodes_y = self%edge_cells_y
  n_edges = 2 * self%edge_cells_x * self%edge_cells_y

  ! Modify properties based on any non-periodic axes.
  if (.not. self%periodic_xy(1)) then
    nodes_x = nodes_x + 1
    n_edges = n_edges + self%edge_cells_y
  end if

  if (.not. self%periodic_xy(2)) then
    nodes_y = nodes_y + 1
    n_edges = n_edges + self%edge_cells_x
  end if

  ! Update properties for the resulting object.
  self%n_nodes = nodes_x * nodes_y
  self%n_edges = n_edges
  self%n_faces = self%edge_cells_x * self%edge_cells_y

  allocate(self%edges_on_cell(self%edges_per_face, self%n_faces))
  allocate(self%verts_on_edge(self%nodes_per_edge, self%n_edges))

  self%edges_on_cell = imdi
  self%verts_on_edge = imdi

  if ( present(target_edge_cells_x) .and. &
       present(target_edge_cells_y) .and. &
       present(target_mesh_names) ) then

    if ( size(target_edge_cells_x) == size(target_mesh_names) .and. &
         size(target_edge_cells_y) == size(target_mesh_names) ) then

      self%nmaps = size(target_mesh_names)
      allocate(self%target_mesh_names(self%nmaps))
      allocate(self%target_edge_cells_x(self%nmaps))
      allocate(self%target_edge_cells_y(self%nmaps))


      self%target_mesh_names(:)   = target_mesh_names(:)
      self%target_edge_cells_x(:) = target_edge_cells_x(:)
      self%target_edge_cells_y(:) = target_edge_cells_y(:)

      do i=1, self%nmaps

        ! Check that mesh is not being mapped to itself.
        if (self%edge_cells_x == target_edge_cells_x(i) .and. &
            self%edge_cells_y == target_edge_cells_y(i)) then
          write(log_scratch_space, '(A)') &
               'Invalid target while attempting to map mesh to itself'
          call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
        end if

        ! for x-axis.
        if (target_edge_cells_x(i) < self%edge_cells_x) then
          remainder = mod(self%edge_cells_x, target_edge_cells_x(i))
          write(log_scratch_space,'(2(A,I0),A)')             &
              'Target edge_cells_x[',target_edge_cells_x(i), &
              '] must be a factor of source edge_cells_x[',  &
              self%edge_cells_x, ']'

        else if (target_edge_cells_x(i) > self%edge_cells_x) then
          remainder = mod(target_edge_cells_x(i), self%edge_cells_x)
          write(log_scratch_space,'(2(A,I0),A)')              &
               'Source edge_cells_x[',target_edge_cells_x(i), &
               '] must be a factor of target edge_cells_x[',  &
               self%edge_cells_x, ']'
        end if

        if (remainder == 0_i_def) then
          self%target_edge_cells_x(i) = target_edge_cells_x(i)
        else
          call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
        end if

        ! for y-axis.
        if (target_edge_cells_y(i) < self%edge_cells_y) then
          remainder = mod(self%edge_cells_y, target_edge_cells_y(i))
          write(log_scratch_space,'(2(A,I0),A)')               &
               'Target edge_cells_y[',target_edge_cells_y(i),  &
               '] must be a factor of source edge_cells_y[',   &
               self%edge_cells_y, ']'
        else if (target_edge_cells_y(i) > self%edge_cells_y) then
          remainder = mod(target_edge_cells_y(i), self%edge_cells_y)
          write(log_scratch_space,'(2(A,I0),A)')               &
               'Source edge_cells_y[',target_edge_cells_y(i),  &
               '] must be a factor of target edge_cells_y[',   &
               self%edge_cells_y, ']'
        end if

        if (remainder == 0_i_def) then
          self%target_edge_cells_y(i) = target_edge_cells_y(i)
        else
          call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
        end if

      end do

    else

      write(log_scratch_space,'(A)') &
           'All optional array inputs for target meshes must be of the same length.'
      call log_event(trim(log_scratch_space),  LOG_LEVEL_ERROR)

    end if

    write(target_mesh_names_str, '(A)')   "'"//trim(adjustl(target_mesh_names(1)))//"'"
    write(target_edge_cells_x_str,'(I0)') target_edge_cells_x(1)
    write(target_edge_cells_y_str,'(I0)') target_edge_cells_y(1)
    if (size(target_mesh_names) > 1) then
      do i=2, self%nmaps
        write(target_mesh_names_str,'(A)')                  &
            trim(adjustl(target_mesh_names_str)) // ",'" // &
            trim(adjustl(target_mesh_names(i)))  // "'"
        write(target_edge_cells_x_str,'(A,I0)')             &
            trim(adjustl(target_edge_cells_x_str)) // ',',  &
            target_edge_cells_x(i)
        write(target_edge_cells_y_str,'(A,I0)')             &
            trim(adjustl(target_edge_cells_y_str)) // ',',  &
            target_edge_cells_y(i)
      end do
    end if


    write(target_mesh_names_str,'(A)')    &
        'target_mesh_names=['//trim(adjustl(target_mesh_names_str))//']'
    write(target_edge_cells_x_str,'(A)')  &
        'target_edge_cells_x=['//trim(adjustl(target_edge_cells_x_str))//']'
    write(target_edge_cells_y_str,'(A)')  &
        'target_edge_cells_y=['//trim(adjustl(target_edge_cells_y_str))//']'


    write(self%constructor_inputs,'(A)')        &
        trim(self%constructor_inputs) // ';' // &
        trim(target_mesh_names_str)   // ';' // &
        trim(target_edge_cells_x_str) // ';' // &
        trim(target_edge_cells_y_str)

  end if

  return
end function gen_planar_constructor


!-------------------------------------------------------------------------------
!> @brief   Calculates mesh cell adjacency. (private subroutine).
!> @details Allocates and populates the instance's cell_next(:,:) array
!>          with the id of each cell to which the index cell is adjacent.
!>
!-------------------------------------------------------------------------------
subroutine calc_adjacency(self)

  implicit none

  class(gen_planar_type), intent(inout) :: self

  integer(i_def) :: ncells_x, ncells_y, ncells, nedges_cell
  integer(i_def) :: cell, astat

  integer(i_def) :: border_types(nborders)

  ncells_x    = self%edge_cells_x
  ncells_y    = self%edge_cells_y
  ncells      = ncells_x * ncells_y
  nedges_cell = self%edges_per_face

  allocate( self%cell_next(nedges_cell, ncells ), stat=astat)
  if (astat /= 0_i_def) then
    call log_event( PREFIX//"Failure to allocate cell_next.", &
                    LOG_LEVEL_ERROR )
  end if

  if (self%periodic_xy(1)) then
    border_types(W) = periodic
    border_types(E) = periodic
  else
    border_types(W) = closed
    border_types(E) = closed
  end if

  if (self%periodic_xy(2)) then
    border_types(S) = periodic
    border_types(N) = periodic
  else
    border_types(S) = closed
    border_types(N) = closed
  end if

  self%cell_next = imdi

  !======================================
  allocate( self%north_cells (ncells_x) )
  allocate( self%south_cells (ncells_x) )
  allocate( self%east_cells  (ncells_y) )
  allocate( self%west_cells  (ncells_y) )

  ! Capture cells on edge of domain.
  do cell=1, ncells_x
    self%north_cells(cell) = cell
    self%south_cells(cell) = cell + (ncells_y-1)*ncells_x
  end do

  do cell=1, ncells_y
    self%east_cells(cell) = cell*ncells_x
    self%west_cells(cell) = (cell-1)*ncells_x + 1
  end do

  !======================================
  do cell=1, ncells
    ! Cell default values.
    self%cell_next(W, cell) = cell - 1
    self%cell_next(S, cell) = cell + ncells_x
    self%cell_next(E, cell) = cell + 1
    self%cell_next(N, cell) = cell - ncells_x
  end do

  call set_border_adjacency( self, W, self%west_cells,  border_types(W) )
  call set_border_adjacency( self, S, self%south_cells, border_types(S) )
  call set_border_adjacency( self, E, self%east_cells,  border_types(E) )
  call set_border_adjacency( self, N, self%north_cells, border_types(N) )

end subroutine calc_adjacency


!-------------------------------------------------------------------------------
!> @brief   Calculates the adjacency of border cells.
!>          (private subroutine).
!> @details Special attention is required for adjacency on border cells dependng
!>          on whether the domain is periodic in the x and y directions.
!>
!> @param[in]  edge_index    The index of the edge to set on the border cells.
!> @param[in]  border_cells  Array of cell ids of those cells on the given border.
!> @param[in]  border_type   Specifies if the border is periodic of not using the
!>                           PERIODIC or CLOSED parameters.
!-------------------------------------------------------------------------------
subroutine set_border_adjacency( self, edge_index, border_cells, border_type )

  implicit none

  type(gen_planar_type), intent(inout) :: self

  integer(i_def), intent(in) :: edge_index
  integer(i_def), intent(in) :: border_cells(:)
  integer(i_def), intent(in) :: border_type

  integer(i_def) :: i
  integer(i_def) :: ncells_x
  integer(i_def) :: ncells_y
  integer(i_def) :: ncells_border
  integer(i_def) :: cell_id

  ncells_x = self%edge_cells_x
  ncells_y = self%edge_cells_y
  ncells_border = size(border_cells)

  select case(border_type)
  case(CLOSED)
    do i=1, ncells_border
      self%cell_next(edge_index, border_cells(i)) = VOID_ID
    end do

  case(PERIODIC)
    select case(edge_index)
    case(N)
      do i=1, ncells_border
        cell_id = border_cells(i)
        self%cell_next(edge_index, cell_id) = (ncells_y-1)*ncells_x + cell_id
      end do

    case(E)
      do i=1, ncells_border
        cell_id = border_cells(i)
        self%cell_next(edge_index, cell_id) = border_cells(i) - (ncells_x-1)
      end do

    case(S)
      do i=1, ncells_border
        cell_id = border_cells(i)
        self%cell_next(edge_index, cell_id) = cell_id - (ncells_y-1)*ncells_x
      end do

    case(W)
      do i=1, ncells_border
        cell_id = border_cells(i)
        self%cell_next(edge_index, cell_id) = cell_id + ncells_x-1
      end do
    end select

  case default
    ! (Non-Periodic), As a namelist input, absent logicals are
    ! default .false., so it is consistence to make the default
    ! border type as non-periodic.
    do i=1, ncells_border
      self%cell_next(edge_index, border_cells(i)) = VOID_ID
    end do

  end select

  return
end subroutine set_border_adjacency


!-------------------------------------------------------------------------------
!> @brief   For each cell, calculates the four vertices whih comprise it.
!>          (private subsroutine)
!> @details Allocates and populates the instance's mesh(:,:) array with
!>          the vertices which form each cell.
!>
!-------------------------------------------------------------------------------
subroutine calc_face_to_vert(self)

  implicit none

  class(gen_planar_type), intent(inout) :: self

  integer(i_def), allocatable :: verts_on_cell(:,:)

  integer(i_def) :: edge_cells_x, edge_cells_y, ncells
  integer(i_def) :: y, vert, cell, astat
  integer(i_def) :: nverts_cell
  integer(i_def) :: node_id

  nverts_cell  = self%nodes_per_face
  edge_cells_x = self%edge_cells_x
  edge_cells_y = self%edge_cells_y
  ncells       = self%n_faces

  allocate(verts_on_cell(nverts_cell, ncells), stat=astat)
  if (astat /= 0) then
    call log_event( PREFIX//"Failure to allocate verts_on_cell array.", &
                    LOG_LEVEL_ERROR )
  end if

  verts_on_cell(:,:) = imdi

  !=====================================================================
  ! FIRST ROW of panel.
  !=====================================================================

  ! First cell of first row.
  !------------------------
  cell     = 1
  node_id  = 1

  do vert = 1, nverts_cell
    verts_on_cell(vert, cell) = node_id
    node_id = node_id+1
  end do


  ! East neighbour.
  if (self%cell_next(E, cell) /= VOID_ID ) then
    verts_on_cell(NW , self%cell_next(E, cell)) = verts_on_cell(NE, cell)
    verts_on_cell(SW , self%cell_next(E, cell)) = verts_on_cell(SE, cell)
  end if

  ! South neighbour.
  if (self%cell_next(S, cell) /= VOID_ID ) then
    verts_on_cell(NW , self%cell_next(S, cell)) = verts_on_cell(SW, cell)
    verts_on_cell(NE , self%cell_next(S, cell)) = verts_on_cell(SE, cell)
  end if


  ! Inner cells of first row.
  !--------------------------
  if (edge_cells_x > 2) then
    do cell = 2, edge_cells_x-1
      verts_on_cell(SE, cell) = node_id
      verts_on_cell(NE, cell) = node_id+1
      node_id = node_id + 2

      ! East neighbour.
      verts_on_cell(NW , self%cell_next(E, cell)) = verts_on_cell(NE, cell)
      verts_on_cell(SW , self%cell_next(E, cell)) = verts_on_cell(SE, cell)

      ! South neighbour.
      if (self%cell_next(S, cell) /= VOID_ID ) then
        verts_on_cell(NW , self%cell_next(S, cell)) = verts_on_cell(SW, cell)
        verts_on_cell(NE , self%cell_next(S, cell)) = verts_on_cell(SE, cell)
      end if
    end do
  end if

  ! Last cell of first row
  !-------------------------
  if (edge_cells_x > 1) then
    cell = edge_cells_x
    if (self%periodic_xy(1)) then
      ! Copy node id from left edge of domain.
      verts_on_cell(SE, cell) = verts_on_cell(SW, cell-edge_cells_x+1)
      verts_on_cell(NE, cell) = verts_on_cell(NW, cell-edge_cells_x+1)
    else
      verts_on_cell(SE, cell) = node_id
      verts_on_cell(NE, cell) = node_id + 1
      node_id = node_id + 2
    end if

    ! South neighbour.
    if (self%cell_next(S, cell) /= VOID_ID ) then
      verts_on_cell(NW , self%cell_next(S, cell)) = verts_on_cell(SW, cell)
      verts_on_cell(NE , self%cell_next(S, cell)) = verts_on_cell(SE, cell)
    end if
  end if
  ! END FIRST ROW of panel.
  !=====================================================================


  !=====================================================================
  ! INNER ROWS of panel.
  !
  if (edge_cells_y > 2) then

    do y = 1, edge_cells_y-2
      ! First cell in row.
      !-------------------
      cell = (y*edge_cells_x) + 1
      verts_on_cell(SW, cell) = node_id
      verts_on_cell(SE, cell) = node_id+1
      node_id = node_id+2

      ! South neighbour.
      verts_on_cell(NW , self%cell_next(S, cell)) = verts_on_cell(SW, cell)
      verts_on_cell(NE , self%cell_next(S, cell)) = verts_on_cell(SE, cell)

      ! East neighbour.
      if (edge_cells_x > 1) then
        verts_on_cell(SW , self%cell_next(E, cell)) = verts_on_cell(SE, cell)
      end if

      ! Inner cells of row.
      !--------------------
      if (edge_cells_x > 2) then
        do cell = y*edge_cells_x+2, (y+1)*edge_cells_x-1
          verts_on_cell(SE, cell) = node_id
          node_id = node_id+1

          ! South neighbour.
          verts_on_cell(NW, self%cell_next(S, cell)) = verts_on_cell(SW, cell)
          verts_on_cell(NE, self%cell_next(S, cell)) = verts_on_cell(SE, cell)

          ! East neighbour.
          verts_on_cell(SW , self%cell_next(E, cell)) = verts_on_cell(SE, cell)
        end do
      end if

      ! Last cell of row.
      !-------------------
      if (edge_cells_x > 1) then
        cell = (y+1)*edge_cells_x
        if (self%periodic_xy(1)) then
          verts_on_cell(SE, cell) = verts_on_cell(SW, self%cell_next(E, cell))
        else
          verts_on_cell(SE, cell) = node_id
          node_id = node_id+1
        end if

        ! South neighbour.
        verts_on_cell(NW, self%cell_next(S, cell)) = verts_on_cell(SW, cell)
        verts_on_cell(NE, self%cell_next(S, cell)) = verts_on_cell(SE, cell)
      end if

    end do
  end if

  !========================
  ! BOTTOM EDGE of Panel.
  !========================
  if (self%periodic_xy(2)) then

    ! Copy from top edge row.
    do cell = ncells-edge_cells_x+1, ncells
      verts_on_cell(SW, cell) = verts_on_cell(NW, self%cell_next(S, cell) )
      verts_on_cell(SE, cell) = verts_on_cell(NE, self%cell_next(S, cell) )
    end do

    node_id = node_id - 1

  else

    if (edge_cells_y > 1) then
      ! Always do first cell in bottom row.
      ! First cell in bottom row.
      cell = ncells-edge_cells_x+1
      verts_on_cell(SW, cell) = node_id
      verts_on_cell(SE, cell) = node_id+1
      node_id = node_id+2

      ! Inner cells in bottom row.
      do cell = ncells-edge_cells_x+2, ncells-1
        verts_on_cell(SW, cell) = verts_on_cell(SE, self%cell_next(W, cell))
        verts_on_cell(SE, cell) = node_id
        node_id = node_id+1
      end do

      ! Last cell in bottom row.
      cell = ncells
      if (edge_cells_x > 1) then
        verts_on_cell(SW, cell) = verts_on_cell(SE, self%cell_next(W, cell))
        if (self%periodic_xy(1)) then
          verts_on_cell(SE, cell) = verts_on_cell(SW, self%cell_next(E, cell))
          node_id = node_id - 1
        else
          verts_on_cell(SE, cell) = node_id
        end if
      end if
    end if
  end if ! self%periodic_xy(2)

  call move_alloc(verts_on_cell, self%verts_on_cell)

  return
end subroutine calc_face_to_vert


!-------------------------------------------------------------------------------
!> @brief   Calculates the edges which are found on each cell and the
!>          pair of vertices which are found on each edge.(private subroutine)
!> @details Allocates and populates both the edges_on_cell and
!>          verts_on_edge arrays for the instance.
!>
!-------------------------------------------------------------------------------
subroutine calc_edges(self)

  implicit none

  class(gen_planar_type), intent(inout)  :: self

  integer(i_def) :: edge_cells_x, edge_cells_y, ncells
  integer(i_def) :: cell
  integer(i_def) :: edge_id

  edge_cells_x = self%edge_cells_x
  edge_cells_y = self%edge_cells_y
  ncells = self%edge_cells_x * self%edge_cells_y

  cell     = 1
  edge_id  = 1

  ! The assumed working layout for the folowing code with respect to
  ! cell/panel is
  !
  ! ID Origin     Top (North)
  !    (NW) O--------------------O (NE)
  !         |                    +
  !         |   edge/vertex id   +
  !         |  + anti-clockwise  +
  !   Left  |   \   Numbering    + Right
  !  (West) |    \               + (East)
  !         |     +------->      +
  !         |                    +
  !         |                    +
  !    (SE) O+-+-+-+-+-+-+-+-+-+-O (SE)
  !             Bottom (South)
  !
  ! This is a 2D grid, and compass points are used to orientate
  ! around the panel/cell. (top, left, right, bottom) could
  ! have been used, though top and bottom could have been
  ! confused with 3D cells. However, for readability,
  ! the terms (left, bottom, top, right) are used interchangably
  ! with the cardinal compass directions in the above diagram.
  ! It should also be noted that the cardinal compass directions
  ! are only used as an alternative way to reference the elements
  ! on the the panel/cells, i.e. the Top (North) edge of a cell
  ! may not actually be aligned in the same direction as geographic
  ! north.
  !
  ! When considered as a panel, the Eastern and Southern edges
  ! are taken as the `ghost` edges/vertices in the case of any
  ! periodicity.
  !
  ! So that vertices on edges are consistent, listed vertices
  ! connected to edges will go from N-S, or W-E.

  !--------------------------------------------
  ! Top panel row, cell@left panel edge (ID:1)
  !--------------------------------------------
  ! Cell western edge.
  self%edges_on_cell(W, cell)      = edge_id
  self%verts_on_edge(1, edge_id)   = self%verts_on_cell(NW, cell)
  self%verts_on_edge(2, edge_id)   = self%verts_on_cell(SW, cell)

  ! Cell southern edge.
  self%edges_on_cell(S, cell)  = edge_id+1
  self%verts_on_edge(1, edge_id+1) = self%verts_on_cell(SW, cell)
  self%verts_on_edge(2, edge_id+1) = self%verts_on_cell(SE, cell)

  ! Cell eastern edge.
  self%edges_on_cell(E, cell)  = edge_id+2
  self%verts_on_edge(1, edge_id+2) = self%verts_on_cell(NE, cell)
  self%verts_on_edge(2, edge_id+2) = self%verts_on_cell(SE, cell)

  ! Cell northern edge.
  self%edges_on_cell(N, cell)  = edge_id+3
  self%verts_on_edge(1, edge_id+3) = self%verts_on_cell(NW, cell)
  self%verts_on_edge(2, edge_id+3) = self%verts_on_cell(NE, cell)

  edge_id = edge_id+4

  !-----------------------------------------------------------
  ! Top panel row, remaining cells, i.e. IDs = 2:edge_cells_x.
  !-----------------------------------------------------------
  if (edge_cells_x > 1) then
    do cell = 2, edge_cells_x

      ! Cell western edge.
      self%edges_on_cell(W, cell) = self%edges_on_cell(E,self%cell_next(W,cell))

      ! Cell southern edge.
      self%edges_on_cell(S, cell) = edge_id
      self%verts_on_edge(1, edge_id)  = self%verts_on_cell(SW, cell)
      self%verts_on_edge(2, edge_id)  = self%verts_on_cell(SE, cell)
      edge_id = edge_id + 1

      !-------------------------
      if (cell == edge_cells_x) then
        ! This cell is on the right-hand edge of the panel.

        ! For a planar panel, if the cell is on the right-hand edge
        ! of the panel then the cell's eastern neighbour is actually
        ! the left-most cell on this row (because of the periodicity).
        ! At this point, left-most cell on this row has already
        ! had its edge ids assigned.

        ! In other words, the eastern edge of this cell exists on
        ! another cell where it has already been assigned an id.

        ! Cell eastern edge.
        if (self%periodic_xy(1)) then
          self%edges_on_cell(E, cell) = self%edges_on_cell(W, self%cell_next(E,cell))

        else
          self%edges_on_cell(E, cell) = edge_id
          self%verts_on_edge(1, edge_id) = self%verts_on_cell(NE, cell)
          self%verts_on_edge(2, edge_id) = self%verts_on_cell(SE, cell)
          edge_id = edge_id + 1
        end if

        ! Cell northern edge.
        self%edges_on_cell(N, cell)  = edge_id
        self%verts_on_edge(1, edge_id) = self%verts_on_cell(NW, cell)
        self%verts_on_edge(2, edge_id) = self%verts_on_cell(NE, cell)
        edge_id = edge_id+1

      else

        ! Cell eastern edge.
        self%edges_on_cell(E, cell)  = edge_id
        self%verts_on_edge(1, edge_id) = self%verts_on_cell(NE, cell)
        self%verts_on_edge(2, edge_id) = self%verts_on_cell(SE, cell)
        edge_id = edge_id+1

        ! Cell northern edge.
        self%edges_on_cell(N, cell)  = edge_id
        self%verts_on_edge(1, edge_id) = self%verts_on_cell(NW, cell)
        self%verts_on_edge(2, edge_id) = self%verts_on_cell(NE, cell)
        edge_id = edge_id+1

      end if

    end do
  end if

  !-----------------------------------------
  ! Internal panel rows.
  !-----------------------------------------
  if (edge_cells_y > 2) then
    do cell = edge_cells_x+1, ncells-edge_cells_x

      if (mod(cell,edge_cells_x) == 1) then
        ! This cell is on the left-hand panel edge.
        !-------------------------------------------

        ! Cell western edge.
        self%edges_on_cell(W, cell)  = edge_id
        self%verts_on_edge(1, edge_id)   = self%verts_on_cell(NW, cell)
        self%verts_on_edge(2, edge_id)   = self%verts_on_cell(SW, cell)

        ! Cell southern edge.
        self%edges_on_cell(S, cell)  = edge_id+1
        self%verts_on_edge(1, edge_id+1) = self%verts_on_cell(SW, cell)
        self%verts_on_edge(2, edge_id+1) = self%verts_on_cell(SE, cell)

        ! Cell eastern edge.
        self%edges_on_cell(E, cell)  = edge_id+2
        self%verts_on_edge(1, edge_id+2) = self%verts_on_cell(NE, cell)
        self%verts_on_edge(2, edge_id+2) = self%verts_on_cell(SE, cell)

        edge_id = edge_id+3

      else if (mod(cell,edge_cells_x) == 0) then
        ! This cell is on the right-hand panel edge.
        !--------------------------------------------
        if (edge_cells_x==1) then
          self%edges_on_cell(W, cell)  = edge_id
          self%verts_on_edge(1, edge_id)   = self%verts_on_cell(NW, cell)
          self%verts_on_edge(2, edge_id)   = self%verts_on_cell(SW, cell)
          edge_id = edge_id + 1
        else
          ! Cell western edge.
          self%edges_on_cell(W, cell) = self%edges_on_cell(E,self%cell_next(W,cell))
        end if

        ! Cell southern edge.
        self%edges_on_cell(S, cell) = edge_id
        self%verts_on_edge(1, edge_id)  = self%verts_on_cell(SW, cell)
        self%verts_on_edge(2, edge_id)  = self%verts_on_cell(SE, cell)

        edge_id = edge_id+1

        ! Cell eastern edge.
        if (self%periodic_xy(1)) then
          self%edges_on_cell(E, cell) = self%edges_on_cell(W, self%cell_next(E, cell))
        else

          self%edges_on_cell(E, cell) = edge_id
          self%verts_on_edge(1, edge_id)  = self%verts_on_cell(NE, cell)
          self%verts_on_edge(2, edge_id)  = self%verts_on_cell(SE, cell)
          edge_id = edge_id+1
        end if

      else
        ! This cell is internal to the panel.
        !-------------------------------------------

        ! Cell western edge.
        self%edges_on_cell(W, cell) = self%edges_on_cell(E,self%cell_next(W,cell))

        ! Cell southern edge.
        self%edges_on_cell(S, cell) = edge_id
        self%verts_on_edge(1, edge_id)  = self%verts_on_cell(SW, cell)
        self%verts_on_edge(2, edge_id)  = self%verts_on_cell(SE, cell)

        ! Cell eastern edge.
        self%edges_on_cell(E, cell)  = edge_id+1
        self%verts_on_edge(1, edge_id+1) = self%verts_on_cell(NE, cell)
        self%verts_on_edge(2, edge_id+1) = self%verts_on_cell(SE, cell)

        edge_id=edge_id+2

      end if

      ! Cell northern edge.
      ! The northern edges on these cells exist on the cells in
      ! the row above/previous to this one. Those edges have
      ! already been assigned ids.
      self%edges_on_cell(N, cell) = self%edges_on_cell(S,self%cell_next(N,cell))

    end do ! Panel inner rows
  end if

  if (edge_cells_y > 1) then
    ! Panel bottom row.
    do cell = ncells-edge_cells_x+1, ncells

      if (mod(cell,edge_cells_x) == 1 ) then

        ! This cell is on the left-hand panel edge.
        !-------------------------------------------

        ! Cell western edge.
        self%edges_on_cell(W, cell) = edge_id
        self%verts_on_edge(1, edge_id)  = self%verts_on_cell(NW, cell)
        self%verts_on_edge(2, edge_id)  = self%verts_on_cell(SW, cell)
        edge_id = edge_id+1

        ! Cell southern edge.
        if (self%periodic_xy(2)) then
          ! For a planar panel, if the cell is on the bottom edge
          ! of the panel then the cell's southern neighbour is actually
          ! the on the top-edge of the panel (because of the periodicity).
          ! At this point, all cells on the top-edge of the panel have
          ! had their edge ids assigned.

          ! In other words, the southern edge of this cell exists on
          ! another cell where it has already been assigned an id.
          self%edges_on_cell(S, cell) = self%edges_on_cell(N, self%cell_next(S, cell))
        else
          self%edges_on_cell(S, cell) = edge_id
          self%verts_on_edge(1, edge_id)  = self%verts_on_cell(SW, cell)
          self%verts_on_edge(2, edge_id)  = self%verts_on_cell(SE, cell)
          edge_id = edge_id+1
        end if

        ! Cell eastern edge.
        self%edges_on_cell(E, cell)  = edge_id
        self%verts_on_edge(1, edge_id) = self%verts_on_cell(NE, cell)
        self%verts_on_edge(2, edge_id) = self%verts_on_cell(SE, cell)
        edge_id = edge_id+1

      else if (mod(cell,edge_cells_x) == 0) then
        ! This cell is on the right-hand panel edge.
        !--------------------------------------------
        ! Cell western edge.
        if (edge_cells_x == 1) then
          self%edges_on_cell(W, cell) = edge_id
          self%verts_on_edge(1, edge_id)  = self%verts_on_cell(NW, cell)
          self%verts_on_edge(2, edge_id)  = self%verts_on_cell(SW, cell)
          edge_id = edge_id+1
        else if (edge_cells_x > 1) then
          self%edges_on_cell(W, cell) = self%edges_on_cell(E,self%cell_next(W,cell))
        end if

        ! Cell southern edge.
        if (self%periodic_xy(2)) then
          ! For a planar panel, if the cell is on the bottom edge
          ! of the panel then the cell's southern neighbour is actually
          ! the on the top-edge of the panel (because of the periodicity).
          ! At this point, all cells on the top-edge of the panel have
          ! had their edge ids assigned.

          ! In other words, the southern edge of this cell exists on
          ! another cell where it has already been assigned an id.
          self%edges_on_cell(S, cell) = self%edges_on_cell(N, self%cell_next(S, cell))
        else
          self%edges_on_cell(S, cell) = edge_id
          self%verts_on_edge(1, edge_id)  = self%verts_on_cell(SW, cell)
          self%verts_on_edge(2, edge_id)  = self%verts_on_cell(SE, cell)
          edge_id = edge_id+1
        end if

        ! Cell eastern edge.
        if (self%periodic_xy(1)) then
          ! For a planar panel, if the cell is on the right-hand edge
          ! of the panel then the cell's eastern neighbour is actually
          ! the left-most cell on this row (because of the periodicity).
          ! At this point, left-most cell on this row has already
          ! had its edge ids assigned.

          ! In other words, the eastern edge of this cell exists on
          ! another cell where it has already been assigned an id.
          self%edges_on_cell(E, cell) = self%edges_on_cell(W, self%cell_next(E,cell))
        else
          self%edges_on_cell(E, cell) = edge_id
          self%verts_on_edge(1, edge_id)  = self%verts_on_cell(NE, cell)
          self%verts_on_edge(2, edge_id)  = self%verts_on_cell(SE, cell)
          edge_id = edge_id+1
        end if

      else

        ! This cell is an inner cell on the bottom row.
        !---------------------------------------------
        ! Cell western edge.
        self%edges_on_cell(W, cell) = self%edges_on_cell(E, self%cell_next(W,cell))

        ! Cell southern edge.
        if (self%periodic_xy(2)) then
          ! For a planar panel, if the cell is on the bottom edge
          ! of the panel then the cell's southern neighbour is actually
          ! the on the top-edge of the panel (because of the periodicity).
          ! At this point, all cells on the top-edge of the panel have
          ! had their edge ids assigned.

          ! In other words, the southern edge of this cell exists on
          ! another cell where it has already been assigned an id.
          self%edges_on_cell(S, cell) = self%edges_on_cell(N, self%cell_next(S, cell))

        else

          self%edges_on_cell(S, cell) = edge_id
          self%verts_on_edge(1, edge_id)  = self%verts_on_cell(SW, cell)
          self%verts_on_edge(2, edge_id)  = self%verts_on_cell(SE, cell)
          edge_id = edge_id+1

        end if

        ! Cell eastern edge.
        self%edges_on_cell(E, cell) = edge_id
        self%verts_on_edge(1, edge_id) = self%verts_on_cell(NE, cell)
        self%verts_on_edge(2, edge_id) = self%verts_on_cell(SE, cell)
        edge_id = edge_id+1

      end if

      ! Cell north edge.
      ! The northern edges on these cells exist on the cells in
      ! the row above/previous to this one. Those edges have
      ! already been assigned ids.
      self%edges_on_cell(N, cell) = self%edges_on_cell(S,self%cell_next(N,cell))

    end do
  end if

  return
end subroutine calc_edges

!-------------------------------------------------------------------------------
!> @brief  Define the cooordinates using a stretching function.
!> @details Assigns an (x,y) coordinate to vertices of the mesh. If the mesh is
!>         a support mesh then the coordinate calculation is based on the
!>         transform_mesh (i.e. For multigrid, the coordinates of the coarse
!>         level mesh are defined using the fine-level mesh).
!>
!-------------------------------------------------------------------------------
subroutine assign_stretched_mesh_coords(self)

  implicit none

  class(gen_planar_type), intent(inout)  :: self

  real(r_def), allocatable :: vert_coords(:,:)

  integer(i_def) :: cell, astat, row, column
  integer(i_def) :: x_ratio, y_ratio

  real(r_def) :: node_coord, stretch_coord
  real(r_def) :: reverse_node_coord
  real(r_def) :: reverse_node_coord_bottom, reverse_node_coord_top

  ! Stretch transform settings
  integer(i_def) :: axis_direction
  integer(i_def) :: n_inner, n_stretch
  real(r_def) :: outer_ends_l, outer_ends_r, inner_ends_l, inner_ends_r
  real(r_def) :: inflation

  allocate(vert_coords(2, self%n_nodes), stat=astat)
  if (astat /= 0) then
    call log_event( PREFIX//"Failure to allocate vert_coords.", &
                    LOG_LEVEL_ERROR )
  end if

  self%domain_extents(:,1) = [           0.0_r_def, -1.0_r_def*self%domain_size(2) ]
  self%domain_extents(:,2) = [ self%domain_size(1), -1.0_r_def*self%domain_size(2) ]
  self%domain_extents(:,3) = [ self%domain_size(1),            0.0_r_def ]
  self%domain_extents(:,4) = [           0.0_r_def,            0.0_r_def ]

  ! Set the ratio between the edge cells for the fine level mesh
  ! and the current multigrid mesh. This takes the value 1 if the
  ! current mesh is the fine level mesh.

  x_ratio = self%fine_mesh_edge_cells_x/self%edge_cells_x
  y_ratio = self%fine_mesh_edge_cells_y/self%edge_cells_y

  ! Loop over cells and assign coordinates based on the stretching
  ! function. The cells begin numbering in rows from NW corner of panel.
  ! Extra coordinates are assigned if it the first or last row/column to
  ! cover both the East and West (North and South) sides of the cell.

  ! --------- Longitude coordinates node_coord--------------------------------

  axis_direction = 1 ! Longitude

  call calculate_settings(axis_direction, self%fine_mesh_edge_cells_x, &
                          n_inner, n_stretch, inflation,               &
                          outer_ends_l, outer_ends_r,                  &
                          inner_ends_l, inner_ends_r)

  cell=1

  ! Loop over the fine resolution mesh rows
  do row=1, self%fine_mesh_edge_cells_y

    ! Set the longitude associated with NW corner (column = 0)
    call stretch_transform(0, node_coord, stretch_coord,  &
                           axis_direction,                &
                           n_inner, n_stretch, inflation, &
                           outer_ends_l, outer_ends_r,    &
                           inner_ends_l, inner_ends_r)

    ! Assign coordinates if this coincides with the multigrid mesh
    if (mod(row, y_ratio) == 0) then
      vert_coords(1, self%verts_on_cell(NW, cell)) &
        = node_coord

      if (row == self%fine_mesh_edge_cells_y) then

        vert_coords(1, self%verts_on_cell(SW, cell)) &
          = node_coord

      end if

      ! Loop over the fine resolution mesh columns
      do column=1, self%fine_mesh_edge_cells_x

        call stretch_transform(column,                        &
                               node_coord, stretch_coord,     &
                               axis_direction,                &
                               n_inner, n_stretch, inflation, &
                               outer_ends_l, outer_ends_r,    &
                               inner_ends_l, inner_ends_r)

        if (mod(column, x_ratio) == 0) then

          vert_coords(1, self%verts_on_cell(NE, cell)) &
            = node_coord

          if (row == self%fine_mesh_edge_cells_y) then

            vert_coords(1, self%verts_on_cell(SE, cell)) &
              = node_coord

          end if

          ! Increment the cell counter

          cell = cell + 1

        end if

      end do

    end if
  end do

  ! --------- Latitudes coordinates ---------------------------------------
  !
  ! Latitudes are a little more complicated than longitude, because the
  ! stretching function increments from bottom to top, but we need
  ! to define from top to bottom. Therefore, we take advantage of the fact
  ! that it is a symmetrical function and define the coordinates in reverse.

  axis_direction = 2 ! Latitude

  call calculate_settings(axis_direction, self%fine_mesh_edge_cells_y, &
                          n_inner, n_stretch, inflation,               &
                          outer_ends_l, outer_ends_r,                  &
                          inner_ends_l, inner_ends_r)

  ! Set the latitude of the NW corner of the last cell (column = edge_cells_y)
  call stretch_transform(self%fine_mesh_edge_cells_y,           &
                         reverse_node_coord_top, stretch_coord, &
                         axis_direction,                        &
                         n_inner, n_stretch, inflation,         &
                         outer_ends_l, outer_ends_r,            &
                         inner_ends_l, inner_ends_r)

  ! Set the latitude of the SW corner (column = 0)
  call stretch_transform(0, reverse_node_coord_bottom,  &
                         stretch_coord,                 &
                         axis_direction,                &
                         n_inner, n_stretch, inflation, &
                         outer_ends_l, outer_ends_r,    &
                         inner_ends_l, inner_ends_r)

  cell=1

  !Start the stretching function at the bottom (South)
  reverse_node_coord = reverse_node_coord_bottom

  ! Loop over the fine resolution mesh rows
  do row=1, self%fine_mesh_edge_cells_y

    ! Assign coordinates if this coincides with the multigrid mesh
    if (mod((row-1), y_ratio) == 0) then

      vert_coords(2, self%verts_on_cell(NW, cell) ) &
        = (reverse_node_coord_top &
        + (reverse_node_coord_bottom - reverse_node_coord))

      if ((row-1) == (self%fine_mesh_edge_cells_y - y_ratio)) then

        vert_coords(2, self%verts_on_cell(SW, cell)) &
          = (reverse_node_coord_top &
          + (reverse_node_coord_bottom - reverse_node_coord_top))

      end if

      ! Loop over the fine resolution mesh columns
      do column=1, self%fine_mesh_edge_cells_x

        if (mod((column-1), x_ratio) == 0) then

          vert_coords(2, self%verts_on_cell(NE, cell)) &
            = (reverse_node_coord_top &
            + (reverse_node_coord_bottom - reverse_node_coord))

          if ((row-1) == (self%fine_mesh_edge_cells_y - y_ratio)) then

            vert_coords(2, self%verts_on_cell(SE, cell)) &
              = (reverse_node_coord_top &
              + (reverse_node_coord_bottom - reverse_node_coord_top))

          end if

          ! Increment the cell counter

          cell = cell + 1

        end if

      end do

    end if

    call stretch_transform(row,                               &
                           reverse_node_coord, stretch_coord, &
                           axis_direction,                    &
                           n_inner, n_stretch, inflation,     &
                           outer_ends_l, outer_ends_r,        &
                           inner_ends_l, inner_ends_r)
  end do

  select case (self%coord_sys)

  case(coord_sys_xyz)
    self%coord_units_x = 'm'
    self%coord_units_y = 'm'

  case(coord_sys_ll)
    self%coord_units_x = 'radians'
    self%coord_units_y = 'radians'

    vert_coords = vert_coords * degrees_to_radians

  case default
    write(log_scratch_space,'(A,I0)') &
        'Unset coordinate system enumeration: ', self%coord_sys
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )

  end select

  call move_alloc(vert_coords, self%vert_coords)

  return
end subroutine assign_stretched_mesh_coords


!-------------------------------------------------------------------------------
!> @brief   Calculates the coordinates of vertices in the mesh.(private subroutine)
!> @details Assigns an (x,y) coordinate in units of dx and dy to each mesh
!>          vertex according to its Cartesian position in the mesh.
!>
!-------------------------------------------------------------------------------
subroutine calc_coords(self)

  implicit none

  class(gen_planar_type), intent(inout)  :: self

  real(r_def), allocatable :: vert_coords(:,:)

  integer(i_def) :: ncells, edge_cells_x, edge_cells_y
  integer(i_def) :: cell, astat, row, column

  real(r_def) :: x_coord
  real(r_def) :: y_coord
  real(r_def) :: offset(2)

  edge_cells_x = self%edge_cells_x
  edge_cells_y = self%edge_cells_y
  ncells = edge_cells_x*edge_cells_y

  allocate(vert_coords(2, self%n_nodes), stat=astat)
  if (astat /= 0) then
    call log_event( PREFIX//"Failure to allocate vert_coords.", &
                    LOG_LEVEL_ERROR )
  end if

!==========================================================
! Domain coordinates initial layout
! After which an offset is applied to place the
! target domain translated on an unrotated frame
! of reference
!
!     |
!   4 |                 3
! ----+---------------+---
!     |               |
!     |               |
!   1 |               | 2
!     +---------------+
!
!==========================================================

  self%domain_extents(:,1) = [           0.0_r_def, -1.0_r_def*self%domain_size(2) ]
  self%domain_extents(:,2) = [ self%domain_size(1), -1.0_r_def*self%domain_size(2) ]
  self%domain_extents(:,3) = [ self%domain_size(1),            0.0_r_def ]
  self%domain_extents(:,4) = [           0.0_r_def,            0.0_r_def ]

  offset(1) = self%domain_centre(1) - (self%domain_size(1)*0.5_r_def)
  offset(2) = self%domain_centre(2) + (self%domain_size(2)*0.5_r_def)

  ! The cells begin numbering in rows from NW corner of panel.
  cell=1
  do row=1, self%edge_cells_y
    do column=1, self%edge_cells_x
      vert_coords(1, self%verts_on_cell(NW, cell)) = (column-1) * self%dx
      vert_coords(2, self%verts_on_cell(NW, cell)) = (row-1)    * self%dy * (-1.0_r_def)
      cell=cell+1
    end do
  end do

  if (.not. self%periodic_xy(1)) then
    ! Vertices on east edge of panel.
    x_coord = edge_cells_x*self%dx
    y_coord = 0.0_r_def
    do cell=1, size(self%east_cells)
      vert_coords(1, self%verts_on_cell(NE, self%east_cells(cell))) = x_coord
      vert_coords(2, self%verts_on_cell(NE, self%east_cells(cell))) = y_coord - (self%dy*(cell-1))
    end do
  end if

  ! Vertices on south edge of panel.
  if (.not. self%periodic_xy(2)) then
    x_coord = 0.0_r_def
    y_coord = -1.0_r_def*edge_cells_y*self%dy

    do cell=1, size(self%south_cells)
      vert_coords(1, self%verts_on_cell(SW, self%south_cells(cell))) = x_coord + (self%dx*(cell-1))
      vert_coords(2, self%verts_on_cell(SW, self%south_cells(cell))) = y_coord
    end do
  end if

  ! Coords of SE panel vertex.
  if (.not. self%periodic_xy(1) .and. .not. self%periodic_xy(2)) then
    cell=ncells
    vert_coords(1, self%verts_on_cell(SE, cell)) = self%dx * self%edge_cells_x
    vert_coords(2, self%verts_on_cell(SE, cell)) = self%dy * self%edge_cells_y * (-1.0_r_def)
  end if

  vert_coords(1,:)    = vert_coords(1,:) + offset(1)
  vert_coords(2,:)    = vert_coords(2,:) + offset(2)

  self%domain_extents(1,:) = self%domain_extents(1,:) + offset(1)
  self%domain_extents(2,:) = self%domain_extents(2,:) + offset(2)

  select case (self%coord_sys)

  case(coord_sys_xyz)
    self%coord_units_x = 'm'
    self%coord_units_y = 'm'

  case(coord_sys_ll)
    self%coord_units_x = 'radians'
    self%coord_units_y = 'radians'

  case default
    write(log_scratch_space,'(A,I0)') &
        'Unset coordinate system enumeration: ', self%coord_sys
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )

  end select

  call move_alloc(vert_coords, self%vert_coords)

  return
end subroutine calc_coords


!-------------------------------------------------------------------------------
!> @brief    Get the global cell id at a specified corner of the planar domain.
!> @details  Returns the global cell id of the cell at a given corner of the
!>           domain.
!>
!> @param[in] corner      [NW|NE|SW|SE] enumerations from the reference element.
!>                        Identifies which corner cell is being requested.
!> @return    corner_gid  Global cell id of cell at requested domain corner.
!>
!-------------------------------------------------------------------------------
function get_corner_gid(self, corner) result(corner_gid)

  implicit none

  class(gen_planar_type), intent(in) :: self
  integer(i_def),         intent(in) :: corner

  integer(i_def) :: corner_gid

  corner_gid = imdi

  select case (corner)

  case(NW)
    corner_gid = self%north_cells(1)

  case(NE)
    corner_gid = self%north_cells(self%edge_cells_x)

  case(SW)
    corner_gid = self%south_cells(1)

  case(SE)
    corner_gid = self%south_cells(self%edge_cells_x)

  case default
    write(log_scratch_space,'(A,I0)') &
        'Unrecognised corner enumeration, use (NW|NE|SW|SE)'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )

  end select

  return
end function get_corner_gid

!-------------------------------------------------------------------------------
!> @brief   Calculates the mesh cell centres.(private subroutine)
!> @details The face centres for the mesh are calculated based on the current
!>          node coordinates of the node on the NW corner of the cell.
!>          The node_cordinates are assumed to be in [m] in the
!>          x,y plane. Resulting face centre coordinates are in [m].
!>
!-------------------------------------------------------------------------------
subroutine calc_cell_centres(self)

  implicit none

  class(gen_planar_type), intent(inout) :: self

  integer(i_def) :: ncells

  ! Counters.
  integer(i_def) :: cell, base_vert,i,j
  integer(i_def), parameter :: NVERTS_PER_CELL = 4
  integer(i_def) :: cell_verts(NVERTS_PER_CELL)

  ncells = NPANELS*self%edge_cells_x*self%edge_cells_y

  ! 1.0 Initialise the face centres.
  if ( .not. allocated(self%cell_coords) ) allocate( self%cell_coords(2,ncells) )
  self%cell_coords(:,:) = 0.0_r_def

  if ( self%topology == topology_non_periodic)then
    ! 2.1 for non_periodic domains, we use the standard approach of taking the
    ! mean of the coordinates at the vertices. This will give a value for the
    ! centre in the target coordinates (e.g. ll or xyz).

    do cell=1, ncells
      self%cell_coords(:,cell) = 0.0_r_def
      cell_verts(:) = self%verts_on_cell(1:NVERTS_PER_CELL, cell)

      do i=1,NVERTS_PER_CELL
        j=cell_verts(i)
        self%cell_coords(:,cell) = self%cell_coords(:,cell) + self%vert_coords(:, j)
      end do

      self%cell_coords(:,cell) = self%cell_coords(:,cell)/NVERTS_PER_CELL

    end do

  else
    ! 2.2 Open cells have `ghost` nodes/edges and are located
    ! on the Eastern/Southern edges of the panel. Identify the open cells
    ! assuming that the numbering is along rows beginning from the NW corner
    ! of the panel. All cells, including the open cells have a unique NW
    ! vertex. Use this NW and self%dx, self%dy to calculate the face centre,
    ! assuming the cells are parallel in the x and y directions.

    do cell=1, ncells
      base_vert = self%verts_on_cell(NW, cell)
      self%cell_coords(1,cell) = self%vert_coords(1, base_vert) + self%dx/2.0_r_def
      self%cell_coords(2,cell) = self%vert_coords(2, base_vert) - self%dy/2.0_r_def
    end do

  end if

end subroutine calc_cell_centres


!-------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------
subroutine get_dimensions( self,                                       &
                           num_nodes, num_edges, num_faces,            &
                           num_nodes_per_face, num_edges_per_face,     &
                           num_nodes_per_edge, max_num_faces_per_node )
  implicit none

  class(gen_planar_type),    intent(in) :: self

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


!-------------------------------------------------------------------------------
!> @brief   Populates the argument array with the coordinates of the
!>          mesh's vertices.
!> @details Exposes the instance's vert_coords array to the caller.
!>
!> @param[out]  node_coordinates  The argument to receive the vert_coords data.
!> @param[out]  cell_coordinates  The argument to receive cell centre coordinates.
!> @param[out]  domain_extents    Principal coordinates describing the domain.
!> @param[out]  coord_units_x     Units for x-coordinates.
!> @param[out]  coord_units_y     Units for y-coordinates.
!-------------------------------------------------------------------------------
subroutine get_coordinates(self, node_coordinates, &
                                 cell_coordinates, &
                                 domain_extents,   &
                                 coord_units_x,    &
                                 coord_units_y)

  implicit none

  class(gen_planar_type), intent(in)  :: self
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
!> @param[out]  face_node_connectivity  Face-node connectivity.
!> @param[out]  face_edge_connectivity  Face-edge connectivity.
!> @param[out]  face_face_connectivity  Face-face connectivity.
!> @param[out]  edge_node_connectivity  Edge-node connectivity.
!-------------------------------------------------------------------------------
subroutine get_connectivity( self,                   &
                             face_node_connectivity, &
                             face_edge_connectivity, &
                             face_face_connectivity, &
                             edge_node_connectivity )

  implicit none

  class(gen_planar_type), intent(in) :: self

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
!> @brief   Generates the mesh and connectivity.
!> @details Calls each of the instance methods which calculate the
!>          specified mesh and populate the arrays.
!>
!-------------------------------------------------------------------------------
subroutine generate(self)

  implicit none

  class(gen_planar_type), intent(inout) :: self

  call calc_adjacency(self)
  call calc_face_to_vert(self)
  call calc_edges(self)

  if (self%nmaps > 0) call calc_global_mesh_maps(self)

  if (apply_stretch_transform) then
    call assign_stretched_mesh_coords(self)
  else
    call calc_coords(self)
  end if

  ! NOTE that due to the way cell centres are calculated for periodic meshes
  ! this calculation must be done before rotation.
  call calc_cell_centres(self)

  if (self%rotate_mesh)then
    call rotate_mesh_coords(self%vert_coords,    self%north_pole)
    call rotate_mesh_coords(self%cell_coords,    self%north_pole)
    call rotate_mesh_coords(self%domain_extents, self%north_pole)
  end if

  ! Convert coordinate units to degrees to be CF compliant.
  if (trim(self%coord_units_x) == 'radians') then
    self%vert_coords(1,:)    = self%vert_coords(1,:)    * radians_to_degrees
    self%cell_coords(1,:)    = self%cell_coords(1,:)    * radians_to_degrees
    self%domain_extents(1,:) = self%domain_extents(1,:) * radians_to_degrees
    self%coord_units_x = 'degrees_east'
  end if

  if (trim(self%coord_units_y) == 'radians') then
    self%vert_coords(2,:)    = self%vert_coords(2,:)    * radians_to_degrees
    self%cell_coords(2,:)    = self%cell_coords(2,:)    * radians_to_degrees
    self%domain_extents(2,:) = self%domain_extents(2,:) * radians_to_degrees
    self%coord_units_y = 'degrees_north'
  end if

  if (DEBUG) call write_mesh(self)

  self%generated = .true.

  return
end subroutine generate


!-------------------------------------------------------------------------------
!> @brief   Generates requested global mesh maps.(private subroutine)
!> @details A map is generated for each requested target mesh based on this
!>          mesh objects mesh details, and those of the requested target
!>          meshes using calc_global_cell_map.
!>
!-------------------------------------------------------------------------------
subroutine calc_global_mesh_maps(self)

  implicit none

  class(gen_planar_type), intent(inout) :: self

  integer(i_def) :: source_id, source_cpp, source_ncells, target_ncells,  &
                    target_edge_cells_x, target_edge_cells_y, target_cpp, &
                    ntarget_per_source_x, ntarget_per_source_y, i
  integer(i_def), allocatable :: cell_map(:,:,:)

  if (.not. allocated( self%global_mesh_maps )) then
    allocate( self%global_mesh_maps, source=global_mesh_map_collection_type() )
  end if

  source_id  = 1
  source_cpp = self%edge_cells_x*self%edge_cells_y
  source_ncells = source_cpp*NPANELS

  do i=1, size(self%target_mesh_names)

    target_edge_cells_x  = self%target_edge_cells_x(i)
    target_edge_cells_y  = self%target_edge_cells_y(i)
    target_cpp           = target_edge_cells_x*target_edge_cells_y
    target_ncells        = target_cpp*NPANELS
    ntarget_per_source_x = max(1,target_edge_cells_x/self%edge_cells_x)
    ntarget_per_source_y = max(1,target_edge_cells_y/self%edge_cells_y)
    allocate(cell_map(ntarget_per_source_x,ntarget_per_source_y,source_ncells))

    call calc_global_cell_map( self,                &
                               target_edge_cells_x, &
                               target_edge_cells_y, &
                               cell_map )

    call self%global_mesh_maps%add_global_mesh_map( source_id, &
                                                    i+1,       &
                                                    cell_map )

    deallocate(cell_map)

  end do

  return
end subroutine calc_global_mesh_maps

!-----------------------------------------------------------------------------
!> @brief Returns the number of panels in the mesh topology.
!> @description Panels are a subset of cells in the mesh domain which may
!>              exhibit common properties.
!> @return answer Number of panels resulting from this generation strategy.
!-----------------------------------------------------------------------------
function get_number_of_panels( self ) result( answer )

  implicit none

  class(gen_planar_type), intent(in) :: self

  integer(i_def) :: answer

  answer = self%npanels

end function get_number_of_panels

!-----------------------------------------------------------------------------
!> @brief Returns mesh metadata information.
!> @details This subroutine is provided as a means to request specific metadata
!>          from the current mesh configuration.
!>
!> @param[out, optional]  mesh_name           Name of mesh instance to generate.
!> @param[out, optional]  geometry            Mesh domain surface type.
!> @param[out, optional]  topology            Mesh boundary/connectivity type.
!> @param[out, optional]  coord_sys           Coordinate system to position nodes.
!> @param[out, optional]  periodic_x          Periodic in E-W direction.
!> @param[out, optional]  periodic_y          Periodic in N-S direction.
!> @param[out, optional]  npanels             Number of panels use to describe mesh.
!> @param[out, optional]  edge_cells_x        Number of panel edge cells (x-axis).
!> @param[out, optional]  edge_cells_y        Number of panel edge cells (y-axis).
!> @param[out, optional]  constructor_inputs  Inputs used to create this mesh from
!>                                            the mesh_generator.
!> @param[out, optional]  nmaps               Number of maps to create with this mesh
!>                                            as source mesh.
!> @param[out, optional]  rim_depth           Rim depth of LBC mesh (LAMs).
!> @param[out, optional]  target_mesh_names   Mesh names of the target meshes that
!>                                            this mesh has maps for.
!> @param[out, optional]  maps_edge_cells_x   Number of panel edge cells (x-axis) of
!>                                            target mesh(es) to create map(s) for.
!> @param[out, optional]  maps_edge_cells_y   Number of panel edge cells (y-axis) of
!>                                            target mesh(es) to create map(s) for.
!> @param[out, optional]  north_pole          [Longitude, Latitude] of north pole
!>                                            used for domain orientation (degrees).
!> @param[out, optional]  null_island         [Longitude, Latitude] of null island
!>                                            used for domain orientation (degrees).
!> @param[out, optional]  equatorial_latitude Latitude of the equator of the mesh
!>                                            implying a stretching towards one
!>                                            of the poles (degrees)
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

  class(gen_planar_type),        intent(in)  :: self
  character(str_def),  optional, intent(out) :: mesh_name
  character(str_def),  optional, intent(out) :: geometry
  character(str_def),  optional, intent(out) :: topology
  character(str_def),  optional, intent(out) :: coord_sys
  logical(l_def),      optional, intent(out) :: periodic_xy(2)

  integer(i_def),      optional, intent(out) :: edge_cells_x
  integer(i_def),      optional, intent(out) :: edge_cells_y
  integer(i_def),      optional, intent(out) :: nmaps
  integer(i_def),      optional, intent(out) :: rim_depth
  integer(i_def),      optional, intent(out) :: void_cell

  character(str_longlong), optional, intent(out) :: constructor_inputs

  character(str_def), optional, allocatable, intent(out) :: target_mesh_names(:)
  integer(i_def),     optional, allocatable, intent(out) :: maps_edge_cells_x(:)
  integer(i_def),     optional, allocatable, intent(out) :: maps_edge_cells_y(:)

  real(r_def), optional, intent(out) :: north_pole(2)
  real(r_def), optional, intent(out) :: null_island(2)
  real(r_def), optional, intent(out) :: equatorial_latitude

  real(r_def) :: factor

  if (present(mesh_name))    mesh_name      = self%mesh_name
  if (present(geometry))     geometry       = key_from_geometry(self%geometry)
  if (present(topology))     topology       = key_from_topology(self%topology)
  if (present(coord_sys))    coord_sys      = key_from_coord_sys(self%coord_sys)
  if (present(periodic_xy))  periodic_xy(:) = self%periodic_xy(:)
  if (present(edge_cells_x)) edge_cells_x   = self%edge_cells_x
  if (present(edge_cells_y)) edge_cells_y   = self%edge_cells_y
  if (present(nmaps))        nmaps          = self%nmaps
  if (present(rim_depth))    rim_depth      = imdi
  if (present(void_cell))    void_cell    = VOID_ID

  if (present(constructor_inputs)) constructor_inputs = trim(self%constructor_inputs)

  if (self%nmaps > 0) then
    if (present(target_mesh_names)) target_mesh_names = self%target_mesh_names
    if (present(maps_edge_cells_x)) maps_edge_cells_x = self%target_edge_cells_x
    if (present(maps_edge_cells_y)) maps_edge_cells_y = self%target_edge_cells_y
  end if

  ! Convert to degrees for cf-compliance if required
  if ( self%coord_sys == coord_sys_ll .and. &
       self%geometry  == geometry_spherical ) then
    factor = radians_to_degrees
  else
    factor = 1.0_r_def
  end if

  if (present(north_pole))  north_pole(:)  = factor * self%north_pole(:)
  if (present(null_island)) null_island(:) = factor * self%null_island(:)
  if (present(equatorial_latitude)) equatorial_latitude = factor * self%equatorial_latitude

  return
end subroutine get_metadata


!-------------------------------------------------------------------------------
!> @brief  Gets the global mesh map collection which uses
!>         this mesh as the source mesh.
!>
!> @return global_mesh_maps global_mesh_map_collection_type.
!-------------------------------------------------------------------------------
function get_global_mesh_maps(self) result(global_mesh_maps)

  implicit none

  class(gen_planar_type), target, intent(in) :: self

  type(global_mesh_map_collection_type), pointer :: global_mesh_maps

  global_mesh_maps => self%global_mesh_maps

  return
end function get_global_mesh_maps


!-------------------------------------------------------------------------------
!> @brief Subroutine to manually deallocate any memory used by the object.
!-------------------------------------------------------------------------------
subroutine clear(self)

    implicit none

    class(gen_planar_type), intent(inout) :: self

    if (allocated(self%north_cells))   deallocate( self%north_cells )
    if (allocated(self%east_cells))    deallocate( self%east_cells  )
    if (allocated(self%south_cells))   deallocate( self%south_cells )
    if (allocated(self%west_cells))    deallocate( self%west_cells  )

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

!-------------------------------------------------------------------------------
!> @brief Finaliser for the planar ugrid mesh generator (gen_planar_type)
!-------------------------------------------------------------------------------
subroutine gen_planar_final(self)

    implicit none

    type (gen_planar_type), intent(inout) :: self

    call self%clear()

    return
end subroutine gen_planar_final


!-------------------------------------------------------------------------------
!> @brief Writes out the mesh and connectivity for debugging purposes.
!>
!-------------------------------------------------------------------------------
subroutine write_mesh(self)

  use global_mesh_map_collection_mod, only: global_mesh_map_collection_type
  use global_mesh_map_mod, only: global_mesh_map_type
  use, intrinsic :: iso_fortran_env, only: stdout => output_unit

  implicit none

  class(gen_planar_type), intent(in) :: self

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

  ncells = self%edge_cells_x * self%edge_cells_y

  write(stdout,'(A)')    "====DEBUG INFO===="
  write(stdout,'(A)')    "Mesh name: "// trim(self%mesh_name)
  write(stdout,'(A)')    "Geometry:  "// trim(key_from_geometry(self%geometry))
  write(stdout,'(A)')    "Topology:  "// trim(key_from_topology(self%topology))
  write(stdout,'(A,L1)') "Periodic in x-axis: ", self%periodic_xy(1)
  write(stdout,'(A,L1)') "Periodic in y-axis: ", self%periodic_xy(2)
  write(stdout,'(A,I0)') "Panels:    ", NPANELS
  write(stdout,'(A,I0)') "Panel edge cells (x): ", self%edge_cells_x
  write(stdout,'(A,I0)') "Panel edge cells (y): ", self%edge_cells_y
  write(stdout,'(A,I0)') 'Number of nodes: ', self%n_nodes
  write(stdout,'(A,I0)') 'Number of edges: ', self%n_edges
  write(stdout,'(A,I0)') 'Number of cells: ', self%n_faces
  write(stdout,'(A)')    "Coord_sys:  "// trim(key_from_coord_sys(self%coord_sys))
  write(stdout,'(A)')    "Co-ord (x) units: "// trim(self%coord_units_x)
  write(stdout,'(A)')    "Co-ord (y) units: "// trim(self%coord_units_y)
  write(stdout,'(A,I0)') "Mappings to other meshes: ", self%nmaps
  do i=1, self%nmaps
    write(stdout,'(T4,A,2(I0,A))') trim(self%target_mesh_names(i))//     &
                                   '(', self%target_edge_cells_x(i),',', &
                                        self%target_edge_cells_y(i),')'
  end do

  write(stdout,'(A)') ''
  write(stdout,'(A)') '=========================='
  write(stdout,'(A)') ' Cell adjacency (W,S,E,N)'
  write(stdout,'(A)') '=========================='
  do cell=1, self%n_faces
    tmp_str=''
    write(tmp_str,'(I07,A,4(I07,"  "))') cell,' => ', self%cell_next(:,cell)
    write(stdout,('(A)')) trim(tmp_str)
  end do


  write(stdout,'(A)') ''
  write(stdout,'(A)') '================================='
  write(stdout,'(A)') ' Vertices on cells (SW,SE,NW,NE)'
  write(stdout,'(A)') '================================='
  do cell=1, self%n_faces
    tmp_str=''
    write(tmp_str,'(I07,A,4(I07,"  "))') cell,' => ', self%verts_on_cell(:,cell)
    write(stdout,('(A)')) trim(tmp_str)
  end do


  write(stdout,'(A)') ''
  write(stdout,'(A)') '================================='
  write(stdout,'(A)') ' Edges on cells (W,S,E,N)'
  write(stdout,'(A)') '================================='
  do cell=1, self%n_faces
    tmp_str=''
    write(tmp_str,'(I07,A,4(I07,"  "))') cell, ' => ', self%edges_on_cell(:,cell)
    write(stdout,('(A)')) trim(tmp_str)
  end do


  write(stdout,'(A)') ''
  write(stdout,'(A)') '================================='
  write(stdout,'(A)') ' Vertices on edges, N-S/W-E'
  write(stdout,'(A)') '================================='
  do edge=1, self%n_edges
    tmp_str=''
    write(stdout,'(I07,A,I07,A,I07)')             &
        edge, ' => ', self%verts_on_edge(1,edge), &
        ' -- ', self%verts_on_edge(2,edge)
  end do


  write(stdout,'(A)')  ''
  write(stdout,'(A)') '================================='
  write(stdout,'(A)') ' Node Coordinates (x,y)'
  write(stdout,'(A)') '================================='
  do vert=1, self%n_nodes
    tmp_str=''
    write(tmp_str,'(I07,A,F10.4,A,F10.4,A)')     &
        vert,' => ( ', self%vert_coords(1,vert), &
        ', ', self%vert_coords(2,vert), ' )'
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
    write(stdout,'(4(A,I0),A)')                                                        &
        trim(self%mesh_name)//'(', self%edge_cells_x, ',', self%edge_cells_y,') => '// &
        trim(self%target_mesh_names(i))//'(', self%target_edge_cells_x(i), ',',        &
        self%target_edge_cells_y(i), '):'
    do j=1, nsource
      write(stdout,'(I7,A,10(I0," "))') j,' => ' , cell_map(:, j)
    end do
  end do

  write(stdout,'(A)')    "====END DEBUG INFO===="

  return
end subroutine write_mesh

!-------------------------------------------------------------------------------
!> @brief    Returns whether the strategy data has been generated.
!> @details  On instantiation, object data such as connectivities, coordianates
!>           etc are not calculated until the object has been "generated".
!>           Objects such as LBC strategy require these data for themselves to
!>           be generated. This function provides a means to inquire about
!>           this requirement.
!> @return   answer  Has this strategy be generated?, <<logical>>
!-------------------------------------------------------------------------------
function is_generated(self) result(answer)

  implicit none
  class(gen_planar_type), intent(in) :: self

  logical(l_def) :: answer

  answer = self%generated

  return
end function is_generated

!>==============================================================================
!> @brief Sets common partition parameters to be applied to global meshes
!>        of this type.
!>
!> @param[out]  decomposition     Object containing decomposition parameters and
!>                                method
!> @param[out]  partitioner_ptr   Mesh partitioning strategy.
!>==============================================================================
subroutine set_partition_parameters( decomposition, partitioner_ptr )

  use panel_decomposition_mod, only: panel_decomposition_type,           &
                                     auto_decomposition_type,            &
                                     row_decomposition_type,             &
                                     column_decomposition_type,          &
                                     custom_decomposition_type,          &
                                     auto_nonuniform_decomposition_type, &
                                     guided_nonuniform_decomposition_type

  use partition_mod, only: partitioner_interface, &
                           partitioner_planar

  ! Configuration modules.
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

  class(panel_decomposition_type), intent(out), allocatable :: decomposition

  procedure(partitioner_interface), &
                  intent(out), pointer :: partitioner_ptr

  partitioner_ptr => null()

  partitioner_ptr => partitioner_planar
  call log_event( "Using planar partitioner", LOG_LEVEL_INFO )

  select case(panel_decomposition)
    case( panel_decomposition_auto )
      decomposition = auto_decomposition_type()

    case( panel_decomposition_row )
      decomposition = row_decomposition_type()

    case( panel_decomposition_column )
      decomposition = column_decomposition_type()

    case( panel_decomposition_custom )
      ! use the values provided from the partitions namelist.
      decomposition = custom_decomposition_type( panel_xproc, panel_yproc )

    case ( panel_decomposition_auto_nonuniform )
      decomposition = auto_nonuniform_decomposition_type()

    case ( panel_decomposition_guided_nonuniform )
      decomposition = guided_nonuniform_decomposition_type( panel_xproc )

    case default
      call log_event( "Missing entry for panel decomposition, "// &
                    "specify 'auto' if unsure.", LOG_LEVEL_ERROR )

  end select

  call log_event( log_scratch_space, LOG_LEVEL_INFO )

end subroutine set_partition_parameters

end module gen_planar_mod
