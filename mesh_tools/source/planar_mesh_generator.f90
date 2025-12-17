!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @mainpage Planar mesh generator
!>
!> @brief   Utility to generate a planar surface mesh and write to a file
!>          which conforms to the UGRID format convention.
!> @details Usage:
!>
!>          planar_mesh_generator <filename>
!>          filename - Controlling namelist file
!>
!-----------------------------------------------------------------------------
program planar_mesh_generator

  use cli_mod,             only: get_initial_filename
  use constants_mod,       only: i_def, l_def, r_def, str_def, &
                                 cmdi, imdi, emdi, str_max_filename
  use configuration_mod,   only: read_configuration, final_configuration
  use coord_transform_mod, only: rebase_longitude_range
  use gen_lbc_mod,         only: gen_lbc_type
  use gen_planar_mod,      only: gen_planar_type, &
                                 set_partition_parameters

  use generate_global_objects_mod, only: generate_global_objects
  use generate_local_objects_mod,  only: generate_local_objects
  use global_mesh_collection_mod,  only: global_mesh_collection, &
                                         global_mesh_collection_type
  use halo_comms_mod,              only: initialise_halo_comms, &
                                         finalise_halo_comms
  use io_utility_mod,              only: open_file, close_file
  use lfric_mpi_mod,               only: global_mpi, create_comm, &
                                         destroy_comm, lfric_comm_type
  use local_mesh_collection_mod,   only: local_mesh_collection_type

  use log_mod,       only: initialise_logging, finalise_logging, &
                           log_event, log_set_level,             &
                           log_scratch_space, LOG_LEVEL_INFO,    &
                           LOG_LEVEL_ERROR

  use namelist_collection_mod, only: namelist_collection_type
  use namelist_mod,            only: namelist_type

  use ncdf_quad_mod, only: ncdf_quad_type
  use omp_lib,       only: omp_get_thread_num
  use partition_mod, only: partition_type, partitioner_interface

  use reference_element_mod,  only: reference_element_type, &
                                    reference_cube_type
  use remove_duplicates_mod,  only: any_duplicates
  use rotation_mod,           only: get_target_north_pole,  &
                                    get_target_null_island
  use timer_mod,              only: timer, init_timer, output_timer
  use ugrid_2d_mod,           only: ugrid_2d_type
  use ugrid_file_mod,         only: ugrid_file_type
  use write_local_meshes_mod, only: write_local_meshes

  use panel_decomposition_mod, only: panel_decomposition_type

  ! Configuration modules.
  use mesh_config_mod,     only: COORD_SYS_LL,          &
                                 COORD_SYS_XYZ,         &
                                 key_from_coord_sys,    &
                                 TOPOLOGY_PERIODIC,     &
                                 TOPOLOGY_NON_PERIODIC, &
                                 TOPOLOGY_CHANNEL,      &
                                 key_from_topology,     &
                                 GEOMETRY_PLANAR,       &
                                 GEOMETRY_SPHERICAL,    &
                                 key_from_geometry
  use rotation_config_mod, only: ROTATION_TARGET_NULL_ISLAND, &
                                 ROTATION_TARGET_NORTH_POLE

  implicit none

  type(lfric_comm_type) :: communicator
  integer(i_def) :: total_ranks, local_rank

  character(:), allocatable :: filename

  type(reference_cube_type) :: cube_element


  class(ugrid_file_type), allocatable :: ugrid_file
  type(gen_planar_type),  allocatable :: mesh_gen(:)
  type(ugrid_2d_type),    allocatable :: ugrid_2d(:)
  type(gen_lbc_type)  :: lbc_mesh_gen
  type(ugrid_2d_type) :: lbc_ugrid_2d

  class(panel_decomposition_type), allocatable :: decomposition

  integer(i_def) :: fsize
  integer(i_def) :: n_mesh_maps = 0
  integer(i_def) :: n_targets

  character(str_def), allocatable :: target_mesh_names(:)
  integer(i_def),     allocatable :: target_edge_cells_x(:)
  integer(i_def),     allocatable :: target_edge_cells_y(:)

  integer(i_def),     allocatable :: target_edge_cells_x_tmp(:)
  integer(i_def),     allocatable :: target_edge_cells_y_tmp(:)
  character(str_def), allocatable :: target_mesh_names_tmp(:)

  ! Partition variables.
  procedure(partitioner_interface), pointer :: partitioner_ptr
  integer(i_def) :: start_partition
  integer(i_def) :: end_partition
  integer(i_def) :: partition_id

  ! Switches.
  logical(l_def) :: l_found = .false.
  logical(l_def) :: lbc_generated = .false.
  logical(l_def) :: any_duplicate_names = .false.

  ! Temporary variables.
  character(str_def), allocatable :: requested_mesh_maps(:)
  character(str_def) :: first_mesh
  character(str_def) :: second_mesh
  character(str_def) :: tmp_str
  character(str_def) :: check_mesh(2)
  integer(i_def)     :: fine_mesh_edge_cells_x, fine_mesh_edge_cells_y
  integer(i_def)     :: first_mesh_edge_cells_x, first_mesh_edge_cells_y
  integer(i_def)     :: second_mesh_edge_cells_x, second_mesh_edge_cells_y
  real(r_def)        :: set_north_pole(2)
  real(r_def)        :: set_null_island(2)
  character(str_def) :: lon_str, lat_str
  character(str_def) :: x_str, y_str
  integer(i_def)     :: log_level

  character(str_max_filename) :: output_basename
  character(str_max_filename) :: output_file

  character(str_def) :: name

  ! Counters.
  integer(i_def) :: i, j, k, l, n_voids

  ! Configuration variables
  type(namelist_collection_type) :: configuration
  type(namelist_type), pointer   :: nml_obj

  character(str_max_filename) :: mesh_file_prefix

  logical :: partition_mesh
  logical :: rotate_mesh

  integer(i_def) :: n_meshes
  character(str_def), allocatable :: mesh_names(:)
  character(str_def), allocatable :: mesh_maps(:)

  integer(i_def) :: coord_sys
  integer(i_def) :: topology
  integer(i_def) :: geometry

  logical(l_def) :: generate_inner_halos
  integer(i_def) :: max_stencil_depth
  integer(i_def) :: n_partitions
  integer(i_def), allocatable :: partition_range(:)
  real(r_def),    allocatable :: domain_size(:)
  real(r_def),    allocatable :: domain_centre(:)

  integer(i_def) :: rotation_target
  real(r_def), allocatable :: target_north_pole(:)
  real(r_def), allocatable :: target_null_island(:)

  integer(i_def), allocatable :: edge_cells_x(:)
  integer(i_def), allocatable :: edge_cells_y(:)

  logical :: periodic_x
  logical :: periodic_y
  logical :: create_lbc_mesh
  integer(i_def) :: lbc_rim_depth

  character(str_def) :: lbc_parent_mesh

  logical :: apply_stretch_transform

  character(str_def) :: transform_mesh

  type(local_mesh_collection_type) :: local_mesh_collection

  integer(i_def) :: thread_id
  character(9), parameter :: timer_file = 'timer.txt'

  nullify(partitioner_ptr)
  nullify(nml_obj)

  !===================================================================
  ! Set the logging level for the run, should really be able
  ! to set it from the command line as an option.
  !===================================================================
  call log_set_level( LOG_LEVEL_INFO )

  !===================================================================
  ! Start up.
  !===================================================================
  cube_element = reference_cube_type()

  call create_comm(communicator)
  call global_mpi%initialise(communicator)

  ! Initialise halo functionality.
  call initialise_halo_comms( communicator )

  total_ranks = global_mpi%get_comm_size()
  local_rank  = global_mpi%get_comm_rank()
  call initialise_logging( communicator%get_comm_mpi_val(), "PlanarGen" )


  !===================================================================
  ! Read in the control namelists from file.
  !===================================================================
  call get_initial_filename( filename )
  call configuration%initialise( 'PlanarGen', table_len=10 )
  call read_configuration( filename, configuration )

  deallocate( filename )
  if (configuration%namelist_exists('mesh')) then
    nml_obj => configuration%get_namelist('mesh')
    call nml_obj%get_value( 'mesh_file_prefix', mesh_file_prefix )
    call nml_obj%get_value( 'n_meshes',         n_meshes )
    call nml_obj%get_value( 'mesh_names',       mesh_names )
    call nml_obj%get_value( 'mesh_maps',        mesh_maps )
    call nml_obj%get_value( 'partition_mesh',   partition_mesh )
    call nml_obj%get_value( 'rotate_mesh',      rotate_mesh )
    call nml_obj%get_value( 'coord_sys',        coord_sys )
    call nml_obj%get_value( 'topology',         topology )
    call nml_obj%get_value( 'geometry',         geometry )
  end if

  if (configuration%namelist_exists('partitions')) then
    nml_obj => configuration%get_namelist('partitions')
    call nml_obj%get_value( 'max_stencil_depth', max_stencil_depth )
    call nml_obj%get_value( 'n_partitions', n_partitions )
    call nml_obj%get_value( 'partition_range', partition_range )
    call nml_obj%get_value( 'partition_range', partition_range )
    call nml_obj%get_value( 'generate_inner_halos', generate_inner_halos )
  end if

  if (configuration%namelist_exists('rotation')) then
    nml_obj => configuration%get_namelist('rotation')
    call nml_obj%get_value( 'rotation_target', rotation_target )
    call nml_obj%get_value( 'target_north_pole', target_north_pole )
    call nml_obj%get_value( 'target_null_island', target_null_island )
  end if

  if (configuration%namelist_exists('planar_mesh')) then
    nml_obj => configuration%get_namelist('planar_mesh')
    call nml_obj%get_value( 'edge_cells_x', edge_cells_x )
    call nml_obj%get_value( 'edge_cells_y', edge_cells_y )
    call nml_obj%get_value( 'periodic_x', periodic_x )
    call nml_obj%get_value( 'periodic_y', periodic_y )
    call nml_obj%get_value( 'domain_size', domain_size )
    call nml_obj%get_value( 'domain_centre', domain_centre )
    call nml_obj%get_value( 'create_lbc_mesh', create_lbc_mesh )
    call nml_obj%get_value( 'lbc_rim_depth', lbc_rim_depth )
    call nml_obj%get_value( 'lbc_parent_mesh', lbc_parent_mesh )
    call nml_obj%get_value( 'apply_stretch_transform', apply_stretch_transform )
  end if

  if (configuration%namelist_exists('stretch_transform')) then
    nml_obj => configuration%get_namelist('stretch_transform')
    call nml_obj%get_value( 'transform_mesh', transform_mesh )
  end if

  call init_timer(timer_file)
  call timer('Setup-Checks')

  ! The number of mesh maps in the namelist array is unbounded
  ! and so may contain unset/empty array elements. Remove
  ! these from the initial count of mesh-maps.
  n_voids = count(cmdi == mesh_maps) + count('' == mesh_maps)
  if ( n_voids == 0 ) then
    n_mesh_maps = size(mesh_maps)
  else
    n_mesh_maps = size(mesh_maps) - n_voids
  end if

  !===================================================================
  ! Perform some error checks on the namelist inputs.
  !===================================================================
  ! Check the namelist file enumeration: geometry.
  log_level = LOG_LEVEL_ERROR
  select case (geometry)

  case (geometry_spherical, geometry_planar)

  case (emdi)
    write( log_scratch_space,'(A)' ) &
        'Enumeration key for geometry has not been set.'
    call log_event( log_scratch_space, log_level )
  case default
    write( log_scratch_space,'(A)' )                   &
        'Unrecognised enumeration key for geometry:'// &
        trim(key_from_geometry(geometry))//'.'
    call log_event( log_scratch_space, log_level )
  end select


  ! Check the namelist file enumeration: topology.
  log_level = LOG_LEVEL_ERROR
  select case (topology)

  case ( topology_periodic, &
         topology_channel,  &
         topology_non_periodic )

    select case (topology)

    case (topology_periodic)
      if (geometry == geometry_planar) then
        if (.not. (periodic_x .and. periodic_y)) then
          write( log_scratch_space,'(A)' )                &
              'A periodic planar regional mesh should '// &
              'have all boundaries as periodic.'
          call log_event( log_scratch_space, log_level )
        end if
      else
         write( log_scratch_space,'(A)' )               &
             'A periodic spherical regional mesh is '// &
             'unsupported by this generator.'
          call log_event( log_scratch_space, log_level )
      end if

    case (topology_channel)
      if ( periodic_x .eqv. periodic_y ) then
        write( log_scratch_space,'(A)' )                  &
            'A channel regional mesh should only have '// &
            'a single periodic axis.'
        call log_event( log_scratch_space, log_level )
      end if

    case (topology_non_periodic)
      if ( periodic_x .or. periodic_y ) then
        write( log_scratch_space,'(A)' )                 &
            'A non-periodic regional mesh should not '// &
            'have any periodic boundaries.'
        call log_event( log_scratch_space, log_level )
      end if

    end select

  case (emdi)
    write( log_scratch_space,'(A)' ) &
        'Enumeration key for topology has not been set.'
    call log_event( log_scratch_space, log_level )

  case default
    write( log_scratch_space,'(A)' )                    &
        'Unrecognised enumeration key for topology: '// &
        key_from_topology(topology)//'.'
    call log_event( log_scratch_space, log_level )
  end select


  ! Check the namelist file enumeration: coord_sys.
  log_level = LOG_LEVEL_ERROR
  select case (coord_sys)

  case (coord_sys_ll, coord_sys_xyz)

  case (emdi)
    write( log_scratch_space,'(A)' ) &
        'Enumeration key for coord_sys has not been set.'
    call log_event( log_scratch_space, log_level )

  case default
    write( log_scratch_space,'(A)' )                    &
        'Unrecognised enumeration key for coord_sys:'// &
        trim(key_from_coord_sys(coord_sys))//'.'
    call log_event( log_scratch_space, log_level )

  end select


  ! 4.2 Check the number of meshes requested.
  if (n_meshes < 1) then
    write( log_scratch_space,'(A,I0,A)' ) &
        'Invalid number of meshes requested, (',n_meshes,').'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! Check that there are enough entries of edge cells
  ! to match the number of meshes requested.
  if ( size(edge_cells_x) < n_meshes .or. &
       size(edge_cells_y) < n_meshes ) then
    write( log_scratch_space,'(A,I0,A)' )                    &
        'Not enough data in edge_cells_x/edge_cells_y for ', &
        n_meshes,' meshe(s).'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! Check for missing data.
  if ( any(edge_cells_x == imdi) .or. &
       any(edge_cells_y == imdi) ) then
    write( log_scratch_space,'(A)' ) &
        'Missing data in namelist variable, edge_cells_x/edge_cells_y.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! Check that all meshes requested have unique names.
  any_duplicate_names = any_duplicates(mesh_names)
  if (any_duplicate_names)  then
    write( log_scratch_space,'(A)' )     &
        'Duplicate mesh names found, '// &
        'all requested meshes must have unique names.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! Check that all mesh map requests are unique.
  any_duplicate_names = any_duplicates(mesh_maps)
  if (any_duplicate_names)  then
    write( log_scratch_space,'(A)' )        &
        'Duplicate mesh requests found, '// &
        'please remove duplicate requests.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  if (apply_stretch_transform) then
    ! This enables support meshes to be created with a variable
    ! resolution stretching function.
    do j=1, n_meshes
      if (trim(mesh_names(j)) == trim(transform_mesh)) then
        fine_mesh_edge_cells_x = edge_cells_x(j)
        fine_mesh_edge_cells_y = edge_cells_y(j)
      end if
    end do
  end if

  ! Perform a number of checks related to mesh map
  ! requests.
  if (n_mesh_maps > 0) then
    do i=1, n_mesh_maps

      tmp_str = mesh_maps(i)
      check_mesh(1) = tmp_str(:index(tmp_str,':')-1)
      check_mesh(2) = tmp_str(index(tmp_str,':')+1:)
      first_mesh    = check_mesh(1)
      second_mesh   = check_mesh(2)

      first_mesh_edge_cells_x  = imdi
      first_mesh_edge_cells_y  = imdi
      second_mesh_edge_cells_x = imdi
      second_mesh_edge_cells_y = imdi

      do j=1, n_meshes
        if (trim(mesh_names(j)) == trim(first_mesh)) then
          first_mesh_edge_cells_x = edge_cells_x(j)
          first_mesh_edge_cells_y = edge_cells_y(j)
        end if

        if (trim(mesh_names(j)) == trim(second_mesh)) then
          second_mesh_edge_cells_x = edge_cells_x(j)
          second_mesh_edge_cells_y = edge_cells_y(j)
        end if
      end do

      ! Check that mesh names in the map request exist.
      do j=1, size(check_mesh)

        l_found = .false.
        do k=1, n_meshes
          if (trim(check_mesh(j)) == trim(mesh_names(k))) then
            l_found = .true.
          end if
        end do

        if ( .not. l_found ) then
          write( log_scratch_space,'(A)' )    &
              'Mesh "'//trim(check_mesh(j))// &
              '" not configured for this file.'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end do

      ! Check the map request is not mapping at mesh
      ! to itself.
      if (trim(first_mesh) == trim(second_mesh)) then
        write( log_scratch_space,'(A)' )              &
            'Found identical adjacent mesh names "'// &
           trim(mesh_maps(i))//'", requested for mapping.'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      ! Check that the number of edge cells of the meshes
      ! are not the same.
      if ( (first_mesh_edge_cells_x == second_mesh_edge_cells_x ) .and. &
            first_mesh_edge_cells_y == second_mesh_edge_cells_y ) then
        write( log_scratch_space,'(A,I0,A)' )                    &
            'Found identical adjacent mesh edge cells, (',       &
            first_mesh_edge_cells_x,',',first_mesh_edge_cells_y, &
            '), requested for mapping "'//                       &
            trim(first_mesh)//'"-"'//trim(second_mesh)//'".'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end do
  end if

  ! Check the requested LBC parent LAM exists.
  if (create_lbc_mesh) then
    l_found = .false.
    do i=1, n_meshes
      if ( trim(mesh_names(i)) == trim(lbc_parent_mesh) ) then
        l_found=.true.
        exit
      end if
    end do
    if ( .not. l_found ) then
      write( log_scratch_space,'(A)' )                  &
          'The parent mesh, '// trim(lbc_parent_mesh)// &
          ' specified for LBC mesh generation does not exist.'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
  end if

  ! Checks related to partitioning.
  if (partition_mesh) then

    ! Checks for valid output partition range.
    start_partition = partition_range(1)
    end_partition   = partition_range(2)

    do i=1, 2
      if ( partition_range(i) <  0 .or. &
           partition_range(i) >= n_partitions ) then
        write( log_scratch_space,'(A,I0)' )         &
            'Invalid partition ID range bound, ' // &
            'valid IDs are 0:', n_partitions-1
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
    end do

    if ( partition_range(1) > partition_range(2) ) then
      write( log_scratch_space,'(A,I0)' )                     &
          'Invalid start/end partitions, start partition ' // &
          'ID should be less than end partition ID.'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    ! Check for valid number of partitions on global mesh.
    if ( n_partitions <  1 ) then
      write( log_scratch_space,'(A,I0)' ) &
          'At least 1 partition must be requested.'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

  end if ! partition_mesh


  !===================================================================
  ! Create unique list of Requested Mesh maps.
  ! Each map request will create two maps, one in each direction.
  !===================================================================
  if (n_mesh_maps > 0) then
    allocate(requested_mesh_maps(n_mesh_maps*2))
    j=1
    do i=1, n_mesh_maps
      tmp_str = mesh_maps(i)
      first_mesh  = tmp_str(:index(tmp_str,':')-1)
      second_mesh = tmp_str(index(tmp_str,':')+1:)
      write( requested_mesh_maps(j),  '(A)' ) &
          trim(first_mesh)//':'//trim(second_mesh)
      write( requested_mesh_maps(j+1),'(A)' ) &
          trim(second_mesh)//':'//trim(first_mesh)
      j=j+2
    end do
  end if


  !===================================================================
  ! Report/Check what the code thinks is requested by user.
  !===================================================================
  log_level=LOG_LEVEL_INFO

  write( log_scratch_space,'(A)' )    &
      '===================================================================='
  call log_event( log_scratch_space, log_level )

  write( log_scratch_space,'(A)' )    &
      'Mesh geometry: ' // trim(key_from_geometry(geometry))
  call log_event( log_scratch_space, log_level )

  write( log_scratch_space,'(A)' )    &
      'Mesh topology: ' // trim(key_from_topology(topology))
  call log_event( log_scratch_space, log_level )

  write( log_scratch_space,'(A)' )    &
      'Co-ordinate system: '// trim(key_from_coord_sys(coord_sys))
  call log_event( log_scratch_space, log_level )

  write( log_scratch_space,'(A,L1)' ) &
      'Periodic in x-axis: ', periodic_x
  call log_event( log_scratch_space, log_level )

  write( log_scratch_space,'(A,L1)' ) &
      'Periodic in y-axis: ', periodic_y
  call log_event( log_scratch_space, log_level )

  if (coord_sys == coord_sys_ll) then
    write( log_scratch_space,'(A)') &
        'Domain centre [lon,lat]: ['
  else
    write( log_scratch_space,'(A)') &
        'Domain centre [x,y]: ['
  end if
  write( x_str,'(F10.2)' ) domain_centre(1)
  write( y_str,'(F10.2)' ) domain_centre(2)

  write( log_scratch_space, '(A)' ) trim(log_scratch_space) // &
      trim(adjustl(x_str)) // ',' // trim(adjustl(y_str)) // ']'
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  if (coord_sys == coord_sys_ll) then
    write( log_scratch_space,'(A)') &
        'Domain size [degrees]: ['
  else
    write( log_scratch_space,'(A)') &
        'Domain size [m]: ['
  end if
  write( x_str,'(F10.2)' ) domain_size(1)
  write( y_str,'(F10.2)' ) domain_size(2)

  write( log_scratch_space, '(A)' ) trim(log_scratch_space) // &
      trim(adjustl(x_str)) // ',' // trim(adjustl(y_str)) // ']'
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  write( log_scratch_space,'(A)' )    &
      '===================================================================='
  call log_event( log_scratch_space, log_level )
  call log_event( "Generating mesh(es):", log_level )

  ! Generate objects which know how to generate each requested
  ! unique mesh.
  allocate( mesh_gen (n_meshes) )
  allocate( ugrid_2d (n_meshes) )

  ! Assign temporary arrays for target meshes in requested maps.
  if (n_mesh_maps > 0) then
    if (allocated(target_mesh_names_tmp))   deallocate(target_mesh_names_tmp)
    if (allocated(target_edge_cells_x_tmp)) deallocate(target_edge_cells_x_tmp)
    if (allocated(target_edge_cells_y_tmp)) deallocate(target_edge_cells_y_tmp)
    allocate( target_mesh_names_tmp(n_mesh_maps*2) )
    allocate( target_edge_cells_x_tmp(n_mesh_maps*2) )
    allocate( target_edge_cells_y_tmp(n_mesh_maps*2) )
  end if

  if ( geometry == geometry_spherical .and. &
       coord_sys == coord_sys_ll ) then
    if (rotate_mesh) then

      write( log_scratch_space,'(A)' ) &
         '  Rotation of mesh requested with: '
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      select case ( rotation_target )
      case ( ROTATION_TARGET_NULL_ISLAND )
        ! Use the domain_centre (Null Island) rather than pole as input.
        set_north_pole(:) = get_target_north_pole(target_null_island)
        set_null_island(:) = target_null_island
        write( log_scratch_space,'(A)' ) &
           '    Target pole will be derived from Null Island.'
        call log_event( log_scratch_space, LOG_LEVEL_INFO )

        write( lon_str,'(F10.2)' ) set_null_island(1)
        write( lat_str,'(F10.2)' ) set_null_island(2)
        write( log_scratch_space,'(A)' )     &
           '    Null Island [lon,lat]: [' // &
           trim(adjustl(lon_str)) // ',' //  &
           trim(adjustl(lat_str)) // ']'
        call log_event( log_scratch_space, LOG_LEVEL_INFO )
      case ( ROTATION_TARGET_NORTH_POLE )
        set_north_pole(:)  = target_north_pole(:)
        set_null_island(:) = get_target_null_island(target_north_pole)
      end select

      ! Ensure the requested target longitudes are in the range -180,180.
      set_null_island(1) = rebase_longitude_range( set_null_island(1), -180.0_r_def)
      set_north_pole(1)  = rebase_longitude_range( set_north_pole(1), -180.0_r_def)

      write( lon_str,'(F10.2)' ) set_north_pole(1)
      write( lat_str,'(F10.2)' ) set_north_pole(2)
      write( log_scratch_space,'(A)' )     &
         '    Target pole [lon,lat]: [' // &
         trim(adjustl(lon_str)) // ',' //  &
         trim(adjustl(lat_str)) // ']'
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

    end if
  end if

  call timer('Setup-Checks')

  write( log_scratch_space,'(A)' ) &
      '===================================================================='
  call log_event( log_scratch_space, log_level_info )
  call log_event( "Generating mesh(es):", log_level )

  call timer('Global mesh generation')

  do i=1, n_meshes

    write( log_scratch_space,'(A,2(I0,A))' )           &
       '  Creating Mesh: '// trim(mesh_names(i))//'(', &
                            edge_cells_x(i), ',',      &
                            edge_cells_y(i), ')'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )


    ! Get any target mappings requested for this mesh.
    n_targets = 0
    if (n_mesh_maps > 0) then

      target_mesh_names_tmp = cmdi
      target_edge_cells_x_tmp = imdi
      target_edge_cells_y_tmp = imdi
      l=1
      do j=1, size(requested_mesh_maps)
        tmp_str= requested_mesh_maps(j)
        if (tmp_str( :index(tmp_str,':')-1) == trim(mesh_names(i))) then
          do k=1, n_meshes
            if ( trim(tmp_str( index(tmp_str,':')+1:)) ==  &
                 trim(mesh_names(k)) ) then
              target_mesh_names_tmp(l)   = trim(mesh_names(k))
              target_edge_cells_x_tmp(l) = edge_cells_x(k)
              target_edge_cells_y_tmp(l) = edge_cells_y(k)
              l=l+1
            end if
          end do
        end if
      end do
      n_targets=l-1
    end if  ! n_mesh_maps > 0

    ! Call generation strategy.
    if (n_targets == 0 .or. n_meshes == 1 ) then

      mesh_gen(i) = gen_planar_type(                                 &
                        reference_element  = cube_element,           &
                        mesh_name          = mesh_names(i),          &
                        geometry           = geometry,               &
                        topology           = topology,               &
                        coord_sys          = coord_sys,              &
                        edge_cells_x       = edge_cells_x(i),        &
                        edge_cells_y       = edge_cells_y(i),        &
                        fine_mesh_edge_cells_x                       &
                                           = fine_mesh_edge_cells_x, &
                        fine_mesh_edge_cells_y                       &
                                           = fine_mesh_edge_cells_y, &
                        periodic_x         = periodic_x,             &
                        periodic_y         = periodic_y,             &
                        domain_size        = domain_size,            &
                        domain_centre      = domain_centre,          &
                        rotate_mesh        = rotate_mesh,            &
                        target_north_pole  = set_north_pole,         &
                        target_null_island = set_null_island )

    else if (n_meshes > 1) then

      if (allocated(target_mesh_names))   deallocate(target_mesh_names)
      if (allocated(target_edge_cells_x)) deallocate(target_edge_cells_x)
      if (allocated(target_edge_cells_y)) deallocate(target_edge_cells_y)

      allocate( target_mesh_names(n_targets)   )
      allocate( target_edge_cells_x(n_targets) )
      allocate( target_edge_cells_y(n_targets) )
      target_mesh_names(:)   = target_mesh_names_tmp(:n_targets)
      target_edge_cells_x(:) = target_edge_cells_x_tmp(:n_targets)
      target_edge_cells_y(:) = target_edge_cells_y_tmp(:n_targets)

      write( log_scratch_space,'(A,I0)' ) '    Maps to:'
      do j=1, n_targets
        write( log_scratch_space,'(2(A,I0),A)' ) &
            trim(log_scratch_space)//' '//       &
            trim(target_mesh_names(j))//         &
            '(',target_edge_cells_x(j),',',      &
            target_edge_cells_y(j),')'
      end do

      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      mesh_gen(i) = gen_planar_type(                               &
                        cube_element, mesh_names(i),               &
                        geometry, topology, coord_sys,             &
                        edge_cells_x(i), edge_cells_y(i),          &
                        fine_mesh_edge_cells_x,                    &
                        fine_mesh_edge_cells_y,                    &
                        periodic_x, periodic_y,                    &
                        domain_size, domain_centre,                &
                        target_mesh_names   = target_mesh_names,   &
                        target_edge_cells_x = target_edge_cells_x, &
                        target_edge_cells_y = target_edge_cells_y, &
                        rotate_mesh         = rotate_mesh,         &
                        target_north_pole   = set_north_pole,      &
                        target_null_island  = set_null_island )

    else
      write( log_scratch_space,'(A,I0,A)' ) &
         '  Number of meshes is negative [', n_meshes,']'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    ! Pass the generation object to the ugrid file writer.
    call ugrid_2d(i)%set_by_generator(mesh_gen(i))

    ! Generate the LBC mesh generation strategy.
    if ( create_lbc_mesh .and. &
       ( trim(mesh_names(i)) == trim(lbc_parent_mesh) ) ) then
      lbc_mesh_gen = gen_lbc_type(mesh_gen(i), lbc_rim_depth)
      call lbc_ugrid_2d%set_by_generator(lbc_mesh_gen)
    end if

    if (allocated(target_mesh_names))   deallocate(target_mesh_names)
    if (allocated(target_edge_cells_x)) deallocate(target_edge_cells_x)
    if (allocated(target_edge_cells_y)) deallocate(target_edge_cells_y)

  end do ! n_meshes

  call timer('Global mesh generation')

  call log_event( "...generation complete.", LOG_LEVEL_INFO )
  write( log_scratch_space,'(A)' ) &
      '===================================================================='
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  output_basename = trim(mesh_file_prefix)

  if (partition_mesh) then

    call timer('Partition meshes')

    !=================================================================
    ! Mesh partitioning.
    !=================================================================

    !-----------------------------------------------------------------
    ! Create global meshes.
    !-----------------------------------------------------------------
    allocate( global_mesh_collection, source=global_mesh_collection_type() )
    if ( create_lbc_mesh ) then
      call generate_global_objects( ugrid_2d, global_mesh_collection, &
                                    lbc_ugrid_2d )
    else
      call generate_global_objects( ugrid_2d, global_mesh_collection )
    end if

    !---------------------------------------------------------------
    ! Get partitioning parameters.
    !---------------------------------------------------------------
    call set_partition_parameters( decomposition, partitioner_ptr )




    !---------------------------------------------------------------
    ! Create local meshes for partitions.
    !---------------------------------------------------------------
    start_partition = partition_range(1)
    end_partition   = partition_range(2)

    thread_id = 0

!$omp parallel default(firstprivate) shared(global_mesh_collection)
!$omp do schedule(static)

    do partition_id=start_partition, end_partition

      thread_id = omp_get_thread_num()

      if ( create_lbc_mesh ) then
        call generate_local_objects(                   &
                      local_mesh_collection,           &
                      global_mesh_collection,          &
                      mesh_names, partition_id,        &
                      n_partitions, max_stencil_depth, &
                      generate_inner_halos,            &
                      decomposition, partitioner_ptr,  &
                      lbc_parent_mesh )
      else
        call generate_local_objects(                   &
                      local_mesh_collection,           &
                      global_mesh_collection,          &
                      mesh_names, partition_id,        &
                      n_partitions, max_stencil_depth, &
                      generate_inner_halos,            &
                      decomposition, partitioner_ptr )
      end if

      !---------------------------------------------------------------
      ! Output local meshes to UGRID file
      !---------------------------------------------------------------
      call write_local_meshes( partition_id,           &
                               global_mesh_collection, &
                               local_mesh_collection,  &
                               output_basename )

      call local_mesh_collection%clear()

    end do ! partition_id

!$omp end do
!$omp end parallel

    call timer('Partition meshes')

  else

    !=================================================================
    ! Write out global meshes to UGRID file.
    !=================================================================
    write(output_file,'(2(A,I0),A)') trim(output_basename)//'.nc'

    do i=1, n_meshes

      if (.not. allocated(ugrid_file)) allocate(ncdf_quad_type::ugrid_file)

      call ugrid_2d(i)%set_file_handler(ugrid_file)

      if (i==1) then
        call ugrid_2d(i)%write_to_file( trim(output_file) )
      else
        call ugrid_2d(i)%append_to_file( trim(output_file) )
      end if

      inquire(file=output_file, size=fsize)
      write( log_scratch_space, '(A,I0,A)')               &
          'Adding mesh (' // trim(mesh_names(i)) //       &
          ') to ' // trim(adjustl(output_file)) // ' - ', &
          fsize, ' bytes written.'

      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      if (allocated(ugrid_file)) deallocate(ugrid_file)

    end do ! n_meshes

    !===================================================================
    ! Now create/output LBC mesh
    !===================================================================
    ! A LBC mesh is created from a parent planar mesh strategy that has
    ! been generated. The name of the resulting LBC mesh will be:
    !
    !    <parent mesh name>-lbc
    !
    if (create_lbc_mesh) then

      lbc_generated = .false.

      do i=1, size(mesh_gen)

        if (lbc_generated) exit

        call mesh_gen(i)%get_metadata(mesh_name=name)
        if (trim(name) == trim(lbc_parent_mesh)) then

          if (.not. allocated(ugrid_file)) allocate(ncdf_quad_type::ugrid_file)

          call lbc_ugrid_2d%set_file_handler(ugrid_file)
          call lbc_ugrid_2d%append_to_file( trim(output_file) )

          inquire(file=output_file, size=fsize)
          write( log_scratch_space,'(A,I0,A)' )                &
              'Adding lbc mesh for ' // trim(mesh_names(i)) // &
              ' to ' // trim(adjustl(output_file)) // ' - ',   &
              fsize, ' bytes written.'

          call log_event( log_scratch_space, LOG_LEVEL_INFO )
          if (allocated(ugrid_file)) deallocate(ugrid_file)

          lbc_generated = .true.

        end if
      end do

    end if ! create_lbc_mesh

  end if ! partition_mesh

  call output_timer()

  !===================================================================
  ! Clean up and Finalise.
  !===================================================================

  if ( allocated( mesh_gen ) ) deallocate (mesh_gen)

  if ( allocated( requested_mesh_maps     ) ) deallocate (requested_mesh_maps)
  if ( allocated( target_mesh_names       ) ) deallocate (target_mesh_names)
  if ( allocated( target_edge_cells_x     ) ) deallocate (target_edge_cells_x)
  if ( allocated( target_edge_cells_y     ) ) deallocate (target_edge_cells_y)
  if ( allocated( target_mesh_names_tmp   ) ) deallocate (target_mesh_names_tmp)
  if ( allocated( target_edge_cells_x_tmp ) ) deallocate (target_edge_cells_x_tmp)
  if ( allocated( target_edge_cells_y_tmp ) ) deallocate (target_edge_cells_y_tmp)

  if ( allocated( mesh_names         ) ) deallocate(mesh_names)
  if ( allocated( partition_range    ) ) deallocate(partition_range)
  if ( allocated( target_north_pole  ) ) deallocate(target_north_pole)
  if ( allocated( target_null_island ) ) deallocate(target_null_island)
  if ( allocated( edge_cells_x       ) ) deallocate(edge_cells_x)
  if ( allocated( edge_cells_y       ) ) deallocate(edge_cells_y)
  if ( allocated( domain_size        ) ) deallocate(domain_size)
  if ( allocated( domain_centre      ) ) deallocate(domain_centre)

  call finalise_halo_comms()

  call global_mpi%finalise()
  call destroy_comm()

  call finalise_logging()

  call final_configuration()

end program planar_mesh_generator
