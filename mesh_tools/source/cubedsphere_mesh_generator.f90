!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @mainpage Cubedsphere mesh generator
!> @brief   Utility to generate a cubedsphere surface mesh and write to a file
!>          which conforms to the UGRID format convention.
!> @details Usage:
!>
!>          cubedsphere_mesh_generator <filename>
!>          filename - Controlling namelist file
!>
!-----------------------------------------------------------------------------
program cubedsphere_mesh_generator

  use cli_mod,             only: get_initial_filename
  use constants_mod,       only: i_def, l_def, r_def, str_def, &
                                 cmdi, imdi, emdi, str_max_filename
  use configuration_mod,   only: read_configuration, final_configuration
  use coord_transform_mod, only: rebase_longitude_range
  use gencube_ps_mod,      only: gencube_ps_type, &
                                 set_partition_parameters

  use generate_op_global_objects_mod, only: generate_op_global_objects
  use generate_op_local_objects_mod,  only: generate_op_local_objects
  use global_mesh_collection_mod,     only: global_mesh_collection, &
                                            global_mesh_collection_type
  use halo_comms_mod,                 only: initialise_halo_comms, &
                                            finalise_halo_comms
  use io_utility_mod,                 only: open_file, close_file
  use namelist_collection_mod,        only: namelist_collection_type
  use local_mesh_collection_mod,      only: local_mesh_collection, &
                                            local_mesh_collection_type

  use log_mod,       only: initialise_logging, finalise_logging, &
                           log_event, log_set_level,             &
                           log_scratch_space, LOG_LEVEL_INFO,    &
                           LOG_LEVEL_ERROR, LOG_LEVEL_WARNING
  use mpi_mod,       only: global_mpi, create_comm, destroy_comm
  use ncdf_quad_mod, only: ncdf_quad_type
  use partition_mod, only: partition_type, partitioner_interface

  use remove_duplicates_mod,  only: any_duplicates
  use rotation_mod,           only: get_target_north_pole, &
                                    get_target_null_island
  use ugrid_2d_mod,           only: ugrid_2d_type
  use ugrid_file_mod,         only: ugrid_file_type
  use write_local_meshes_mod, only: write_local_meshes

  ! Configuration modules.
  use cubedsphere_mesh_config_mod, only: edge_cells, smooth_passes,  &
                                         stretch_factor
  use mesh_config_mod,             only: mesh_file_prefix, rotate_mesh, &
                                         n_meshes, mesh_names,          &
                                         mesh_maps, partition_mesh,     &
                                         coord_sys, coord_sys_ll,       &
                                         coord_sys_xyz,                 &
                                         key_from_coord_sys,            &
                                         topology,                      &
                                         topology_periodic,             &
                                         topology_non_periodic,         &
                                         topology_channel,              &
                                         key_from_topology,             &
                                         geometry, geometry_planar,     &
                                         geometry_spherical,            &
                                         key_from_geometry
  use partitions_config_mod,       only: max_stencil_depth, &
                                         n_partitions,      &
                                         partition_range
  use rotation_config_mod,         only: target_north_pole,           &
                                         target_null_island,          &
                                         rotation_target,             &
                                         ROTATION_TARGET_NULL_ISLAND, &
                                         ROTATION_TARGET_NORTH_POLE

  implicit none

  integer(i_def) :: communicator = -999
  integer(i_def) :: total_ranks, local_rank

  character(:), allocatable :: filename

  integer(i_def),         allocatable :: ncells(:)
  integer(i_def),         allocatable :: cpp(:)
  type(gencube_ps_type),  allocatable :: mesh_gen(:)
  type(ugrid_2d_type),    allocatable :: ugrid_2d(:)
  class(ugrid_file_type), allocatable :: ugrid_file

  integer(i_def) :: fsize
  integer(i_def) :: max_res

  integer(i_def) :: xproc, yproc
  integer(i_def) :: n_mesh_maps = 0
  integer(i_def) :: n_targets

  character(len=1),   allocatable :: map_spec(:)

  integer(i_def),     allocatable :: target_edge_cells(:)
  character(str_def), allocatable :: target_mesh_names(:)

  integer(i_def),     allocatable :: target_edge_cells_tmp(:)
  character(str_def), allocatable :: target_mesh_names_tmp(:)

  character(str_def) :: map_entry

  ! Switches.
  logical(l_def) :: l_found = .false.
  logical(l_def) :: any_duplicate_names = .false.

  ! Partition variables.
  procedure(partitioner_interface), pointer :: partitioner_ptr => null()
  integer(i_def) :: start_partition
  integer(i_def) :: end_partition

  ! Temporary variables.
  character(str_def), allocatable :: requested_mesh_maps(:)
  character(str_def) :: first_mesh
  character(str_def) :: second_mesh
  character(str_def) :: source_mesh
  character(str_def) :: target_mesh
  character(str_def) :: tmp_str
  character(str_def) :: check_mesh(2)
  integer(i_def)     :: first_mesh_edge_cells
  integer(i_def)     :: second_mesh_edge_cells
  real(r_def)        :: set_north_pole(2)
  real(r_def)        :: set_null_island(2)
  character(str_def) :: lon_str
  character(str_def) :: lat_str
  integer(i_def)     :: log_level

  ! Variables for output filenames.
  character(str_max_filename) :: output_basename
  character(str_max_filename) :: output_file

  ! Counters.
  integer(i_def) :: i, j, k, l, n_voids

  type(namelist_collection_type), save :: nml_bank

  !===================================================================
  ! 1.0 Set the logging level for the run, should really be able
  !     to set it from the command line as an option.
  !===================================================================
  call log_set_level(LOG_LEVEL_INFO )

  !===================================================================
  ! 2.0 Start up.
  !===================================================================
  call create_comm(communicator)
  call global_mpi%initialise(communicator)

  ! Initialise halo functionality.
  call initialise_halo_comms( communicator )

  total_ranks = global_mpi%get_comm_size()
  local_rank  = global_mpi%get_comm_rank()
  call initialise_logging( communicator, 'CubeGen' )
  call nml_bank%initialise( 'CubeGen', table_len=10 )

  !===================================================================
  ! 3.0 Read in the control namelists from file.
  !===================================================================
  call get_initial_filename( filename )
  call read_configuration( filename, nml_bank )
  deallocate( filename )

  max_res = maxval(edge_cells(:n_meshes))

  n_voids = count(cmdi == mesh_maps)
  if ( n_voids == 0 ) then
    n_mesh_maps = size(mesh_maps)
  else
    n_mesh_maps = size(mesh_maps) - n_voids
  end if

  !===================================================================
  ! 4.0 Perform some error checks on the namelist inputs.
  !===================================================================
  ! 4.1a Check the namelist file enumeration: geometry.
  log_level = LOG_LEVEL_ERROR
  select case (geometry)

  case (geometry_spherical)

  case (geometry_planar)
    write( log_scratch_space,'(A)' ) &
        'Invalid domain geometry, Cubed-sphere geometry is spherical '
    call log_event( log_scratch_space, log_level )

  case (emdi)
    write( log_scratch_space,'(A)' ) &
        'Enumeration key for geometry has not been set.'
    call log_event( log_scratch_space, log_level )

  case default
    write( log_scratch_space,'(A)' )                    &
        'Unrecognised enumeration key for geometry:' // &
        trim(key_from_geometry(geometry))
    call log_event( log_scratch_space, log_level )

  end select



  ! 4.1b Check the namelist file enumeration: topology.
  log_level = LOG_LEVEL_ERROR
  select case (topology)

  case (topology_periodic)

  case (topology_channel, topology_non_periodic)
    write( log_scratch_space,'(A)' ) &
        'Cubed-sphere topology should be periodic.'
    call log_event( log_scratch_space, log_level )

  case (emdi)
    write( log_scratch_space,'(A)' ) &
        'Enumeration key for topology has not been set.'
    call log_event( log_scratch_space, log_level )

  case default
    write( log_scratch_space,'(A)' )                    &
        'Unrecognised enumeration key for topology:' // &
        trim(key_from_topology(topology))
    call log_event( log_scratch_space, log_level )

  end select


  ! 4.1c Check the namelist file enumeration: coord_sys.
  log_level = LOG_LEVEL_ERROR
  select case (coord_sys)

  case (coord_sys_ll)

  case (coord_sys_xyz)
    write( log_scratch_space,'(A)' )        &
        'Cartesian co-ordinate space ' //   &
        'currently unsupported for cubed-sphere'
    call log_event( log_scratch_space, log_level )

  case (emdi)
    write( log_scratch_space,'(A)' )        &
        'Enumeration key for coord_sys ' // &
        'has not been set.'
    call log_event( log_scratch_space, log_level )

  case default
    write( log_scratch_space,'(A)' )        &
        'Unrecognised enumeration key ' //  &
        'for coord_sys:' //                 &
        trim(key_from_coord_sys(coord_sys))
    call log_event( log_scratch_space, log_level )

  end select


  ! 4.2 Check the number of meshes requested.
  if (n_meshes < 1) then
    write( log_scratch_space,'(A,I0,A)' )                  &
       'Invalid number of meshes requested, (n_meshes = ', &
       n_meshes,')'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  if ( any(mesh_names == cmdi) .or. &
       any(edge_cells == imdi) ) then
    write( log_scratch_space,'(A,I0,A)' )                   &
       'Invalid number of meshes/edge_cells requested, ' // &
       '(n_meshes = ',n_meshes,')'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! 4.3 Check for missing data.
  if (ANY(edge_cells == imdi)) then
    write( log_scratch_space,'(A)' ) &
       'Missing data in namelist variable, edge_cells'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! 4.4 Check that all meshes requested have unique names.
  any_duplicate_names = any_duplicates(mesh_names)
  if (any_duplicate_names)  then
    write( log_scratch_space,'(A)' )     &
       'Duplicate mesh names found, ' // &
       'all requested meshes must have unique names.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! 4.5 Check that all mesh map requests are unique.
  any_duplicate_names = any_duplicates(mesh_maps)
  if (any_duplicate_names)  then
    write( log_scratch_space,'(A)' )        &
       'Duplicate mesh requests found, ' // &
       'please remove duplicate requests.'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! 4.6 Perform a number of checks related to mesh map
  !     requests.
  if (n_mesh_maps > 0) then
    do i=1, n_mesh_maps

      ! Process entry for checks.
      map_entry = mesh_maps(i)
      allocate(map_spec(len(map_entry)))
      do j=1,len(map_entry)
        map_spec(j) = map_entry(j:j)
      end do

      check_mesh(1) = trim(adjustl(map_entry(:index(map_entry,':')-1)))
      check_mesh(2) = trim(adjustl(map_entry(index(map_entry,':')+1:)))
      first_mesh    = check_mesh(1)
      second_mesh   = check_mesh(2)

      do j=1, n_meshes
        if (trim(mesh_names(j)) == trim(first_mesh)) then
          first_mesh_edge_cells = edge_cells(j)
        end if

        if (trim(mesh_names(j)) == trim(second_mesh)) then
          second_mesh_edge_cells = edge_cells(j)
        end if
      end do

      ! 4.6a Check that there is a ':' in the map specification.
      if (count(':' == map_spec ) /= 1) then
        write( log_scratch_space,'(A)' )                              &
           '['//trim(adjustl(map_entry))//'] ' //                     &
           'Incorrectly specified map entry must contain one ":" ' // &
           'separating the mesh identifiers'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      deallocate(map_spec)

      ! 4.6b Check that mesh names in the map request exist.
      do j=1, size(check_mesh)

        l_found = .false.
        do k=1, n_meshes
          if (trim(check_mesh(j)) == trim(mesh_names(k))) then
            l_found = .true.
          end if
        end do

        if ( .not. l_found ) then
          write( log_scratch_space,'(A)' )       &
           '['//trim(adjustl(map_entry))//'] '// &
             'Mesh "'//trim(check_mesh(j))//     &
             '" not configured for this file.'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end do

      ! 4.6c Check the map request is not mapping at mesh
      !      to itself.
      if (trim(first_mesh) == trim(second_mesh)) then
        write( log_scratch_space,'(A)' )             &
           '['//trim(adjustl(map_entry))//'] '//     &
           'Found identical adjacent mesh names "'// &
           trim(mesh_maps(i))//'", requested for mapping.'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      ! 4.6d Check that the number of edge cells of the meshes
      !      are not the same.
      if (first_mesh_edge_cells == second_mesh_edge_cells) then
        write( log_scratch_space,'(A,I0,A)' )                  &
           '['//trim(adjustl(map_entry))//'] '//               &
           'Found identical adjacent mesh edge cells,',        &
           first_mesh_edge_cells,', requested for mapping "'// &
           trim(first_mesh)//'"-"'//trim(second_mesh)//'".'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      ! 4.6e Check that the number of edge cells for one mesh is a
      !      factor of the number of edges cells on the other mesh.
      if ( mod(first_mesh_edge_cells, second_mesh_edge_cells) /= 0 .and. &
           mod(second_mesh_edge_cells, first_mesh_edge_cells) /= 0 ) then
        write( log_scratch_space,'(2(A,I0))' )                           &
          '['//trim(adjustl(map_entry))//'] '//                          &
          'Edges cells of one mesh must be a factor of the other. '//    &
          trim(first_mesh)//' edge cells=',first_mesh_edge_cells,', '//  &
          trim(second_mesh)//' edge cells=',second_mesh_edge_cells
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end do  ! n_mesh_maps
  end if  ! n_mesh_maps > 0

  ! 4.7 Checks related to partitioning.
  if (partition_mesh) then

    ! 4.7a Checks on valid output partition range.
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

    ! 4.7b Checks on valid number of partitions.
    if ( n_partitions <  1 ) then
      write( log_scratch_space,'(A,I0)' ) &
          'At least 1 partition must be requested.'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

  end if ! partition_mesh


  !===================================================================
  ! 5.0 Create unique list of requested mesh maps.
  !     Each map request will create two maps, one in each direction.
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
  ! 6.0 Generate each mesh and any maps associated with it
  !     where the mesh is that source.
  !===================================================================
  log_level=LOG_LEVEL_INFO
  write( log_scratch_space,'(A)' ) &
      '===================================================================='
  call log_event( log_scratch_space, log_level )

  write( log_scratch_space,'(A)' ) &
      'Mesh geometry: ' // trim(key_from_geometry(geometry))
  call log_event( log_scratch_space, log_level )

  write( log_scratch_space,'(A)' ) &
      'Mesh topology: ' // trim(key_from_topology(topology))
  call log_event( log_scratch_space, log_level )

  write( log_scratch_space,'(A)' ) &
      'Co-ordinate system: '// trim(key_from_coord_sys(coord_sys))
  call log_event( log_scratch_space, log_level )

  write( log_scratch_space,'(A)' ) &
      '===================================================================='
  call log_event( log_scratch_space, log_level )
  call log_event( 'Additional transformations', log_level )

  write( log_scratch_space,'(A,I0)' ) &
      '  Smoothing passes: ', smooth_passes
  call log_event( log_scratch_space, log_level )

  if (smooth_passes > 0) then
    write( log_scratch_space,'(A,I0)' )               &
        '  Smoothing on meshes for intergrid maps '// &
        'may result in non-colocated nodes.'
    call log_event( log_scratch_space, LOG_LEVEL_WARNING )
  end if

  if (stretch_factor >= 0.0_r_def) then
    write( log_scratch_space,'(A,F4.2)' ) &
      '  Stretch factor:   ', stretch_factor
    call log_event( log_scratch_space, log_level )
  end if

  ! 6.1 Generate objects which know how to generate each requested
  !     unique mesh.
  allocate( cpp      (n_meshes) )
  allocate( mesh_gen (n_meshes) )
  allocate( ncells   (n_meshes) )
  allocate( ugrid_2d (n_meshes) )


  ! 6.2 Assign temporary arrays for target meshes in requested maps.
  if (n_mesh_maps > 0) then
    if (allocated( target_mesh_names_tmp)) deallocate( target_mesh_names_tmp )
    if (allocated( target_edge_cells_tmp)) deallocate( target_edge_cells_tmp )
    allocate( target_mesh_names_tmp(n_mesh_maps*2) )
    allocate( target_edge_cells_tmp(n_mesh_maps*2) )
  end if

  if (rotate_mesh) then

    write(log_scratch_space,'(A)') &
        '  Rotation of mesh requested with: '
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    select case( rotation_target )
    case ( ROTATION_TARGET_NULL_ISLAND )
      ! Use the domain_centre (Null Island) rather than pole as input.
      set_north_pole(:)  = get_target_north_pole(target_null_island)
      set_null_island(:) = target_null_island(:)
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
    write( log_scratch_space,'(A)' )    &
       '    North pole [lon,lat]: [' // &
       trim(adjustl(lon_str)) // ',' // &
       trim(adjustl(lat_str)) // ']'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    write( lon_str,'(F10.2)' ) set_null_island(1)
    write( lat_str,'(F10.2)' ) set_null_island(2)
    write( log_scratch_space,'(A)' )     &
       '    Null island [lon,lat]: [' // &
       trim(adjustl(lon_str)) // ','  // &
       trim(adjustl(lat_str)) // ']'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end if

  write( log_scratch_space,'(A)' ) &
      '===================================================================='
  call log_event( log_scratch_space, LOG_LEVEL_INFO )
  call log_event( "Generating mesh(es):", log_level )
  do i=1, n_meshes
    cpp(i)    = edge_cells(i)*edge_cells(i)
    ncells(i) = cpp(i)*mesh_gen(i)%get_number_of_panels()

    write( log_scratch_space,'(2(A,I0),A)' )         &
        '  Creating Mesh: '// trim(mesh_names(i)) // &
        '(',edge_cells(i),',',edge_cells(i),')'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )


    ! 6.3 Get any target mappings requested for this mesh.
    n_targets = 0
    if (n_mesh_maps > 0) then

      target_mesh_names_tmp = cmdi
      target_edge_cells_tmp = imdi
      l=1
      do j=1, size(requested_mesh_maps)
        tmp_str= requested_mesh_maps(j)
        source_mesh = tmp_str( :index(tmp_str,':')-1)
        target_mesh = tmp_str( index(tmp_str,':')+1:)
        if (trim(source_mesh) == trim(mesh_names(i))) then
          do k=1, n_meshes
            if ( trim(target_mesh) == trim(mesh_names(k)) ) then
              target_mesh_names_tmp(l) = trim(mesh_names(k))
              target_edge_cells_tmp(l) = edge_cells(k)
              l=l+1
            end if
          end do
        end if
      end do
      n_targets=l-1
    end if

    ! 6.4 Call generation strategy.
    if (n_targets == 0 .or. n_meshes == 1 ) then

      ! 6.4a No mesh maps required for this mesh.
      mesh_gen(i) = gencube_ps_type( mesh_name=mesh_names(i),            &
                                     edge_cells=edge_cells(i),           &
                                     nsmooth=smooth_passes,              &
                                     coord_sys=coord_sys,                &
                                     rotate_mesh=rotate_mesh,            &
                                     target_north_pole=set_north_pole,   &
                                     target_null_island=set_null_island, &
                                     stretch_factor=stretch_factor )

    else if (n_meshes > 1) then

      ! 6.4b Mesh maps may be required for this mesh.
      if (allocated(target_mesh_names)) deallocate(target_mesh_names)
      if (allocated(target_edge_cells)) deallocate(target_edge_cells)

      allocate(target_mesh_names(n_targets))
      allocate(target_edge_cells(n_targets))
      target_mesh_names(:) = target_mesh_names_tmp(:n_targets)
      target_edge_cells(:) = target_edge_cells_tmp(:n_targets)

      write( log_scratch_space,'(A,I0)' ) '    Maps to: '
      do j=1, n_targets
        write( log_scratch_space,'(2(A,I0),A)' ) &
            trim(log_scratch_space)//' '//       &
            trim(target_mesh_names(j))//         &
            '(',target_edge_cells(j), ',',       &
            target_edge_cells(j),')'
      end do
      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      mesh_gen(i) = gencube_ps_type(                         &
                        mesh_name=mesh_names(i),             &
                        edge_cells=edge_cells(i),            &
                        nsmooth=smooth_passes,               &
                        coord_sys=coord_sys,                 &
                        rotate_mesh=rotate_mesh,             &
                        target_north_pole=set_north_pole,    &
                        target_null_island=set_null_island,  &
                        target_mesh_names=target_mesh_names, &
                        target_edge_cells=target_edge_cells, &
                        stretch_factor=stretch_factor )

    else
      write( log_scratch_space,'(A,I0,A)' ) &
          '  Number of unique meshes is negative [', n_meshes,']'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    ! Pass the cubesphere generation object to the UGRID file writer.
    call ugrid_2d(i)%set_by_generator(mesh_gen(i))

    if (allocated(target_mesh_names)) deallocate(target_mesh_names)

  end do  ! n_meshes

  if ( allocated(target_edge_cells) ) deallocate(target_edge_cells)
  if ( allocated(ncells)            ) deallocate(ncells)

  call log_event( "...generation complete.", LOG_LEVEL_INFO )
  write( log_scratch_space,'(A)' ) &
      '===================================================================='
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  output_basename = trim(mesh_file_prefix)

  if (partition_mesh) then
    !===============================================================
    ! 7.0 Mesh partitioning.
    !===============================================================

    !---------------------------------------------------------------
    ! 7.1 Create global meshes.
    !---------------------------------------------------------------
    allocate( global_mesh_collection, source=global_mesh_collection_type() )
    call generate_op_global_objects( ugrid_2d, global_mesh_collection )

    !---------------------------------------------------------------
    ! 7.2 Get partitioning parameters.
    !---------------------------------------------------------------
    call set_partition_parameters( xproc, yproc, partitioner_ptr )

    !---------------------------------------------------------------
    ! 7.3 Create local meshes for partitions.
    !---------------------------------------------------------------
    allocate( local_mesh_collection, source=local_mesh_collection_type() )
    call generate_op_local_objects( local_mesh_collection,              &
                                    mesh_names, global_mesh_collection, &
                                    n_partitions, partition_range,      &
                                    max_stencil_depth,                  &
                                    xproc, yproc, partitioner_ptr )

    !---------------------------------------------------------------
    ! 7.4 Output local meshes to UGRID file.
    !---------------------------------------------------------------
    call write_local_meshes( global_mesh_collection, &
                             local_mesh_collection,  &
                             output_basename )

  else

    !=================================================================
    ! 8.0 Write out global meshes to UGRID file.
    !=================================================================
    write( output_file,'(2(A,I0),A)' ) trim(output_basename)//'.nc'

    do i=1, n_meshes

      if (.not. allocated(ugrid_file)) allocate(ncdf_quad_type::ugrid_file)

      call ugrid_2d(i)%set_file_handler(ugrid_file)

      if ( i==1 ) then
        call ugrid_2d(i)%write_to_file( trim(output_file) )
      else
        call ugrid_2d(i)%append_to_file( trim(output_file) )
      end if

      inquire(file=trim(output_file), size=fsize)
      write( log_scratch_space,'(A,I0,A)' )               &
          'Adding mesh (' // trim(mesh_names(i)) //       &
          ') to ' // trim(adjustl(output_file)) // ' - ', &
          fsize, ' bytes written.'

      call log_event( log_scratch_space, LOG_LEVEL_INFO )

      if (allocated(ugrid_file)) deallocate(ugrid_file)

    end do  ! n_meshes

  end if ! partition_mesh

  !===================================================================
  ! 9.0 Clean up and Finalise.
  !===================================================================
  if ( allocated( ncells   ) ) deallocate (ncells)
  if ( allocated( cpp      ) ) deallocate (cpp)
  if ( allocated( mesh_gen ) ) deallocate (mesh_gen)

  if ( allocated( requested_mesh_maps   ) ) deallocate (requested_mesh_maps)
  if ( allocated( target_edge_cells     ) ) deallocate (target_edge_cells)
  if ( allocated( target_mesh_names     ) ) deallocate (target_mesh_names)
  if ( allocated( target_edge_cells_tmp ) ) deallocate (target_edge_cells_tmp)
  if ( allocated( target_mesh_names_tmp ) ) deallocate (target_mesh_names_tmp)

  call finalise_halo_comms()

  call global_mpi%finalise()
  call destroy_comm()

  call finalise_logging()

  call final_configuration()

end program cubedsphere_mesh_generator
