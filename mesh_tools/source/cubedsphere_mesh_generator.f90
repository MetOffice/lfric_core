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

  use cli_mod,           only: get_initial_filename
  use constants_mod,     only: i_def, imdi, l_def, str_def, str_long, cmdi
  use configuration_mod, only: read_configuration, final_configuration
  use mesh_config_mod,   only: n_meshes, mesh_names, mesh_filename, mesh_maps
  use cubedsphere_mesh_config_mod,                        &
                         only: edge_cells, smooth_passes, &
                               do_rotate, lat_north,      &
                               lon_north, rotate_angle,   &
                               stretch_factor
  use mpi_mod,           only: initialise_comm, store_comm, finalise_comm, &
                               get_comm_size, get_comm_rank
  use gencube_ps_mod,    only: gencube_ps_type
  use io_utility_mod,    only: open_file, close_file
  use log_mod,           only: initialise_logging, finalise_logging, &
                               log_event, log_set_level, &
                               log_scratch_space, &
                               LOG_LEVEL_INFO, LOG_LEVEL_ERROR
  use ncdf_quad_mod,     only: ncdf_quad_type

  use remove_duplicates_mod, only: any_duplicates
  use ugrid_2d_mod,      only: ugrid_2d_type
  use ugrid_file_mod,    only: ugrid_file_type

  implicit none

  integer(i_def) :: comm, total_ranks, local_rank

  character(:), allocatable   :: filename

  integer(i_def),         allocatable :: ncells(:)
  integer(i_def),         allocatable :: cpp(:)
  type(gencube_ps_type),  allocatable :: mesh_gen(:)
  type(ugrid_2d_type),    allocatable :: ugrid_2d(:)
  class(ugrid_file_type), allocatable :: ugrid_file

  integer(i_def) :: fsize
  integer(i_def) :: max_res

  integer(i_def) :: nsmooth
  integer(i_def) :: n_mesh_maps = 0
  integer(i_def) :: n_targets

  character(len=1),   allocatable :: map_spec(:)

  integer(i_def),     allocatable :: target_edge_cells(:)
  character(str_def), allocatable :: target_mesh_names(:)

  integer(i_def),     allocatable :: target_edge_cells_tmp(:)
  character(str_def), allocatable :: target_mesh_names_tmp(:)

  character(str_def) :: map_entry

  ! Switches
  logical(l_def) :: l_found = .false.
  logical(l_def) :: any_duplicate_names = .false.

  ! Parametes
  integer(i_def), parameter :: npanels = 6
  integer(i_def), parameter :: max_n_targets = 6

  ! Temporary variables
  character(str_def), allocatable :: requested_mesh_maps(:)
  character(str_def) :: first_mesh
  character(str_def) :: second_mesh
  character(str_def) :: source_mesh
  character(str_def) :: target_mesh
  character(str_def) :: tmp_str
  character(str_def) :: check_mesh(2)
  integer(i_def)     :: first_mesh_edge_cells
  integer(i_def)     :: second_mesh_edge_cells

  ! Counters
  integer(i_def) :: i, j, k, l, n_voids

  !===================================================================
  ! 1.0 Set the logging level for the run, should really be able
  !     to set it from the command line as an option
  !===================================================================
  call log_set_level(LOG_LEVEL_INFO)

  !===================================================================
  ! 2.0 Start up
  !===================================================================
  call initialise_comm(comm)
  call store_comm(comm)
  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()
  call initialise_logging(local_rank, total_ranks, "cubedsphere")

  !===================================================================
  ! 3.0 Read in the control namelists from file
  !===================================================================
  call get_initial_filename( filename )
  call read_configuration( filename )
  deallocate( filename )

  max_res = maxval(edge_cells(:n_meshes))

  n_voids = count(cmdi == mesh_maps)
  if ( n_voids == 0 ) then
    n_mesh_maps = size(mesh_maps)
  else
    n_mesh_maps = size(mesh_maps) - n_voids
  end if

  !===================================================================
  ! 4.0 Perform some error checks on the namelist inputs
  !===================================================================
  ! 4.1 Check the number of meshes requested.
  if (n_meshes < 1) then
    write(log_scratch_space,'(A,I0,A)') &
       'Invalid number of meshes requested, (n_meshes = ',&
       n_meshes,')'
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if

  if ( any(mesh_names == cmdi) .or. &
       any(edge_cells == imdi) ) then
    write(log_scratch_space,'(A,I0,A)') &
       'Invalid number of meshes/edge_cells requested, '//&
       '(n_meshes = ',n_meshes,')'
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if

  ! 4.2 Check for missing data.
  if (ANY(edge_cells == imdi)) then
    write(log_scratch_space,'(A)') &
       'Missing data in namelist variable, edge_cells'
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if

  ! 4.3 Check that all meshes requested have unique names.
  any_duplicate_names = any_duplicates(mesh_names)
  if (any_duplicate_names)  then
    write(log_scratch_space,'(A)')          &
       'Duplicate mesh names found, '// &
       'all requested meshes must have unique names.'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  ! 4.4 Check that all mesh map requests are unique.
  any_duplicate_names = any_duplicates(mesh_maps)
  if (any_duplicate_names)  then
    write(log_scratch_space,'(A)')          &
       'Duplicate mesh requests found, '//  &
       'please remove duplicate requests.'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  ! Perform a number of checks related to mesh map
  ! requests.
  if (n_mesh_maps > 0) then
    do i=1, n_mesh_maps

      ! Process entry for checks
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

      ! 4.4 Check that there is a ':' in the map specification.
      if (count(':' == map_spec ) /= 1) then
        write(log_scratch_space,'(A)')                               &
           '['//trim(adjustl(map_entry))//'] '//                     &
           'Incorrectly specified map entry must contain one ":" '// &
           'separating the mesh identifiers'
        call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
      end if
      deallocate(map_spec)

      ! 4.4 Check that mesh names in the map request exist.
      do j=1, size(check_mesh)

        l_found = .false.
        do k=1, n_meshes
          if (trim(check_mesh(j)) == trim(mesh_names(k))) then
            l_found = .true.
          end if
        end do

        if ( .not. l_found ) then
          write(log_scratch_space,'(A)')         &
           '['//trim(adjustl(map_entry))//'] '// &
             'Mesh "'//trim(check_mesh(j))//     &
             '" not configured for this file.'
          call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
        end if
      end do

      ! 4.5 Check the map request is not mapping at mesh
      !     to itself.
      if (trim(first_mesh) == trim(second_mesh)) then
        write(log_scratch_space,'(A)')               &
           '['//trim(adjustl(map_entry))//'] '//     &
           'Found identical adjacent mesh names "'// &
           trim(mesh_maps(i))//'", requested for mapping.'
        call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
      end if

      ! 4.6 Check that the number of edge cells of the meshes
      !     are not the same.
      if (first_mesh_edge_cells == second_mesh_edge_cells) then
        write(log_scratch_space,'(A,I0,A)')                    &
           '['//trim(adjustl(map_entry))//'] '//               &
           'Found identical adjacent mesh edge cells,',        &
           first_mesh_edge_cells,', requested for mapping "'// &
           trim(first_mesh)//'"-"'//trim(second_mesh)//'".'
        call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
      end if

      ! 4.7 Check that the number of edge cells for one mesh is a
      !     factor of the number of edges cells on the other mesh.
      if ( mod(first_mesh_edge_cells, second_mesh_edge_cells) /= 0 .and. &
           mod(second_mesh_edge_cells, first_mesh_edge_cells) /= 0 ) then
        write(log_scratch_space, '(2(A,I0))')                            &
          '['//trim(adjustl(map_entry))//'] '//                          &
          'Edges cells of one mesh must be a factor of the other. '//    &
          trim(first_mesh)//' edge cells=',first_mesh_edge_cells,', '//  &
          trim(second_mesh)//' edge cells=',second_mesh_edge_cells
        call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
      end if

    end do  ! n_mesh_maps
  end if  ! n_mesh_maps > 0


  !===================================================================
  ! 5.0 Create unique list of Requested Mesh maps
  !     Each map request will create two maps, one in each direction
  !===================================================================
  if (n_mesh_maps > 0) then
    allocate(requested_mesh_maps(n_mesh_maps*2))
    j=1
    do i=1, n_mesh_maps
      tmp_str = mesh_maps(i)
      first_mesh  = tmp_str(:index(tmp_str,':')-1)
      second_mesh = tmp_str(index(tmp_str,':')+1:)
      write(requested_mesh_maps(j),   '(A)') &
          trim(first_mesh)//':'//trim(second_mesh)
      write(requested_mesh_maps(j+1), '(A)') &
          trim(second_mesh)//':'//trim(first_mesh)
      j=j+2
    end do
  end if


  !===================================================================
  ! 6.0 Generate each mesh and any maps associated with it
  !     where the mesh is that source.
  !===================================================================
  call log_event( "Generating cubed-sphere mesh(es):", LOG_LEVEL_INFO )


  ! 6.1 Generate objects which know how to generate each requested
  !     unique mesh.
  allocate( cpp      (n_meshes) )
  allocate( mesh_gen (n_meshes) )
  allocate( ncells   (n_meshes) )
  allocate( ugrid_2d (n_meshes) )


  ! 6.2 Assign temporary arrays for target meshes in requested maps
  if (n_mesh_maps > 0) then
    if (allocated( target_mesh_names_tmp)) deallocate( target_mesh_names_tmp )
    if (allocated( target_edge_cells_tmp)) deallocate( target_edge_cells_tmp )
    allocate( target_mesh_names_tmp(n_mesh_maps*2) )
    allocate( target_edge_cells_tmp(n_mesh_maps*2) )
  end if


  if (do_rotate) then
    write(log_scratch_space, '(A)') &
       '  Rotation of mesh requested with: '
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )
    write(log_scratch_space, '(A,F6.1)') &
       '  New North lat: ', lat_north
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )
    write(log_scratch_space, '(A,F6.1)') &
       '  New North lon: ', lon_north
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )
    write(log_scratch_space, '(A,F6.1)') &
       '  Rotation about north: ', rotate_angle
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )
  end if


  do i=1, n_meshes
    cpp(i)    = edge_cells(i)*edge_cells(i)
    ncells(i) = cpp(i)*npanels

    ! Only smooth this mesh if it is the mesh with the highest
    ! number of cells in the chain.

    ! NOTE: Do not consider smoothing unless all requested
    !       mesh edge cells are a factor of the maximum requested
    !       mesh edge cell.
    if (edge_cells(i) == max_res) then
      nsmooth = smooth_passes
    else
      nsmooth = 0
    end if

    write(log_scratch_space,'(2(A,I0),A)')           &
        '  Creating Mesh: '// trim(mesh_names(i)) // &
        '(',edge_cells(i),',',edge_cells(i),')'
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO)

    write(log_scratch_space,'(A,I0)') &
        '    Smoothing passes: ', nsmooth
    call log_event( trim(log_scratch_space), LOG_LEVEL_INFO)



    ! 6.3 Get any target mappings requested for this mesh
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

    ! 6.4 Call generation stratedgy
    if (n_targets == 0 .or. n_meshes == 1 ) then
      ! No mesh maps required for this mesh
      mesh_gen(i) = gencube_ps_type( mesh_name      = mesh_names(i), &
                                     edge_cells     = edge_cells(i), &
                                     nsmooth        = nsmooth,       &
                                     do_rotate      = do_rotate,     &
                                     lat_north      = lat_north,     &
                                     lon_north      = lon_north,     &
                                     rotate_angle   = rotate_angle,  &
                                     stretch_factor = stretch_factor )

    else if (n_meshes > 1) then

      if (allocated(target_mesh_names)) deallocate(target_mesh_names)
      if (allocated(target_edge_cells)) deallocate(target_edge_cells)

      allocate(target_mesh_names(n_targets))
      allocate(target_edge_cells(n_targets))
      target_mesh_names(:) = target_mesh_names_tmp(:n_targets)
      target_edge_cells(:) = target_edge_cells_tmp(:n_targets)

      write(log_scratch_space,'(A,I0)') '    Maps to: '
      do j=1, n_targets
        write(log_scratch_space,'(2(A,I0),A)') trim(log_scratch_space)//' '//&
              trim(target_mesh_names(j))//&
              '(',target_edge_cells(j),',',target_edge_cells(j),')'
      end do
      call log_event( trim(log_scratch_space), LOG_LEVEL_INFO)

      mesh_gen(i) = gencube_ps_type(                         &
                        mesh_name=mesh_names(i),             &
                        edge_cells=edge_cells(i),            &
                        target_mesh_names=target_mesh_names, &
                        target_edge_cells=target_edge_cells, &
                        nsmooth=nsmooth,                     &
                        do_rotate=do_rotate,                 &
                        lat_north=lat_north,                 &
                        lon_north=lon_north,                 &
                        rotate_angle=rotate_angle,           &
                        stretch_factor=stretch_factor )
    else
      write(log_scratch_space, "(A,I0,A)") &
           '  Number of unique meshes is negative [', n_meshes,']'
      call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR)
    end if

    ! Pass the cubesphere generation object to the ugrid file writer
    call ugrid_2d(i)%set_by_generator(mesh_gen(i))
    if (allocated(target_mesh_names)) deallocate(target_mesh_names)

  end do  ! n_meshes

  if ( allocated(target_edge_cells) ) deallocate(target_edge_cells)
  if ( allocated(ncells)            ) deallocate(ncells)

  call log_event( "...generation complete.", LOG_LEVEL_INFO )


  !===================================================================
  ! 9.0 Write out to ugrid file
  !===================================================================
  do i=1, n_meshes

    if (.not. allocated(ugrid_file)) allocate(ncdf_quad_type::ugrid_file)

    call ugrid_2d(i)%set_file_handler(ugrid_file)

    if ( i==1 ) then
      call ugrid_2d(i)%write_to_file( trim(mesh_filename) )
    else
      call ugrid_2d(i)%append_to_file( trim(mesh_filename) )
    end if

    inquire(file=trim(mesh_filename), size=fsize)
    write( log_scratch_space, '(A,I0,A)')                 &
        'Adding mesh (' // trim(mesh_names(i)) //         &
        ') to ' // trim(adjustl(mesh_filename)) // ' - ', &
        fsize, ' bytes written.'

    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    if (allocated(ugrid_file)) deallocate(ugrid_file)

  end do  ! n_meshes


  !===================================================================
  ! 10.0 Clean up and Finalise
  !===================================================================
  if ( allocated( ncells   ) ) deallocate (ncells)
  if ( allocated( cpp      ) ) deallocate (cpp)
  if ( allocated( mesh_gen ) ) deallocate (mesh_gen)

  if ( allocated( requested_mesh_maps   ) ) deallocate (requested_mesh_maps)
  if ( allocated( target_edge_cells     ) ) deallocate (target_edge_cells)
  if ( allocated( target_mesh_names     ) ) deallocate (target_mesh_names)
  if ( allocated( target_edge_cells_tmp ) ) deallocate (target_edge_cells_tmp)
  if ( allocated( target_mesh_names_tmp ) ) deallocate (target_mesh_names_tmp)

  call finalise_comm()

  call finalise_logging()

  call final_configuration()

end program cubedsphere_mesh_generator
