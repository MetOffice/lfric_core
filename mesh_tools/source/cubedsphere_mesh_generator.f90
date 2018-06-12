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
  use constants_mod,     only: i_def, imdi, l_def, str_def, str_long
  use cubedsphere_mesh_generator_config_mod,                                    &
                         only: read_cubedsphere_mesh_generator_namelist,        &
                               postprocess_cubedsphere_mesh_generator_namelist, &
                               edge_cells, smooth_passes, nmeshes, mesh_names,  &
                               mesh_filename
  use mpi_mod,           only: initialise_comm, store_comm, finalise_comm, &
                               get_comm_size, get_comm_rank
  use gencube_ps_mod,    only: gencube_ps_type
  use io_utility_mod,    only: open_file, close_file
  use log_mod,           only: initialise_logging, finalise_logging, &
                               log_event, log_set_level, &
                               log_scratch_space, &
                               LOG_LEVEL_INFO, LOG_LEVEL_ERROR
  use ncdf_quad_mod,     only: ncdf_quad_type
  use remove_duplicates_mod, &
                         only: remove_duplicates
  use ugrid_2d_mod,      only: ugrid_2d_type
  use ugrid_file_mod,    only: ugrid_file_type

  implicit none

  integer(i_def) :: comm, total_ranks, local_rank

  character(:), allocatable   :: filename
  integer(i_def)              :: namelist_unit

  integer(i_def),         allocatable :: ncells(:)
  integer(i_def),         allocatable :: cpp(:)
  type(gencube_ps_type),  allocatable :: csgen(:)
  type(ugrid_2d_type),    allocatable :: ugrid_2d(:)
  class(ugrid_file_type), allocatable :: ugrid_file

  integer(i_def) :: fsize
  integer(i_def) :: max_res

  integer(i_def) :: target
  integer(i_def) :: nsmooth
  integer(i_def) :: targets
  integer(i_def) :: n_unique_meshes
  integer(i_def) :: test_int

  integer(i_def),     allocatable :: unique_edge_cells(:)
  integer(i_def),     allocatable :: target_edge_cells(:)
  character(str_def), allocatable :: target_mesh_names(:)

  integer(i_def),     pointer :: unique_target_edge_cells(:) => null()
  character(str_def), pointer :: unique_mesh_names(:)        => null()
  
  ! Switches
  logical(l_def) :: l_found = .false.

  ! Parametes
  integer(i_def), parameter :: npanels = 6
  integer(i_def), parameter :: max_n_targets = 6

  ! Temporary variables
  character(str_long) :: tmp_str1
  character(str_long) :: tmp_str2

  ! Counters
  integer(i_def) :: i, j, k


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
  namelist_unit = open_file( filename )
  call read_cubedsphere_mesh_generator_namelist( namelist_unit, 0 )
  call postprocess_cubedsphere_mesh_generator_namelist()
  call close_file( namelist_unit )
  deallocate( filename )

  !===================================================================
  ! 4.0 Perform some error checks on the namelist inputs
  !===================================================================
  ! 4.1 Check the number of meshes requested.
  if (nmeshes < 1) then
    write(log_scratch_space,'(A,I0,A)') &
       'Invalid number of meshes requested, (',nmeshes,')'
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if

  ! 4.2 Check for missing data.
  if (ANY(edge_cells == imdi)) then
    write(log_scratch_space,'(A)') &
       'Missing data in namelist variable, edge_cells'
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if


  !===================================================================
  ! 5.0 Get the unique mesh_names list as meshes could appear more than
  !     once the chain
  !===================================================================
  unique_mesh_names => remove_duplicates(mesh_names)
  n_unique_meshes = size(unique_mesh_names)

  allocate(unique_edge_cells(n_unique_meshes))

  do i=1, n_unique_meshes
    l_found=.false.
    test_int = imdi
    unique_edge_cells(i) = imdi

    do j=1, nmeshes
      if (trim(mesh_names(j)) == trim(unique_mesh_names(i))) then
        test_int = edge_cells(j)
        if (l_found) then
          if (test_int /= unique_edge_cells(i)) then
            write(log_scratch_space,'(A)')        &
                'All instances of a mesh tag "'// &
                trim(mesh_names(j))//             &
                '" must have the same mesh specification.'
            call log_event(log_scratch_space, LOG_LEVEL_ERROR)
          end if
        else
          unique_edge_cells(i) = test_int
          l_found = .true.
        end if

      end if
    end do
  end do

  max_res = maxval(unique_edge_cells)

  !===================================================================
  ! 5.1  Check that all meshes have edge_cells which are a factor of the
  !      highest edge_cells value.
  !===================================================================
  do i=1, n_unique_meshes
    if (mod(max_res, unique_edge_cells(i)) /= 0) then
      write(log_scratch_space, '(A,I0,A)')                             &
          '  All mesh edge cell values must be a factor of the ' //    &
          'maximum edge cell value[', max_res,']'
      call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
    end if
  end do

  !===================================================================
  ! 6.0 Report/Check what the code thinks is requested by user
  !===================================================================
  call log_event( "Generating ordered cubed-sphere mesh(es):", &
                  LOG_LEVEL_INFO )
  tmp_str1=''
  tmp_str2=''
  do i=1, nmeshes
    write(tmp_str1,'(2(A,I0),A)')            &
        trim(adjustl(mesh_names(i))) // '(', &
        edge_cells(i), ',', edge_cells(i), ')'
    if (i==1) then
      tmp_str2 = trim(adjustl(tmp_str1))
    else 
      tmp_str2 = trim(adjustl(tmp_str2))//'-'//trim(adjustl(tmp_str1))
    end if
  end do
  write(log_scratch_space, '(A)') &
      '  Names(edge_cells): '//trim(tmp_str2)
  call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )
  write(log_scratch_space, '(A,I0)') &
      '  Smoothing passes for maximum resolution: ', smooth_passes
  call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )

  !===================================================================
  ! 7.0 Generate objects which know how to generate each requested
  !     unique mesh.
  !===================================================================
  allocate( cpp      (n_unique_meshes) )
  allocate( csgen    (n_unique_meshes) )
  allocate( ncells   (n_unique_meshes) )
  allocate( ugrid_2d (n_unique_meshes) )



  !===================================================================
  ! 8.0 Determine which targets meshes are required for each
  !     unique mesh.
  !===================================================================
  allocate( target_edge_cells(max_n_targets) )

  do i=1, n_unique_meshes
    cpp(i)    = unique_edge_cells(i)*unique_edge_cells(i)
    ncells(i) = cpp(i)*npanels

    ! Only smooth this mesh if it is the mesh with the highest 
    ! number of cells in the chain
    if (unique_edge_cells(i) == max_res) then
      nsmooth = smooth_passes
    else
      nsmooth = 0
    end if

    if (n_unique_meshes > 1) then

      ! From the requested chain, get all the edge_cell values for
      ! all the other meshes this mesh need to map to
      target_edge_cells = imdi
      target            = 1

      do j=1, nmeshes
        if (unique_edge_cells(i) == edge_cells(j)) then
          if (j==1) then
            target_edge_cells(target) = edge_cells(j+1)

            target=target+1
          else if (j==nmeshes) then
            target_edge_cells(target) = edge_cells(j-1)

            target=target+1
          else
            target_edge_cells(target)   = edge_cells(j+1)
            target_edge_cells(target+1) = edge_cells(j-1)
            target=target+2
          end if
        end if
      end do


      unique_target_edge_cells => remove_duplicates(target_edge_cells)
      targets =size(unique_target_edge_cells)
      allocate(target_mesh_names(targets))
      target_mesh_names=''

      do k=1, targets
        do j=1,nmeshes
          if (unique_target_edge_cells(k) == edge_cells(j) ) then
            target_mesh_names(k) =trim(mesh_names(j))
          end if
        end do
      end do


      tmp_str1=''
      tmp_str2=''
      do j=1, targets
        write(tmp_str1,'(A,I0,A)') &
          trim(adjustl(target_mesh_names(j))) // '(' , unique_target_edge_cells(j), ')'
        if (j==1) then
          tmp_str2 = trim(adjustl(tmp_str1))
        else 
          tmp_str2 = trim(adjustl(tmp_str2))//', '//trim(adjustl(tmp_str1))
        end if
      end do

      write(log_scratch_space,'(2(A,I0),A)')               &
          '  Creating Mesh: '// trim(unique_mesh_names(i)) &
                             //'(',unique_edge_cells(i),',',unique_edge_cells(i),')'
      call log_event( trim(log_scratch_space), LOG_LEVEL_INFO)

      write(log_scratch_space,'(A,I0)') '    Smoothing passes: ', nsmooth
      call log_event( trim(log_scratch_space), LOG_LEVEL_INFO)

      csgen(i) = gencube_ps_type( mesh_name=unique_mesh_names(i),             & 
                                  edge_cells=unique_edge_cells(i),            &
                                  target_mesh_names=target_mesh_names,        &
                                  target_edge_cells=unique_target_edge_cells, &
                                  nsmooth=nsmooth )

    else if ( n_unique_meshes == 1 ) then

      ! Only 1 mesh requested, so it must be the prime mesh
      ! and so no optional target_ndivs required
      csgen(i) = gencube_ps_type( mesh_name  = unique_mesh_names(i), & 
                                  edge_cells = unique_edge_cells(i), &
                                  nsmooth    = nsmooth )

    else
      write(log_scratch_space, "(A,I0,A)") &
           '  Number of unique meshes is negative [', n_unique_meshes,']'
      call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR)
    end if

    ! Pass the cubesphere generation object to the ugrid file writer
    call ugrid_2d(i)%set_by_generator(csgen(i))
    if (allocated(target_mesh_names)) deallocate(target_mesh_names)
  end do

  if ( allocated(target_edge_cells) ) deallocate(target_edge_cells)
  if ( allocated(ncells)            ) deallocate(ncells)

  call log_event( "...generation complete.", LOG_LEVEL_INFO )


  !===================================================================
  ! 9.0 Write out to ugrid file
  !===================================================================
  do i=1, n_unique_meshes
    if (.not. allocated(ugrid_file)) allocate(ncdf_quad_type::ugrid_file)

    call ugrid_2d(i)%set_file_handler(ugrid_file)

    if ( i==1 ) then
      call ugrid_2d(i)%write_to_file( trim(mesh_filename) )
    else
      call ugrid_2d(i)%append_to_file( trim(mesh_filename) )
    end if

    inquire(file=trim(mesh_filename), size=fsize)
    write( log_scratch_space, '(A,I0,A)')                 &
        'Adding mesh (' // trim(unique_mesh_names(i)) //  &
        ') to ' // trim(adjustl(mesh_filename)) // ' - ', &
        fsize, ' bytes written.'

    call log_event( log_scratch_space, LOG_LEVEL_INFO )


    if (allocated(ugrid_file)) deallocate(ugrid_file)
  end do

  call finalise_comm()

  if ( allocated( ncells ) )   deallocate (ncells)
  if ( allocated( cpp ) )      deallocate (cpp)
  if ( allocated( csgen ) )    deallocate (csgen)

  ! Finalise the logging system
  call finalise_logging()

end program cubedsphere_mesh_generator
