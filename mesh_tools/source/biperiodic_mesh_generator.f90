!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @mainpage Biperiodic mesh generator
!>
!> @brief   Utility to generate a biperiodic surface mesh and write to a file
!>          which conforms to the UGRID format convention.
!> @details Usage:
!>
!>          biperiodic_mesh_generator <filename>
!>          filename - Controlling namelist file
!>
!-----------------------------------------------------------------------------
program biperiodic_mesh_generator

  use biperiodic_mesh_generator_config_mod,                                    &
                         only: read_biperiodic_mesh_generator_namelist,        &
                               postprocess_biperiodic_mesh_generator_namelist, &
                               edge_cells_x, edge_cells_y, domain_x, domain_y, &
                               nmeshes, mesh_names, mesh_filename

  use cli_mod,           only: get_initial_filename
  use constants_mod,     only: i_def, str_def, str_long, l_def, imdi
  use mpi_mod,           only: initialise_comm, store_comm, finalise_comm, &
                               get_comm_size, get_comm_rank
  use genbiperiodic_mod, only: genbiperiodic_type
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

  character(:), allocatable :: filename
  integer(i_def)            :: namelist_unit

  type(genbiperiodic_type), allocatable :: bpgen(:)
  type(ugrid_2d_type),      allocatable :: ugrid_2d(:)
  class(ugrid_file_type),   allocatable :: ugrid_file

  integer(i_def) :: fsize
  integer(i_def) :: target
  integer(i_def) :: targets
  integer(i_def) :: n_unique_meshes
  integer(i_def), allocatable :: unique_edge_cells_x(:)
  integer(i_def), allocatable :: unique_edge_cells_y(:)
  integer(i_def), allocatable :: unique_target_edge_cells_x(:)
  integer(i_def), allocatable :: unique_target_edge_cells_y(:)

  character(str_def), allocatable :: target_mesh_names(:)
  character(str_def), pointer :: unique_target_mesh_names(:) => null()
  character(str_def), pointer :: unique_mesh_names(:)        => null()
  character(str_def) :: test_str, ref_str

  ! Switches
  logical(l_def) :: l_found = .false.

  ! Parametes
  integer(i_def), parameter :: npanels = 1
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
  call initialise_logging(local_rank, total_ranks, "biperiodic")

  !===================================================================
  ! 3.0 Read in the control namelists from file
  !===================================================================
  call get_initial_filename( filename )
  namelist_unit = open_file( filename )
  call read_biperiodic_mesh_generator_namelist( namelist_unit, 0 )
  call postprocess_biperiodic_mesh_generator_namelist( )
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
  if ( ANY(edge_cells_x == imdi) .OR. &
       ANY(edge_cells_y == imdi) ) then
    write(log_scratch_space,'(A)') &
       'Missing data in namelist variable, edge_cells_x/edge_cells_y'
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if

  

  !===================================================================
  ! 5.0 Get the unique edge_cells list as meshes could appear more than
  !     once the chain
  !===================================================================
  unique_mesh_names => remove_duplicates(mesh_names)
  n_unique_meshes = size(unique_mesh_names)

  allocate(unique_edge_cells_x(n_unique_meshes))
  allocate(unique_edge_cells_y(n_unique_meshes))

  do i=1, n_unique_meshes
    l_found=.false.
    test_str=''
    ref_str=''
    unique_edge_cells_x(i) = imdi
    unique_edge_cells_y(i) = imdi

    do j=1, nmeshes
      if (trim(mesh_names(j)) == trim(unique_mesh_names(i))) then
        write(test_str, '(I0,A,I0)') &
             edge_cells_x(j),'x',edge_cells_y(j)
        if (l_found) then
          if (trim(ref_str) /= trim(test_str)) then
            write(log_scratch_space,'(A)')        &
                'All instances of a mesh tag "'// &
                trim(mesh_names(j))//             &
                '" must have the same mesh specification.'
            call log_event(log_scratch_space, LOG_LEVEL_ERROR)
          end if
        else
          ref_str = test_str
          unique_edge_cells_x(i) = edge_cells_x(j)
          unique_edge_cells_y(i) = edge_cells_y(j)
          l_found=.true.
        end if

      end if
    end do
  end do


  !===================================================================
  ! 6.0 Report/Check what the code thinks is requested by user
  !===================================================================
  call log_event( "Generating ordered bi-periodic mesh(es):", &
                  LOG_LEVEL_INFO )
  tmp_str1=''
  tmp_str2=''
  do i=1, nmeshes
    write(tmp_str1,'(A,2(I0,A))')          &
        trim(adjustl(mesh_names(i)))//'(', &
        edge_cells_x(i), ',', edge_cells_y(i), ')'
    if (i==1) then
      tmp_str2 = trim(adjustl(tmp_str1))
    else 
      tmp_str2 = trim(adjustl(tmp_str2))//'-'//trim(adjustl(tmp_str1))
    end if
  end do
  write(log_scratch_space, '(A)') &
      '  Names(edge_cells): '//trim(tmp_str2)
  call log_event( trim(log_scratch_space), LOG_LEVEL_INFO )


  ! Create objects to manipulate UGRID conforming NetCDF file
  allocate( bpgen    (n_unique_meshes) )
  allocate( ugrid_2d (n_unique_meshes) )


  !===================================================================
  ! 7.0 Determine which targets meshes are required for each
  !     unique mesh.
  !===================================================================
  allocate( target_mesh_names(max_n_targets) )

  do i=1, n_unique_meshes

    target_mesh_names(:) = ''
    if (n_unique_meshes > 1) then

      ! From the requested chain, get all the edge_cell values for
      ! all the other meshes this mesh need to map to
      target = 1
      do j=1, nmeshes
        if (unique_mesh_names(i) == mesh_names(j)) then
          if (j==1) then
            target_mesh_names(target) = trim(mesh_names(j+1))
            target=target+1
          else if (j==nmeshes) then
            target_mesh_names(target) = trim(mesh_names(j-1))
            target=target+1
          else
            target_mesh_names(target)   = trim(mesh_names(j+1))
            target_mesh_names(target+1) = trim(mesh_names(j-1))
            target=target+2
          end if
        end if
      end do


      unique_target_mesh_names => remove_duplicates(target_mesh_names)
      targets =size(unique_target_mesh_names)

      allocate(unique_target_edge_cells_x(targets))
      allocate(unique_target_edge_cells_y(targets))
      target_mesh_names=''

      do k=1, targets
        do j=1,nmeshes
          if (unique_target_mesh_names(k) == mesh_names(j) ) then
            unique_target_edge_cells_x(k) = edge_cells_x(j)
            unique_target_edge_cells_y(k) = edge_cells_y(j)
          end if
        end do
      end do


      tmp_str1=''
      tmp_str2=''
      do j=1, targets
        write(tmp_str1,'(A,2(I0,A))')                      &
          trim(adjustl(unique_target_mesh_names(j)))//'(', &
          unique_target_edge_cells_x(j), ',',              &
          unique_target_edge_cells_y(j), ')'
        if (j==1) then
          tmp_str2 = trim(adjustl(tmp_str1))
        else 
          tmp_str2 = trim(adjustl(tmp_str2))//', '//trim(adjustl(tmp_str1))
        end if
      end do

      write(log_scratch_space,'(A,2(I0,A))')                     &
          '  Creating Mesh: '// trim(unique_mesh_names(i))//'(', &
                                unique_edge_cells_x(i), ',',     &
                                unique_edge_cells_y(i), ')'
      call log_event( trim(log_scratch_space), LOG_LEVEL_INFO)

      bpgen(i) = genbiperiodic_type( mesh_name=unique_mesh_names(i),                 & 
                                     edge_cells_x=unique_edge_cells_x(i),            &
                                     edge_cells_y=unique_edge_cells_y(i),            &
                                     target_mesh_names=unique_target_mesh_names,     &
                                     target_edge_cells_x=unique_target_edge_cells_x, &
                                     target_edge_cells_y=unique_target_edge_cells_y, &
                                     domain_x=domain_x,                              &
                                     domain_y=domain_y )




    else if ( n_unique_meshes == 1 ) then

      ! Only 1 mesh requested, so it must be the prime mesh
      ! and so no optional target_ndivs required
      bpgen(i) = genbiperiodic_type( mesh_name=unique_mesh_names(i),      & 
                                     edge_cells_x=unique_edge_cells_x(i), &
                                     edge_cells_y=unique_edge_cells_y(i), &
                                     domain_x=domain_x,                   &
                                     domain_y=domain_y )



    else
      write(log_scratch_space, "(A,I0,A)") &
           '  Number of unique meshes is negative [', n_unique_meshes,']'
      call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR)
    end if

    ! Pass the cubesphere generation object to the ugrid file writer
    call ugrid_2d(i)%set_by_generator(bpgen(i))

    if (allocated(unique_target_edge_cells_x)) deallocate(unique_target_edge_cells_x)
    if (allocated(unique_target_edge_cells_y)) deallocate(unique_target_edge_cells_y)
  end do

  if (allocated(target_mesh_names)) deallocate(target_mesh_names)

  call log_event( "...generation complete.", LOG_LEVEL_INFO )


  !===================================================================
  ! 8.0 Now the write out mesh to the NetCDF file
  !===================================================================
  do i=1, n_unique_meshes

    if (.not. allocated(ugrid_file)) allocate(ncdf_quad_type::ugrid_file)

    call ugrid_2d(i)%set_file_handler(ugrid_file)

    if (i==1) then
      call ugrid_2d(i)%write_to_file( trim(mesh_filename) )
    else
      call ugrid_2d(i)%append_to_file( trim(mesh_filename) )
    end if

    inquire(file=mesh_filename, size=fsize)
    write( log_scratch_space, '(A,I0,A)')                 &
        'Adding mesh (' // trim(unique_mesh_names(i)) //  &
        ') to ' // trim(adjustl(mesh_filename)) // ' - ', &
        fsize, ' bytes written.'

    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    if (allocated(ugrid_file)) deallocate(ugrid_file)

  end do

  call finalise_comm()

  if ( allocated( bpgen ) ) deallocate (bpgen)

  ! Finalise the logging system
  call finalise_logging()

end program biperiodic_mesh_generator
