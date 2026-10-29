!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Support routine to write out local mesh objects to UGRID files
!> @details All local mesh objects across a given partition are written
!>          to a single UGRID file. File names are the requested basename
!>          appended with the relevant partition number.
module write_local_meshes_mod

  use calc_cell_centres_mod, only: calc_cell_centres_global_model, &
                                   calc_cell_centres_regional_model

  use constants_mod, only: i_def, r_def, str_def, &
                           str_max_filename,      &
                           radians_to_degrees,    &
                           degrees_to_radians
  use log_mod,       only: log_event, log_scratch_space, &
                           log_level_info
  use omp_lib,       only: omp_get_thread_num
  use sci_query_mod, only: valid_for_global_model

  ! Derived types
  use global_mesh_mod,                only: global_mesh_type
  use global_mesh_collection_mod,     only: global_mesh_collection_type
  use global_mesh_map_collection_mod, only: global_mesh_map_collection_type
  use global_mesh_map_mod,            only: global_mesh_map_type
  use local_mesh_mod,                 only: local_mesh_type
  use local_mesh_collection_mod,      only: local_mesh_collection_type
  use ncdf_quad_mod,                  only: ncdf_quad_type
  use ugrid_2d_mod,                   only: ugrid_2d_type
  use ugrid_file_mod,                 only: ugrid_file_type

  ! Configuration modules
  use mesh_config_mod,        only: coord_sys, coord_sys_ll, n_meshes,   &
                                    mesh_names
  use partitions_config_mod,  only: n_partitions, partition_range
  use planar_mesh_config_mod, only: edge_cells_x, edge_cells_y,   &
                                    domain_size, create_lbc_mesh, &
                                    lbc_parent_mesh

  implicit none

  private
  public ::  write_local_meshes

contains

!-----------------------------------------------------------------------------
!> @brief   Generates the local mesh objects by based on partitions
!>          of specified global mesh objects.
!> @details Partitioned global mesh objects are returned as local mesh objects
!>          (1 per partition). Additionally, for planar meshes, local LBC mesh
!>          objects are produce if an lbc_parent_name is specified.
!>
!> @param[in] partition_id       Partition number being written.
!> @param[in] global_mesh_bank   Collection of generated global meshes.
!> @param[in] local_mesh_bank    Collection of generated local meshes.
!> @param[in] output_basename    UGRID output file basename.
!-----------------------------------------------------------------------------
subroutine write_local_meshes( partition_id,     &
                               global_mesh_bank, &
                               local_mesh_bank,  &
                               output_basename )

  implicit none

  integer(i_def),                    intent(in) :: partition_id
  type(global_mesh_collection_type), intent(in) :: global_mesh_bank
  type(local_mesh_collection_type),  intent(in) :: local_mesh_bank
  character(str_max_filename),       intent(in) :: output_basename

  ! Local variables
  !-----------------
  class(ugrid_file_type), allocatable :: ugrid_file

  type(ugrid_2d_type) :: ugrid_2d

  type(global_mesh_type), pointer :: global_mesh_ptr
  type(global_mesh_type), pointer :: global_lbc_mesh_ptr
  type(local_mesh_type),  pointer :: local_mesh_ptr
  type(local_mesh_type),  pointer :: local_lbc_mesh_ptr

  type(global_mesh_map_collection_type), pointer :: global_lbc_mesh_maps_ptr
  type(global_mesh_map_type),            pointer :: global_lbc_mesh_map_ptr


  integer(i_def), allocatable :: global_ids(:)
  integer(i_def), allocatable :: local_gids(:)
  integer(i_def), allocatable :: verts_on_cells(:,:)
  integer(i_def), allocatable :: global_lbc_lam_cell_map(:,:,:)

  real(r_def), allocatable :: vert_coords(:,:)
  real(r_def), allocatable :: local_cell_coords(:,:)
  real(r_def), allocatable :: global_cell_coords(:,:)

  character(str_max_filename) :: output_file

  character(str_def) :: partition_of
  character(str_def) :: units_xy(2)
  character(str_def) :: source_name
  character(str_def) :: lbc_name
  character(str_def) :: name
  integer(i_def)     :: fsize, cell, thread_id
  real(r_def)        :: dx
  real(r_def)        :: dy
  real(r_def)        :: factor

  integer(i_def)     :: n_digit
  character(str_def) :: fmt_str, number_str

  ! Counters
  integer(i_def) :: i

  nullify(global_lbc_mesh_maps_ptr)
  nullify(global_lbc_mesh_map_ptr)

  !===================================================================

  if ( allocated(ugrid_file) ) deallocate(ugrid_file)
  call ugrid_2d%clear()

  allocate(ncdf_quad_type::ugrid_file)
  call ugrid_2d%set_file_handler(ugrid_file)

  do i=1, n_meshes

    source_name  = mesh_names(i)
    partition_of = trim(mesh_names(i))//'_partition'
    write(name,'(A,I0)') trim(source_name)//'_', partition_id
    local_mesh_ptr => local_mesh_bank%get_local_mesh(name)

    !---------------------------------------------
    ! Populate the ugrid_2d object with
    ! local mesh object. The name to be
    ! written will need to refer to the
    ! base mesh name.
    !---------------------------------------------
    call local_mesh_ptr%as_ugrid_2d(ugrid_2d)
    call ugrid_2d%set_metadata( mesh_name    = mesh_names(i), &
                                partition_of = partition_of )

    ! Create cell coords for output file
    !---------------------------------------------
    ! Extract corresponding cell GIDs for
    ! cells on this local mesh.
    !---------------------------------------------
    local_gids = local_mesh_ptr%get_all_gid()

    ! Calculate cell centre coordinates from
    ! the global mesh node coordinates.
    !---------------------------------------------
    global_mesh_ptr => global_mesh_bank%get_global_mesh(source_name)
    verts_on_cells  =  global_mesh_ptr%get_vert_on_all_cells()
    call global_mesh_ptr%get_vert_coords( vert_coords )

    if (allocated(global_cell_coords)) deallocate(global_cell_coords)
    allocate(global_cell_coords(2,size(verts_on_cells,2)))

    if ( .not. valid_for_global_model(global_mesh_ptr) ) then
      units_xy = global_mesh_ptr%get_coord_units()

      ! This is only required to calculated the cell centres for
      ! planar meshes. In future cell centres should be calculated
      ! using the nodes around the cell. Then numbers of edges cells
      ! will not be needed.
      dx = domain_size(1) / edge_cells_x(i)
      dy = domain_size(2) / edge_cells_y(i)

      if ( (coord_sys == coord_sys_ll)      .and. &
           (trim(units_xy(1)) == 'radians') .and. &
           (trim(units_xy(2)) == 'radians') ) then

        ! All values should be in as radians
        ! for lon-lat.
        dx = dx * degrees_to_radians
        dy = dy * degrees_to_radians
      end if

      call calc_cell_centres_regional_model( dx, dy,         &
                                             verts_on_cells, &
                                             vert_coords,    &
                                             global_cell_coords )
    else
      call calc_cell_centres_global_model( verts_on_cells, &
                                           vert_coords,    &
                                           global_cell_coords )
    end if

    ! Add cell centre coordinates to
    ! the ugrid_2d object.
    !---------------------------------------------
    if ( allocated(local_cell_coords) ) deallocate( local_cell_coords )
    allocate( local_cell_coords(2,size(local_gids)) )

    do cell=1, size(local_gids)
      local_cell_coords(:,cell) = global_cell_coords(:,local_gids(cell))
    end do

    if (coord_sys == coord_sys_ll) then
      ! Cell coords returned as radians, need to
      ! convert to degrees for output to file.
      factor = radians_to_degrees
    else
      factor = 1.0_r_def
    end if
    call ugrid_2d%set_coords( face_coords=factor*local_cell_coords )


    !---------------------------------------------
    ! Add ugrid_2d mesh to the output file
    ! for this partition.
    !---------------------------------------------
    n_digit = int(log10(real(n_partitions))) + 1
    write(fmt_str, '(A,I0,A,I0,A)') "(I", n_digit, ".", n_digit, ")"
    write(number_str, fmt_str) partition_id
    write(output_file, '(A, "_", A, "-", I0, ".nc")') &
        trim(output_basename), trim(number_str), n_partitions

!$omp critical
    thread_id = omp_get_thread_num()

    if (i==1) then
      call ugrid_2d%write_to_file( trim(output_file) )
    else
      call ugrid_2d%append_to_file( trim(output_file) )
    end if

    inquire(file=output_file, size=fsize)
    write( log_scratch_space, '(2(A,I0),A)' )     &
        'Thread ', thread_id, ' adding mesh (' // &
        trim(source_name) // ') to ' //           &
        trim(output_file) // ' - ', fsize,        &
        ' bytes written.'
    call log_event( log_scratch_space, log_level_info )
!$omp end critical

    !---------------------------------------------
    ! Add the lbc-mesh to the output file for
    ! this partition if required.
    !---------------------------------------------
    if ( create_lbc_mesh .and. &
       ( trim(mesh_names(i)) == trim(lbc_parent_mesh) ) ) then

      ! Add the local lbc mesh to the output file.
      ! Extract lbc mesh from the local mesh collection
      ! and write to file.
      write(name,'(A,I0,A)') trim(source_name)//'_', partition_id,'-lbc'

      local_lbc_mesh_ptr => local_mesh_bank%get_local_mesh(name)

      if (associated(local_lbc_mesh_ptr)) then

        if (.not. allocated(ugrid_file)) then
          allocate(ncdf_quad_type::ugrid_file)
        end if
        call ugrid_2d%set_file_handler( ugrid_file )
        call local_lbc_mesh_ptr%as_ugrid_2d( ugrid_2d )

        name = trim(source_name)//'-lbc'
        partition_of = trim(source_name)//'-lbc_partition'
        call ugrid_2d%set_metadata( mesh_name    = name, &
                                    partition_of = partition_of )


        ! Need to get  global cell coords of LBC local mesh from LBC GLobal mesh
        lbc_name = trim(source_name)//'-lbc'
        global_lbc_mesh_ptr      => global_mesh_bank%get_global_mesh(lbc_name)
        global_lbc_mesh_maps_ptr => global_lbc_mesh_ptr%get_mesh_maps()
        global_lbc_mesh_map_ptr  => global_lbc_mesh_maps_ptr%get_global_mesh_map(1,2)

        if (local_lbc_mesh_ptr%get_last_edge_cell() > 0) then

          call global_lbc_mesh_map_ptr%get_cell_map(global_lbc_lam_cell_map)

          global_ids = local_lbc_mesh_ptr%get_all_gid()

          if (allocated(local_cell_coords)) deallocate(local_cell_coords)
          allocate(local_cell_coords(2, local_lbc_mesh_ptr%get_last_edge_cell()))

          do cell=1, local_lbc_mesh_ptr%get_last_edge_cell()
            local_cell_coords(:,cell) = &
                global_cell_coords(:, global_lbc_lam_cell_map(1,1,global_ids(cell)))
          end do

          call ugrid_2d%set_coords( face_coords=local_cell_coords*factor )

        end if ! LBC on this partition

!$omp critical
        thread_id = omp_get_thread_num()

        call ugrid_2d%append_to_file( trim(output_file) )

        inquire(file=output_file, size=fsize)
        write( log_scratch_space, '(2(A,I0),A)' )         &
            'Thread ', thread_id, ' adding mesh (' //     &
            trim(name) // ') to ' // trim(output_file) // &
            ' - ', fsize, ' bytes written.'
        call log_event( log_scratch_space, log_level_info )
!$omp end critical

      end if ! local mesh pointer associated

    end if ! lbc required

  end do ! i (n_meshes)

end subroutine write_local_meshes

end module write_local_meshes_mod
