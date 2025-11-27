!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Support routine to generate local mesh objects for mesh generators.
module generate_local_objects_mod

  use constants_mod,                  only: i_def, l_def, str_def
  use local_mesh_mod,                 only: local_mesh_type
  use local_mesh_collection_mod,      only: local_mesh_collection_type
  use local_mesh_map_collection_mod,  only: local_mesh_map_collection_type
  use log_mod,                        only: log_event, log_scratch_space, &
                                            log_level_error,              &
                                            log_level_debug
  use global_mesh_mod,                only: global_mesh_type
  use global_mesh_collection_mod,     only: global_mesh_collection_type
  use global_mesh_map_collection_mod, only: global_mesh_map_collection_type
  use global_mesh_map_mod,            only: global_mesh_map_type
  use partition_mod,                  only: partition_type, partitioner_interface
  use panel_decomposition_mod,        only: panel_decomposition_type, &
                                            calc_mapping_factor

  implicit none

  private
  public :: generate_local_objects

contains

!-----------------------------------------------------------------------------
!> @brief   Generates the local mesh objects by based on partitions
!>          of specified global mesh objects.
!> @details Partitioned global mesh objects are returned as local mesh objects
!>          (1 per partition). Additionally, for planar meshes, local LBC mesh
!>          objects are produce if an lbc_parent_name is specified.
!>
!> @param[in, out] local_mesh_bank    Collection for local meshes.
!> @param[in]      global_mesh_bank   Collection of global meshes to partition.
!> @param[in]      mesh_names         Names of meshes to partition.
!> @param[in]      partition_id       Partition number being generated.

!> @param[in]      n_partitions       Total number of partitions for each mesh.
!> @param[in]      max_stencil_depth  Maximum stencil depth that the partitions
!>                                    should support.
!> @param [in]     decomposition      Object containing decomposition parameters
!>                                    and method
!> @param[in]      partitioner        Partitioner to apply on meshes.
!> @param[in]      lbc_parent_name    Optional, Name of mesh to produce
!>                                              corresponding local LBC meshes
!>                                              (Planar meshes only).
!-----------------------------------------------------------------------------
subroutine generate_local_objects( local_mesh_bank,       &
                                   global_mesh_bank,      &
                                   mesh_names,            &
                                   partition_id,          &
                                   n_partitions,          &
                                   max_stencil_depth,     &
                                   generate_inner_halos,  &
                                   decomposition,         &
                                   partitioner,           &
                                   lbc_parent_name )

  implicit none

  type(local_mesh_collection_type),  intent(inout) :: local_mesh_bank
  type(global_mesh_collection_type), intent(in)    :: global_mesh_bank

  character(str_def), intent(in) :: mesh_names(:)
  integer(i_def),     intent(in) :: partition_id
  integer(i_def),     intent(in) :: n_partitions
  integer(i_def),     intent(in) :: max_stencil_depth
  logical(l_def),     intent(in) :: generate_inner_halos

  class(panel_decomposition_type),   intent(in) :: decomposition
  procedure(partitioner_interface),  intent(in), pointer :: partitioner

  character(str_def), optional,      intent(in) :: lbc_parent_name

  ! Local variables
  type(local_mesh_type) :: local_mesh
  type(partition_type)  :: partition

  type(global_mesh_type), pointer :: source_global_mesh_ptr
  type(global_mesh_type), pointer :: target_global_mesh_ptr
  type(local_mesh_type),  pointer :: source_local_mesh_ptr
  type(local_mesh_type),  pointer :: target_local_mesh_ptr

  type(global_mesh_map_collection_type), pointer :: global_mesh_maps_ptr
  type(local_mesh_map_collection_type),  pointer :: local_mesh_maps_ptr

  type(global_mesh_map_type), pointer :: global_mesh_map_ptr

  integer(i_def), allocatable :: cell_map(:,:,:)
  integer(i_def), allocatable :: local_map_gid(:,:,:)
  integer(i_def), allocatable :: local_map_lid(:,:,:)
  integer(i_def), allocatable :: source_local_gids(:)

  integer(i_def) :: source_local_mesh_ncells

  character(str_def), allocatable  :: target_names(:)

  character(str_def) :: source_name, name, lbc_name

  integer(i_def) :: n_meshes, n_maps, map_xcells, map_ycells, local_id
  integer(i_def) :: i, p, q, target, local_cell, cell_global_id, mapping_factor

  ! Local variables for LBC meshes
  type(local_mesh_type)       :: local_lbc_mesh
  integer(i_def), allocatable :: cell_lbc_lam_map(:)

  n_meshes = size(mesh_names)

  !====================================================================
  ! Create requested range of partitions for each mesh and store in
  ! the local mesh bank.
  !====================================================================
  do i=1, n_meshes

    source_name = mesh_names(i)
    source_global_mesh_ptr => global_mesh_bank%get_global_mesh(source_name)

    mapping_factor = calc_mapping_factor(global_mesh_bank, source_global_mesh_ptr)

    call log_event( 'Partitioning mesh:'//trim(source_name), log_level_debug )

    write( name,'(A,I0)' ) trim(adjustl(source_name))//'_', partition_id
    call log_event( 'Initialising local mesh object '//trim(name), &
                    log_level_debug )

    partition = partition_type( source_global_mesh_ptr, &
                                partitioner,            &
                                decomposition,          &
                                max_stencil_depth,      &
                                generate_inner_halos,   &
                                partition_id,           &
                                n_partitions,           &
                                mapping_factor )

    call local_mesh%initialise( source_global_mesh_ptr, &
                                partition, name=name )

    local_id = local_mesh_bank%add_new_local_mesh(local_mesh)

  end do ! i

  !====================================================================
  ! Partition an LBC mesh if an lbc_parent_name is provided.
  !====================================================================
  if ( present(lbc_parent_name) ) then

    ! Source mesh is the LBC mesh.
    ! Target mesh is the LBC parent mesh.

    lbc_name = trim(lbc_parent_name)//'-lbc'
    call log_event( 'Partioning LBC mesh:'//trim(lbc_name), log_level_debug )

    ! Set-up Global mesh pointers.
    !----------------------------------
    source_global_mesh_ptr => global_mesh_bank%get_global_mesh(lbc_name)
    target_global_mesh_ptr => global_mesh_bank%get_global_mesh(lbc_parent_name)

    global_mesh_maps_ptr => source_global_mesh_ptr%get_mesh_maps()
    global_mesh_map_ptr  => global_mesh_maps_ptr%get_global_mesh_map(1,2)

    ! Initialise local LBC mesh objects.
    !----------------------------------------
    write( name,'(A,I0)' ) trim(lbc_parent_name)//'_', partition_id
    call log_event( 'Initialising local LBC mesh object using '// &
                     trim(name), log_level_debug )

    target_local_mesh_ptr => local_mesh_bank%get_local_mesh( name )

    call local_lbc_mesh%initialise( target_global_mesh_ptr,    &  ! Global LAM
                                    source_global_mesh_ptr,    &  ! Global LBC
                                    target_local_mesh_ptr,     &  ! Local LAM
                                    name=name,                 &
                                    cell_map=cell_lbc_lam_map )

    write( log_scratch_space,'(2(I03,A))' )            &
           local_lbc_mesh%get_last_edge_cell(),         &
           ' cells in partition ', partition_id, '('//trim(name)// &
           '), map to lbc-mesh'
    call log_event( log_scratch_space, log_level_debug )

    ! The partition is determined by the parent LAM and so may
    ! occur over a location where there are no corresponding LBC
    ! cells. Only add to collection if the partition contains LBC
    ! cells.
    !
    ! Note: This may be moved later but test here for now
    ! Add the cell LBC-LAM map to the local mesh. Done here as
    ! when initialised from within the local mesh object the map
    ! appears to go out of scope. Should be investigated at later
    ! date.
    if ( local_lbc_mesh%get_last_edge_cell() > 0 ) then

      local_id = local_mesh_bank%add_new_local_mesh(local_lbc_mesh)
      source_local_mesh_ptr => local_mesh_bank%get_local_mesh(local_id)

      if (allocated(cell_map)) deallocate(cell_map)
      allocate(cell_map (1,1,source_local_mesh_ptr%get_last_edge_cell()))
      cell_map = reshape( cell_lbc_lam_map, &
                          [1,1,source_local_mesh_ptr%get_last_edge_cell()] )
      local_mesh_maps_ptr => source_local_mesh_ptr%get_mesh_maps()
      call local_mesh_maps_ptr%add_local_mesh_map(1,2,cell_map)

    end if

  end if  ! lbc_parent_name


  !====================================================================
  ! Create LiD -> LiD intergrid maps for local meshes -> target
  ! meshes on same partition over the requested range of partitions
  ! for each mesh. Assign the resultant integrid maps with the
  ! source local mesh.
  !====================================================================
  do i=1, n_meshes

    source_name = mesh_names(i)
    source_global_mesh_ptr => global_mesh_bank%get_global_mesh(source_name)

    n_maps = source_global_mesh_ptr%get_nmaps()

    if (n_maps > 0) then

      global_mesh_maps_ptr => source_global_mesh_ptr%get_mesh_maps()
      call source_global_mesh_ptr%get_target_mesh_names(target_names)

      do target=1, n_maps

        target_global_mesh_ptr => &
                     global_mesh_bank%get_global_mesh(target_names(target))

        global_mesh_map_ptr =>    &
                     global_mesh_maps_ptr%get_global_mesh_map(1,target+1)
        map_xcells = global_mesh_map_ptr%get_ntarget_cells_per_source_x()
        map_ycells = global_mesh_map_ptr%get_ntarget_cells_per_source_y()

        ! At this point:
        ! global_mesh_map_ptr is the Mesh map object for
        ! Source => Target global meshes.

        ! Extract Source/Target Local meshes.
        write( name,'(A,I0)' ) trim(source_name)//'_', partition_id
        source_local_mesh_ptr => local_mesh_bank%get_local_mesh(name)

        source_local_mesh_ncells = source_local_mesh_ptr%get_num_cells_in_layer()
        source_local_gids        = source_local_mesh_ptr%get_all_gid()

        ! Constructed an intergrid cell map for the partition in
        ! Global IDs(GID) and Local IDs(LID).
        if ( allocated(local_map_gid) ) deallocate (local_map_gid)
        if ( allocated(local_map_lid) ) deallocate (local_map_lid)
        allocate( local_map_gid( map_xcells, map_ycells, &
                                 source_local_mesh_ncells ) )
        allocate( local_map_lid( map_xcells, map_ycells, &
                                 source_local_mesh_ncells ) )

        call global_mesh_map_ptr%get_cell_map(                                  &
                                 source_local_gids(1:source_local_mesh_ncells), &
                                 local_map_gid )

        write( name,'(A,I0)' ) trim(target_names(target))//'_', partition_id
        target_local_mesh_ptr => local_mesh_bank%get_local_mesh(name)

        ! Array local_map_gid should contain the integrid mesh map, where the
        ! index is the local cell id on the source partition. The values are
        ! the cell IDs of the target cells which map to the same location,
        ! these target cell IDs respect to the target mesh, i.e. GIDS.
        !
        ! Convert the contents of local_map_gid from GIDs to LIDs, using the
        ! LID/GID map of the target local meshes to give a local source/target
        ! intergrid map in local IDs of the respective meshes.
        do local_cell=1, source_local_mesh_ncells
          do p=1, map_xcells
            do q=1, map_ycells
              cell_global_id = local_map_gid(p,q,local_cell)
              local_map_lid(p,q,local_cell) = &
                  target_local_mesh_ptr%get_lid_from_gid( cell_global_id )
              if (local_map_lid(p,q,local_cell) == -1) then
                local_map_lid(p,q,local_cell) = source_local_mesh_ptr%get_void_cell()
              end if
            end do
          end do
        end do

        ! Now add the cell map to the local source mesh map collection.
        local_mesh_maps_ptr => source_local_mesh_ptr%get_mesh_maps()
        call local_mesh_maps_ptr%add_local_mesh_map( 1, target+1, local_map_lid )

      end do    ! target
    end if    ! n_maps > 0
  end do    ! n_meshes

  !=====================================================================

  nullify( source_local_mesh_ptr )
  nullify( target_local_mesh_ptr )

  nullify( source_global_mesh_ptr )
  nullify( target_global_mesh_ptr )

  nullify( local_mesh_maps_ptr  )
  nullify( global_mesh_maps_ptr )

  nullify( global_mesh_map_ptr  )

  if (allocated( cell_map      )) deallocate( cell_map      )
  if (allocated( local_map_gid )) deallocate( local_map_gid )
  if (allocated( local_map_lid )) deallocate( local_map_lid )

end subroutine generate_local_objects

end module  generate_local_objects_mod
