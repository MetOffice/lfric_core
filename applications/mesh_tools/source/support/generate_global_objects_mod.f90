!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Support routine to generate global mesh objects for mesh generators.
!> @details Mesh generators were initially based on outputting meshes from
!>          ugrid_2d_types. Longer term, global mesh objects are to be used.
!>          This routine produces the required global mesh objects for
!>          cubedsphere and planar mesh generators in preparation for
!>          partitioning.
module generate_global_objects_mod

  use constants_mod,                  only: i_def, r_def, str_def
  use global_mesh_collection_mod,     only: global_mesh_collection_type
  use global_mesh_map_collection_mod, only: global_mesh_map_collection_type
  use global_mesh_mod,                only: global_mesh_type
  use global_mesh_map_mod,            only: global_mesh_map_type
  use ugrid_2d_mod,                   only: ugrid_2d_type
  use ugrid_mesh_data_mod,            only: ugrid_mesh_data_type

  implicit none

  private
  public :: generate_global_objects

contains

!-----------------------------------------------------------------------------
!> @brief   Generates the required global mesh objects as specified by the
!>          generated ugrid_2d_types.
!> @details Global mesh objects are returned in the global mesh collection
!>          argument. This routine allows for an lbc mesh that may be
!>          requested with planar mesh creation.
!>
!> @param[in]   ugridders         Array of generated ugrid_2d_types
!>                                specifying each mesh
!> @param[in]   global_mesh_bank  Global mesh collection to contain
!>                                the requested global mesh objects
!> @param[in]   ugridder_lbc      Optional, Ugridder object specifying
!>                                the lbc-mesh
!-----------------------------------------------------------------------------

subroutine generate_global_objects( ugridders,        &
                                    global_mesh_bank, &
                                    ugridder_lbc )

  implicit none


  type(ugrid_2d_type),               intent(in)                 :: ugridders(:)
  type(global_mesh_collection_type), intent(inout), allocatable :: global_mesh_bank
  type(ugrid_2d_type),               intent(in), optional       :: ugridder_lbc


  ! Local variables
  type(ugrid_mesh_data_type) :: ugrid_data
  type(global_mesh_type)     :: global_mesh

  type(global_mesh_type), pointer :: source_global_mesh_ptr => null()
  type(global_mesh_type), pointer :: target_global_mesh_ptr => null()

  type(global_mesh_map_collection_type), pointer :: ugrid_mesh_maps_ptr         => null()
  type(global_mesh_map_collection_type), pointer :: source_global_mesh_maps_ptr => null()

  type(global_mesh_map_type), pointer :: ugrid_mesh_map_ptr => null()

  character(str_def) :: source_name, lbc_name

  integer(i_def), allocatable :: cell_map(:,:,:)

  integer(i_def) :: i, target
  integer(i_def) :: n_maps, n_meshes

  ! In order to obtain partititon informtion, the local
  ! mesh objects must be generated before being output
  ! to file.
  !
  ! For consistency global mesh objects will also be able
  ! to be output to file.

  !===========================================
  ! Generate requested global mesh objects
  !===========================================
  ! Global meshes are generated and stored in the global
  ! mesh bank. InterGrid maps are then transferred
  ! to the respective global meshes.

  n_meshes = size( ugridders )

  ! 1.0 Create the requested global mesh objects.
  do i=1, n_meshes

    call ugrid_data%set_by_ugrid_2d( ugridders(i) )
    global_mesh = global_mesh_type( ugrid_data )
    call ugrid_data%clear()
    call global_mesh_bank%add_new_global_mesh(global_mesh)

  end do


  ! 1.1 Create global lbc mesh object (only for planar meshes).
  if (present(ugridder_lbc)) then
    call ugrid_data%set_by_ugrid_2d( ugridder_lbc )
    global_mesh = global_mesh_type( ugrid_data )
    call ugrid_data%clear()
    call global_mesh_bank%add_new_global_mesh(global_mesh)
  end if


  ! 2.0 Add in intergrid maps.
  do i=1, n_meshes

    call ugridders(i)%get_metadata( mesh_name=source_name )
    source_global_mesh_ptr => global_mesh_bank%get_global_mesh( source_name )
    n_maps = source_global_mesh_ptr%get_nmaps()

    if (n_maps > 0) then

      call ugridders(i)%get_global_mesh_maps(ugrid_mesh_maps_ptr)
      source_global_mesh_maps_ptr => source_global_mesh_ptr%get_mesh_maps()

      do target=1, n_maps
        ugrid_mesh_map_ptr => ugrid_mesh_maps_ptr%get_global_mesh_map(1,target+1)
        call ugrid_mesh_map_ptr%get_cell_map(cell_map)
        call source_global_mesh_maps_ptr%add_global_mesh_map(1,target+1,cell_map)
      end do

    end if

  end do ! n_meshes


  ! 2.1 Add lbc-map to parent planar mesh.
  if ( present(ugridder_lbc) ) then

    call ugridder_lbc%get_metadata(mesh_name=lbc_name)
    source_global_mesh_ptr => global_mesh_bank%get_global_mesh(lbc_name)

    call ugridder_lbc%get_global_mesh_maps(ugrid_mesh_maps_ptr)
    source_global_mesh_maps_ptr => source_global_mesh_ptr%get_mesh_maps()

    ugrid_mesh_map_ptr => ugrid_mesh_maps_ptr%get_global_mesh_map(1,2)
    call ugrid_mesh_map_ptr%get_cell_map(cell_map)
    call source_global_mesh_maps_ptr%add_global_mesh_map(1,2,cell_map)

  end if

  nullify( source_global_mesh_ptr )
  nullify( target_global_mesh_ptr )

  nullify( ugrid_mesh_maps_ptr  )
  nullify( source_global_mesh_maps_ptr )
  nullify( ugrid_mesh_map_ptr   )

  if (allocated( cell_map )) deallocate( cell_map )

end subroutine generate_global_objects

end module generate_global_objects_mod
