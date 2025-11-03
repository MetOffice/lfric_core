!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!> @brief    Routines to assign mesh maps to mesh_type objects.
!> @details  Mesh type objects are initially created with only
!!           the target mesh names for their inter-grid maps.
!!           This module setups up the correct source target meshes
!!           to add intergrid maps to, allowing that the mesh names
!!           may differ to the local mesh from which the inter-grid maps
!!           originate.
module add_mesh_map_mod

  use constants_mod, only: i_def, str_def, cmdi
  use log_mod,       only: log_event,         &
                           log_scratch_space, &
                           LOG_LEVEL_ERROR,   &
                           LOG_LEVEL_INFO


  use extrusion_mod,       only: extrusion_type,           &
                                 uniform_extrusion_type,   &
                                 geometric_extrusion_type, &
                                 quadratic_extrusion_type
  use local_mesh_mod,      only: local_mesh_type
  use mesh_mod,            only: mesh_type
  use sci_query_mod,       only: check_lbc
  use ugrid_mesh_data_mod, only: ugrid_mesh_data_type


  use local_mesh_collection_mod, only: local_mesh_collection
  use mesh_collection_mod,       only: mesh_collection

  implicit none

  private
  public :: assign_mesh_maps, add_mesh_map

contains

!> @brief Assign mesh maps to specied mesh_type objects.
!> @param[in] mesh_names Names of meshes in application mesh
!!                       collection to assign intergrid maps.
subroutine assign_mesh_maps( mesh_names )

  implicit none

  character(str_def), intent(in) :: mesh_names(:)

  character(str_def) :: local_mesh_name
  character(str_def) :: mesh_name_A, mesh_name_B

  character(str_def), allocatable :: target_mesh_names(:)
  character(str_def), allocatable :: local_mesh_names(:)

  type(mesh_type),       pointer :: mesh       => null()
  type(local_mesh_type), pointer :: local_mesh => null()

  integer(i_def) :: i, j, k

  if (size(mesh_names) > 1) then

    !============================================================================
    ! 1.0 Acquire the local mesh object used to create mesh object.
    !============================================================================
    ! Names of the local_mesh_type object may differ from
    ! mesh_type object name provided.
    allocate(local_mesh_names(size(mesh_names)))
    do i=1, size(mesh_names)
      mesh => mesh_collection%get_mesh(mesh_names(i))

      if (.not. associated(mesh)) then
        if (check_lbc(mesh_names(i))) then
          cycle
        end if
      end if

      local_mesh => mesh%get_local_mesh()
      local_mesh_names(i) = local_mesh%get_mesh_name()

    end do

    !============================================================================
    ! 2.0 Assign maps to meshes
    !============================================================================
    ! Intergrid maps reference the local mesh objects that they were read/created
    ! from. The names of these local mesh objects may differ from the name
    ! of the mesh_type object. So the correct pairing of mesh_type names needs
    ! to be found from the connected local meshes.

    do i=1, size(mesh_names)
      mesh_name_A = mesh_names(i)

      ! Find all of the target local meshes associated with the local
      ! mesh object that this mesh object was extruded from.
      mesh => mesh_collection%get_mesh(mesh_name_A)

      if ( .not. associated(mesh)) then
        if (check_lbc(mesh_names(i))) then
          cycle
        end if
      end if

      local_mesh => mesh%get_local_mesh()
      local_mesh_name = local_mesh%get_mesh_name()

      call local_mesh%get_target_mesh_names(target_mesh_names)

      ! Find the name of the mesh_type object that was extruded from
      ! the local mesh object that corresponds to the specified
      ! local mesh target.
      if (allocated(target_mesh_names)) then

        do j=1, size(target_mesh_names)
          mesh_name_B = cmdi
          do k=1, size(local_mesh_names)
            if (local_mesh_names(k) == target_mesh_names(j)) then
              mesh_name_B = mesh_names(k)
              exit
            end if
          end do

          if (mesh_name_B /= cmdi) then
            call add_mesh_map( mesh_name_A, mesh_name_B )
          end if
        end do

        deallocate(target_mesh_names)
      end if
    end do

  end if

end subroutine assign_mesh_maps


!> @brief       Creates integrid map between two mesh_type objects.
!> @description The meshes should be contain valid local mesh integrid maps.
!> @param[in] source_mesh_name  Name of source_mesh in the
!!                              application mesh collection.
!> @param[in] target_mesh_name  Name of source_mesh in the
!!                              application mesh collection
subroutine add_mesh_map( source_mesh_name, &
                         target_mesh_name )

  implicit none

  character(len=str_def), intent(in) :: source_mesh_name
  character(len=str_def), intent(in) :: target_mesh_name

  type(mesh_type), pointer :: source_mesh => null()
  type(mesh_type), pointer :: target_mesh => null()


  source_mesh => mesh_collection % get_mesh( source_mesh_name )
  target_mesh => mesh_collection % get_mesh( target_mesh_name )

  if ( associated(source_mesh) .and. &
       associated(target_mesh) ) then

    ! Mesh tag names may be different but "could point to the same mesh
    ! So check the IDs are not the same
    if (source_mesh%get_id() == target_mesh%get_id()) then
      write(log_scratch_space,'(A)')                  &
          'Unable to create intergrid map: Source('// &
          trim(source_mesh_name)//' and target('//    &
          trim(target_mesh_name)//') mesh IDs are the same'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    call source_mesh % add_mesh_map (target_mesh)
    write(log_scratch_space,'(A,I0,A)')     &
        'Adding intergrid map "'//          &
         trim(source_mesh_name)//'"-->"'//  &
         trim(target_mesh_name)//'"'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  else
    write(log_scratch_space,'(A,I0,A)')          &
        'Unable to create mesh map between "'//  &
        trim(source_mesh_name)//'"-"'//          &
        trim(target_mesh_name)//'"'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  nullify(source_mesh)
  nullify(target_mesh)

  return
end subroutine add_mesh_map

end module add_mesh_map_mod
