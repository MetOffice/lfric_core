!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
module check_local_mesh_mod

  use constants_mod,             only: i_def, str_def, &
                                       str_max_filename
  use local_mesh_collection_mod, only: local_mesh_collection
  use local_mesh_mod,            only: local_mesh_type
  use log_mod,                   only: log_event,         &
                                       log_scratch_space, &
                                       log_level_error
  use namelist_collection_mod,   only: namelist_collection_type
  use namelist_mod,              only: namelist_type
  use sci_query_mod,             only: is_lbc, check_lbc

  use base_mesh_config_mod, only: key_from_geometry,       &
                                  key_from_topology,       &
                                  geometry_spherical,      &
                                  geometry_planar,         &
                                  topology_fully_periodic, &
                                  topology_non_periodic
  implicit none

  private
  public :: check_local_mesh

contains

!> @brief Basic validation that local meshes are suitable
!!        for the specified configuration.
!> @param[in]  configuration Configuration object.
!> @param[in]  stencil_depth Stencil depth that local meshes
!>                           need to support.
!> @param[in]  mesh_names    Local meshes held in application
!!                           local mesh collection object.
subroutine check_local_mesh( configuration, &
                             stencil_depth, &
                             mesh_names )

  implicit none

  type(namelist_collection_type), intent(in) :: configuration
  integer(i_def),                 intent(in) :: stencil_depth
  character(str_def),             intent(in) :: mesh_names(:)

  integer(i_def) :: topology
  integer(i_def) :: geometry

  logical :: valid_geometry
  logical :: valid_topology

  type(local_mesh_type), pointer :: local_mesh    => null()
  type(namelist_type),   pointer :: base_mesh_nml => null()

  integer(i_def) :: i
  integer(i_def) :: max_stencil_depth

  base_mesh_nml => configuration%get_namelist('base_mesh')

  call base_mesh_nml%get_value( 'geometry', geometry )
  call base_mesh_nml%get_value( 'topology', topology )

  base_mesh_nml => null()

  do i=1, size(mesh_names)

    local_mesh => local_mesh_collection%get_local_mesh(mesh_names(i))

    if (associated(local_mesh)) then

      ! Check mesh has valid domain geometry
      !=====================================
      valid_geometry = .false.
      select case ( geometry )

      case ( geometry_spherical )
        if ( local_mesh%is_geometry_spherical() ) valid_geometry = .true.

      case ( geometry_planar )
        if ( local_mesh%is_geometry_planar() ) valid_geometry = .true.

      end select

      if ( .not. valid_geometry ) then
        write(log_scratch_space, '(A)')        &
            'Mesh (' // trim(mesh_names(i)) // &
            ') in file is not valid as a ' //  &
            trim(key_from_geometry(geometry)) // ' domain geometry'
        call log_event(log_scratch_space, log_level_error)
      end if

      ! Check mesh has valid domain toplogy
      !=====================================
      valid_topology = .false.
      select case ( topology )

      case ( topology_fully_periodic )
        if ( local_mesh%is_topology_periodic() ) valid_topology = .true.

      case ( topology_non_periodic )
        if ( local_mesh%is_topology_non_periodic() ) valid_topology = .true.

      end select

      if ( .not. valid_topology ) then
        write(log_scratch_space, '(A)')           &
            'Mesh (' // trim(mesh_names(i)) //    &
            ') in file does not have a valid ' // &
            trim(key_from_topology(topology)) // ' topology'
        call log_event(log_scratch_space, log_level_error)
      end if

      if ( .not. is_lbc(local_mesh) ) then

        ! Check if local meshes are able to
        ! support the requested stencil depth.
        !=====================================
        max_stencil_depth = local_mesh%get_max_stencil_depth()

        if ( max_stencil_depth < stencil_depth ) then
          write(log_scratch_space,'(2(A,I0),A)')                     &
              'Insufficient stencil depth, configuration requires ', &
              stencil_depth, '. Mesh "'//trim(mesh_names(i))//       &
              '" supports a maximum stencil depth of ',              &
              max_stencil_depth, '.'
          call log_event(log_scratch_space, log_level_error)
        end if

      end if ! is_lbc

    else

      if ( .not. check_lbc(mesh_names(i)) ) then
        write(log_scratch_space,'(2(A,I0),A)')               &
            'Local mesh "' // trim(mesh_names(i)) // '" ' // &
            'not present in local mesh collection.'
        call log_event(log_scratch_space, log_level_error)
      end if

    end if

  end do

  nullify(local_mesh)

end subroutine check_local_mesh

end module check_local_mesh_mod

