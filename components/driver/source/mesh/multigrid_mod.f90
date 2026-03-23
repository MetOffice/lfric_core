!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module multigrid_mod

  use extrusion_mod, only: extrusion_type, prime_extrusion, &
                           shifted, double_level
  use config_mod,    only: config_type
  use constants_mod, only: i_def, l_def, str_def

  implicit none

  public :: get_multigrid_tile_size


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> @brief
!>
!> @param[in] config
!> @param[in] local_mesh_name
!> @param[in] extrusion
!>
!> @return tile_size
!>
subroutine get_multigrid_tile_size( config, local_mesh_names, extrusion, &
                                    tile_size )

  implicit none

  type(config_type),    intent(in) :: config
  character(str_def),   intent(in) :: local_mesh_names(:)
  type(extrusion_type), intent(in) :: extrusion

  integer(i_def), intent(inout) :: tile_size(:,:)



  integer(i_def) :: multigrid_level
  integer(i_def) :: max_multigrid_level
  logical(l_def) :: coarsen_multigrid_tiles
  logical(l_def) :: set_tile_size

  character(str_def), allocatable :: chain_mesh_tags(:)

  !=========================================================================
  ! This whole section should probably be in gungho science. It allows the
  ! Gungho multigrid scheme to override the tile settings in the
  ! configuration. This should really be written in the gungho science,
  ! though the decision to call it should be made by the application, i.e.
  ! the application may wish to use it's own tileing settings.
  !
  ! In partitioning namelist, should be in multigrid
  ! max_tiled_multigrid_level = config%multigrid%max_tiled_multigrid_level()
  ! coarsen_multigrid_tiles   = config%multigrid%coarsen_multigrid_tiles()
  coarsen_multigrid_tiles = config%multigrid%coarsen_multigrid_tiles()
  max_multigrid_level     = config%multigrid%max_tiled_multigrid_level()
  chain_mesh_tags         = config%multigrid%chain_mesh_tags()

  extrusion_id = extrusion%get_id()

  !=========================================================================
  if (coarsen_multigrid_tiles) then

    select case (extrusion_id)
    case(prime_extrusion, shifted, double_level)

      ! Set coarsest multigrid level that will be tiled;
      ! restrict to the finest grid by default
      if (max_multigrid_level == imdi) then
        call log_event('no max multigrid level set', log_level_error)
      end if

      do i=1, size(local_mesh_names)
        set_tile_size = .false.
        name =local_mesh_names(i)

        ! Multigrid setup - use tiling if multigrid level is allowed, and
        ! if mesh name includes the mesh tag at that level
        do multigrid_level=1, size(chain_mesh_tags)
          if ( index( trim(name),                                  &
                      trim(chain_mesh_tags(multigrid_level)) ) > 0 &
               .and. multigrid_level <= max_multigrid_level ) then
            set_tile_size = .true.
            exit
          end if
        end do

        if (set_tile_size) then
          do multigrid_level=1, size(chain_mesh_tags)
            if ( index( trim(name), &
                        trim(chain_mesh_tags(multigrid_level)) ) > 0 ) then
            exit
            end if
            tile_size(:,i) = max( tile_size(:,i)/2, 1 )
          end do
        end if ! set_tile_size
      end do ! local_mesh_names

    case default
      return
    end select

  end if ! Coarsen multigrid_tiles

end subroutine get_multigrid_tile_size

end module multigrid_mod
