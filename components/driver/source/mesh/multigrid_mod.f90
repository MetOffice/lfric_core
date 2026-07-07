!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module multigrid_mod

  use constants_mod,     only: i_def, l_def, str_def, imdi
  use extrusion_mod,     only: prime_extrusion, shifted, double_level, twod
  use fs_continuity_mod, only: w2, w3, wtheta, w2v, w2h
  use log_mod,           only: log_event, log_scratch_space, &
                               log_level_info, log_level_error

  ! Collections
  use mesh_collection_mod,           only: mesh_collection
  use function_space_collection_mod, only: function_space_collection

  ! Object types
  use config_mod,    only: config_type
  use extrusion_mod, only: extrusion_type
  use mesh_mod,      only: mesh_type

  use function_space_mod,       only: function_space_type
  use function_space_chain_mod, only: function_space_chain_type

  implicit none

  public :: init_multigrid_fs_chain, &
            get_multigrid_tile_size
  public :: single_layer_function_space_chain,     &
            multigrid_function_space_chain,        &
            w2_multigrid_function_space_chain,     &
            wtheta_multigrid_function_space_chain, &
            w2h_multigrid_function_space_chain,    &
            w2v_multigrid_function_space_chain

  !> @name Global variables
  !>
  !> @todo An alternative to global variables will be needed
  !> in order to support multi-instance models.
  !>
  !> @{
  type(function_space_chain_type), allocatable :: &
           single_layer_function_space_chain,     &
           multigrid_function_space_chain,        &
           w2_multigrid_function_space_chain,     &
           wtheta_multigrid_function_space_chain, &
           w2h_multigrid_function_space_chain,    &
           w2v_multigrid_function_space_chain
  !> @}

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> @brief Routine returns tile sizes for supplied mesh names/extrusion where
!>        applicable to the multigrid configuration.
!>
!> @param[in] config           Application namelist configuration object
!> @param[in] local_mesh_names Meshes to set multigrid tile sizes
!> @param[in] extrusion        Extrusion object being applied to meshes.
!>
!> @return tile_size  Updated tile sizes for multigrid, if applicable.
!>                    Missing data indicator is returned for where the local
!>                    mesh tile size is not to be updated for multigrid.
!
function get_multigrid_tile_size( config, local_mesh_names, extrusion) &
                         result ( tile_size )

  implicit none

  type(config_type),     intent(in) :: config
  character(str_def),    intent(in) :: local_mesh_names(:)
  class(extrusion_type), intent(in) :: extrusion

  integer(i_def), allocatable :: tile_size(:,:)

  integer(i_def) :: multigrid_level
  integer(i_def) :: max_multigrid_level
  logical(l_def) :: coarsen_multigrid_tiles
  logical(l_def) :: set_tile_size

  character(str_def), allocatable :: chain_mesh_tags(:)

  integer(i_def)     :: extrusion_id, i
  character(str_def) :: name

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

  !=========================================================================
  if (allocated(tile_size)) deallocate(tile_size)
  allocate(tile_size(2,(size(local_mesh_names))))
  tile_size = imdi

  if (coarsen_multigrid_tiles) then

    extrusion_id = extrusion%get_id()
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

end function get_multigrid_tile_size

!> @brief  Initialises the function space chains used in multigrid.
!> @param[in] multigrid_mesh_names  Names of the multigrid meshes
subroutine init_multigrid_fs_chain(multigrid_mesh_names)

    implicit none

    character(str_def), intent(in) :: multigrid_mesh_names(:)

    type(mesh_type), pointer :: mesh
    type(mesh_type), pointer :: twod_mesh

    type(function_space_type), pointer :: fs

    integer(i_def) :: i

    nullify(mesh, twod_mesh, fs)

    call log_event( 'FEM specifics: creating function space chains...', &
                    log_level_info )

    ! ======================================================================== !
    ! Create function space chains
    ! ======================================================================== !

    multigrid_function_space_chain        = function_space_chain_type()
    w2_multigrid_function_space_chain     = function_space_chain_type()
    w2v_multigrid_function_space_chain    = function_space_chain_type()
    w2h_multigrid_function_space_chain    = function_space_chain_type()
    wtheta_multigrid_function_space_chain = function_space_chain_type()

    write(log_scratch_space,'(A,I1,A)') &
        'Initialising MultiGrid ', size(multigrid_mesh_names), &
        '-level function space chain.'
    call log_event( log_scratch_space, log_level_info )

    do i = 1, size(multigrid_mesh_names)

      mesh => mesh_collection%get_mesh( multigrid_mesh_names(i) )

      ! Make sure this function_space is in the collection
      fs => function_space_collection%get_fs( mesh, 0, 0, w3 )
      call multigrid_function_space_chain%add( fs )

      fs => function_space_collection%get_fs( mesh, 0, 0, w2 )
      call w2_multigrid_function_space_chain%add( fs )

      fs => function_space_collection%get_fs( mesh, 0, 0, w2v )
      call w2v_multigrid_function_space_chain%add( fs )

      fs => function_space_collection%get_fs( mesh, 0, 0, w2h )
      call w2h_multigrid_function_space_chain%add( fs )

      fs => function_space_collection%get_fs( mesh, 0, 0, wtheta )
      call wtheta_multigrid_function_space_chain%add( fs )
    end do

    single_layer_function_space_chain = function_space_chain_type()
    do i = 1, size(multigrid_mesh_names)
      mesh => mesh_collection%get_mesh( multigrid_mesh_names(i) )
      twod_mesh => mesh_collection%get_mesh( mesh, twod )
      fs => function_space_collection%get_fs( twod_mesh, 0, 0, w3 )
      call single_layer_function_space_chain%add( fs )
    end do

    call log_event( 'Function space chains created', log_level_info )

  end subroutine init_multigrid_fs_chain

end module multigrid_mod
