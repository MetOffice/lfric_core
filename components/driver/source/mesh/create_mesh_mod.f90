!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!> @brief Functions/Routines related to creating a <mesh_object_type>
module create_mesh_mod

  use constants_mod, only: i_def, str_def, r_def, l_def, imdi, &
                           str_max_filename
  use log_mod,       only: log_event,         &
                           log_scratch_space, &
                           LOG_LEVEL_DEBUG,   &
                           LOG_LEVEL_ERROR,   &
                           LOG_LEVEL_INFO

  use extrusion_mod,       only: extrusion_type,           &
                                 uniform_extrusion_type,   &
                                 geometric_extrusion_type, &
                                 quadratic_extrusion_type, &
                                 PRIME_EXTRUSION,          &
                                 SHIFTED,                  &
                                 DOUBLE_LEVEL
  use local_mesh_mod,      only: local_mesh_type
  use mesh_mod,            only: mesh_type
  use sci_query_mod,       only: check_lbc
  use ugrid_mesh_data_mod, only: ugrid_mesh_data_type

  use local_mesh_collection_mod,  only: local_mesh_collection
  use mesh_collection_mod,        only: mesh_collection

  ! Configuration modules
  use extrusion_config_mod, only: method_uniform,   &
                                  method_geometric, &
                                  method_quadratic

  use multigrid_config_mod, only: chain_mesh_tags

  use partitioning_config_mod, only: tile_size_x,               &
                                     tile_size_y,               &
                                     inner_halo_tiles,          &
                                     max_tiled_multigrid_level, &
                                     coarsen_multigrid_tiles

  implicit none

  private
  public :: create_extrusion, create_mesh

  interface create_mesh
    module procedure create_mesh_single
    module procedure create_mesh_multiple
  end interface create_mesh

contains


!> @brief  Creates vertical mesh extrusion.
!> @return  Resulting extrusion object
function create_extrusion( extrusion_method, &
                           domain_height,       &
                           domain_bottom,    &
                           n_layers,         &
                           extrusion_id ) result(new)

  implicit none

  class(extrusion_type), allocatable :: new

  integer(i_def), intent(in) :: extrusion_method
  real(r_def),    intent(in) :: domain_height
  real(r_def),    intent(in) :: domain_bottom
  integer(i_def), intent(in) :: n_layers
  integer(i_def), intent(in) :: extrusion_id

  if (allocated(new)) deallocate(new)

  select case (extrusion_method)
    case (method_uniform)
      allocate( new, source=uniform_extrusion_type(        &
                                domain_bottom, domain_height, &
                                n_layers, extrusion_id ) )
    case (method_quadratic)
      allocate( new, source=quadratic_extrusion_type(      &
                                domain_bottom, domain_height, &
                                n_layers, extrusion_id ) )
    case (method_geometric)
      allocate( new, source=geometric_extrusion_type(      &
                                domain_bottom, domain_height, &
                                n_layers, extrusion_id ) )
    case default
      call log_event("Invalid method for simple extrusion", LOG_LEVEL_ERROR)
  end select

end function create_extrusion


!> @brief    Generates mesh object types.
!> @details  Creates multiple mesh objects from local_mesh_type objects
!!           held in the application local_mesh_collection and the
!!           specified extrusion.
!!
!> @param[in]  local_mesh_names  Names of the local_mesh_types to extrude.
!> @param[in]  extrusion         Extrusion to employ.
!> @param[in]  alt_name          Optional, Alternative names for the
!!                               extruded meshes, defaults to local_mesh_names
!!                               if absent.
subroutine create_mesh_multiple( local_mesh_names, &
                                 extrusion,        &
                                 alt_name )
  implicit none

  character(str_def),    intent(in) :: local_mesh_names(:)
  class(extrusion_type), intent(in) :: extrusion
  character(str_def),    intent(in), &
                         optional   :: alt_name(:)

  ! Local variables
  integer(i_def) :: i
  character(str_def), allocatable :: names(:)

  if (present(alt_name)) then

    if ( size(alt_name) /= size(local_mesh_names) ) then
      write(log_scratch_space, '(A)')                          &
          'Number of alternative mesh names does not match '// &
          'number of requested meshes.'
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

    allocate(names, source=alt_name)

  else

    allocate(names, source=local_mesh_names)

  end if

  do i=1, size(local_mesh_names)
    call create_mesh_single( local_mesh_names(i), &
                             extrusion,           &
                             alt_name=names(i) )
  end do

  deallocate(names)

end subroutine create_mesh_multiple

!> @brief    Generates a single mesh_type object.
!> @details  Instantiates a mesh_type object and adds it to the applications
!!           mesh collection. Multiple meshes may be generated in the model
!!           based on the same global mesh but with differing extrusions.
!!
!> @param[in]  local_mesh_name  Name of local_mesh_type object in
!!                              application local_mesh_collection.
!> @param[in]  extrusion        Extrusion to employ for this mesh_type object
!> @param[in]  alt_name         Optional, Alternative name for the
!!                              extruded mesh, defaults to local_mesh_name
!!                              if absent.
subroutine create_mesh_single( local_mesh_name, &
                               extrusion,       &
                               alt_name )

  implicit none

  character(str_def),    intent(in) :: local_mesh_name
  class(extrusion_type), intent(in) :: extrusion
  character(str_def),    intent(in), &
                         optional   :: alt_name

  type(local_mesh_type), pointer :: local_mesh_ptr => null()

  type(mesh_type)        :: mesh
  integer(kind=i_def)    :: mesh_id
  character(len=str_def) :: name

  integer(kind=i_def) :: tile_size(2)
  integer(kind=i_def) :: multigrid_level
  integer(kind=i_def) :: max_multigrid_level
  logical(kind=l_def) :: set_tile_size

  if ( .not. present(alt_name) ) then
    name = local_mesh_name
  else
    name = alt_name
  end if


  ! 1.0 Check if mesh_type already exists.
  !===============================================
  if ( mesh_collection%check_for(name) ) then
    write(log_scratch_space,'(A)')                          &
        'No action taken: Mesh '//trim(name)//' already '// &
        'exists in the program mesh_collection object.'
    call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    return
  end if


  ! 2.0 Extrude the local_mesh_object.
  !===============================================
  local_mesh_ptr => local_mesh_collection%get_local_mesh(local_mesh_name)

  if ( .not. associated(local_mesh_ptr) ) then
    if ( check_lbc(local_mesh_name) ) then
      return
    else
      write(log_scratch_space,'(A)')                                &
          'Specified local mesh object ('//trim(local_mesh_name)//  &
          ') was not found in the program local_mesh_collection '// &
          'object.'
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if
  end if


  ! 3.0 Set up tiling
  !===============================================
  ! Set coarsest multigrid level that will be tiled;
  ! restrict to the finest grid by default
  max_multigrid_level = 1
  if ( max_tiled_multigrid_level /= imdi ) then
    max_multigrid_level = max_tiled_multigrid_level
  end if

  ! The tiling module uses 1x1 tiles (equivalent to colouring) by
  ! default; allow user-specified tile sizes in case of 3D meshes
  ! (PRIME_EXTRUSION, SHIFTED, and DOUBLE_LEVEL extrusions) and up to
  ! the specified multigrid level (count levels until mesh name
  ! includes the chain mesh tag). This relies on mesh name conventions
  ! and a tag order from finest (level 1) to coarsest mesh (level n).
  set_tile_size = .false.
  if ( extrusion%get_id() == PRIME_EXTRUSION .or. &
       extrusion%get_id() == SHIFTED         .or. &
       extrusion%get_id() == DOUBLE_LEVEL ) then
    if ( allocated(chain_mesh_tags) ) then
      ! Multigrid setup - use tiling if multigrid level is allowed, and
      ! if mesh name includes the mesh tag at that level
      do multigrid_level = 1, SIZE(chain_mesh_tags)
        if ( index( trim(name), trim(chain_mesh_tags(multigrid_level)) ) > 0 &
             .and. multigrid_level <= max_multigrid_level ) then
          set_tile_size = .true.
          exit
        end if
      end do
    else
      ! Not a multigrid setup - use tiling
      set_tile_size = .true.
    end if
  end if

  ! Set user-specified tile size if tiling is allowed and adapt it to coarser
  ! multigrid levels if requested and applicable
  tile_size = 1
  if ( set_tile_size ) then
    if ( tile_size_x /= imdi ) tile_size(1) = tile_size_x
    if ( tile_size_y /= imdi ) tile_size(2) = tile_size_y
    if ( coarsen_multigrid_tiles .and. allocated( chain_mesh_tags ) ) then
      do multigrid_level = 1, SIZE(chain_mesh_tags)
        if ( index( trim(name), &
                    trim(chain_mesh_tags(multigrid_level)) ) > 0 ) exit
        tile_size = max( tile_size / 2, 1 )
      end do
    end if
  end if

  mesh = mesh_type( local_mesh_ptr, extrusion, mesh_name=name, &
                    tile_size=tile_size, inner_halo_tiles=inner_halo_tiles )

  mesh_id = mesh_collection%add_new_mesh( mesh )
  call mesh%clear()


  ! 4.0 Report on mesh_type creation.
  !===============================================
  write(log_scratch_space,'(A,I0,A)')                 &
      '   ... "'//trim(name)//'"(id:', mesh_id,') '// &
      'based on mesh "'//trim(local_mesh_name)//'"'

  if (mesh_id /= imdi) then
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
  else
    write(log_scratch_space,'(A,I0,A)') &
        trim(log_scratch_space)//' (FAILED)'
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

end subroutine create_mesh_single

end module create_mesh_mod
