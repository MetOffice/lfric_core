!==============================================================================
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!==============================================================================
module calc_global_cell_map_mod

  use ugrid_generator_mod, only: ugrid_generator_type
  use constants_mod,       only: i_def, l_def
  use log_mod,             only: log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  private

  public :: calc_global_cell_map

contains

!-------------------------------------------------------------------------------
!>  @brief   Returns a mapping from this mesh (source) to a target mesh
!>  @details
!>
!>  @param[in]  self    The gencube_ps_type instance reference.
!>  @param[in]  target_edge_cells_x
!>  @param[in]  target_edge_cells_y
!>
!>  @param[out] cell_map The map from source to target.
!-------------------------------------------------------------------------------
subroutine calc_global_cell_map( source_mesh,         &
                                 target_edge_cells_x, &
                                 target_edge_cells_y, &
                                 cell_map,            &
                                 panel_rotations      )

  use, intrinsic :: iso_fortran_env, only : stdout => output_unit

  implicit none

  class(ugrid_generator_type), intent(in) :: source_mesh

  integer(i_def), allocatable, intent(out) :: cell_map(:,:,:)
  integer(i_def), optional,    intent(in)  :: panel_rotations(:)

  integer(i_def) :: npanels

  integer(i_def) :: fine_cpp
  integer(i_def) :: fine_ncells
  integer(i_def) :: fine_start_id
  integer(i_def) :: fine_end_id
  integer(i_def) :: fine_x
  integer(i_def) :: fine_y
  integer(i_def), allocatable :: fine_ids(:,:)
  integer(i_def), allocatable :: fine_to_coarse_gid_map(:,:,:)

  integer(i_def) :: coarse_cpp
  integer(i_def) :: coarse_ncells
  integer(i_def) :: coarse_start_id
  integer(i_def) :: coarse_end_id
  integer(i_def) :: coarse_x
  integer(i_def) :: coarse_y
  integer(i_def), allocatable :: coarse_ids(:,:)
  integer(i_def), allocatable :: coarse_to_fine_gid_map(:,:,:)

  integer(i_def) :: source_edge_cells_x
  integer(i_def) :: source_edge_cells_y
  integer(i_def) :: target_edge_cells_x
  integer(i_def) :: target_edge_cells_y

  integer(i_def) :: start_x, end_x
  integer(i_def) :: start_y, end_y
  integer(i_def) :: mod_x, mod_y
  integer(i_def) :: edge_factor_x, edge_factor_y

  logical(l_def) :: refining
  logical(l_def) :: refine_x
  logical(l_def) :: refine_y
  logical(l_def) :: coarsen_x
  logical(l_def) :: coarsen_y

  integer(i_def), allocatable :: tmp_panel_ids(:)
  integer(i_def), allocatable :: tmp_map(:,:)
  integer(i_def) :: n,i,j,k, count
  integer(i_def) :: ifine, jfine



  ! Get some information from the source mesh
  call source_mesh%get_metadata( npanels=npanels,                  &
                                 edge_cells_x=source_edge_cells_x, &
                                 edge_cells_y=source_edge_cells_y  )

  ! Set logicals to describe refine/coarsen operations
  refine_x  = ( source_edge_cells_x < target_edge_cells_x )
  refine_y  = ( source_edge_cells_y < target_edge_cells_y )
  coarsen_x = ( source_edge_cells_x > target_edge_cells_x )
  coarsen_y = ( source_edge_cells_y > target_edge_cells_y )

  ! Perform some checks and set refining logical
  if ( (.not. refine_x) .and. (.not. coarsen_x) .and. &
       (.not. refine_y) .and. (.not. coarsen_y) ) then
    write(log_scratch_space,'(A)') &
        'Attempting to map mesh to itself .... pointless'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  if ( (refine_x  .and. coarsen_y) .or. &
       (coarsen_x .and. refine_y) ) then
    write(log_scratch_space,'(A)')                              &
        'Mixed refine/coarsen operations for differing axis '// &
        'is not supported for single mappings'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  else if (refine_x .or. refine_y) then
    refining = .true.
  else if (coarsen_x .or. coarsen_y) then
    refining = .false.
  end if


  ! Now we know we are only performing coarsening or refining operations

  ! Now check that number of target mesh cells are a factor
  ! of the source mesh cells if we are coarsening
  mod_x = mod(source_edge_cells_x, target_edge_cells_x)
  mod_y = mod(source_edge_cells_y, target_edge_cells_y)

  if ((.not. refining) .and. (mod_x /= 0 .or. mod_y /= 0)) then
    write(log_scratch_space,'(A)')                        &
        'When coarsening, number of target axis cells '// &
        'must be factor source axis cells'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  ! Now we cat set which of the inputs meshes is the coarse
  ! and fine mesh
  if (refining) then
    coarse_x   = source_edge_cells_x
    coarse_y   = source_edge_cells_y
    fine_x     = target_edge_cells_x
    fine_y     = target_edge_cells_y
  else
    coarse_x   = target_edge_cells_x
    coarse_y   = target_edge_cells_y
    fine_x     = source_edge_cells_x
    fine_y     = source_edge_cells_y
  end if

  ! Set up values and arrays to do mapping
  fine_cpp   = fine_x * fine_y
  coarse_cpp = coarse_x * coarse_y

  edge_factor_x = max(1,fine_x / coarse_x)
  edge_factor_y = max(1,fine_y / coarse_y)

  fine_ncells   = nPanels * fine_cpp
  coarse_ncells = nPanels * coarse_cpp

  allocate( coarse_ids ( coarse_x, coarse_y ))
  allocate( fine_ids   ( fine_x,   fine_y   ))

  allocate( coarse_to_fine_gid_map( edge_factor_x, edge_factor_y, coarse_ncells ) )

  if (.not. refining) then
    allocate( fine_to_coarse_gid_map( 1, 1, fine_ncells ) )
  end if



  ! Now loop over each panel numbering coarse / fine cell ids
  ! and then populating the cell maps
  do n=1, nPanels

    ! =======================================================
    ! Get the id of the start/end cells for the panel
    ! =======================================================
    coarse_start_id = ((n-1)*coarse_cpp) + 1
    coarse_end_id   = coarse_start_id + coarse_cpp-1
    fine_start_id   = ((n-1)*fine_cpp) + 1
    fine_end_id     = fine_start_id + fine_cpp-1

    ! =======================================================
    ! Populate arrays with ids and reshape to panel layout
    ! =======================================================

    ! For coarse mesh
    allocate(tmp_panel_ids(coarse_cpp))
    count=1
    do i=coarse_start_id, coarse_end_id
      tmp_panel_ids(count) = i
      count=count+1
    end do

    coarse_ids = reshape(tmp_panel_ids,(/coarse_x, coarse_y/))
    deallocate(tmp_panel_ids)

    ! For fine mesh
    allocate(tmp_panel_ids(fine_cpp))
    count=1
    do i=fine_start_id, fine_end_id
      tmp_panel_ids(count) = i
      count=count+1
    end do
    fine_ids = reshape(tmp_panel_ids,(/fine_x, fine_y/))
    deallocate(tmp_panel_ids)

    ! ===============================
    ! Populate Global Id maps
    ! ===============================

    allocate(tmp_map(edge_factor_x, edge_factor_y))

    ! Coarse to Fine
    do j=1, coarse_y
      start_y = ((j-1)*edge_factor_y) + 1
      end_y   = start_y + edge_factor_y - 1
      do i=1, coarse_x
        start_x = ((i-1)*edge_factor_x) + 1
        end_x   = start_x + edge_factor_x - 1

        ! Rotation of map for LFRic cubed sphere panels
        ! The global_mesh_map is a private object so must happen here
        tmp_map(:,:) = reshape( fine_ids( start_x:end_x, start_y:end_y ),  &
                                          (/edge_factor_x, edge_factor_y/) )

        if ( ( edge_factor_x /= edge_factor_y ) .or. &
             ( .not. present(panel_rotations) ) ) then

          coarse_to_fine_gid_map(:, :,coarse_ids(i,j)) = tmp_map(:,:)

        else if ( panel_rotations(n) == 0_i_def ) then
          ! No rotation
          coarse_to_fine_gid_map(:, :,coarse_ids(i,j)) = tmp_map(:,:)

        else if ( panel_rotations(n) == 1_i_def ) then
          ! Rotate left
          do jfine = 1, edge_factor_y
            do ifine = 1, edge_factor_x
              coarse_to_fine_gid_map(ifine, jfine, coarse_ids(i,j)) = &
                tmp_map(jfine, edge_factor_x+1-ifine)
            end do
          end do

        else if ( panel_rotations(n) == -1_i_def ) then
          ! Rotate right
          do jfine = 1, edge_factor_y
            do ifine = 1, edge_factor_x
              coarse_to_fine_gid_map(ifine, jfine, coarse_ids(i,j)) = &
                tmp_map(edge_factor_y+1-jfine, ifine)
            end do
          end do
        else
          call log_event('Unexpected panel rotation', LOG_LEVEL_ERROR)
        end if
      end do
    end do

    deallocate(tmp_map)

    if (.not. refining) then
      ! Fine to Coarse
      do k=coarse_start_id, coarse_end_id
        do j=1, edge_factor_y
          do i=1, edge_factor_x
            fine_to_coarse_gid_map(1, 1, coarse_to_fine_gid_map(i,j,k)) = k
          end do
        end do
      end do
    end if

  end do

  if (refining) then
    cell_map = coarse_to_fine_gid_map(:,:,:)
  else
    cell_map = fine_to_coarse_gid_map(:,:,:)
  end if

  if ( allocated( coarse_ids ) ) deallocate( coarse_ids )
  if ( allocated( fine_ids ) )   deallocate( fine_ids )
  if ( allocated( coarse_to_fine_gid_map ) ) deallocate( coarse_to_fine_gid_map )
  if ( allocated( fine_to_coarse_gid_map ) ) deallocate( fine_to_coarse_gid_map )

  return
end subroutine calc_global_cell_map

end module calc_global_cell_map_mod
