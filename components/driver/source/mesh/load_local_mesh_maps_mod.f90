!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!> @brief Load local mesh maps from file.
module load_local_mesh_maps_mod

  use constants_mod,             only: i_def, str_def, str_max_filename
  use local_mesh_mod,            only: local_mesh_type
  use log_mod,                   only: log_event,         &
                                       log_scratch_space, &
                                       LOG_LEVEL_ERROR
  use ncdf_quad_mod,             only: ncdf_quad_type
    use ugrid_2d_mod,   only: ugrid_2d_type


  use local_mesh_collection_mod, only: local_mesh_collection

  implicit none

  private
  public :: load_local_mesh_maps

  interface load_local_mesh_maps
    module procedure load_local_mesh_maps_single_source
    module procedure load_local_mesh_maps_multiple_source
  end interface load_local_mesh_maps

contains

!> @brief    Loads local intergrid cell maps from a UGRID file.
!> @details  For given meshes, (that have been previosly loaded), the
!>           the intergrid cell maps for target local mesh(es) are
!>           loaded from the <input_mesh_file> and assigned to the
!>           specfied local meshes.
!>
!> @param[in]  input_mesh_file    UGRID file containing data to
!>                                populate <local_mesh_type> object.
!> @param[in]  source_mesh_name   The name of the local source mesh to
!>                                load maps from the <input_mesh_file>.
subroutine load_local_mesh_maps_multiple_source( input_mesh_file, &
                                                 source_mesh_names, ugrid_2d )

  implicit none

  character(str_max_filename), intent(in) :: input_mesh_file
  character(str_def),          intent(in) :: source_mesh_names(:)
  type(ugrid_2d_type), intent(inout) :: ugrid_2d

  integer(i_def) :: i

  do i=1, size(source_mesh_names)
    call load_local_mesh_maps_single_source( input_mesh_file, &
                                             source_mesh_names(i), ugrid_2d )
  end do

end subroutine load_local_mesh_maps_multiple_source


!> @brief    Loads local intergrid cell maps from a UGRID file.
!> @details  For a given mesh, (that has been previosly loaded), the
!>           the intergrid cell maps for target local mesh(es) are
!>           loaded the <input_mesh_file> and assigned to the specfied
!>           local mesh.
!>
!> @param[in]  input_mesh_file    UGRID file containing data to
!>                                populate <local_mesh_type> object.
!> @param[in]  source_mesh_name   The name of the local source mesh to load
!>                                maps from the <input_mesh_file>.
subroutine load_local_mesh_maps_single_source( input_mesh_file, &
                                               source_mesh_name, ugrid_2d )

  implicit none

  character(str_max_filename), intent(in) :: input_mesh_file
  character(str_def),          intent(in) :: source_mesh_name

  character(str_def), allocatable :: target_mesh_names(:)
  integer(i_def),     allocatable :: lid_mesh_map(:,:,:)

  integer(i_def) :: i

  type(local_mesh_type), pointer :: source_mesh => null()
  type(local_mesh_type), pointer :: target_mesh => null()
  type(ugrid_2d_type), intent(inout) :: ugrid_2d

  integer(i_def) :: target_mesh_id


  ! Read in the maps for each local mesh
  !=================================================================
  source_mesh => local_mesh_collection%get_local_mesh( source_mesh_name )
  if (.not. associated(source_mesh)) then
    write(log_scratch_space,'(A)') ' Mesh "'//trim(source_mesh_name)// &
                                   '" not found in collection'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  call source_mesh%get_target_mesh_names(target_mesh_names)

  if (allocated(target_mesh_names)) then

    do i=1, size(target_mesh_names)
      if ( local_mesh_collection%check_for( target_mesh_names(i) ) ) then

        ! Read in the local mesh map.
        call ugrid_2d%file_handler%read_map( source_mesh_name,     &
                                    target_mesh_names(i), &
                                    lid_mesh_map )

        target_mesh &
            => local_mesh_collection%get_local_mesh(target_mesh_names(i))
        target_mesh_id =  target_mesh%get_id()

        ! Assign local mesh map to the local mesh object.
        call source_mesh%add_local_mesh_map( target_mesh_id, &
                                             lid_mesh_map )

        if (allocated( lid_mesh_map )) deallocate( lid_mesh_map )

      end if
    end do

    deallocate( target_mesh_names )

  end if

end subroutine load_local_mesh_maps_single_source

end module load_local_mesh_maps_mod
