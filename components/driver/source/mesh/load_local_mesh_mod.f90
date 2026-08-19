!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!> @brief Load local mesh object data from file.
module load_local_mesh_mod

  use constants_mod,       only: i_def, str_def, &
                                 str_max_filename
  use local_mesh_mod,      only: local_mesh_type
  use log_mod,             only: log_event,         &
                                 log_scratch_space, &
                                 LOG_LEVEL_INFO
  use ugrid_mesh_data_mod, only: ugrid_mesh_data_type
  use sci_query_mod,       only: is_lbc

  use local_mesh_collection_mod, only: local_mesh_collection

  implicit none

  private
  public :: load_local_mesh

  interface load_local_mesh
    module procedure load_local_mesh_single
    module procedure load_local_mesh_multiple
  end interface load_local_mesh

contains


!> @brief Loads multiple local mesh object data from a UGRID file.
!>        and adds them to the application local mesh collection.
!> @param[in] input_mesh_file  UGRID file containing data to
!>                             populate <local_mesh_type> objects.
!> @param[in] mesh_names       The names of the local mesh data to
!>                             load from the <input_mesh_file>.
subroutine load_local_mesh_multiple( input_mesh_file, &
                                     mesh_names )
  implicit none

  character(str_max_filename), intent(in) :: input_mesh_file
  character(str_def),          intent(in) :: mesh_names(:)

  integer(i_def) :: i

  do i=1, size(mesh_names)
    call load_local_mesh_single( input_mesh_file, mesh_names(i) )
  end do

end subroutine load_local_mesh_multiple


!> @brief Loads a single local mesh object data from a UGRID file.
!>        and adds them to the application local mesh collection.
!> @param[in] input_mesh_file  UGRID file containing data to
!>                             populate <local_mesh_type> objects.
!> @param[in] mesh_name        The name of the local mesh data to
!>                             load from the <input_mesh_file>.
subroutine load_local_mesh_single( input_mesh_file, &
                                   mesh_name )

  implicit none

  character(str_max_filename), intent(in) :: input_mesh_file
  character(str_def),          intent(in) :: mesh_name

  type(ugrid_mesh_data_type) :: ugrid_mesh_data
  type(local_mesh_type)      :: local_mesh

  integer(i_def) :: local_mesh_id

  if (.not. local_mesh_collection%check_for(mesh_name)) then

    ! Load mesh data into local_mesh
    call ugrid_mesh_data%read_from_file( input_mesh_file, mesh_name )

    if (ugrid_mesh_data%contains_mesh()) then

      write(log_scratch_space,'(A)') &
          'Local mesh: "'//trim(mesh_name)//'"'//' loaded.'
      call local_mesh%initialise_from_ugrid_data( ugrid_mesh_data )

      ! Assign cell ownership for lookup.
      call local_mesh%init_cell_owner()

      local_mesh_id = local_mesh_collection%add_new_local_mesh( local_mesh )

    else

      write(log_scratch_space,'(A)') &
          'Local mesh: "'//trim(mesh_name)//'"'//' NOT loaded.'

    end if ! ugrid data object contains a mesh

    call log_event(log_scratch_space, LOG_LEVEL_INFO)

  end if ! not already in collection

end subroutine load_local_mesh_single

end module load_local_mesh_mod
