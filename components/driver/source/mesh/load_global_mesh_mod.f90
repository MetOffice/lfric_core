!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
!> @brief Load global mesh object data from file.
module load_global_mesh_mod

  use constants_mod,       only: i_def, str_def, &
                                 str_max_filename
  use global_mesh_mod,     only: global_mesh_type
  use log_mod,             only: log_event,         &
                                 log_scratch_space, &
                                 LOG_LEVEL_INFO
  use ugrid_mesh_data_mod, only: ugrid_mesh_data_type

  use global_mesh_collection_mod, only: global_mesh_collection_type, &
                                        global_mesh_collection

  implicit none

  private
  public :: load_global_mesh

  interface load_global_mesh
    module procedure load_global_mesh_single
    module procedure load_global_mesh_multiple
  end interface load_global_mesh

contains


!> @brief Loads multiple global mesh objetc data from a UGRID file.
!>        and adds them to the global_mesh_collection object.
!> @param[in] input_mesh_file  UGRID file containing data to
!>                             populate <global_mesh_type> object.
!> @param[in] mesh_names       The names of the global meshes to load
!>                             from the <input_mesh_file>.
subroutine load_global_mesh_multiple( input_mesh_file, &
                                      mesh_names,      &
                                      rename_to )


  implicit none

  character(str_max_filename), intent(in) :: input_mesh_file
  character(str_def),          intent(in) :: mesh_names(:)

  character(str_def), optional :: rename_to(:)

  character(str_def), allocatable :: names(:)

  integer(i_def) :: i


  allocate(names, source=mesh_names)

  if ( present(rename_to) ) then
    if (size(rename_to) == size(mesh_names)) then
      deallocate(names)
      allocate(names, source=rename_to)
    end if
  end if

  do i=1, size(mesh_names)
    call load_global_mesh_single( input_mesh_file, &
                                  mesh_names(i), rename_to=names(i) )
  end do

end subroutine load_global_mesh_multiple


!> @brief Loads a single global mesh topology data from a UGRID file.
!>        and adds it to the global_mesh_collection object.
!> @param[in] input_mesh_file  UGRID file containing data to
!>                             populate <global_mesh_type> object.
!> @param[in] mesh_name        The name of the global mesh to load
!>                             from the <input_mesh_file>.
subroutine load_global_mesh_single( input_mesh_file, &
                                    mesh_name, rename_to )

  implicit none

  character(*), intent(in) :: input_mesh_file
  character(*), intent(in) :: mesh_name

  character(*), optional, intent(in) :: rename_to


  type(ugrid_mesh_data_type) :: ugrid_mesh_data
  type(global_mesh_type)     :: global_mesh

  character(str_def) :: name


  if ( present(rename_to) ) then
    name=rename_to
  else
    name=mesh_name
  end if

  if (.not. global_mesh_collection%check_for(name)) then

    write(log_scratch_space,'(A)') &
        'Reading global mesh: "'//trim(mesh_name)//'"'
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    ! Load mesh data into global_mesh
    call ugrid_mesh_data%read_from_file( trim(input_mesh_file), &
                                         mesh_name )

    global_mesh = global_mesh_type( ugrid_mesh_data, rename_to=name )
    call ugrid_mesh_data%clear()

    call global_mesh_collection%add_new_global_mesh ( global_mesh )

  end if

end subroutine load_global_mesh_single

end module load_global_mesh_mod
