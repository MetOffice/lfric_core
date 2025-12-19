!-------------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Setup of files to be read for the multifile IO demo
!> @details Handles the setup of the filelist to later be passed to the XIOS
!>          context describing what files to read, when and what will be read
!>          from them.
module multifile_file_setup_mod

  use constants_mod,         only: i_def, str_def, l_def
  use driver_modeldb_mod,    only: modeldb_type
  use field_collection_mod,  only: field_collection_type
  use file_mod,              only: file_type, FILE_MODE_READ
  use lfric_xios_file_mod,   only: lfric_xios_file_type, OPERATION_ONCE, &
                                   OPERATION_TIMESERIES
  use linked_list_mod,       only: linked_list_type

  implicit none

  private
  public :: init_multifile_files

contains

  !> @details Creates the fields needed for the multifile IO
  !> @param[in,out] modeldb The model database in which to store model data.
  subroutine init_multifile_files(files_list, modeldb, filename)

    implicit none

    type(linked_list_type), intent(out)   :: files_list
    type(modeldb_type),     intent(inout) :: modeldb
    character(str_def),     intent(in)    :: filename

    character(str_def) :: timestep_start, timestep_end

    integer(i_def) :: ts_start, ts_end
    integer(i_def) :: rc
    logical(l_def) :: use_xios_io

    type(field_collection_type), pointer :: multifile_fields

    use_xios_io    = modeldb%config%io%use_xios_io()
    timestep_start = modeldb%config%time%timestep_start()
    timestep_end   = modeldb%config%time%timestep_end()

    multifile_fields => modeldb%fields%get_field_collection("multifile_io_fields")

    ! Get time configuration in integer form
    read(timestep_start,*,iostat=rc) ts_start
    read(timestep_end,  *,iostat=rc) ts_end

    if (use_xios_io) then

      call files_list%insert_item( &
        lfric_xios_file_type( filename, &
        xios_id = "multifile_io_fields", &
        io_mode=FILE_MODE_READ, &
        freq=1, &
        operation=OPERATION_TIMESERIES, &
        fields_in_file=multifile_fields))

    end if ! use_xios_io

  end subroutine init_multifile_files

end module multifile_file_setup_mod
