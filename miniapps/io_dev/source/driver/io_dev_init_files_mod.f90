!-------------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Sets up file configuration fof IO_Dev
!> @details Collects file configuration information and formats it so that it
!>          can be passed to the infrastructure
module io_dev_init_files_mod

  use constants_mod,         only: i_def, &
                                   str_def, str_max_filename
  use file_mod,              only: file_type, FILE_MODE_READ, &
                                   FILE_MODE_WRITE
  use lfric_xios_file_mod,   only: lfric_xios_file_type, OPERATION_ONCE, &
                                   OPERATION_TIMESERIES
  use linked_list_mod,       only: linked_list_type
  use log_mod,               only: log_event, log_level_error
  use io_dev_data_mod,       only: io_dev_data_type
  use driver_model_data_mod, only: model_data_type

  ! Configuration modules
  use files_config_mod,    only: diag_stem_name,            &
                                 checkpoint_stem_name,      &
                                 start_dump_filename,       &
                                 start_dump_directory,      &
                                 time_varying_input_path,   &
                                 time_data_path
  use io_dev_config_mod,   only: field_initialisation,            &
                                 field_initialisation_start_dump, &
                                 time_variation,                  &
                                 time_variation_ancil
  use io_config_mod,       only: use_xios_io,               &
                                 diagnostic_frequency,      &
                                 checkpoint_write,          &
                                 checkpoint_read,           &
                                 write_diag, write_dump
  use time_config_mod,     only: timestep_start,            &
                                 timestep_end

  implicit none

  private
  public :: init_io_dev_files

  contains

  subroutine init_io_dev_files(files_list, model_data)

    implicit none

    type(linked_list_type),                   intent(out) :: files_list
    class(model_data_type), optional, target, intent(in)  :: model_data

    type(io_dev_data_type), pointer  :: io_dev_data => null()
    character(len=str_max_filename)  :: checkpoint_write_fname, &
                                        checkpoint_read_fname,  &
                                        dump_fname,             &
                                        input_fname
    integer(i_def)                   :: ts_start, ts_end
    integer(i_def)                   :: rc

    ! Get time configuration in integer form
    read(timestep_start,*,iostat=rc)  ts_start
    read(timestep_end,*,iostat=rc)  ts_end

    ! Make sure model data is of the correct type
    if (present(model_data)) then
      select type(model_data)
      type is (io_dev_data_type)
        io_dev_data => model_data
      end select
    else
      call log_event( "Model data required for IO_Dev file initialisation", &
                      log_level_error )
    end if

    if ( use_xios_io) then

      call files_list%insert_item( &
          lfric_xios_file_type( "io_dev_initial",         &
                                xios_id="io_dev_initial", &
                                io_mode=FILE_MODE_WRITE,  &
                                freq=1,                   &
                                operation=OPERATION_ONCE, &
                                fields_in_file=io_dev_data%dump_fields ) )

      ! Setup temporary diagnostic file information - we need this information
      ! to be able to process the files: in the future we would wish to obtain
      ! the diagnostic file information via XIOS API
      if ( write_diag) then
        call files_list%insert_item( &
          lfric_xios_file_type( "io_dev_diag",             &
                                xios_id="io_dev_diag",     &
                                io_mode=FILE_MODE_WRITE,   &
                                freq=diagnostic_frequency, &
                                operation=OPERATION_TIMESERIES ) )
      end if

      ! Setup dump-writing context information
      if ( write_dump ) then
        ! Create dump filename from base name and end timestep
        write(dump_fname,'(A,A,I6.6)') &
          trim(start_dump_directory)//'/'//trim(start_dump_filename),"_",ts_end

        call files_list%insert_item( &
          lfric_xios_file_type( dump_fname,                           &
                                xios_id="io_dev_dump_out",            &
                                io_mode=FILE_MODE_WRITE, freq=ts_end, &
                                operation=OPERATION_ONCE,             &
                                fields_in_file=io_dev_data%dump_fields ) )
      end if

      ! Setup dump-reading context information
      if ( field_initialisation == field_initialisation_start_dump ) then
        ! Create dump filename from stem
        write(dump_fname,'(A)') trim(start_dump_directory)//'/'// &
                                trim(start_dump_filename)

        call files_list%insert_item( &
          lfric_xios_file_type( dump_fname,                     &
                                xios_id="io_dev_dump_in",       &
                                io_mode=FILE_MODE_READ, freq=1, &
                                operation=OPERATION_ONCE,       &
                                fields_in_file=io_dev_data%dump_fields ) )
      end if

      ! Setup checkpoint writing context information
      if ( checkpoint_write ) then
        ! Create checkpoint filename from stem and end timestep
        write(checkpoint_write_fname,'(A,A,I6.6)') &
                             trim(checkpoint_stem_name),"_", ts_end

        call files_list%insert_item( &
          lfric_xios_file_type( checkpoint_write_fname,            &
                                xios_id="io_dev_checkpoint_write", &
                                io_mode=FILE_MODE_WRITE,           &
                                freq=ts_end - (ts_start - 1),      &
                                operation=OPERATION_ONCE,          &
                                fields_in_file=io_dev_data%core_fields ) )
      end if

      ! Setup checkpoint reading context information
      if ( checkpoint_read ) then
        ! Create checkpoint filename from stem and (start - 1) timestep
        write(checkpoint_read_fname,'(A,A,I6.6)') &
                     trim(checkpoint_stem_name),"_", (ts_start - 1)

        call files_list%insert_item( &
          lfric_xios_file_type( checkpoint_read_fname,            &
                                xios_id="io_dev_checkpoint_read", &
                                io_mode=FILE_MODE_READ, freq=1,   &
                                operation=OPERATION_ONCE,         &
                                fields_in_file=io_dev_data%core_fields ) )
      end if

      ! Setup time-varying input files
      if ( time_variation == time_variation_ancil ) then
        ! Set time-varying input filename from namelist
        write(input_fname,'(A)') trim(start_dump_directory)//'/'// &
                                 trim(time_varying_input_path)

        call files_list%insert_item( &
          lfric_xios_file_type( input_fname,                         &
                                xios_id="io_dev_time_varying_input", &
                                io_mode=FILE_MODE_READ ) )

        ! Set time data input filename from namelist
        write(input_fname,'(A)') trim(start_dump_directory)//'/'// &
                                 trim(time_data_path)

        call files_list%insert_item( &
          lfric_xios_file_type( input_fname,                    &
                                xios_id="io_dev_times",         &
                                io_mode=FILE_MODE_READ, freq=1, &
                                operation=OPERATION_ONCE ) )

      end if

    end if ! use_xios_io

  end subroutine init_io_dev_files

end module io_dev_init_files_mod
