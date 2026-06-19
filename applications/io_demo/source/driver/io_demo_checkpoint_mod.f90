!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Sets up checkpointing for the IO_demo app
!> @details Creates an IO context for checkpointing and adds the appropriate
!!          files and fields to it
module io_demo_checkpoint_mod

  use constants_mod,          only: i_def, str_max_filename, str_def, r_second
  use driver_modeldb_mod,     only: modeldb_type
  use event_mod,              only: event_action
  use event_actor_mod,        only: event_actor_type
  use field_mod,              only: field_type
  use field_collection_mod,   only: field_collection_type
  use file_mod,               only: FILE_MODE_WRITE, FILE_MODE_READ
  use io_context_mod,         only: io_context_type
  use linked_list_mod,        only: linked_list_type
  use lfric_xios_action_mod,  only: advance
  use lfric_xios_context_mod, only: lfric_xios_context_type
  use lfric_xios_file_mod,    only: lfric_xios_file_type, OPERATION_ONCE
  use log_mod,                only: log_event, log_scratch_space, &
                                    LOG_LEVEL_DEBUG, LOG_LEVEL_ERROR

  implicit none

  private
  public :: setup_checkpoint_io

contains

  !> @details Sets up checkpointing for the IO_demo app
  !> @param[in,out] modeldb   The model database
  !> @param[in]     chi       The co-ordinate field
  !> @param[in]     panel_id  The panel id field
  subroutine setup_checkpoint_io(modeldb, chi, panel_id)

    type(modeldb_type), intent(inout) :: modeldb
    type(field_type),   intent(in)    :: chi(:)
    type(field_type),   intent(in)    :: panel_id

    type(lfric_xios_context_type) :: tmp_io_context
    type(lfric_xios_context_type), pointer :: cp_context, io_context
    type(linked_list_type),        pointer :: file_list
    type(field_collection_type),   pointer :: checkpoint_fields
    real(r_second), allocatable            :: checkpoint_times(:)

    class(event_actor_type), pointer :: event_actor_ptr
    procedure(event_action), pointer :: context_advance

    character(len=str_max_filename) :: checkpoint_write_filename
    character(len=str_max_filename) :: checkpoint_read_filename
    character(len=str_def)          :: checkpoint_id
    integer(i_def) :: ts_start, ts_end, t_cp, freq_ts

    call log_event( 'io_demo: Setting up checkpoint I/O', LOG_LEVEL_DEBUG )

    ts_start = modeldb%calendar%parse_instance(modeldb%config%time%timestep_start())
    ts_end   = modeldb%calendar%parse_instance(modeldb%config%time%timestep_end())
    checkpoint_fields => modeldb%fields%get_field_collection("depository")

    call tmp_io_context%initialise( "checkpoint_context", start=ts_start, stop=ts_end )
    call modeldb%io_contexts%add_context(tmp_io_context)

    ! Get pointer to persistent context
    call modeldb%io_contexts%get_io_context("checkpoint_context", cp_context)
    file_list => cp_context%get_filelist()

    ! Set up file definitions for checkpoint writing
    if (modeldb%config%io%checkpoint_write()) then
      ! End of run checkpoint definition
      if (modeldb%config%io%end_of_run_checkpoint()) then
        write(checkpoint_write_filename, '(A,I0)') &
              trim(modeldb%config%files%checkpoint_stem_name()), ts_end
        call file_list%insert_item( lfric_xios_file_type( checkpoint_write_filename, &
                                          xios_id = "io_demo_checkpoint",            &
                                          io_mode = FILE_MODE_WRITE,                 &
                                          freq = ts_end - ts_start + 1,              &
                                          operation = OPERATION_ONCE,                &
                                          fields_in_file = checkpoint_fields ) )
      end if

      ! Flexible checkpoint definition
      checkpoint_times = modeldb%config%io%checkpoint_times()
      if (size(checkpoint_times) > 0) then
        do t_cp = 1, size(checkpoint_times)

          if (mod(checkpoint_times(t_cp), modeldb%clock%get_seconds_per_step()) /= 0) then
            write(log_scratch_space, '(A,F6.1, A)') "io_demo: Checkpoint time ", &
                                      checkpoint_times(t_cp),                    &
                                      " is not an integer multiple of the model timestep."
            call log_event(log_scratch_space, LOG_LEVEL_ERROR)
          else
              freq_ts = int(checkpoint_times(t_cp) / modeldb%clock%get_seconds_per_step())
              write(checkpoint_write_filename, '(A,I0)') &
                    trim(modeldb%config%files%checkpoint_stem_name()), freq_ts
              write(checkpoint_id, '(A,I0)') "io_demo_checkpoint_", freq_ts
              call file_list%insert_item( &
                    lfric_xios_file_type( trim(checkpoint_write_filename), &
                                          xios_id = trim(checkpoint_id),   &
                                          io_mode = FILE_MODE_WRITE,       &
                                          freq = freq_ts,                  &
                                          operation = OPERATION_ONCE,      &
                                          fields_in_file = checkpoint_fields ) )
          end if
        end do
      end if
    end if

    ! Set up file definitions for checkpoint reading
    if (modeldb%config%io%checkpoint_read()) then
      write(checkpoint_read_filename, '(A,I0)') &
            trim(modeldb%config%files%checkpoint_stem_name()), ts_start - 1
      write(checkpoint_id, '(A,I0)') "io_demo_restart_", ts_start - 1
      call file_list%insert_item( &
            lfric_xios_file_type( checkpoint_read_filename,      &
                                  xios_id = trim(checkpoint_id), &
                                  io_mode = FILE_MODE_READ,      &
                                  freq = 1,                      &
                                  operation = OPERATION_ONCE,    &
                                  fields_in_file = checkpoint_fields ) )
    end if

    ! Add checkpoint context to clock events so that it is advanced at each timestep
    event_actor_ptr => cp_context
    context_advance => advance
    call cp_context%initialise_xios_context( modeldb%mpi%get_comm(), chi, panel_id, &
                                             modeldb%clock, modeldb%calendar )

    call modeldb%clock%add_event(context_advance, event_actor_ptr)
    call cp_context%set_active(.true.)

    ! Set current context back to main
    call modeldb%io_contexts%get_io_context("io_demo", io_context)
    call io_context%set_current()

    nullify(cp_context)
    nullify(file_list)
    nullify(checkpoint_fields)
    nullify(io_context)

  end subroutine setup_checkpoint_io

end module io_demo_checkpoint_mod
