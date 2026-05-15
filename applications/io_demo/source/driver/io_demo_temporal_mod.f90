!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Sets up temporal reading for the IO_demo app
module io_demo_temporal_mod

  use constants_mod,          only: i_def, str_max_filename
  use driver_modeldb_mod,     only: modeldb_type
  use event_mod,              only: event_action
  use event_actor_mod,        only: event_actor_type
  use field_mod,              only: field_type
  use field_collection_mod,   only: field_collection_type
  use field_parent_mod,       only: read_interface, write_interface
  use file_mod,               only: FILE_MODE_WRITE, FILE_MODE_READ
  use function_space_collection_mod,          only : function_space_collection
  use function_space_mod,                     only : function_space_type
  use fs_continuity_mod,      only: Wtheta
  use io_context_mod,         only: io_context_type, callback_clock_arg
  use linked_list_mod,        only: linked_list_type
  use lfric_xios_action_mod,  only: advance
  use lfric_xios_context_mod, only: lfric_xios_context_type
  use lfric_xios_file_mod,    only: lfric_xios_file_type, OPERATION_TIMESERIES
  use lfric_xios_read_mod,    only: read_field_generic
  use lfric_xios_write_mod,   only: write_field_generic
  use log_mod,                only: log_event, LOG_LEVEL_DEBUG
  use mesh_mod,               only: mesh_type

  implicit none

  private
  public :: init_temporal_fields, setup_temporal_io

contains

  !> @details Initialises fields for temporal I/O in the io_demo app
  !!
  !> @param[in]     mesh      The model mesh
  !> @param[in,out] modeldb   The model database
  subroutine init_temporal_fields(mesh, modeldb)

    type(mesh_type),    intent(in), pointer :: mesh
    type(modeldb_type), intent(inout)       :: modeldb

    type(field_type) :: monthly_field
    type(field_collection_type), pointer :: temporal_fields
    type(function_space_type),   pointer :: fs

    procedure(read_interface),  pointer :: read_method
    procedure(write_interface), pointer :: write_method

    ! Create field collection for temporal fields
    call modeldb%fields%add_empty_field_collection("temporal_fields")
    temporal_fields =>modeldb%fields%get_field_collection("temporal_fields")

    fs => function_space_collection%get_fs(mesh, 0, 0, Wtheta)
    call monthly_field%initialise(fs, name="monthly_field")

    ! Set up field to be read from monthly ancil file
    if (modeldb%config%io%use_xios_io()) then
       write_method => write_field_generic
       call monthly_field%set_write_behaviour(write_method)
       read_method => read_field_generic
       call monthly_field%set_read_behaviour(read_method)
    end if

    ! Add field to temporal field collection
    call temporal_fields%add_field(monthly_field)

  end subroutine init_temporal_fields

  !> @details Sets up a temporal reading context for the IO_demo app
  !> @param[in,out] modeldb   The model database
  !> @param[in]     chi       The co-ordinate field
  !> @param[in]     panel_id  The panel id field
  subroutine setup_temporal_io(modeldb, chi, panel_id)

    type(modeldb_type), intent(inout) :: modeldb
    type(field_type),   intent(in)    :: chi(:)
    type(field_type),   intent(in)    :: panel_id

    type(lfric_xios_context_type) :: tmp_io_context
    type(lfric_xios_context_type), pointer :: temporal_context, io_context
    type(linked_list_type),        pointer :: file_list
    type(field_collection_type),   pointer :: temporal_fields

    class(event_actor_type), pointer :: event_actor_ptr
    procedure(event_action), pointer :: context_advance
    procedure(callback_clock_arg), pointer :: before_close

    call log_event( 'io_demo: Setting up temporal I/O', LOG_LEVEL_DEBUG )

    temporal_fields => modeldb%fields%get_field_collection("temporal_fields")

    ! Set up new I/O context for temporal reading
    call tmp_io_context%initialise( "temporal_context", &
                                    start=modeldb%calendar%parse_instance(modeldb%config%time%timestep_start()), &
                                    stop=modeldb%calendar%parse_instance(modeldb%config%time%timestep_end()) )
    ! Add context to modeldb
    call modeldb%io_contexts%add_context(tmp_io_context)

    ! Get pointer to context from modeldb - this context is persistent beyond the scope of routine
    call modeldb%io_contexts%get_io_context("temporal_context", temporal_context)
    file_list => temporal_context%get_filelist()

    ! Set up definition of temporal read file - we use a monthly ancil as an example
    call file_list%insert_item( &
              lfric_xios_file_type( modeldb%config%files%temporal_file_path(), &
                                    xios_id = "monthly_ancil",                 &
                                    io_mode = FILE_MODE_READ,                  &
                                    operation = OPERATION_TIMESERIES,          &
                                    cyclic = .true.,                           &
                                    fields_in_file = temporal_fields ) )

    if (modeldb%config%io%write_diag()) then
       ! Set up definition of temporal write file - this will contain the time series output of the model
       call file_list%insert_item( &
                      lfric_xios_file_type( "io_demo_temporal_diag",          &
                                            xios_id = "temporal_diag",        &
                                            io_mode = FILE_MODE_WRITE,        &
                                            operation = OPERATION_TIMESERIES, &
                                            freq = 1,                         &
                                            fields_in_file = temporal_fields ) )
    end if

    ! Initialise the XIOS context attached to the temporal context object
    before_close => null()
    call temporal_context%initialise_xios_context( modeldb%mpi%get_comm(), chi, panel_id, &
                                             modeldb%clock, modeldb%calendar, before_close )

    ! Add context object to the model clock's event loop, this means that the
    ! temporal context will be advanced at each model time step, and the
    ! appropriate files read/written to
    event_actor_ptr => temporal_context
    context_advance => advance
    call modeldb%clock%add_event(context_advance, event_actor_ptr)
    call temporal_context%set_active(.true.)

    ! Set current context back to main
    call modeldb%io_contexts%get_io_context("io_demo", io_context)
    call io_context%set_current()

    nullify(temporal_context)
    nullify(file_list)
    nullify(temporal_fields)
    nullify(io_context)

  end subroutine setup_temporal_io

end module io_demo_temporal_mod
