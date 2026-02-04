!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Module controlling the initialisation, stepping, and finalisation of
!> the multifile IO.
module multifile_io_mod

  use calendar_mod,            only: calendar_type
  use constants_mod,           only: str_def, i_def
  use driver_model_data_mod,   only: model_data_type
  use driver_modeldb_mod,      only: modeldb_type
  use empty_io_context_mod,    only: empty_io_context_type
  use event_mod,               only: event_action
  use event_actor_mod,         only: event_actor_type
  use field_collection_mod,    only: field_collection_type
  use field_mod,               only: field_type
  use multifile_file_setup_mod,only: init_multifile_files
  use inventory_by_mesh_mod,   only: inventory_by_mesh_type
  use io_context_collection_mod, only: io_context_collection_type
  use io_context_mod,          only: io_context_type, callback_clock_arg
  use log_mod,                 only: log_event, log_level_error, &
                                     log_level_trace, log_level_info, &
                                     log_scratch_space
  use lfric_xios_context_mod,  only: lfric_xios_context_type
  use linked_list_mod,         only: linked_list_type
  use lfric_xios_action_mod,   only: advance_read_only
  use mesh_mod,                only: mesh_type
  use mesh_collection_mod,     only: mesh_collection
  use model_clock_mod,         only: model_clock_type
  use step_calendar_mod,       only: step_calendar_type

  use multifile_io_nml_iterator_mod, only: multifile_io_nml_iterator_type
  use multifile_io_nml_mod,          only: multifile_io_nml_type

  implicit none

  private

  public :: init_multifile_io
  public :: step_multifile_io
  private :: context_init

contains

  !> @brief Initialise the multifile IO
  !> @param[inout] modeldb Modeldb object
  subroutine init_multifile_io(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    type(lfric_xios_context_type), pointer :: io_context

    character(str_def) :: context_name

    integer(i_def)     :: start_timestep
    integer(i_def)     :: stop_timestep
    character(str_def) :: filename

    type(linked_list_type), pointer :: file_list

    type(multifile_io_nml_iterator_type) :: iter
    type(multifile_io_nml_type), pointer :: multifile_nml

    call iter%initialise(modeldb%config%multifile_io)
    do while (iter%has_next())

      multifile_nml => iter%next()

      filename       = multifile_nml%filename()
      start_timestep = multifile_nml%start_timestep()
      stop_timestep  = multifile_nml%stop_timestep()

      context_name = "multifile_context_" // trim(filename)
      call context_init(modeldb, context_name, start_timestep, stop_timestep)

      call modeldb%io_contexts%get_io_context(context_name, io_context)

      file_list => io_context%get_filelist()
      call init_multifile_files(file_list, modeldb, filename)

    end do

  end subroutine init_multifile_io

  !> @brief Step the multifile IO
  !> @param[inout] modeldb            Model database object
  !> @param[in]    chi_inventory      Inventory object, containing all of
  !!                                  the chi fields indexed by mesh
  !> @param[in]    panel_id_inventory Inventory object, containing all of
  !!                                  the fields with the ID of mesh panels
  subroutine step_multifile_io(modeldb, chi_inventory, panel_id_inventory)

    implicit none

    type(modeldb_type),           intent(inout) :: modeldb
    type(inventory_by_mesh_type), intent(in)    :: chi_inventory
    type(inventory_by_mesh_type), intent(in)    :: panel_id_inventory

    type(lfric_xios_context_type), pointer :: io_context

    class(event_actor_type), pointer :: event_actor_ptr

    type(mesh_type),  pointer :: mesh
    type(field_type), pointer :: chi(:)
    type(field_type), pointer :: panel_id

    class(calendar_type), allocatable :: tmp_calendar

    type(multifile_io_nml_iterator_type) :: iter
    type(multifile_io_nml_type), pointer :: multifile_nml

    character(str_def) :: context_name
    character(str_def) :: prime_mesh_name
    character(str_def) :: filename
    character(str_def) :: time_origin
    character(str_def) :: time_start

    procedure(event_action), pointer :: context_advance
    procedure(callback_clock_arg), pointer :: before_close

    nullify(mesh)
    nullify(chi)
    nullify(panel_id)
    nullify(before_close)

    call iter%initialise(modeldb%config%multifile_io)
    do while (iter%has_next())

      multifile_nml => iter%next()

      filename     = multifile_nml%filename()
      context_name = "multifile_context_" // trim(filename)

      call modeldb%io_contexts%get_io_context(context_name, io_context)

      if (modeldb%clock%get_step() == io_context%get_stop_time()) then

        ! Finalise XIOS context
        call io_context%set_current()
        call io_context%set_active(.false.)
        call modeldb%clock%remove_event(context_name)
        call io_context%finalise_xios_context()

      else if (modeldb%clock%get_step() == io_context%get_start_time()) then

        prime_mesh_name = modeldb%config%base_mesh%prime_mesh_name()
        time_origin     = modeldb%config%time%calendar_origin()
        time_start      = modeldb%config%time%calendar_start()

        ! Initialise XIOS context
        mesh => mesh_collection%get_mesh(prime_mesh_name)
        call chi_inventory%get_field_array(mesh, chi)
        call panel_id_inventory%get_field(mesh, panel_id)

        allocate(tmp_calendar, source=step_calendar_type(time_origin, time_start))

        call io_context%initialise_xios_context( modeldb%mpi%get_comm(),      &
                                                 chi, panel_id,               &
                                                 modeldb%clock, tmp_calendar, &
                                                 before_close,                &
                                                 start_at_zero=.true. )

        ! Attach context advancement to the model's clock
        context_advance => advance_read_only
        event_actor_ptr => io_context
        call modeldb%clock%add_event( context_advance, event_actor_ptr )
        call io_context%set_active(.true.)
      end if
    end do

    call modeldb%io_contexts%get_io_context("io_demo", io_context)
    call io_context%set_current()

  end subroutine step_multifile_io

  !> @brief Helper function for initialising the lfric IO context and adding it
  !>        to the io_contexts collection in the modeldb.
  !> @param[inout] modeldb                   Model database object
  !> @param[in]    context_name              The name of the IO context
  !> @param[in]    multifile_start_timestep  The start time step for the IO context
  !> @param[in]    multifile_stop_timestep   The end time step for the IO context
  subroutine context_init(modeldb, &
                          context_name, &
                          multifile_start_timestep, &
                          multifile_stop_timestep)
    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    character(*),   intent(in) :: context_name
    integer(i_def), intent(in) :: multifile_start_timestep
    integer(i_def), intent(in) :: multifile_stop_timestep

    type(lfric_xios_context_type) :: tmp_io_context

    call tmp_io_context%initialise( context_name, start=multifile_start_timestep, &
                                    stop=multifile_stop_timestep )
    call modeldb%io_contexts%add_context(tmp_io_context)

  end subroutine context_init

end module multifile_io_mod
