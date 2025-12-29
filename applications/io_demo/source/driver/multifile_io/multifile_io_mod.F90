!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Module controlling the initialisation, stepping, and finalisation of
!> the multifile IO.
module multifile_io_mod

  use base_mesh_config_mod,    only: prime_mesh_name
  use calendar_mod,            only: calendar_type
  use constants_mod,           only: str_def, i_def
  use driver_modeldb_mod,      only: modeldb_type
  use driver_model_data_mod,   only: model_data_type
  use empty_io_context_mod,    only: empty_io_context_type
  use event_mod,               only: event_action
  use event_actor_mod,         only: event_actor_type
  use field_collection_mod,    only: field_collection_type
  use field_mod,               only: field_type
  use multifile_file_setup_mod,only: init_multifile_files
  use inventory_by_mesh_mod,   only: inventory_by_mesh_type
  use io_context_collection_mod, only: io_context_collection_type
  use io_context_mod,          only: io_context_type, callback_clock_arg
  use io_config_mod,           only: use_xios_io
  use log_mod,                 only: log_event, log_level_error, &
                                     log_level_trace, log_level_info, &
                                     log_scratch_space
  use lfric_xios_context_mod,  only: lfric_xios_context_type
  use linked_list_mod,         only: linked_list_type
  use lfric_xios_action_mod,   only: advance_read_only
  use mesh_mod,                only: mesh_type
  use mesh_collection_mod,     only: mesh_collection
  use model_clock_mod,         only: model_clock_type
  use namelist_mod,            only: namelist_type
  use step_calendar_mod,       only: step_calendar_type

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

    type(lfric_xios_context_type), pointer          :: io_context
    type(namelist_type), pointer :: multifile_nml

    character(str_def) :: context_name

    integer(i_def) :: multifile_start_timestep
    integer(i_def) :: multifile_stop_timestep
    character(str_def) :: filename
    character(str_def), allocatable :: multifile_io_profiles(:)
    integer(i_def) :: i

    type(linked_list_type), pointer :: file_list

    allocate(multifile_io_profiles, source=modeldb%configuration%get_namelist_profiles("multifile_io"))

    do i=1, size(multifile_io_profiles)

      multifile_nml => modeldb%configuration%get_namelist('multifile_io', &
        profile_name=trim(multifile_io_profiles(i)))
      call multifile_nml%get_value('filename', filename)
      call multifile_nml%get_value('start_timestep', multifile_start_timestep)
      call multifile_nml%get_value('stop_timestep', multifile_stop_timestep)
      context_name = "multifile_context_" // trim(filename)
      call context_init(modeldb, context_name, multifile_start_timestep, &
                        multifile_stop_timestep)

      call modeldb%io_contexts%get_io_context(context_name, io_context)

      file_list => io_context%get_filelist()
      call init_multifile_files(file_list, modeldb, filename)

    end do

    deallocate(multifile_io_profiles)

  end subroutine init_multifile_io

  !> @brief Step the multifile IO
  !> @param[inout] modeldb            Model database object
  !> @param[in]    chi_inventory      Inventory object, containing all of
  !!                                  the chi fields indexed by mesh
  !> @param[in]    panel_id_inventory Inventory object, containing all of
  !!                                  the fields with the ID of mesh panels
  subroutine step_multifile_io(modeldb, chi_inventory, panel_id_inventory)
    implicit none
    type(modeldb_type),                  intent(inout) :: modeldb
    type(inventory_by_mesh_type),        intent(in)    :: chi_inventory
    type(inventory_by_mesh_type),        intent(in)    :: panel_id_inventory

    type(lfric_xios_context_type), pointer          :: io_context
    class(event_actor_type), pointer :: event_actor_ptr
    type(mesh_type),                  pointer       :: mesh => null()
    type(field_type),                 pointer       :: chi(:) => null()
    type(field_type),                 pointer       :: panel_id => null()

    class(calendar_type), allocatable :: tmp_calendar

    type(namelist_type), pointer :: multifile_nml
    type(namelist_type), pointer :: time_nml
    character(str_def) :: context_name
    character(str_def) :: filename

    character(str_def) :: time_origin
    character(str_def) :: time_start

    character(str_def), allocatable :: multifile_io_profiles(:)
    integer(i_def) :: i

    procedure(event_action), pointer :: context_advance
    procedure(callback_clock_arg), pointer :: before_close


    nullify(before_close)

    allocate(multifile_io_profiles, source=modeldb%configuration%get_namelist_profiles("multifile_io"))

    do i=1, size(multifile_io_profiles)

      multifile_nml => modeldb%configuration%get_namelist('multifile_io', &
        profile_name=trim(multifile_io_profiles(i)))
      call multifile_nml%get_value('filename', filename)

      context_name = "multifile_context_" // trim(filename)
      call modeldb%io_contexts%get_io_context(context_name, io_context)

      if (modeldb%clock%get_step() == io_context%get_stop_time()) then
        ! Finalise XIOS context
        call io_context%set_current()
        call io_context%set_active(.false.)
        call modeldb%clock%remove_event(context_name)
        call io_context%finalise_xios_context()

      elseif (modeldb%clock%get_step() == io_context%get_start_time()) then
        ! Initialise XIOS context
        mesh => mesh_collection%get_mesh(prime_mesh_name)
        call chi_inventory%get_field_array(mesh, chi)
        call panel_id_inventory%get_field(mesh, panel_id)

        time_nml => modeldb%configuration%get_namelist('time')

        call time_nml%get_value('calendar_origin', time_origin)
        call time_nml%get_value('calendar_start', time_start)

        allocate(tmp_calendar, source=step_calendar_type(time_origin, time_start))

        call io_context%initialise_xios_context( modeldb%mpi%get_comm(),      &
                                                 chi, panel_id,               &
                                                 modeldb%clock, tmp_calendar, &
                                                 before_close, start_at_zero=.true. )

        ! Attach context advancement to the model's clock
        context_advance => advance_read_only
        event_actor_ptr => io_context
        call modeldb%clock%add_event( context_advance, event_actor_ptr )
        call io_context%set_active(.true.)
      end if
    end do

    call modeldb%io_contexts%get_io_context("io_demo", io_context)
    call io_context%set_current()

    deallocate(multifile_io_profiles)

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
    character(*), intent(in) :: context_name
    integer(i_def), intent(in) :: multifile_start_timestep
    integer(i_def), intent(in) :: multifile_stop_timestep

    type(lfric_xios_context_type) :: tmp_io_context

    call tmp_io_context%initialise( context_name, start=multifile_start_timestep, &
                                    stop=multifile_stop_timestep )
    call modeldb%io_contexts%add_context(tmp_io_context)

  end subroutine context_init

end module multifile_io_mod
