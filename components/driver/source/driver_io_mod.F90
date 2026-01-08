!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Module controlling the initialisation and finalisation of IO within the
!! driver layer. This module also contains the model's io_context object and
!! associated getter routines
!>
module driver_io_mod

  use constants_mod,           only: str_def, i_def
  use driver_modeldb_mod,      only: modeldb_type
  use driver_model_data_mod,   only: model_data_type
  use empty_io_context_mod,    only: empty_io_context_type
  use event_mod,               only: event_action
  use event_actor_mod,         only: event_actor_type
  use field_collection_mod,    only: field_collection_type
  use field_mod,               only: field_type
  use inventory_by_mesh_mod,   only: inventory_by_mesh_type
  use io_context_collection_mod, only: io_context_collection_type
  use io_context_mod,          only: io_context_type, callback_clock_arg
  use log_mod,                 only: log_event, log_level_error, &
                                     log_level_trace, log_level_info, &
                                     log_scratch_space
#ifdef USE_XIOS
  use lfric_xios_context_mod,  only: lfric_xios_context_type
#endif
  use linked_list_mod,         only: linked_list_type
  use lfric_xios_action_mod,   only: advance
  use mesh_mod,                only: mesh_type
  use mesh_collection_mod,     only: mesh_collection
  use model_clock_mod,         only: model_clock_type
  use namelist_mod,            only: namelist_type

  implicit none

  private
  public :: init_io, final_io, &
            filelist_populator

  abstract interface
    subroutine filelist_populator(files_list, modeldb)
      import linked_list_type, modeldb_type
      type(linked_list_type), intent(out) :: files_list
      type(modeldb_type), optional, intent(inout) :: modeldb
    end subroutine filelist_populator
  end interface

contains

  !> @brief  Initialises the model I/O
  !>
  !> @param[in] context_name             Unique context identifier
  !> @param[in] mesh_name                Mesh used by this context
  !> @param[inout] modeldb               Model state
  !> @param[in] chi_inventory            Contains the model's coordinate fields
  !> @param[in] panel_id_inventory       Contains the model's panel ID fields
  !> @param[in] populate_filelist        Optional procedure for creating a list of
  !!                                     file descriptions used by the model I/O

  !> @param[in] alt_mesh_names           Optional array of names for other meshes
  !!                                     to initialise I/O for
  !> @param[in] before_close             Optional routine to be called before
  !!                                     context closes
  subroutine init_io( context_name,       &
                      mesh_name,          &
                      modeldb,            &
                      chi_inventory,      &
                      panel_id_inventory, &
                      populate_filelist,  &
                      alt_mesh_names,     &
                      before_close )

    implicit none

    character(*),                     intent(in)    :: context_name
    character(*),                     intent(in)    :: mesh_name
    class(modeldb_type),              intent(inout) :: modeldb
    type(inventory_by_mesh_type),     intent(in)    :: chi_inventory
    type(inventory_by_mesh_type),     intent(in)    :: panel_id_inventory
    procedure(filelist_populator), &
                   pointer, optional, intent(in)    :: populate_filelist
    character(len=str_def), optional, intent(in)    :: alt_mesh_names(:)
    procedure(callback_clock_arg), optional         :: before_close

    procedure(callback_clock_arg), pointer :: before_close_ptr

    type(namelist_type), pointer :: io_nml

    logical :: use_xios_io

    io_nml => modeldb%configuration%get_namelist('io')
    call io_nml%get_value( 'use_xios_io', use_xios_io )

    ! Allocate IO context type based on model configuration
    if ( use_xios_io ) then
#ifdef USE_XIOS
      if (present(before_close)) then
        before_close_ptr => before_close
      else
        before_close_ptr => null()
      end if

      call init_xios_io_context( context_name,       &
                                 mesh_name,          &
                                 modeldb,            &
                                 chi_inventory,      &
                                 panel_id_inventory, &
                                 before_close_ptr,   &
                                 populate_filelist,  &
                                 alt_mesh_names )
#else
      call log_event( "Cannot use XIOS I/O: model has not been built with " // &
                      "XIOS enabled", log_level_error )
#endif
    else
      call init_empty_io_context(context_name, modeldb%io_contexts)
    end if

  end subroutine init_io

  !> @brief Finalises the model I/O
  !> @param[inout] modeldb The modeldb in which the contexts are stored
  subroutine final_io(modeldb)
    implicit none

    class(modeldb_type), intent(inout) :: modeldb

    call modeldb%io_contexts%clear()

  end subroutine final_io


  !> @brief Initialises an empty I/O context and puts it in and I/O context collection
  !> @param[in]    context_name  Unique context identifier
  !> @param[inout] modeldb       Model state
  subroutine init_empty_io_context(context_name, io_context_collection)

    implicit none
    character(*),                     intent(in)    :: context_name
    type(io_context_collection_type), intent(inout) :: io_context_collection
    type(empty_io_context_type) :: tmp_io_context

    call tmp_io_context%initialise(context_name)
    call io_context_collection%add_context(tmp_io_context)
    write(log_scratch_space, "(A25)")"Adding empty IO context: " // context_name
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

  end subroutine init_empty_io_context

#ifdef USE_XIOS
  !> @brief  Initialises an xios I/O context based on user input
  !>
  !> @param[in] context_name        Unique context identifier
  !> @param[in] mesh_name           Mesh used by this context
  !> @param[in] modeldb             Model state
  !> @param[in] chi_inventory       Contains the model's coordinate fields
  !> @param[in] panel_id_inventory  Contains the model's panel ID fields
  !> @param[in] before_close        Routine to be called before context closes
  !> @param[in] populate_filelist   Optional procedure for creating a list of
  !!                                file descriptions used by the model I/O
  !> @param[in] alt_mesh_names      Optional array of names for other meshes
  !!                                to initialise I/O for
  subroutine init_xios_io_context( context_name,       &
                                   mesh_name,          &
                                   modeldb,            &
                                   chi_inventory,      &
                                   panel_id_inventory, &
                                   before_close,       &
                                   populate_filelist,  &
                                   alt_mesh_names )

    implicit none

    character(*),                        intent(in)    :: context_name
    character(*),                        intent(in)    :: mesh_name
    class(modeldb_type),                 intent(inout) :: modeldb
    type(inventory_by_mesh_type),        intent(in)    :: chi_inventory
    type(inventory_by_mesh_type),        intent(in)    :: panel_id_inventory
    procedure(callback_clock_arg), pointer, intent(in) :: before_close
    procedure(filelist_populator), &
                   pointer, optional,    intent(in)    :: populate_filelist
    character(len=str_def), optional,    intent(in)    :: alt_mesh_names(:)

    type(mesh_type),  pointer     :: mesh
    type(field_type), pointer     :: chi(:)
    type(field_type), pointer     :: panel_id
    type(field_type), pointer     :: alt_chi_ptr(:)
    type(field_type), pointer     :: alt_panel_id_ptr
    type(field_type), allocatable :: alt_coords(:,:)
    type(field_type), allocatable :: alt_panel_ids(:)

    type(lfric_xios_context_type)          :: tmp_io_context
    type(lfric_xios_context_type), pointer :: io_context
    class(event_actor_type),       pointer :: event_actor_ptr

    type(linked_list_type), pointer :: file_list
    procedure(event_action), pointer :: context_advance

    integer(i_def) :: num_meshes, i, j

    type(namelist_type), pointer :: io_nml
    logical :: subroutine_timers

    io_nml => modeldb%configuration%get_namelist('io')
    call io_nml%get_value( 'subroutine_timers', subroutine_timers )

    subroutine_timers = .false.

    mesh             => null()
    chi              => null()
    panel_id         => null()
    alt_chi_ptr      => null()
    alt_panel_id_ptr => null()

    mesh => mesh_collection%get_mesh(mesh_name)

    !==============================================================

    call tmp_io_context%initialise(context_name)
    call modeldb%io_contexts%add_context(tmp_io_context)
    call modeldb%io_contexts%get_io_context(context_name, io_context)

    ! Populate list of I/O files if procedure passed through
    if (present(populate_filelist)) then
      ! Get the filelist in the io_context and pass it to the populating
      ! procedure which was passed to this routine.
      file_list => io_context%get_filelist()
      call populate_filelist(file_list, modeldb)
    end if
    call io_context%set_timer_flag(subroutine_timers)

    ! ===============================
    ! Check that a mesh exists
    ! ===============================

    if (associated(mesh)) then
      call chi_inventory%get_field_array(mesh, chi)
      call panel_id_inventory%get_field(mesh, panel_id)
    else
      write(log_scratch_space,'(A)') trim(mesh_name) // &
          'mesh not associated, skipping xios io ' //   &
          'initialisation of this mesh for this partition.'
      call log_event(log_scratch_space, log_level_info)
      return
    end if

    ! Unpack alternative meshes and get their coordinates to pass to I/O
    if (present(alt_mesh_names)) then
      num_meshes = SIZE(alt_mesh_names)
      allocate(alt_coords(num_meshes,3))
      allocate(alt_panel_ids(num_meshes))

      do i = 1, num_meshes
        mesh => mesh_collection%get_mesh(alt_mesh_names(i))
        call chi_inventory%get_field_array(mesh, alt_chi_ptr)
        call panel_id_inventory%get_field(mesh, alt_panel_id_ptr)
        ! Copy into array to give to I/O
        do j =1,3
          call alt_chi_ptr(j)%copy_field_serial(alt_coords(i,j))
        end do
        call alt_panel_id_ptr%copy_field_serial(alt_panel_ids(i))
      end do

      call io_context%initialise_xios_context( modeldb%mpi%get_comm(), &
                                               chi, panel_id,          &
                                               modeldb%clock,          &
                                               modeldb%calendar,       &
                                               before_close,           &
                                               alt_coords,             &
                                               alt_panel_ids )
      deallocate(alt_coords)
      deallocate(alt_panel_ids)
    else
      call io_context%initialise_xios_context( modeldb%mpi%get_comm(), &
                                               chi, panel_id,          &
                                               modeldb%clock,          &
                                               modeldb%calendar,       &
                                               before_close )
    end if

    ! Attach context advancement to the model's clock
    context_advance => advance
    event_actor_ptr => io_context
    call modeldb%clock%add_event( context_advance, event_actor_ptr )

    nullify(mesh, chi, panel_id, alt_chi_ptr, alt_panel_id_ptr)

  end subroutine init_xios_io_context
#endif

end module driver_io_mod
