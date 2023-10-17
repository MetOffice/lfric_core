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

  use base_mesh_config_mod,    only: prime_mesh_name
  use calendar_mod,            only: calendar_type
  use constants_mod,           only: str_def, i_def
  use driver_model_data_mod,   only: model_data_type
  use empty_io_context_mod,    only: empty_io_context_type
  use field_mod,               only: field_type
  use inventory_by_mesh_mod,   only: inventory_by_mesh_type
  use io_context_mod,          only: io_context_type, callback_clock_arg
  use io_config_mod,           only: use_xios_io, subroutine_timers
  use log_mod,                 only: log_event, log_level_error, &
                                     log_level_trace
  use linked_list_mod,         only: linked_list_type
  use mesh_mod,                only: mesh_type
  use mesh_collection_mod,     only: mesh_collection
  use model_clock_mod,         only: model_clock_type
#ifdef USE_XIOS
  use lfric_xios_context_mod,  only: lfric_xios_context_type
#endif

  implicit none

  private
  public :: init_io, final_io, &
            get_io_context,    &
            filelist_populator

  class(io_context_type), allocatable, target :: model_context

  abstract interface
    subroutine filelist_populator(files_list, model_data)
      import linked_list_type, model_data_type
      type(linked_list_type), intent(out) :: files_list
      class(model_data_type), optional, target, intent(in) :: model_data
    end subroutine filelist_populator
  end interface

contains

  !> @brief  Initialises the model I/O
  !>
  !> @param[in] id                  A string identifier for the model
  !> @param[in] communicator        The ID for the model MPI communicator
  !> @param[in] chi_inventory       Contains the model's coordinate fields
  !> @param[in] panel_id_inventory  Contains the model's panel ID fields
  !> @param[in] model_clock         The model clock
  !> @param[in] calendar            The model calendar
  !> @param[in] populate_filelist   Optional procedure for creating a list of
  !!                                file descriptions used by the model I/O
  !> @param[in] model_data          Optional Model data object
  !> @param[in] alt_mesh_names      Optional array of names for other meshes
  !!                                to initialise I/O for
  !> @param[in] before_close        Optional routine to be called before
  !!                                context closes
  subroutine init_io( id, communicator,      &
                      chi_inventory,         &
                      panel_id_inventory,    &
                      model_clock, calendar, &
                      populate_filelist,     &
                      model_data,            &
                      alt_mesh_names,        &
                      before_close           )

    implicit none

    character(*),                     intent(in)    :: id
    integer(i_def),                   intent(in)    :: communicator
    type(inventory_by_mesh_type),     intent(in)    :: chi_inventory
    type(inventory_by_mesh_type),     intent(in)    :: panel_id_inventory
    type(model_clock_type),           intent(inout) :: model_clock
    class(calendar_type),             intent(in)    :: calendar
    procedure(filelist_populator), &
                   pointer, optional, intent(in)    :: populate_filelist
    class(model_data_type), optional, intent(in)    :: model_data
    character(len=str_def), optional, intent(in)    :: alt_mesh_names(:)
    procedure(callback_clock_arg), optional         :: before_close

    integer(i_def) :: rc
    procedure(callback_clock_arg), pointer :: before_close_ptr => null()

    ! Allocate IO context type based on model configuration
    if ( use_xios_io ) then
#ifdef USE_XIOS
      allocate( lfric_xios_context_type::model_context, stat=rc )
      if (rc /= 0) then
        call log_event( "Unable to allocate LFRic-XIOS context object", &
                        log_level_error )
      end if
      if (present(before_close)) then
        before_close_ptr => before_close
      end if
      call init_xios_io_context( model_context,      &
                                 id, communicator,   &
                                 chi_inventory,      &
                                 panel_id_inventory, &
                                 model_clock,        &
                                 calendar,           &
                                 before_close_ptr,   &
                                 populate_filelist,  &
                                 model_data,         &
                                 alt_mesh_names )
#else
      call log_event( "Cannot use XIOS I/O: model has not been built with " // &
                      "enabled", log_level_error )
#endif
    else
      allocate( empty_io_context_type::model_context, stat=rc )
      if (rc /= 0) then
        call log_event( "Unable to allocate empty context object", &
                        log_level_error )
      end if
    end if

  end subroutine init_io

#ifdef USE_XIOS
  !> @brief  Initialises an xios I/O context based on user input
  !>
  !> @param[in] io_context          The I/O context to be set up
  !> @param[in] id                  A string identifier for the model
  !> @param[in] communicator        The ID for the model MPI communicator
  !> @param[in] chi_inventory       Contains the model's coordinate fields
  !> @param[in] panel_id_inventory  Contains the model's panel ID fields
  !> @param[in] model_clock         The model clock
  !> @param[in] calendar            The model calendar
  !> @param[in] before_close        Routine to be called before context closes
  !> @param[in] populate_filelist   Optional procedure for creating a list of
  !!                                file descriptions used by the model I/O
  !> @param[in] model_data          Optional Model data object
  !> @param[in] alt_mesh_names      Optional array of names for other meshes
  !!                                to initialise I/O for
  subroutine init_xios_io_context( io_context,            &
                                   id, communicator,      &
                                   chi_inventory,         &
                                   panel_id_inventory,    &
                                   model_clock, calendar, &
                                   before_close,          &
                                   populate_filelist,     &
                                   model_data,            &
                                   alt_mesh_names )

    implicit none

    class(io_context_type), allocatable, intent(inout) :: io_context
    character(*),                        intent(in)    :: id
    integer(i_def),                      intent(in)    :: communicator
    type(inventory_by_mesh_type),        intent(in)    :: chi_inventory
    type(inventory_by_mesh_type),        intent(in)    :: panel_id_inventory
    type(model_clock_type),              intent(inout) :: model_clock
    class(calendar_type),                intent(in)    :: calendar
    procedure(callback_clock_arg), pointer, intent(in) :: before_close
    procedure(filelist_populator), &
                   pointer, optional,    intent(in)    :: populate_filelist
    class(model_data_type), optional,    intent(in)    :: model_data
    character(len=str_def), optional,    intent(in)    :: alt_mesh_names(:)

    type(mesh_type),                  pointer       :: mesh => null()
    type(field_type),                 pointer       :: chi(:) => null()
    type(field_type),                 pointer       :: panel_id => null()
    type(field_type),                 pointer       :: alt_chi_ptr(:) => null()
    type(field_type),                 pointer       :: alt_panel_id_ptr => null()
    type(field_type),                 allocatable   :: alt_coords(:,:)
    type(field_type),                 allocatable   :: alt_panel_ids(:)
    integer(kind=i_def)                             :: num_meshes, i, j

    type(linked_list_type), pointer :: file_list

    select type(io_context)
    type is (lfric_xios_context_type)
      if ( .not. allocated(model_context) ) then
        call log_event( "Cannot initialise unallocated xios context", log_level_error )
      end if

      ! Populate list of I/O files if procedure passed through
      if (present(populate_filelist)) then
        file_list => io_context%get_filelist()
        call populate_filelist(file_list, model_data)
      end if
      call io_context%set_timer_flag(subroutine_timers)

      ! Get coordinate fields for prime mesh
      mesh => mesh_collection%get_mesh(prime_mesh_name)
      call chi_inventory%get_field_array(mesh, chi)
      call panel_id_inventory%get_field(mesh, panel_id)

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

        call io_context%initialise( id, communicator,      &
                                    chi, panel_id,         &
                                    model_clock, calendar, &
                                    before_close,          &
                                    alt_coords, alt_panel_ids )
        deallocate(alt_coords)
        deallocate(alt_panel_ids)
      else
        call io_context%initialise( id, communicator,      &
                                    chi, panel_id,         &
                                    model_clock, calendar, &
                                    before_close )
      end if

      nullify(mesh, chi, panel_id, alt_chi_ptr, alt_panel_id_ptr)
    class default
      call log_event( "Passed a non lfric xios type context.", &
                      log_level_error )
    end select

  end subroutine init_xios_io_context
#endif

  !> @brief  Finalises the model I/O
  subroutine final_io()

    implicit none

    if (allocated(model_context)) deallocate(model_context)

  end subroutine final_io

  !> @brief  Returns the model io context.
  function get_io_context() result(context_ptr)

    implicit none

    class(io_context_type), pointer :: context_ptr

    context_ptr => model_context

  end function get_io_context

end module driver_io_mod
