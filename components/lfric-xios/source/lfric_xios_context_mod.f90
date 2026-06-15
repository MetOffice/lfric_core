!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Wrap the XIOS context in an object for easier management and cleaner code.
!>
module lfric_xios_context_mod

  use calendar_mod,         only : calendar_type
  use clock_mod,            only : clock_type
  use constants_mod,        only : i_def, &
                                   r_second, i_timestep, &
                                   l_def
  use field_mod,            only : field_type
  use file_mod,             only : file_type
  use io_context_mod,       only : io_context_type, callback_clock_arg
  use io_config_mod,        only : file_convention,       &
                                   file_convention_ugrid
  use lfric_xios_file_mod,  only : lfric_xios_file_type
  use lfric_mpi_mod,        only : lfric_comm_type
  use log_mod,              only : log_event, log_scratch_space, &
                                   log_level_error, log_level_debug
  use lfric_xios_setup_mod, only : init_xios_calendar,   &
                                   init_xios_dimensions, &
                                   setup_xios_files
  use lfric_xios_file_mod,  only : lfric_xios_file_type
  use lfric_xios_process_output_mod, only: process_output_file
  use linked_list_mod,      only : linked_list_type, linked_list_item_type
  use mesh_mod,             only : mesh_type
  use model_clock_mod,      only : model_clock_type
  use timing_mod,           only : start_timing, stop_timing, &
                                   tik, LPROF
  use xios,                 only : xios_context,                  &
                                   xios_context_initialize,       &
                                   xios_close_context_definition, &
                                   xios_context_finalize,         &
                                   xios_get_handle,               &
                                   xios_set_current_context
  use mod_wait,             only : init_wait

  implicit none

  private

  !> Contains an instance of an XIOS context and manages interactions between
  !> the model and the context.
  type, public, extends(io_context_type) :: lfric_xios_context_type
    private

    type(xios_context)     :: handle
    type(linked_list_type) :: filelist
    integer(i_def)         :: context_clock_step = 1_i_def

    logical :: uses_timer = .false.
    logical :: xios_context_initialised = .false.
    !> Flag denoting if this file is a UGRID Planar mesh file with
    !> projected coordinates that have been scaled
    logical :: ugrid_scaled_projected_coordinates = .false.

  contains
    private
    procedure, public :: initialise => initialise_lfric_xios_context
    procedure, public :: initialise_xios_context
    procedure, public :: get_filelist
    procedure, public :: set_current
    procedure, public :: tick_context_clock
    procedure, public :: get_context_clock_step
    procedure, public :: finalise_xios_context
    final :: finalise
  end type lfric_xios_context_type

contains

  !> @brief Set up an LFRic-XIOS context object.
  !>
  !> @param [in] name Unique identifying string.
  subroutine initialise_lfric_xios_context(this, name, start, stop)
    class(lfric_xios_context_type), intent(inout) :: this
    character(*), intent(in) :: name
    integer(i_timestep), optional, intent(in) :: start
    integer(i_timestep), optional, intent(in) :: stop

    ! Initialise the parent
    call this%initialise_io_context(name, start, stop)

  end subroutine initialise_lfric_xios_context

  !> @brief Set up an XIOS context.
  !>
  !> @param [in]     communicator      MPI communicator used by context.
  !> @param [in]     chi               Array of coordinate fields
  !> @param [in]     panel_id          Panel ID field
  !> @param [in]     model_clock       The model clock.
  !> @param [in]     calendar          The model calendar.
  !> @param [in]     before_close      Routine to be called before context closes
  !> @param [in]     alt_coords        Array of coordinate fields for alternative meshes
  !> @param [in]     alt_panel_ids     Panel ID fields for alternative meshes
  subroutine initialise_xios_context( this, communicator,    &
                                      chi, panel_id,         &
                                      model_clock, calendar, &
                                      before_close,          &
                                      alt_coords,            &
                                      alt_panel_ids,         &
                                      start_at_zero )

    implicit none

    class(lfric_xios_context_type), intent(inout) :: this
    type(lfric_comm_type),          intent(in)    :: communicator
    type(field_type),               intent(in)    :: chi(:)
    type(field_type),               intent(in)    :: panel_id
    type(model_clock_type),         intent(inout) :: model_clock
    class(calendar_type),           intent(in)    :: calendar
    procedure(callback_clock_arg), pointer, &
                                    intent(in)    :: before_close
    type(field_type),     optional, intent(in)    :: alt_coords(:,:)
    type(field_type),     optional, intent(in)    :: alt_panel_ids(:)
    logical,              optional, intent(in)    :: start_at_zero

    type(mesh_type), pointer             :: mesh => null()
    type(linked_list_item_type), pointer :: loop => null()
    type(lfric_xios_file_type),  pointer :: file => null()
    logical :: zero_start
    integer(tik) :: timing_idlx, timing_idxc

    write(log_scratch_space, "(A)") &
        "Initialising XIOS context: " // this%get_context_name()
    call log_event(log_scratch_space, log_level_debug)
    if ( LPROF ) call start_timing(timing_idlx, 'lfric_xios.init_context')

    if (present(start_at_zero)) then
      zero_start = start_at_zero
    else
      zero_start = .false.
    end if

    call xios_context_initialize( this%get_context_name(), &
                                  communicator%get_comm_mpi_val() )
    call xios_get_handle( this%get_context_name(), this%handle )
    call xios_set_current_context( this%handle )

    ! Run XIOS setup routines
    call init_xios_calendar(model_clock, calendar, zero_start, this%context_clock_step)

    call init_xios_dimensions(chi, panel_id, alt_coords, alt_panel_ids)
    ! Obtain information on whether the mesh is ugrid and planar here?
    ! This is to inform decisions on file post processing work around code path.
    mesh => chi(1)%get_mesh()
    if ( mesh%is_geometry_planar() .and. &
                    file_convention == file_convention_ugrid ) then
      this%ugrid_scaled_projected_coordinates = .true.
    end if
    if (this%filelist%get_length() > 0) call setup_xios_files(this%filelist)

    if (associated(before_close)) call before_close(model_clock)

    ! Close the context definition - no more I/O configuration operations
    ! can be defined after this point
    if ( LPROF ) call start_timing(timing_idxc, 'xios.close_context_definition')
    call log_event('XIOS context definition closing', log_level_debug)
    call xios_close_context_definition()
    if ( LPROF ) call stop_timing(timing_idxc, 'xios.close_context_definition')
    call log_event('XIOS context definition closed', log_level_debug)

    this%xios_context_initialised = .true.

    ! Read all files that need to be read from
    if (this%filelist%get_length() > 0) then
      loop => this%filelist%get_head()
      do while (associated(loop))
        select type(list_item => loop%payload)
          type is (lfric_xios_file_type)
            file => list_item
            if (file%mode_is_read()) call file%recv_fields()
        end select
        loop => loop%next
      end do
    end if
    if ( LPROF ) call stop_timing(timing_idlx, 'lfric_xios.init_context')

  end subroutine initialise_xios_context

  subroutine finalise( this )
    implicit none

    type(lfric_xios_context_type), intent(inout) :: this

    call this%finalise_xios_context()

  end subroutine finalise

  !> Finaliser for lfric_xios_context object.
  subroutine finalise_xios_context( this )

    implicit none

    class(lfric_xios_context_type), intent(inout) :: this

    type(linked_list_item_type), pointer :: loop => null()
    type(lfric_xios_file_type),  pointer :: file => null()
    integer(tik) :: timing_idlx, timing_idxc


    if (this%xios_context_initialised) then
      if ( LPROF ) call start_timing(timing_idlx, 'lfric_xios.finalise_context')
      call log_event( 'Finalising XIOS context: ' // this%get_context_name(), LOG_LEVEL_DEBUG )
      call this%set_current()

      ! Perform final write
      if (this%filelist%get_length() > 0) then
        loop => this%filelist%get_head()
        do while (associated(loop))
          select type( list_item => loop%payload )
            type is (lfric_xios_file_type)
              file => list_item
              if (file%mode_is_write()) call file%send_fields()
          end select
          loop => loop%next
        end do
      end if

      ! Finalise the XIOS context - all data will be written to disk and files
      ! will be closed.
      write(log_scratch_space, "(A)") "Finalising XIOS context: " // this%get_context_name()
      call log_event(log_scratch_space, log_level_debug)
      if ( LPROF ) call start_timing(timing_idxc, 'xios.context_finalize')
      call xios_context_finalize()
      if ( LPROF ) call stop_timing(timing_idxc, 'xios.context_finalize')

      ! Only take action if this is a regional model with UGRID Projected
      ! coordinates, as these are awaiting XIOS feature development
      if ( this%ugrid_scaled_projected_coordinates ) then
        call log_event("Closing file for post processing.", LOG_LEVEL_DEBUG)
        ! We have closed the context on our end, but we need to make sure that XIOS
        ! has closed the files for all servers before we process them.
        call init_wait()

        ! Process and close all files in list
        if (this%filelist%get_length() > 0) then
          loop => this%filelist%get_head()
          do while (associated(loop))
            select type( list_item => loop%payload )
              type is (lfric_xios_file_type)
                file => list_item
                if (file%mode_is_write()) call process_output_file(file)
                call file%file_close()
            end select
            loop => loop%next
          end do
        end if
      end if

      this%xios_context_initialised = .false.
      if ( LPROF ) call stop_timing(timing_idlx, 'lfric_xios.finalise_context')
    end if
    nullify(loop)
    nullify(file)

  end subroutine finalise_xios_context

  !> Gets the file list associated with this context.
  !>
  !> @return Linked list of file objects
  function get_filelist( this ) result(filelist)

    implicit none

    class(lfric_xios_context_type), intent(in), target :: this
    type(linked_list_type), pointer :: filelist

    filelist => this%filelist

  end function get_filelist

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets this context as the model's current I/O context
  !>
  subroutine set_current( this )

    implicit none

    class(lfric_xios_context_type), intent(inout) :: this

    call xios_set_current_context( this%handle )

  end subroutine set_current

  subroutine tick_context_clock(this)
    implicit none
    class(lfric_xios_context_type), intent(inout) :: this

    this%context_clock_step = this%context_clock_step + 1_i_def

  end subroutine tick_context_clock

  function get_context_clock_step(this) result(step)
    implicit none
    class(lfric_xios_context_type), intent(inout) :: this
    integer(i_def) :: step

    step = this%context_clock_step

  end function get_context_clock_step

end module lfric_xios_context_mod
