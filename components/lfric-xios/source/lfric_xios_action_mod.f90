!-------------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
module lfric_xios_action_mod

  use constants_mod, only : str_def
  use timing_mod,    only : start_timing, stop_timing, tik, LPROF

  implicit none

  private

  public :: advance
  public :: advance_read_only

contains

  !> Advances the XIOS context forward in time, performing all I/O operations
  !> expected by XIOS at the end and beginning of the current and subsequent
  !> timesteps.
  !>
  !> @param[in] context     The IO context to be advanced
  !> @param[in] model_clock The model's clock
  subroutine advance(context, model_clock)
    use clock_mod,              only : clock_type
    use event_actor_mod,        only : event_actor_type
    use linked_list_mod,        only : linked_list_item_type, linked_list_type
    use lfric_xios_file_mod,    only : lfric_xios_file_type
    use lfric_xios_context_mod, only : lfric_xios_context_type
    use log_mod,                only : log_event, log_level_error, &
                                       log_level_info

    !> TODO Remove icontext, see ticket #4313.
    !> Use icontext is needed here as the revision of xios used by lfric_coupled
    !> is old enough to not have xios_get_current_context forwarded through the
    !> xios module.
    use icontext,             only : xios_get_current_context
    use xios,                 only : xios_context,                  &
                                     xios_set_current_context,      &
                                     xios_update_calendar

    implicit none

    class(event_actor_type), intent(inout) :: context
    class(clock_type),       intent(in)    :: model_clock

    type(linked_list_item_type), pointer :: loop => null()
    type(lfric_xios_file_type),  pointer :: file => null()
    type(xios_context)                   :: xios_context_handle
    type(linked_list_type), pointer      :: filelist
    integer(tik)                         :: id
    logical                              :: profiling

    ! Get the handle of the current context (Not necessarily the one passed to this routine).
    ! This is used to reset the context on return.
    call xios_get_current_context(xios_context_handle)

    select type(context)
    type is (lfric_xios_context_type)
      ! Write all files that need to be written to
      call context%set_current()

      if (model_clock%get_step() > model_clock%get_first_step()) then
        filelist => context%get_filelist()
        if (filelist%get_length() > 0) then
          loop => filelist%get_head()
          do while (associated(loop))
            select type(list_item => loop%payload)
            type is (lfric_xios_file_type)
              file => list_item
              if (file%mode_is_write()) call file%send_fields()
            end select
            loop => loop%next
          end do
        end if
      end if

      ! Update XIOS calendar
      profiling = (context%get_timer_flag() .and. LPROF )
      if ( profiling ) call start_timing( id, 'xios_update_calendar' )
      call xios_update_calendar( model_clock%get_step() - model_clock%get_first_step() + 1 )
      if ( profiling ) call stop_timing( id, 'xios_update_calendar' )

      ! Read all files that need to be read from
      filelist => context%get_filelist()
      if (filelist%get_length() > 0) then
        loop => filelist%get_head()
        do while (associated(loop))
          select type(list_item => loop%payload)
          type is (lfric_xios_file_type)
            file => list_item
            if (file%mode_is_read()) call file%recv_fields()
          end select
          loop => loop%next
        end do
      end if

    class default
      call log_event( "Can not advance a non lfric xios type context.", &
        log_level_error )
    end select

    nullify(loop)
    nullify(file)

    ! Reset the xios context to what it was before this subroutine was called.
    call xios_set_current_context(xios_context_handle)

  end subroutine advance

  !> Advances the XIOS context forward in time, performing only read operations
  !> with XIOS at the beginning of the current and timesteps.
  !>
  !> @param[in] context     The IO context to be advanced
  !> @param[in] model_clock The model's clock
  subroutine advance_read_only(context, model_clock)
    use clock_mod,              only : clock_type
    use event_actor_mod,        only : event_actor_type
    use linked_list_mod,        only : linked_list_item_type, linked_list_type
    use lfric_xios_file_mod,    only : lfric_xios_file_type
    use lfric_xios_context_mod, only : lfric_xios_context_type
    use log_mod,                only : log_event, log_level_error, &
                                       log_level_info, log_level_debug

    !> TODO Remove icontext, see ticket #4313.
    !> Use icontext is needed here as the revision of xios used by lfric_coupled
    !> is old enough to not have xios_get_current_context forwarded through the
    !> xios module.
    use icontext,             only : xios_get_current_context
    use xios,                 only : xios_context,                  &
                                     xios_date,                     &
                                     xios_set_current_context,      &
                                     xios_update_calendar,          &
                                     xios_get_current_date,         &
                                     xios_date_convert_to_string

    implicit none

    class(event_actor_type), intent(inout) :: context
    class(clock_type),       intent(in)    :: model_clock

    type(linked_list_item_type), pointer :: loop => null()
    type(lfric_xios_file_type),  pointer :: file => null()
    type(xios_context)                   :: xios_context_handle
    type(linked_list_type), pointer      :: filelist
    integer(tik)                         :: id
    logical                              :: profiling

    ! Get the handle of the current context (Not necessarily the one passed to this routine).
    ! This is used to reset the context on return.
    call xios_get_current_context(xios_context_handle)

    select type(context)
    type is (lfric_xios_context_type)
      call context%set_current()
      call context%tick_context_clock()
      ! Update XIOS calendar
      profiling = ( context%get_timer_flag() .and. LPROF )
      if ( profiling ) call start_timing( id, 'xios_update_calendar' )
      call xios_update_calendar( context%get_context_clock_step() )
      if ( profiling ) call stop_timing( id, 'xios_update_calendar' )

      ! Read all files that need to be read from
      filelist => context%get_filelist()
      if (filelist%get_length() > 0) then
        loop => filelist%get_head()
        do while (associated(loop))
          select type(list_item => loop%payload)
          type is (lfric_xios_file_type)
            file => list_item
            if (file%mode_is_read()) call file%recv_fields()
          end select
          loop => loop%next
        end do
      end if

    class default
      call log_event( "Can not advance a non lfric xios type context.", &
        log_level_error )
    end select

    nullify(loop)
    nullify(file)

    ! Reset the xios context to what it was before this subroutine was called.
    call xios_set_current_context(xios_context_handle)

  end subroutine advance_read_only

end module lfric_xios_action_mod
