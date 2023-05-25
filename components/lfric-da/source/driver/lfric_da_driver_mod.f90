!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief    Module containing routines needed to drive LFRic-JEDI interface
!>
module lfric_da_driver_mod

  use calendar_mod,          only: calendar_type
  use constants_mod,         only: i_native
  use driver_io_mod,         only: init_io_context, get_io_context
  use driver_model_data_mod, only: model_data_type
  use field_mod,             only: field_type
  use io_context_mod,        only: io_context_type
  use model_clock_mod,       only: model_clock_type

  implicit none

  class(io_context_type), allocatable :: da_io_context

  private
  public :: init_da, final_da

contains

  !> @brief Initialises the DA component
  !>
  !> @param[in] communicator The ID of the models MPI communicator
  !> @param[in] chi          The model coordinate field
  !> @param[in] panel_id     The model panel_id field
  !> @param[in] model_clock  The model clock
  !> @param[in] calendar     The model calendar
  subroutine init_da( communicator, chi, panel_id, model_clock, calendar )

    implicit none

    integer(i_native),      intent(in)    :: communicator
    class(field_type),      intent(in)    :: chi(:)
    class(field_type),      intent(in)    :: panel_id
    type(model_clock_type), intent(inout) :: model_clock
    class(calendar_type),   intent(in)    :: calendar

    class(io_context_type), pointer :: model_io_context => null()

    ! Create an I/O context for DA
    call init_io_context( da_io_context, "da_context", communicator, chi, &
                          panel_id, model_clock, calendar )

    ! Set the primary I/O context back to the main model context
    model_io_context => get_io_context()
    call model_io_context%set_current()

  end subroutine init_da

  !> @brief Finalises the DA component
  subroutine final_da()

    implicit none

    if (allocated(da_io_context)) deallocate(da_io_context)

  end subroutine final_da

end module lfric_da_driver_mod