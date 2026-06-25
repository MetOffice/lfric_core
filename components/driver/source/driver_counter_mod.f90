!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Lifecycle management of the simple counter profiling system.
!>
module driver_counter_mod

  use config_mod,    only: config_type
  use constants_mod, only: str_max_filename, l_def
  use count_mod,     only: count_type, halo_calls
  use timer_mod,     only: timer, output_timer, init_timer

  implicit none

  private
  public :: init_counters, final_counters

contains

  !> Initialises counters from namelists.
  !>
  !> As well as initialising the system a "top level" counter is set up
  !? for tracking halo calls.
  !>
  !> @param[in] config     Application namelist configuration object
  !> @param[in] identifier Top level halo name.
  !>
  subroutine init_counters(config, identifier)

    implicit none

    type(config_type), intent(in) :: config
    character(*),      intent(in) :: identifier

    if ( config%io%subroutine_counters() ) then
      allocate( halo_calls, source=count_type('halo_calls') )
      call halo_calls%counter( identifier )
    end if

  end subroutine init_counters

  !> Shuts down counters.
  !>
  !> The identifier specified when shutting down should be the same as the one
  !> given on initialisation. There is a chance to mismatch the identifiers
  !> which will cause problems but it avoids the use of a global variable.
  !>
  !> @todo Reconsider the existance of the simple counter system once the
  !>       profiler is integrated.
  !>
  !> @param[in] config     Application namelist configuration object
  !> @param[in] identifier Top level counter name.
  !>
  subroutine final_counters(config, identifier)

    implicit none

    type(config_type), intent(in) :: config
    character(*),      intent(in) :: identifier

    if ( config%io%subroutine_counters() ) then
      call halo_calls%counter( identifier )
      call halo_calls%output_counters( config%io%counter_output_suffix() )
    end if

  end subroutine final_counters

end module driver_counter_mod
