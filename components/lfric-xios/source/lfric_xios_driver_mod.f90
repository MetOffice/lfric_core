!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief    Module containing routines needed to drive XIOS executable
!>
module lfric_xios_driver_mod

  use constants_mod, only: i_def
  use mod_wait,      only: init_wait
  use xios,          only: xios_initialize, xios_finalize

  implicit none

  public :: lfric_xios_initialise, lfric_xios_finalise

contains

  subroutine lfric_xios_initialise( model_name,         &
                                    model_communicator, &
                                    comm_has_been_split )

    implicit none

    character(len=*), intent(in)    :: model_name
    integer(i_def),   intent(inout) :: model_communicator
    logical,          intent(in)    :: comm_has_been_split

    if (comm_has_been_split) then
      call xios_initialize( model_name, local_comm=model_communicator )
    else
      call init_wait()
      call xios_initialize( model_name, return_comm=model_communicator )
    end if

  end subroutine lfric_xios_initialise

  subroutine lfric_xios_finalise()

    implicit none

    call xios_finalize()

  end subroutine lfric_xios_finalise

end module lfric_xios_driver_mod