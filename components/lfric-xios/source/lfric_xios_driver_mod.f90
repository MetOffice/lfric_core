!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief    Module containing routines needed to drive XIOS executable
!>
module lfric_xios_driver_mod

  use constants_mod, only: i_def
  use lfric_mpi_mod, only: lfric_comm_type
  use mod_wait,      only: init_wait
  use xios,          only: xios_initialize, xios_finalize

  implicit none

  public :: lfric_xios_initialise, lfric_xios_finalise

contains

  subroutine lfric_xios_initialise( model_name,         &
                                    model_communicator, &
                                    comm_has_been_split )

    implicit none

    character(len=*),      intent(in)    :: model_name
    type(lfric_comm_type), intent(inout) :: model_communicator
    logical,               intent(in)    :: comm_has_been_split

    integer(i_def) :: comm

    if (comm_has_been_split) then
      ! The model hame has been hardwired as lfric in this case, to match
      ! the conventions used throughout the Met Office coupled model
      ! configurations
      call xios_initialize( "lfric", &
                            local_comm=model_communicator%get_comm_mpi_val() )
    else
      call init_wait()
      call xios_initialize( model_name, return_comm=comm )
      call model_communicator%set_comm_mpi_val(comm)
    end if

  end subroutine lfric_xios_initialise

  subroutine lfric_xios_finalise()

    implicit none

    call xios_finalize()

  end subroutine lfric_xios_finalise

end module lfric_xios_driver_mod
