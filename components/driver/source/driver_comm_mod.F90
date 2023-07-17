!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Controls the initialisation and finalisation of the model
!>        communicator

!> @details This controls how model communications are initialised and
!>          finalised. Two modes of operation are supported:
!>           i) The model can be provided with an MPI communicator to run in.
!>          ii) This module can initialise MPI to create its own communicator
!>              to run in.
!>          Ideally, Oasis and XIOS would then sub-divide that communicator
!>          as they see fit, but XIOS doesn't currently support this.
!>
module driver_comm_mod

  use constants_mod,         only: i_native, l_def
  use halo_comms_mod,        only: initialise_halo_comms, &
                                   finalise_halo_comms
  use mpi_mod,               only: mpi_type, create_comm, destroy_comm

! MCT flag used for models coupled to OASIS-MCT ocean model
#ifdef MCT
  use coupler_mod,           only: cpl_initialize, cpl_finalize
#endif
! USE_XIOS flag used for models using the XIOS I/O server
#ifdef USE_XIOS
  use lfric_xios_driver_mod, only: lfric_xios_initialise, lfric_xios_finalise
#endif

  implicit none

  public :: init_comm, final_comm
  private

  ! MPI can only be initialised once per executable, so the following is a
  ! genuinely global variable to describe if MPI has been initialised from here
  logical(l_def) :: comm_created = .false.

contains

  !> @brief  Initialises the model communicator
  !>
  !> @param[in]    program_name  The model name
  !> @param[inout] mpi           The object that holds MPI information
  !> @param[in]    input_comm    An optional argument that can be supplied if
  !>                             mpi has been initialised outside the model.
  !>                             In that case, this provides the communicator
  !>                             that should be used
  subroutine init_comm( program_name, mpi, input_comm )

    implicit none

    character(len=*),                 intent(in)    :: program_name
    type(mpi_type),                   intent(inout) :: mpi
    integer(kind=i_native), optional, intent(in)    :: input_comm

    integer(kind=i_native) :: start_communicator = -999
    integer(kind=i_native) :: model_communicator = -999

    logical :: comm_is_split

    ! Comm has not been split yet
    comm_is_split = .false.

    ! Get the communicator for the whole system (from which
    ! we can start splitting, if we need to)
    if (present(input_comm)) then
      ! Start by using the communicator that we've been given
      start_communicator = input_comm
    else
      ! Initialise mpi and use MPI_COMM_WORLD as the starting communicator
      if(.not. comm_created)then
       call create_comm( start_communicator )
       comm_created = .true.
      endif
    endif

    ! Call the initialisations for Oasis and XIOS as required. These will
    ! spilt the communicator and return a communicator for the model to run in.

#ifdef MCT
    ! Initialise OASIS coupling and get back the split communicator
    call cpl_initialize( model_communicator, start_communicator )
    comm_is_split = .true.
#endif

#ifdef USE_XIOS
    ! Initialise XIOS and get back the split communicator
    ! (At the moment, XIOS2 can only cope with starting from a split
    ! communicator if it has been split by OASIS. In all other cases it just
    ! splits MPI_COMM_WORLD)
    call lfric_xios_initialise( program_name, model_communicator, comm_is_split )
    comm_is_split = .true.
#endif

    ! If neither OASIS nor XIOS has split the communicator, set the model's
    ! communicator to the starting one created (or input) above
    if (.not. comm_is_split) model_communicator = start_communicator

    !Store the MPI communicator for later use
    call mpi%initialise( model_communicator )

    ! Initialise halo functionality
    call initialise_halo_comms( model_communicator )

  end subroutine init_comm

  !> @brief  Finalises the model communicator
  !> @param[inout] mpi The object that holds MPI information
  subroutine final_comm(mpi)

    implicit none

    type(mpi_type), intent(inout) :: mpi

#ifdef USE_XIOS
    ! Finalise XIOS
    call lfric_xios_finalise()
#endif

#ifdef MCT
    ! FInalise OASIS coupling
    call cpl_finalize()
#endif

    ! Finalise halo exchange functionality
    call finalise_halo_comms()

    ! Finalise the mpi object
    call mpi%finalise()
    ! Release the communicator if it is ours to release. If a communicator has
    ! been provided to LFRic, then that is someone else's responsibility
    if(comm_created)call destroy_comm()

  end subroutine final_comm

end module driver_comm_mod
