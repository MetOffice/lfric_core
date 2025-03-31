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

  use constants_mod,         only: i_def, l_def, str_def
  use driver_modeldb_mod,    only: modeldb_type
  use halo_comms_mod,        only: initialise_halo_comms, &
                                   finalise_halo_comms
  use mpi_mod,               only: create_comm, destroy_comm, &
                                   lfric_comm_type

#ifdef MCT
  use coupling_mod,          only: coupling_type, &
                                   get_coupling_from_collection
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
  !> @param[in]     program_name  The model name
  !> @param[in,out] modeldb       The structure that holds model state
  !> @param[in]     input_comm    An optional argument that can be supplied if
  !>                              mpi has been initialised outside the model.
  !>                              In that case, this provides the communicator
  !>                              that should be used
  subroutine init_comm( program_name, modeldb, input_comm )

    implicit none

    character(len=*),                intent(in)    :: program_name
    class(modeldb_type),             intent(inout) :: modeldb
    type(lfric_comm_type), optional, intent(in)    :: input_comm

    type(lfric_comm_type) :: start_communicator
    type(lfric_comm_type) :: model_communicator

    logical :: comm_is_split

#ifdef MCT
    type(coupling_type)          :: coupling
    type(coupling_type), pointer :: coupling_ptr
    character(str_def),  pointer :: cpl_name
#endif

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
    ! Add a place to store the coupling object in modeldb
    call modeldb%values%add_key_value('coupling', coupling)
    ! Extract the version that was actaully placed in the collection
    coupling_ptr => get_coupling_from_collection(modeldb%values, "coupling" )
    ! Get the coupling component name
    call modeldb%values%get_value("cpl_name", cpl_name)
    ! Initialise OASIS coupling and get back the split communicator
    call coupling_ptr%initialise( trim(cpl_name),     &
                                  model_communicator, &
                                  start_communicator, &
                                  comm_is_split )

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
    call modeldb%mpi%initialise( model_communicator )

    ! Initialise halo functionality
    call initialise_halo_comms( model_communicator )

  end subroutine init_comm

  !> @brief  Finalises the model communicator
  !> @param[in,out] modeldb       The structure that holds model state
  subroutine final_comm(modeldb)

    implicit none

    class(modeldb_type), intent(inout) :: modeldb

#ifdef MCT
    type(coupling_type), pointer :: coupling_ptr
#endif

#ifdef USE_XIOS
    ! Finalise XIOS
    call lfric_xios_finalise()
#endif

#ifdef MCT
    ! Extract the coupling object from the modeldb key-value pair collection
    coupling_ptr => get_coupling_from_collection(modeldb%values, "coupling" )
    ! FInalise OASIS coupling
    call coupling_ptr%finalise()
#endif

    ! Finalise halo exchange functionality
    call finalise_halo_comms()

    ! Finalise the mpi object
    call modeldb%mpi%finalise()
    ! Release the communicator if it is ours to release. If a communicator has
    ! been provided to LFRic, then that is someone else's responsibility
    if(comm_created)call destroy_comm()

  end subroutine final_comm

end module driver_comm_mod
