!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Container for the working data set of the IO_Dev model run, including
!> methods to initialise and finalise the data set
!>
!> This module provides a type to hold all the model fields and methods to
!> initialise (create and read) and finalise (write and destroy) the
!> data contained within the type.
!>
module io_dev_data_mod

  ! Infrastructure
  use clock_mod,                        only : clock_type
  use constants_mod,                    only : i_def, l_def
  use field_mod,                        only : field_type
  use field_parent_mod,                 only : write_interface
  use field_collection_mod,             only : field_collection_type
  use log_mod,                          only : log_event,      &
                                               LOG_LEVEL_INFO, &
                                               LOG_LEVEL_ERROR
  ! Configuration
  use io_config_mod,                    only : checkpoint_read,  &
                                               checkpoint_write, &
                                               write_dump
  ! I/O methods
  use read_methods_mod,                 only : read_state
  use write_methods_mod,                only : write_state
  ! IO_Dev driver modules
  use io_dev_init_mod,                  only : create_io_dev_fields, &
                                               init_io_dev_fields


  implicit none

  private

  !> @brief Holds the working data set for an IO_Dev model run.
  !>
  type :: io_dev_data_type

    private

    !> Stores all the fields used by the model - the remaining collections store
    !> pointers to the data in the depository.
    type( field_collection_type ), public :: depository

    !> @name Collections controlling input/output of fields
    !> @{

    !> Field collection holding fields to be read into
    type( field_collection_type ), public   :: input_fields
    !> Field collection holding fields to be written
    type( field_collection_type ), public   :: output_fields

    !> @}

  end type io_dev_data_type

  public io_dev_data_type, create_model_data, finalise_model_data, &
         initialise_model_data, output_model_data

contains

  !> @brief Create the fields contained in model_data
  !> @param[inout] model_data   The working data set for a model run
  !> @param[in]    mesh_id      The identifier given to the current 3d mesh
  !> @param[in]    twod_mesh_id The identifier given to the current 2d mesh
  !> @param[in]    clock        The model clock object
  subroutine create_model_data( model_data, &
                                mesh_id,    &
                                twod_mesh_id, &
                                clock )

    implicit none

    type( io_dev_data_type ), intent(inout) :: model_data
    integer( i_def ),         intent(in)    :: mesh_id
    integer( i_def ),         intent(in)    :: twod_mesh_id
    ! Clock will be integrated with future I/O testing functionality
    class( clock_type ),      intent(in)    :: clock

    ! Create real fields for input and output
    call create_io_dev_fields( mesh_id,                 &
                               twod_mesh_id,            &
                               model_data%depository,   &
                               model_data%input_fields, &
                               model_data%output_fields )

  end subroutine create_model_data

  !> @brief Initialises the working data set dependent of namelist configuration
  !> @param[inout] model_data The working data set for a model run
  !> @param[in]    clock      The model clock object
  subroutine initialise_model_data( model_data, clock )

    implicit none

    type( io_dev_data_type ), intent(inout) :: model_data
    ! Clock will be integrated with future I/O testing functionality
    class( clock_type ),      intent(in)    :: clock

    ! Initialise all the model fields here.

    ! Call to read routines depending on model configuration - read from ugrid
    ! file, restart, dump etc.

    !call read_state( model_data%input_fields )

    call init_io_dev_fields( model_data%input_fields, &
                             model_data%output_fields )


  end subroutine initialise_model_data

  !> @brief Writes out a checkpoint and dump file dependent on namelist options
  !> @param[inout] model_data The working data set for the model run
  !> @param[in]    clock      The model clock object.
  subroutine output_model_data( model_data, &
                                clock )

    implicit none

    type( io_dev_data_type ), intent(inout), target :: model_data
    ! Clock will be integrated with future I/O testing functionality
    class( clock_type ),      intent(in)            :: clock


    !===================== Write fields to dump ======================!
    ! Functionality to be added

    !=================== Write fields to checkpoint files ====================!
    ! Functionality to be added

    !=================== Write fields to diagnostic files ====================!
    call write_state( model_data%output_fields )


  end subroutine output_model_data

  !> @brief Routine to destroy all the field collections in the working data set
  !> @param[inout] model_data The working data set for a model run
  subroutine finalise_model_data( model_data )

    implicit none

      type(io_dev_data_type), intent(inout) :: model_data

      ! Clear all the fields in each field collection
      call model_data%depository%clear()
      call model_data%input_fields%clear()
      call model_data%output_fields%clear()

      call log_event( 'finalise_model_data: all fields have been cleared', &
                       LOG_LEVEL_INFO )

  end subroutine finalise_model_data

end module io_dev_data_mod
