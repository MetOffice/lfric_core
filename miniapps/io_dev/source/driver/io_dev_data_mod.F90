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
  use constants_mod,                    only : i_def
  use field_mod,                        only : field_type
  use field_collection_mod,             only : field_collection_type
  use linked_list_mod,                  only : linked_list_type
  use log_mod,                          only : log_event,      &
                                               LOG_LEVEL_INFO, &
                                               LOG_LEVEL_ERROR
  use variable_fields_mod,              only : init_variable_fields
  ! Configuration
  use io_config_mod,                    only : write_diag,                &
                                               write_dump
  use initialization_config_mod,        only : init_option,               &
                                               init_option_fd_start_dump, &
                                               ancil_option,              &
                                               ancil_option_basic_gal
  ! I/O methods
  use lfric_xios_read_mod,              only : read_state
  use lfric_xios_write_mod,             only : write_state
  ! IO_Dev modules
  use io_dev_init_mod,                  only : setup_io_dev_fields
  use io_dev_init_fields_alg_mod,       only : io_dev_init_fields_alg
  use io_dev_checksum_alg_mod,          only : io_dev_checksum_alg

  implicit none

  private

  !> @brief Holds the working data set for an IO_Dev model run.
  !>
  type :: io_dev_data_type

    private

    !> Stores all the fields used by the model - the remaining collections store
    !> pointers to the data in the core_fields.
    type( field_collection_type ), public :: core_fields

    !> Field collection holding fields for dumps
    type( field_collection_type ), public :: dump_fields

    !> Linked list of time_axis objects for variable fields
    type( linked_list_type ),      public :: variable_field_times

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
                                twod_mesh_id )

    implicit none

    type( io_dev_data_type ), intent(inout) :: model_data
    integer( i_def ),         intent(in)    :: mesh_id
    integer( i_def ),         intent(in)    :: twod_mesh_id

    ! Create model data fields
    call setup_io_dev_fields( mesh_id,                        &
                              twod_mesh_id,                   &
                              model_data%core_fields,         &
                              model_data%dump_fields,         &
                              model_data%variable_field_times )

  end subroutine create_model_data

  !> @brief Initialises the working data set dependent of namelist configuration
  !> @param[inout] model_data The working data set for a model run
  !> @param[in]    chi_xyz    A size 3 array of fields holding the mesh
  !>                          (X,Y,Z) coordinates
  !> @param[in]    clock      The model clock object
  subroutine initialise_model_data( model_data, chi_xyz, clock )

    implicit none

    type( io_dev_data_type ), intent(inout) :: model_data
    type( field_type ),       intent(in)    :: chi_xyz(3)
    ! Clock will be integrated with future I/O testing functionality
    class( clock_type ),      intent(in)    :: clock

    ! Initialise all the model fields here analytically - setting data value
    ! equal to product of x, y and z coordinates
    call io_dev_init_fields_alg( model_data%core_fields, chi_xyz )

    !---------------------------------------------------------------
    ! Now we make separate init calls based on model configuration
    !---------------------------------------------------------------

    ! If fields need to be read from dump file, read them
    select case( init_option )
      case ( init_option_fd_start_dump )
        call read_state( model_data%dump_fields, prefix='input_' )
    end select

    ! If testing initialisation of time-varying I/O
    select case( ancil_option )
    case ( ancil_option_basic_gal )
      call init_variable_fields( model_data%variable_field_times, &
                                 clock, model_data%core_fields )
  end select


  end subroutine initialise_model_data

  !> @brief Writes out a checkpoint and dump file dependent on namelist options
  !> @param[inout] model_data The working data set for the model run
  !> @param[in]    clock      The model clock object.
  subroutine output_model_data( model_data )

    implicit none

    type( io_dev_data_type ), intent(inout), target :: model_data

    !===================== Write fields to dump ======================!
    if ( write_dump ) then
      call write_state( model_data%dump_fields, prefix='output_' )
    end if

    !=================== Write fields to checkpoint files ====================!
    ! Functionality to be added

    !=================== Write fields to diagnostic files ====================!
    if ( write_diag ) then
      call write_state( model_data%core_fields )
    end if


  end subroutine output_model_data

  !> @brief Routine to destroy all the field collections in the working data set
  !> @param[inout] model_data The working data set for a model run
  subroutine finalise_model_data( model_data )

    implicit none

      type(io_dev_data_type), intent(inout) :: model_data

      !==================== Write checksum output =====================
      ! Remove multi-data fields from core fields as they cannot currently be passed
      ! to Psyclone to create checksum
      call model_data%core_fields%remove_field( "multi_data_field" )

      call io_dev_checksum_alg( model_data%core_fields )

      ! Clear all the fields in each field collection
      call model_data%core_fields%clear()
      call model_data%dump_fields%clear()

      call log_event( 'finalise_model_data: all fields have been cleared', &
                       LOG_LEVEL_INFO )

  end subroutine finalise_model_data

end module io_dev_data_mod
