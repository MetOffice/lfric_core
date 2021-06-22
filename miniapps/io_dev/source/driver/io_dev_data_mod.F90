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
  use timer_mod,                        only : timer
  use variable_fields_mod,              only : init_variable_fields, &
                                               update_variable_fields
  ! Configuration
  use io_config_mod,                    only : write_diag, write_dump, &
                                               checkpoint_read,        &
                                               checkpoint_write,       &
                                               subroutine_timers
  use io_dev_config_mod,                only : field_initialisation,            &
                                               field_initialisation_start_dump, &
                                               time_variation,                  &
                                               time_variation_analytic,         &
                                               time_variation_ancil,            &
                                               time_variation_none
  ! I/O methods
  use lfric_xios_read_mod,              only : read_state, read_checkpoint
  use lfric_xios_write_mod,             only : write_state, write_checkpoint
  ! IO_Dev modules
  use io_dev_init_mod,                  only : setup_io_dev_fields
  use io_dev_init_fields_alg_mod,       only : io_dev_init_fields_alg
  use io_dev_checksum_alg_mod,          only : io_dev_checksum_alg
  use io_dev_timestep_alg_mod,          only : io_dev_timestep_alg

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

    !> Field collection holding fields which can be processed by PSyClone
    type( field_collection_type ), public :: alg_fields

    !> Linked list of time_axis objects for variable fields
    type( linked_list_type ),      public :: variable_field_times

  end type io_dev_data_type

  public io_dev_data_type, create_model_data, initialise_model_data, &
         update_model_data, finalise_model_data, output_model_data

contains

  !> @brief Create the fields contained in model_data
  !> @param[in,out] model_data   The working data set for a model run
  !> @param[in]     mesh_id      The identifier given to the current 3d mesh
  !> @param[in]     twod_mesh_id The identifier given to the current 2d mesh
  !> @param[in]     clock        The model clock object
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
                              model_data%alg_fields,          &
                              model_data%variable_field_times )

  end subroutine create_model_data

  !> @brief Initialises the working data set dependent of namelist configuration
  !> @param[in,out] model_data The working data set for a model run
  !> @param[in]     chi_sph    A size 3 array of fields holding the mesh coordinates
  !> @param[in]     panel_id   A field with the IDs of mesh panels
  !> @param[in]     clock      The model clock object
  subroutine initialise_model_data( model_data, chi_sph, panel_id, clock )

    implicit none

    type( io_dev_data_type ), intent(inout) :: model_data
    type( field_type ),       intent(in)    :: chi_sph(3)
    type( field_type ),       intent(in)    :: panel_id
    class( clock_type ),      intent(in)    :: clock

    ! Initialise all the model fields here analytically
    call io_dev_init_fields_alg( model_data%core_fields, chi_sph, panel_id )


    !---------------------------------------------------------------
    ! Now we make separate init calls based on model configuration
    !---------------------------------------------------------------
    if ( checkpoint_read ) then
      if ( subroutine_timers ) call timer('read_checkpoint')
      call read_checkpoint( model_data%core_fields, &
                            clock%get_first_step() - 1 )
      if ( subroutine_timers ) call timer('read_checkpoint')

    else
      ! If fields need to be read from dump file, read them
      if ( field_initialisation == field_initialisation_start_dump ) then
          if ( subroutine_timers ) call timer('read_state: dump')
          call read_state( model_data%dump_fields, prefix='input_' )
          if ( subroutine_timers ) call timer('read_state: dump')
      end if

      ! If testing initialisation of time-varying I/O
      if ( time_variation == time_variation_ancil ) then
        if ( subroutine_timers ) call timer('init_variable_fields')
        call init_variable_fields( model_data%variable_field_times, &
                                   clock, model_data%core_fields )
        if ( subroutine_timers ) call timer('init_variable_fields')
      end if

    end if


  end subroutine initialise_model_data

  !> @brief Updates the working data set dependent of namelist configuration
  !> @param[in,out] model_data The working data set for a model run
  !> @param[in]     clock      The model clock object
  subroutine update_model_data( model_data, clock )

    implicit none

    type( io_dev_data_type ), intent(inout) :: model_data
    class( clock_type ),      intent(in)    :: clock

    !---------------------------------------------------------------
    ! Separate update calls are made based on model configuration
    !---------------------------------------------------------------
    select case ( time_variation )

    case ( time_variation_analytic )
      call log_event( "IO_Dev: Updating fields analytically", LOG_LEVEL_INFO )
      if (model_data%alg_fields%get_length() /= 0) then
        call io_dev_timestep_alg( model_data%alg_fields, clock )
      end if

    case ( time_variation_ancil )
      call log_event( "IO_Dev: Updating fields from time_varying ancillary", LOG_LEVEL_INFO )
      if ( subroutine_timers ) call timer('update_variable_fields')
      call update_variable_fields( model_data%variable_field_times, &
                                   clock, model_data%core_fields )
      if ( subroutine_timers ) call timer('update_variable_fields')

    case ( time_variation_none )
      call log_event( "IO_Dev: No time variation for this run", LOG_LEVEL_INFO )

    case default
      call log_event( "IO_Dev: Invalid choice for time-variation namelist", LOG_LEVEL_ERROR )

    end select

  end subroutine update_model_data


  !> @brief Writes out a checkpoint and dump file dependent on namelist options
  !> @param[in,out] model_data The working data set for the model run
  subroutine output_model_data( model_data )

    implicit none

    type( io_dev_data_type ), intent(inout), target :: model_data

    !===================== Write fields to dump ======================!
    if ( write_dump ) then
      if ( subroutine_timers ) call timer('write_state: dump')
      call write_state( model_data%dump_fields, prefix='output_' )
      if ( subroutine_timers ) call timer('write_state: dump')
    end if

    !=================== Write fields to diagnostic files ====================!
    if ( write_diag ) then
      if ( subroutine_timers ) call timer('write_state: diagnostic')
      call write_state( model_data%core_fields )
      if ( subroutine_timers ) call timer('write_state: diagnostic')
    end if


  end subroutine output_model_data

  !> @brief Routine to destroy all the field collections in the working data set
  !> @param[in,out] model_data The working data set for a model run
  !> @param[in]     clock      The model clock object
  subroutine finalise_model_data( model_data, clock )

    implicit none

      type(io_dev_data_type), intent(inout) :: model_data
      class(clock_type),      intent(in)    :: clock

      !=================== Write fields to checkpoint files ====================
      if ( checkpoint_write ) then
        if ( subroutine_timers ) call timer('write_checkpoint')
        call write_checkpoint( model_data%core_fields, clock )
        if ( subroutine_timers ) call timer('write_checkpoint')
      end if

      !======================== Write checksum output ==========================
      if (model_data%alg_fields%get_length() /= 0) then
        call io_dev_checksum_alg( model_data%alg_fields )
      end if

      ! Clear all the fields in each field collection
      call model_data%core_fields%clear()
      call model_data%dump_fields%clear()

      call log_event( 'finalise_model_data: all fields have been cleared', &
                       LOG_LEVEL_INFO )

  end subroutine finalise_model_data

end module io_dev_data_mod
