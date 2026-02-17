!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

! Tests the LFRic-XIOS temporal reading functionality using iodef file configuration.
! Correct behaviour is to read only the minimal required time-entries from
! input file at the correct times. The validity of the data written from this
! test is checked against the input data in the python part of the test.
program lfric_xios_temporal_iodef_test

  use constants_mod,          only: r_def
  use event_mod,              only: event_action
  use event_actor_mod,        only: event_actor_type
  use field_mod,              only: field_type, field_proxy_type
  use file_mod,               only: FILE_MODE_READ, FILE_MODE_WRITE
  use io_context_mod,         only: callback_clock_arg
  use lfric_xios_action_mod,  only: advance
  use lfric_xios_context_mod, only: lfric_xios_context_type
  use lfric_xios_driver_mod,  only: lfric_xios_initialise, lfric_xios_finalise
  use lfric_xios_file_mod,    only: lfric_xios_file_type, OPERATION_TIMESERIES
  use linked_list_mod,        only: linked_list_type
  use log_mod,                only: log_event, log_level_info
  use test_db_mod,            only: test_db_type

  implicit none

  type(test_db_type)                                 :: test_db
  type(lfric_xios_context_type), target, allocatable :: io_context

  procedure(callback_clock_arg), pointer :: before_close
  type(linked_list_type),        pointer :: file_list
  class(event_actor_type),       pointer :: context_actor
  procedure(event_action),       pointer :: context_advance
  type(field_type),              pointer :: rfield
  type(field_proxy_type)                 :: rproxy

  call test_db%initialise()
  call lfric_xios_initialise( "test", test_db%comm, .false. )

  ! =============================== Start test ================================

  allocate(io_context)
  call io_context%initialise( "test_io_context", 1, 10 )

  file_list => io_context%get_filelist()
  call file_list%insert_item( lfric_xios_file_type( "lfric_xios_temporal_input",         &
                                                    xios_id="lfric_xios_temporal_input", &
                                                    io_mode=FILE_MODE_READ,              &
                                                    operation=OPERATION_TIMESERIES,      &
                                                    fields_in_file=test_db%temporal_fields ) )
  call file_list%insert_item( lfric_xios_file_type( "lfric_xios_temporal_output",         &
                                                    xios_id="lfric_xios_temporal_output", &
                                                    io_mode=FILE_MODE_WRITE,              &
                                                    operation=OPERATION_TIMESERIES,       &
                                                    freq=1,                               &
                                                    fields_in_file=test_db%temporal_fields ) )

  before_close => null()
  call io_context%initialise_xios_context( test_db%comm,                    &
                                           test_db%chi,  test_db%panel_id,  &
                                           test_db%clock, test_db%calendar, &
                                           before_close )


  context_advance => advance
  context_actor => io_context
  call test_db%clock%add_event( context_advance, context_actor )
  call io_context%set_active(.true.)

  do while (test_db%clock%tick())
    call test_db%temporal_fields%get_field("temporal_field", rfield)
    rproxy = rfield%get_proxy()
    call log_event("Valid data for this TS:", log_level_info)
    print*,rproxy%data(1)
  end do

  deallocate(io_context)

  ! ============================== Finish test =================================

  call lfric_xios_finalise()
  call test_db%finalise()

end program lfric_xios_temporal_iodef_test
