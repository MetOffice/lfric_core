!-----------------------------------------------------------------------------
! (C) Crown copyright 2024-2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

! Simple test which initialises and destroys an XIOS context
!
program lfric_xios_context_test

  use field_mod,              only: field_type
  use io_context_mod,         only: callback_clock_arg
  use lfric_xios_context_mod, only: lfric_xios_context_type
  use lfric_xios_driver_mod,  only: lfric_xios_initialise, lfric_xios_finalise
  use log_mod,                only: log_event, log_level_info
  use test_db_mod,            only: test_db_type

  use local_mesh_mod, only: local_mesh_type
  use mesh_mod, only: mesh_type

  implicit none

  type(test_db_type)                         :: test_db
  type(lfric_xios_context_type), allocatable :: io_context

  procedure(callback_clock_arg), pointer :: before_close => null()

  call test_db%initialise()
  call lfric_xios_initialise( "test", test_db%comm, .false. )

  ! =============================== Start test ================================

  allocate(io_context)
  call io_context%initialise( "test_io_context", 1, 10 )
  call io_context%initialise_xios_context( test_db%comm,                    &
                                           test_db%chi,  test_db%panel_id,  &
                                           test_db%clock, test_db%calendar, &
                                           before_close )
  deallocate(io_context)

  ! ============================== Finish test =================================

  call lfric_xios_finalise()
  call test_db%finalise()

end program lfric_xios_context_test
