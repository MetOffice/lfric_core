!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

! Test for reading time data from a NetCDF file
!
program lfric_xios_time_read_test

    use io_context_mod,         only: callback_clock_arg
    use lfric_xios_context_mod, only: lfric_xios_context_type
    use lfric_xios_driver_mod,  only: lfric_xios_initialise, lfric_xios_finalise
    use lfric_xios_utils_mod,   only: read_time_data
    use log_mod,                only: log_event, log_level_error, log_scratch_space
    use test_db_mod,            only: test_db_type
    use xios,                   only: xios_date, operator(/=)

    implicit none

    type(test_db_type)                         :: test_db
    type(lfric_xios_context_type), allocatable :: io_context
    type(xios_date), allocatable :: result(:), check(:)
    integer :: t

    procedure(callback_clock_arg), pointer :: before_close => null()

    call test_db%initialise()
    call lfric_xios_initialise( "test", test_db%comm, .false. )

    allocate(io_context)

    call io_context%initialise( "test_io_context", 1, 10 )
    call io_context%initialise_xios_context( test_db%comm,                    &
                                             test_db%chi,  test_db%panel_id,  &
                                             test_db%clock, test_db%calendar, &
                                             before_close )

    allocate(check(10))
    check = [ xios_date(2024, 1, 1, 15, 1, 0), &
              xios_date(2024, 1, 1, 15, 2, 0), &
              xios_date(2024, 1, 1, 15, 3, 0), &
              xios_date(2024, 1, 1, 15, 4, 0), &
              xios_date(2024, 1, 1, 15, 5, 0), &
              xios_date(2024, 1, 1, 15, 6, 0), &
              xios_date(2024, 1, 1, 15, 7, 0), &
              xios_date(2024, 1, 1, 15, 8, 0), &
              xios_date(2024, 1, 1, 15, 9, 0), &
              xios_date(2024, 1, 1, 15, 10, 0) ]
    result = read_time_data("lfric_xios_time_read_data")

    do t = 1, size(result)
        if (result(t) /= check(t)) then
            print*, "Expected time:", check(t)
            print*, "but got:      ", result(t)
            call log_event("Time data read in incorrectly", log_level_error)
        end if
    end do

    deallocate(io_context)
    deallocate(result)
    deallocate(check)

    call lfric_xios_finalise()
    call test_db%finalise()

  end program lfric_xios_time_read_test
