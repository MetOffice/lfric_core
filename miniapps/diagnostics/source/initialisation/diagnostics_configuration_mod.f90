!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> Diagnostics miniapp program support functions.
!> Mostly generic boiler plate!
!> Of interest:
!>     - required configuration
!>
module diagnostics_configuration_mod

  use log_mod, only : log_event,         &
                      log_scratch_space, &
                      LOG_LEVEL_ALWAYS,  &
                      LOG_LEVEL_ERROR, &
                      LOG_LEVEL_INFO


  implicit none

  private
  public :: load_configuration, program_name

  character(*), parameter :: program_name = "diagnostics"

contains

  !> Loads run-time configuration and ensures everything is ship-shape.
  !>
  subroutine load_configuration( filename )

    use configuration_mod, only : read_configuration, &
                                  ensure_configuration

    implicit none

    character(*), intent(in) :: filename

    character(*), parameter ::                                &
        required_configuration(9) =  [ 'base_mesh          ', & ! global space setup
                                       'extrusion          ', & ! ""
                                       'finite_element     ', & ! ""
                                       'partitioning       ', & ! ""
                                       'planet             ', & ! ""
                                       'io                 ', & ! ""
                                       'logging            ', & ! ""
                                       'time               ', & ! Run Configuration
                                       'diagnostics_miniapp'] !   ""

    logical              :: okay
    logical, allocatable :: success_map(:)
    integer              :: i

    allocate( success_map(size(required_configuration)) )

    call log_event( 'Loading '//program_name//' configuration ...', &
                    LOG_LEVEL_ALWAYS )

    call read_configuration( filename )

    call log_event( program_name//' configuration read', &
                    LOG_LEVEL_INFO )
    okay = ensure_configuration( required_configuration, success_map )
    call log_event( program_name//' configuration validated', &
                    LOG_LEVEL_INFO )
    if (.not. okay) then
      write( log_scratch_space, '(A)' ) &
                             'The following required namelists were not loaded:'
      do i = 1,size(required_configuration)
        if (.not. success_map(i)) &
          log_scratch_space = trim(log_scratch_space) // ' ' &
                              // required_configuration(i)
      end do
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if


    call log_event( program_name//' configuration loaded', &
                    LOG_LEVEL_INFO )
    deallocate( success_map )

  end subroutine load_configuration
end module diagnostics_configuration_mod
