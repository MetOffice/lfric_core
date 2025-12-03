!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Handles the loading of namelists.
!
module content_mod

  use constants_mod, only : i_def, l_def, str_def, str_max_filename
  use lfric_mpi_mod, only : global_mpi
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR

  use namelist_collection_mod, only: namelist_collection_type
  use namelist_mod,            only: namelist_type
  use config_mod,              only: config_type

  use foo_nml_mod, only: foo_nml_type

  use foo_config_mod, only : read_foo_namelist, &
                             postprocess_foo_namelist, &
                             foo_is_loadable, &
                             foo_is_loaded, &
                             foo_reset_load_status, &
                             foo_final, &
                             get_foo_nml, &
                             get_new_foo_nml

  implicit none

  private
  public :: read_configuration, ensure_configuration, final_configuration

contains

  ! Reads configuration namelists from a file.
  !
  ! [in] filename File holding the namelists.
  !
  ! TODO: Support "namelist file" namelists which recursively call this
  !       procedure to load other namelist files.
  !
  subroutine read_configuration( filename, configuration, config )

    use io_utility_mod, only : open_file, close_file

    implicit none

    character(*), intent(in) :: filename

    type(namelist_collection_type), optional, intent(inout) :: configuration
    type(config_type),              optional, intent(inout) :: config

    integer(i_def) :: local_rank

    character(str_def), allocatable :: namelists(:)
    integer(i_def) :: unit

    if (.not. present(configuration) .and. .not. present(config)) then
      write(log_scratch_space,'(A)')                               &
          'At least one optional argument must be provided for ' //&
          'read_configuration.'
      call log_event(log_scratch_space, log_level_error)
    end if

    local_rank = global_mpi%get_comm_rank()

    unit = -1
    if (local_rank == 0) unit = open_file( filename )

    call get_namelist_names( unit, local_rank, namelists )

    if (present(configuration) .and. present(config)) then
      ! TODO Transition, remove when all code ported to config
      ! access pattern
      call read_configuration_namelists( unit, local_rank,       &
                                         namelists, filename,    &
                                         nml_bank=configuration, &
                                         config=config )

    else if (present(configuration) .and. .not. present(config)) then
      ! TODO Deprecated, remove when all code ported to config
      ! access pattern
      call read_configuration_namelists( unit, local_rank,    &
                                         namelists, filename, &
                                         nml_bank=configuration )

    else if (.not. present(configuration) .and. present(config)) then
      call read_configuration_namelists( unit, local_rank,    &
                                         namelists, filename, &
                                         config=config )

    end if

    if (local_rank == 0) call close_file( unit )

  end subroutine read_configuration

  ! Finds names of all namelists present in file.
  !
  ! [in] unit File holding namelists.
  ! [out] names of namelist in file (in order).
  !
  subroutine get_namelist_names( unit, local_rank, names )

    use io_utility_mod, only : read_line

    implicit none

    integer(i_def),     intent(in)                 :: unit
    integer(i_def),     intent(in)                 :: local_rank
    character(str_def), intent(inout), allocatable :: names(:)

    character(str_def), allocatable :: names_temp(:)
    ! TODO: Buffer is large enough for a fair sized string and a filename.
    !       Ideally it should be dynamically sized for the length of the
    !       incoming data but I'm not sure how best to achieve that at the
    !       moment. #1752
    character(str_def + str_max_filename) :: buffer
    logical(l_def)     :: continue_read
    ! Number of names
    integer(i_def)  :: namecount

    namecount = 0
    if (local_rank == 0) then
      text_line_loop: do

        continue_read = read_line( unit, buffer )
        if ( .not. continue_read ) exit text_line_loop

        buffer = adjustl(buffer)
        if (buffer(1:1) == '&') then
          namecount = namecount + 1
          allocate(names_temp(namecount))
          if (namecount > 1) then
            names_temp(1:namecount-1) = names
          end if
          names_temp(namecount) = trim(buffer(2:))
          call move_alloc(names_temp, names)
        end if
      end do text_line_loop
      rewind(unit)
    end if

    call global_mpi%broadcast( namecount, 0 )

    if (local_rank /= 0) then
      allocate(names(namecount))
    end if

    call global_mpi%broadcast( names, namecount*str_def, 0 )

  end subroutine get_namelist_names

  ! Checks that the requested namelists have been loaded.
  !
  ! [in]  names List of namelists.
  ! [out] success_mask Marks corresponding namelists as having failed.
  !
  ! [return] Overall success.
  !
  function ensure_configuration( names, success_mask )

    implicit none

    character(*),             intent(in)  :: names(:)
    logical(l_def), optional, intent(out) :: success_mask(:)
    logical(l_def)                        :: ensure_configuration

    integer(i_def)    :: i
    logical           :: configuration_found = .True.

    if (present(success_mask) &
        .and. (size(success_mask, 1) /= size(names, 1))) then
      call log_event( 'Arguments "names" and "success_mask" to function' &
                      // '"ensure_configuration" are different shapes',  &
                      LOG_LEVEL_ERROR )
    end if

    ensure_configuration = .True.

    name_loop: do i = 1, size(names)
      select case(trim( names(i) ))
      case ('foo')
        configuration_found = foo_is_loaded()

      case default
        write( log_scratch_space, '(A)' )               &
            'Tried to ensure unrecognised namelist "'// &
            trim(names(i))//'" was loaded.'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select

      ensure_configuration = ensure_configuration .and. configuration_found

      if (present(success_mask)) success_mask(i) = configuration_found

    end do name_loop

  end function ensure_configuration

  subroutine read_configuration_namelists( unit, local_rank,    &
                                           namelists, filename, &
                                           nml_bank, config )
    implicit none

    integer(i_def),     intent(in) :: unit
    integer(i_def),     intent(in) :: local_rank
    character(str_def), intent(in) :: namelists(:)
    character(*),       intent(in) :: filename

    type(namelist_collection_type), optional, intent(inout) :: nml_bank
    type(config_type),              optional, intent(inout) :: config

    type(namelist_type) :: nml_obj
    type(foo_nml_type) :: foo_nml_obj

    integer(i_def) :: i, j

    logical :: scan

    ! Reset load status from any previous file reads
    call foo_reset_load_status()

    ! Read the namelists
    do j=1, 2

      select case(j)
      case(1)
        scan = .true.
      case(2)
        scan = .false.
      end select

      do i=1, size(namelists)

        select case (trim(namelists(i)))

        case ('foo')
          if (foo_is_loadable()) then
            call read_foo_namelist( unit, local_rank, scan )
            if (.not. scan) then
              call postprocess_foo_namelist()

              if (present(nml_bank)) then
                nml_obj = get_foo_nml()
                call nml_bank%add_namelist(nml_obj)
              end if

              if (present(config)) then
                foo_nml_obj = get_new_foo_nml()
                call config%add_namelist(foo_nml_obj)
              end if

            end if
          else
            write( log_scratch_space, '(A)' )      &
                'Namelist "'//trim(namelists(i))// &
                '" can not be read. Too many instances?'
            call log_event( log_scratch_space, LOG_LEVEL_ERROR )
          end if

        case default
          write( log_scratch_space, '(A)' )                   &
              'Unrecognised namelist "'//trim(namelists(i))// &
              '" found in file '//trim(filename)//'.'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select

      end do ! Namelists

      if ( local_rank == 0 ) then
        rewind( unit )
      end if

    end do ! Reading passes

  end subroutine read_configuration_namelists

  subroutine final_configuration()

    implicit none

    call foo_final()

    return
  end subroutine final_configuration

end module content_mod
