!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the test namelist.
!>
module test_config_mod

  use constants_mod, only: i_def
  use lfric_mpi_mod, only: global_mpi
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_DEBUG, LOG_LEVEL_INFO

  use namelist_mod,      only: namelist_type
  use namelist_item_mod, only: namelist_item_type

  use constants_mod, only: cmdi, emdi, imdi, rmdi, str_def, unset_key

  implicit none

  private
  public :: read_test_namelist, postprocess_test_namelist, &
            test_is_loadable, test_is_loaded, &
            test_reset_load_status, &
            test_multiples_allowed, test_final, &
            get_test_nml, get_new_test_nml

  integer(i_def), public, protected :: foo = imdi

  character(*), parameter :: listname = 'test'
  character(str_def) :: profile_name = cmdi

  logical, parameter :: multiples_allowed = .false.

  logical :: nml_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !> @param [in] scan .true. if reading namelist to acquire scalar
  !>                  values which may possbly be required for
  !>                  array sizing during postprocessing.
  !>
  subroutine read_test_namelist( file_unit, local_rank, scan )

    use constants_mod, only: i_def

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan

    call read_namelist( file_unit, local_rank, scan )

  end subroutine read_test_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, scan )

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan

    integer(i_def) :: buffer_integer_i_def(1)

    namelist /test/ foo

    integer(i_def) :: condition

    foo = imdi

    if (local_rank == 0) then

      read( file_unit, nml=test, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_integer_i_def(1) = foo

    call global_mpi%broadcast( buffer_integer_i_def, 1, 0 )

    foo = buffer_integer_i_def(1)

    if (scan) then
      nml_loaded = .false.
    else
      nml_loaded = .true.
    end if

  end subroutine read_namelist


  !> @brief Returns a <<namelist_type>> object populated with the
  !>        current contents of this configuration module.
  !> @return namelist_obj <<namelist_type>> with current namelist contents.
  function get_test_nml() result(namelist_obj)

    implicit none

    type(namelist_type)      :: namelist_obj
    type(namelist_item_type) :: members(1)

    call members(1)%initialise( &
                'foo', foo )

    if (trim(profile_name) /= trim(cmdi) ) then
      call namelist_obj%initialise( trim(listname), &
                                    members, &
                                    profile_name = profile_name )
    else
      call namelist_obj%initialise( trim(listname), &
                                    members )
    end if

  end function get_test_nml

  !> @brief Returns a <<test_nml_type>> object populated with the
  !>        current contents of this configuration module.
  !> @return namelist_obj <<test_nml_type>> with current namelist contents.
  function get_new_test_nml() result(namelist_obj)

    use test_nml_mod, only: test_nml_type

    implicit none

    type(test_nml_type) :: namelist_obj
    type(namelist_item_type) :: members(1)

    call members(1)%initialise( &
                'foo', foo )

    if (trim(profile_name) /= trim(cmdi) ) then
      call namelist_obj%initialise( trim(listname), &
                                    members, &
                                    profile_name = profile_name )
    else
      call namelist_obj%initialise( trim(listname), &
                                    members )
    end if

  end function get_new_test_nml


  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_test_namelist()

    use constants_mod, only: i_def

    implicit none

    ! Computed fields are resolved after everything has been loaded since they
    ! can refer to fields in other namelists.
    !

  end subroutine postprocess_test_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function test_is_loadable()

    implicit none

    logical :: test_is_loadable

    if ( multiples_allowed .or. .not. nml_loaded ) then
      test_is_loadable = .true.
    else
      test_is_loadable = .false.
    end if

  end function test_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function test_is_loaded()

    implicit none

    logical :: test_is_loaded

    test_is_loaded = nml_loaded

  end function test_is_loaded

  !> Are multiple test namelists allowed to be read?
  !>
  !> @return True If multiple test namelists are
  !>              permitted.
  !>
  function test_multiples_allowed()

    implicit none

    logical :: test_multiples_allowed

    test_multiples_allowed = multiples_allowed

  end function test_multiples_allowed

  !> Resets the load status to allow
  !> test namelist to be read.
  !>
  subroutine test_reset_load_status()

    implicit none

    nml_loaded = .false.

  end subroutine test_reset_load_status

  !> Clear out any allocated memory
  !>
  subroutine test_final()

    implicit none

    foo = imdi

    return
  end subroutine test_final


end module test_config_mod
