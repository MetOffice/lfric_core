!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the twoenum namelist.
!>
module twoenum_config_mod

  use constants_mod, only: i_def, &
                           str_def
  use lfric_mpi_mod, only: global_mpi
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_DEBUG, LOG_LEVEL_INFO

  use namelist_mod,      only: namelist_type
  use namelist_item_mod, only: namelist_item_type

  use constants_mod, only: cmdi, emdi, imdi, rmdi, str_def, unset_key

  implicit none

  private
  public :: first_from_key, key_from_first, &
            second_from_key, key_from_second, &
            read_twoenum_namelist, postprocess_twoenum_namelist, &
            twoenum_is_loadable, twoenum_is_loaded, &
            twoenum_reset_load_status, &
            twoenum_multiples_allowed, twoenum_final, &
            get_twoenum_nml, get_new_twoenum_nml

  integer(i_def), public, parameter :: first_one = 1952457118
  integer(i_def), public, parameter :: first_three = 1813125082
  integer(i_def), public, parameter :: first_two = 533081353
  integer(i_def), public, parameter :: second_ay = 1248446338
  integer(i_def), public, parameter :: second_bee = 144118421
  integer(i_def), public, parameter :: second_see = 359914450

  integer(i_def), public, protected :: first = emdi
  integer(i_def), public, protected :: second = emdi

  character(*), parameter :: listname = 'twoenum'
  character(str_def) :: profile_name = cmdi

  logical, parameter :: multiples_allowed = .false.

  logical :: nml_loaded = .false.

  character(str_def), parameter :: first_key(3) &
          = [character(len=str_def) :: 'one', &
                                       'three', &
                                       'two']
  character(str_def), parameter :: second_key(3) &
          = [character(len=str_def) :: 'ay', &
                                       'bee', &
                                       'see']

  integer(i_def), parameter :: first_value(3) &
          = [1952457118_i_def, &
             1813125082_i_def, &
             533081353_i_def]
  integer(i_def), parameter :: second_value(3) &
          = [1248446338_i_def, &
             144118421_i_def, &
             359914450_i_def]

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_def) function first_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_def) :: key_index

    if (key == unset_key) then
      write( log_scratch_space, '(A)') &
          'Missing key for first enumeration in twoenum namelist.'
      first_from_key = emdi
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      return
    end if

    key_index = 1
    do
      if (trim(first_key(key_index)) == trim(key)) then
        first_from_key = first_value(key_index)
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(first_key, 1)) then
          write( log_scratch_space, &
              '("Key ''", A, "'' not recognised for twoenum first")' ) &
              trim(adjustl(key))
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function first_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> @param[in] value Enumeration value.
  !>
  character(str_def) function key_from_first( value )

    implicit none

    integer(i_def), intent(in) :: value

    integer(i_def) :: value_index

    value_index = 1
    do
      if (first_value(value_index) == emdi) then
        key_from_first = unset_key
        return
      else if (first_value(value_index) == value) then
        key_from_first = first_key(value_index)
        return
      else
        value_index = value_index + 1
        if (value_index > ubound(first_key, 1)) then
          write( log_scratch_space, &
                 '("Value ", I0, " is not in twoenum first")' ) value
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function key_from_first

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_def) function second_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_def) :: key_index

    if (key == unset_key) then
      write( log_scratch_space, '(A)') &
          'Missing key for second enumeration in twoenum namelist.'
      second_from_key = emdi
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      return
    end if

    key_index = 1
    do
      if (trim(second_key(key_index)) == trim(key)) then
        second_from_key = second_value(key_index)
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(second_key, 1)) then
          write( log_scratch_space, &
              '("Key ''", A, "'' not recognised for twoenum second")' ) &
              trim(adjustl(key))
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function second_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> @param[in] value Enumeration value.
  !>
  character(str_def) function key_from_second( value )

    implicit none

    integer(i_def), intent(in) :: value

    integer(i_def) :: value_index

    value_index = 1
    do
      if (second_value(value_index) == emdi) then
        key_from_second = unset_key
        return
      else if (second_value(value_index) == value) then
        key_from_second = second_key(value_index)
        return
      else
        value_index = value_index + 1
        if (value_index > ubound(second_key, 1)) then
          write( log_scratch_space, &
                 '("Value ", I0, " is not in twoenum second")' ) value
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function key_from_second

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
  subroutine read_twoenum_namelist( file_unit, local_rank, scan )

    use constants_mod, only: i_def

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan

    call read_namelist( file_unit, local_rank, scan, &
                        first, &
                        second )

  end subroutine read_twoenum_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, scan, &
                            dummy_first, &
                            dummy_second )

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan
    integer(i_def), intent(out) :: dummy_first
    integer(i_def), intent(out) :: dummy_second

    integer(i_def) :: buffer_integer_i_def(2)

    character(str_def) :: first
    character(str_def) :: second

    namelist /twoenum/ first, &
                       second

    integer(i_def) :: condition

    first = unset_key
    second = unset_key

    if (local_rank == 0) then

      read( file_unit, nml=twoenum, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      dummy_first = first_from_key( first )
      dummy_second = second_from_key( second )

    end if

    buffer_integer_i_def(1) = dummy_first
    buffer_integer_i_def(2) = dummy_second

    call global_mpi%broadcast( buffer_integer_i_def, 2, 0 )

    dummy_first = buffer_integer_i_def(1)
    dummy_second = buffer_integer_i_def(2)

    if (scan) then
      nml_loaded = .false.
    else
      nml_loaded = .true.
    end if

  end subroutine read_namelist


  !> @brief Returns a <<namelist_type>> object populated with the
  !>        current contents of this configuration module.
  !> @return namelist_obj <<namelist_type>> with current namelist contents.
  function get_twoenum_nml() result(namelist_obj)

    implicit none

    type(namelist_type)      :: namelist_obj
    type(namelist_item_type) :: members(2)

    call members(1)%initialise( &
                'first', first )

    call members(2)%initialise( &
                'second', second )

    if (trim(profile_name) /= trim(cmdi) ) then
      call namelist_obj%initialise( trim(listname), &
                                    members, &
                                    profile_name = profile_name )
    else
      call namelist_obj%initialise( trim(listname), &
                                    members )
    end if

  end function get_twoenum_nml

  !> @brief Returns a <<twoenum_nml_type>> object populated with the
  !>        current contents of this configuration module.
  !> @return namelist_obj <<twoenum_nml_type>> with current namelist contents.
  function get_new_twoenum_nml() result(namelist_obj)

    use twoenum_nml_mod, only: twoenum_nml_type

    implicit none

    type(twoenum_nml_type) :: namelist_obj
    type(namelist_item_type) :: members(2)

    call members(1)%initialise( &
                'first', first )

    call members(2)%initialise( &
                'second', second )

    if (trim(profile_name) /= trim(cmdi) ) then
      call namelist_obj%initialise( trim(listname), &
                                    members, &
                                    profile_name = profile_name )
    else
      call namelist_obj%initialise( trim(listname), &
                                    members )
    end if

  end function get_new_twoenum_nml


  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_twoenum_namelist()

    use constants_mod, only: i_def

    implicit none

    ! Computed fields are resolved after everything has been loaded since they
    ! can refer to fields in other namelists.
    !

  end subroutine postprocess_twoenum_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function twoenum_is_loadable()

    implicit none

    logical :: twoenum_is_loadable

    if ( multiples_allowed .or. .not. nml_loaded ) then
      twoenum_is_loadable = .true.
    else
      twoenum_is_loadable = .false.
    end if

  end function twoenum_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function twoenum_is_loaded()

    implicit none

    logical :: twoenum_is_loaded

    twoenum_is_loaded = nml_loaded

  end function twoenum_is_loaded

  !> Are multiple twoenum namelists allowed to be read?
  !>
  !> @return True If multiple twoenum namelists are
  !>              permitted.
  !>
  function twoenum_multiples_allowed()

    implicit none

    logical :: twoenum_multiples_allowed

    twoenum_multiples_allowed = multiples_allowed

  end function twoenum_multiples_allowed

  !> Resets the load status to allow
  !> twoenum namelist to be read.
  !>
  subroutine twoenum_reset_load_status()

    implicit none

    nml_loaded = .false.

  end subroutine twoenum_reset_load_status

  !> Clear out any allocated memory
  !>
  subroutine twoenum_final()

    implicit none

    first = emdi
    second = emdi

    return
  end subroutine twoenum_final


end module twoenum_config_mod
