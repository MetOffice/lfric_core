!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the enum namelist.
!>
module enum_config_mod

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
  public :: value_from_key, key_from_value, &
            read_enum_namelist, postprocess_enum_namelist, &
            enum_is_loadable, enum_is_loaded, &
            enum_reset_load_status, &
            enum_multiples_allowed, enum_final, &
            get_enum_nml, get_new_enum_nml

  integer(i_def), public, parameter :: value_one = 1695414371
  integer(i_def), public, parameter :: value_three = 839906103
  integer(i_def), public, parameter :: value_two = 246150388

  integer(i_def), public, protected :: value = emdi

  character(*), parameter :: listname = 'enum'
  character(str_def) :: profile_name = cmdi

  logical, parameter :: multiples_allowed = .false.

  logical :: nml_loaded = .false.

  character(str_def), parameter :: value_key(3) &
          = [character(len=str_def) :: 'one', &
                                       'three', &
                                       'two']

  integer(i_def), parameter :: value_value(3) &
          = [1695414371_i_def, &
             839906103_i_def, &
             246150388_i_def]

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_def) function value_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_def) :: key_index

    if (key == unset_key) then
      write( log_scratch_space, '(A)') &
          'Missing key for value enumeration in enum namelist.'
      value_from_key = emdi
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      return
    end if

    key_index = 1
    do
      if (trim(value_key(key_index)) == trim(key)) then
        value_from_key = value_value(key_index)
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(value_key, 1)) then
          write( log_scratch_space, &
              '("Key ''", A, "'' not recognised for enum value")' ) &
              trim(adjustl(key))
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function value_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> @param[in] value Enumeration value.
  !>
  character(str_def) function key_from_value( value )

    implicit none

    integer(i_def), intent(in) :: value

    integer(i_def) :: value_index

    value_index = 1
    do
      if (value_value(value_index) == emdi) then
        key_from_value = unset_key
        return
      else if (value_value(value_index) == value) then
        key_from_value = value_key(value_index)
        return
      else
        value_index = value_index + 1
        if (value_index > ubound(value_key, 1)) then
          write( log_scratch_space, &
                 '("Value ", I0, " is not in enum value")' ) value
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function key_from_value

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
  subroutine read_enum_namelist( file_unit, local_rank, scan )

    use constants_mod, only: i_def

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan

    call read_namelist( file_unit, local_rank, scan, &
                        value )

  end subroutine read_enum_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, scan, &
                            dummy_value )

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan
    integer(i_def), intent(out) :: dummy_value

    integer(i_def) :: buffer_integer_i_def(1)

    character(str_def) :: value

    namelist /enum/ value

    integer(i_def) :: condition

    value = unset_key

    if (local_rank == 0) then

      read( file_unit, nml=enum, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      dummy_value = value_from_key( value )

    end if

    buffer_integer_i_def(1) = dummy_value

    call global_mpi%broadcast( buffer_integer_i_def, 1, 0 )

    dummy_value = buffer_integer_i_def(1)

    if (scan) then
      nml_loaded = .false.
    else
      nml_loaded = .true.
    end if

  end subroutine read_namelist


  !> @brief Returns a <<namelist_type>> object populated with the
  !>        current contents of this configuration module.
  !> @return namelist_obj <<namelist_type>> with current namelist contents.
  function get_enum_nml() result(namelist_obj)

    implicit none

    type(namelist_type)      :: namelist_obj
    type(namelist_item_type) :: members(1)

    call members(1)%initialise( &
                'value', value )

    if (trim(profile_name) /= trim(cmdi) ) then
      call namelist_obj%initialise( trim(listname), &
                                    members, &
                                    profile_name = profile_name )
    else
      call namelist_obj%initialise( trim(listname), &
                                    members )
    end if

  end function get_enum_nml

  !> @brief Returns a <<enum_nml_type>> object populated with the
  !>        current contents of this configuration module.
  !> @return namelist_obj <<enum_nml_type>> with current namelist contents.
  function get_new_enum_nml() result(namelist_obj)

    use enum_nml_mod, only: enum_nml_type

    implicit none

    type(enum_nml_type) :: namelist_obj
    type(namelist_item_type) :: members(1)

    call members(1)%initialise( &
                'value', value )

    if (trim(profile_name) /= trim(cmdi) ) then
      call namelist_obj%initialise( trim(listname), &
                                    members, &
                                    profile_name = profile_name )
    else
      call namelist_obj%initialise( trim(listname), &
                                    members )
    end if

  end function get_new_enum_nml


  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_enum_namelist()

    use constants_mod, only: i_def

    implicit none

    ! Computed fields are resolved after everything has been loaded since they
    ! can refer to fields in other namelists.
    !

  end subroutine postprocess_enum_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function enum_is_loadable()

    implicit none

    logical :: enum_is_loadable

    if ( multiples_allowed .or. .not. nml_loaded ) then
      enum_is_loadable = .true.
    else
      enum_is_loadable = .false.
    end if

  end function enum_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function enum_is_loaded()

    implicit none

    logical :: enum_is_loaded

    enum_is_loaded = nml_loaded

  end function enum_is_loaded

  !> Are multiple enum namelists allowed to be read?
  !>
  !> @return True If multiple enum namelists are
  !>              permitted.
  !>
  function enum_multiples_allowed()

    implicit none

    logical :: enum_multiples_allowed

    enum_multiples_allowed = multiples_allowed

  end function enum_multiples_allowed

  !> Resets the load status to allow
  !> enum namelist to be read.
  !>
  subroutine enum_reset_load_status()

    implicit none

    nml_loaded = .false.

  end subroutine enum_reset_load_status

  !> Clear out any allocated memory
  !>
  subroutine enum_final()

    implicit none

    value = emdi

    return
  end subroutine enum_final


end module enum_config_mod
