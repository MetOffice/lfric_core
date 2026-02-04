!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Manages the test namelist.
!>
module test_config_mod

  use constants_mod, only: i_def, &
                           i_long, &
                           i_short, &
                           l_def, &
                           r_def, &
                           r_double, &
                           r_second, &
                           r_single, &
                           str_def, &
                           str_max_filename
  use lfric_mpi_mod, only: global_mpi
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_DEBUG, LOG_LEVEL_INFO

  use namelist_mod,      only: namelist_type
  use namelist_item_mod, only: namelist_item_type

  use constants_mod, only: cmdi, emdi, imdi, rmdi, str_def, unset_key

  implicit none

  private
  public :: enum_from_key, key_from_enum, &
            read_test_namelist, postprocess_test_namelist, &
            test_is_loadable, test_is_loaded, &
            test_reset_load_status, &
            test_multiples_allowed, test_final, &
            get_test_nml, get_new_test_nml

  integer(i_def), public, parameter :: enum_one = 189779348
  integer(i_def), public, parameter :: enum_three = 1061269036
  integer(i_def), public, parameter :: enum_two = 1625932035

  integer(i_def), public, protected :: dint = imdi
  logical(l_def), public, protected :: dlog = .false.
  real(r_def), public, protected :: dreal = rmdi
  character(str_def), public, protected :: dstr = cmdi
  integer(i_def), public, protected :: enum = emdi
  character(str_max_filename), public, protected :: fstr = cmdi
  integer(i_long), public, protected :: lint = imdi
  real(r_double), public, protected :: lreal = rmdi
  integer(i_short), public, protected :: sint = imdi
  real(r_single), public, protected :: sreal = rmdi
  real(r_second), public, protected :: treal = rmdi
  integer(i_def), public, protected :: vint = imdi
  real(r_def), public, protected :: vreal = rmdi
  character(str_def), public, protected :: vstr = cmdi

  character(*), parameter :: listname = 'test'
  character(str_def) :: profile_name = cmdi

  logical, parameter :: multiples_allowed = .false.

  logical :: nml_loaded = .false.

  character(str_def), parameter :: enum_key(3) &
          = [character(len=str_def) :: 'one', &
                                       'three', &
                                       'two']

  integer(i_def), parameter :: enum_value(3) &
          = [189779348_i_def, &
             1061269036_i_def, &
             1625932035_i_def]

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_def) function enum_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_def) :: key_index

    if (key == unset_key) then
      write( log_scratch_space, '(A)') &
          'Missing key for enum enumeration in test namelist.'
      enum_from_key = emdi
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      return
    end if

    key_index = 1
    do
      if (trim(enum_key(key_index)) == trim(key)) then
        enum_from_key = enum_value(key_index)
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(enum_key, 1)) then
          write( log_scratch_space, &
              '("Key ''", A, "'' not recognised for test enum")' ) &
              trim(adjustl(key))
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function enum_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> @param[in] value Enumeration value.
  !>
  character(str_def) function key_from_enum( value )

    implicit none

    integer(i_def), intent(in) :: value

    integer(i_def) :: value_index

    value_index = 1
    do
      if (enum_value(value_index) == emdi) then
        key_from_enum = unset_key
        return
      else if (enum_value(value_index) == value) then
        key_from_enum = enum_key(value_index)
        return
      else
        value_index = value_index + 1
        if (value_index > ubound(enum_key, 1)) then
          write( log_scratch_space, &
                 '("Value ", I0, " is not in test enum")' ) value
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function key_from_enum

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

    call read_namelist( file_unit, local_rank, scan, &
                        enum )

  end subroutine read_test_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, scan, &
                            dummy_enum )

    implicit none

    integer(i_def), intent(in) :: file_unit
    integer(i_def), intent(in) :: local_rank
    logical,        intent(in) :: scan
    integer(i_def), intent(out) :: dummy_enum

    character(str_def) :: buffer_character_str_def(2)
    character(str_max_filename) :: buffer_character_str_max_filename(1)
    integer(i_def) :: buffer_integer_i_def(3)
    integer(i_long) :: buffer_integer_i_long(1)
    integer(i_short) :: buffer_integer_i_short(1)
    integer(i_def) :: buffer_logical_l_def(1)
    real(r_def) :: buffer_real_r_def(2)
    real(r_double) :: buffer_real_r_double(1)
    real(r_second) :: buffer_real_r_second(1)
    real(r_single) :: buffer_real_r_single(1)

    character(str_def) :: enum

    namelist /test/ dint, &
                    dlog, &
                    dreal, &
                    dstr, &
                    enum, &
                    fstr, &
                    lint, &
                    lreal, &
                    sint, &
                    sreal, &
                    treal, &
                    vint, &
                    vreal, &
                    vstr

    integer(i_def) :: condition

    dint = imdi
    dlog = .false.
    dreal = rmdi
    dstr = cmdi
    enum = unset_key
    fstr = cmdi
    lint = imdi
    lreal = rmdi
    sint = imdi
    sreal = rmdi
    treal = rmdi
    vint = imdi
    vreal = rmdi
    vstr = cmdi

    if (local_rank == 0) then

      read( file_unit, nml=test, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      dummy_enum = enum_from_key( enum )

    end if

    buffer_integer_i_def(2) = dint
    buffer_logical_l_def(1) = merge( 1, 0, dlog )
    buffer_real_r_def(2) = dreal
    buffer_character_str_def(2) = dstr
    buffer_integer_i_def(3) = dummy_enum
    buffer_character_str_max_filename(1) = fstr
    buffer_integer_i_long(1) = lint
    buffer_real_r_double(1) = lreal
    buffer_integer_i_short(1) = sint
    buffer_real_r_single(1) = sreal
    buffer_real_r_second(1) = treal
    buffer_integer_i_def(1) = vint
    buffer_real_r_def(1) = vreal
    buffer_character_str_def(1) = vstr

    call global_mpi%broadcast( buffer_character_str_def, 2*str_def, 0 )
    call global_mpi%broadcast( buffer_character_str_max_filename, 1*str_max_filename, 0 )
    call global_mpi%broadcast( buffer_integer_i_def, 3, 0 )
    call global_mpi%broadcast( buffer_integer_i_long, 1, 0 )
    call global_mpi%broadcast( buffer_integer_i_short, 1, 0 )
    call global_mpi%broadcast( buffer_logical_l_def, 1, 0 )
    call global_mpi%broadcast( buffer_real_r_def, 2, 0 )
    call global_mpi%broadcast( buffer_real_r_double, 1, 0 )
    call global_mpi%broadcast( buffer_real_r_second, 1, 0 )
    call global_mpi%broadcast( buffer_real_r_single, 1, 0 )

    dint = buffer_integer_i_def(2)
    dlog = buffer_logical_l_def(1) /= 0
    dreal = buffer_real_r_def(2)
    dstr = buffer_character_str_def(2)
    dummy_enum = buffer_integer_i_def(3)
    fstr = buffer_character_str_max_filename(1)
    lint = buffer_integer_i_long(1)
    lreal = buffer_real_r_double(1)
    sint = buffer_integer_i_short(1)
    sreal = buffer_real_r_single(1)
    treal = buffer_real_r_second(1)
    vint = buffer_integer_i_def(1)
    vreal = buffer_real_r_def(1)
    vstr = buffer_character_str_def(1)

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
    type(namelist_item_type) :: members(14)

    call members(1)%initialise( &
                'dint', dint )

    call members(2)%initialise( &
                'dlog', dlog )

    call members(3)%initialise( &
                'dreal', dreal )

    call members(4)%initialise( &
                'dstr', dstr )

    call members(5)%initialise( &
                'enum', enum )

    call members(6)%initialise( &
                'fstr', fstr )

    call members(7)%initialise( &
                'lint', lint )

    call members(8)%initialise( &
                'lreal', lreal )

    call members(9)%initialise( &
                'sint', sint )

    call members(10)%initialise( &
                'sreal', sreal )

    call members(11)%initialise( &
                'treal', treal )

    call members(12)%initialise( &
                'vint', vint )

    call members(13)%initialise( &
                'vreal', vreal )

    call members(14)%initialise( &
                'vstr', vstr )

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
    type(namelist_item_type) :: members(14)

    call members(1)%initialise( &
                'dint', dint )

    call members(2)%initialise( &
                'dlog', dlog )

    call members(3)%initialise( &
                'dreal', dreal )

    call members(4)%initialise( &
                'dstr', dstr )

    call members(5)%initialise( &
                'enum', enum )

    call members(6)%initialise( &
                'fstr', fstr )

    call members(7)%initialise( &
                'lint', lint )

    call members(8)%initialise( &
                'lreal', lreal )

    call members(9)%initialise( &
                'sint', sint )

    call members(10)%initialise( &
                'sreal', sreal )

    call members(11)%initialise( &
                'treal', treal )

    call members(12)%initialise( &
                'vint', vint )

    call members(13)%initialise( &
                'vreal', vreal )

    call members(14)%initialise( &
                'vstr', vstr )

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

    dint = imdi
    dlog = .false.
    dreal = real(rmdi,r_def)
    dstr = cmdi
    enum = emdi
    fstr = cmdi
    lint = imdi
    lreal = real(rmdi,r_double)
    sint = imdi
    sreal = real(rmdi,r_single)
    treal = real(rmdi,r_second)
    vint = imdi
    vreal = real(rmdi,r_def)
    vstr = cmdi

    return
  end subroutine test_final


end module test_config_mod
