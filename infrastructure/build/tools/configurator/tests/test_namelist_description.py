#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
import io
import random
import os
from pathlib import Path
import pytest
from subprocess import Popen
from textwrap import dedent

import configurator.namelistdescription as description

PICKER_EXE = 'rose_picker'


class TestNamelistMeta():
    def test_module_write_empty(self):
        output_file = io.StringIO()

        uut = description.NamelistDescription('test')
        with pytest.raises(description.NamelistDescriptionException):
            uut.write_module(output_file)

    def test_module_write_one_of_each(self):
        expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the test namelist.
!>
module test_config_mod

  use constants_mod, only: i_def, &
                           i_long, &
                           i_native, &
                           i_short, &
                           l_def, &
                           r_def, &
                           r_double, &
                           r_single, &
                           str_def, &
                           str_max_filename
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key

  implicit none

  private
  public :: enum_from_key, key_from_enum, &
            read_test_namelist, postprocess_test_namelist, &
            test_is_loadable, test_is_loaded, test_final

  integer(i_native), public, parameter :: enum_one = 138
  integer(i_native), public, parameter :: enum_three = 870
  integer(i_native), public, parameter :: enum_two = 584

  integer(i_def), public, protected :: dint = imdi
  logical(l_def), public, protected :: dlog = .false.
  real(r_def), public, protected :: dreal = rmdi
  character(str_def), public, protected :: dstr = cmdi
  integer(i_native), public, protected :: enum = emdi
  character(str_max_filename), public, protected :: fstr = cmdi
  integer(i_long), public, protected :: lint = imdi
  real(r_double), public, protected :: lreal = rmdi
  integer(i_short), public, protected :: sint = imdi
  real(r_single), public, protected :: sreal = rmdi
  integer(i_def), public, protected :: vint = imdi
  real(r_def), public, protected :: vreal = rmdi
  character(str_def), public, protected :: vstr = cmdi

  logical :: namelist_loaded = .false.

  character(str_def), parameter :: enum_key(3) &
          = [character(len=str_def) :: 'one', &
                                       'three', &
                                       'two']

  integer(i_native), parameter :: enum_value(3) &
          = [138_i_native, &
             870_i_native, &
             584_i_native]

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_native) function enum_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    if (key == unset_key) then
      write( log_scratch_space, '(A)') &
          'Missing key for enum enumeration in test namelist.'
      enum_from_key = emdi
      call log_event( log_scratch_space, LOG_LEVEL_WARNING )
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

    integer(i_native), intent(in) :: value

    integer(i_native) :: value_index

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
  !>
  subroutine read_test_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank, &
                        enum )

  end subroutine read_test_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, &
                            dummy_enum )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_native), intent(out) :: dummy_enum

    character(str_def) :: buffer_character_str_def(2)
    character(str_max_filename) :: buffer_character_str_max_filename(1)
    integer(i_def) :: buffer_integer_i_def(2)
    integer(i_long) :: buffer_integer_i_long(1)
    integer(i_native) :: buffer_integer_i_native(1)
    integer(i_short) :: buffer_integer_i_short(1)
    integer(i_native) :: buffer_logical_l_def(1)
    real(r_def) :: buffer_real_r_def(2)
    real(r_double) :: buffer_real_r_double(1)
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
                    vint, &
                    vreal, &
                    vstr

    integer(i_native) :: condition

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
    buffer_integer_i_native(1) = dummy_enum
    buffer_character_str_max_filename(1) = fstr
    buffer_integer_i_long(1) = lint
    buffer_real_r_double(1) = lreal
    buffer_integer_i_short(1) = sint
    buffer_real_r_single(1) = sreal
    buffer_integer_i_def(1) = vint
    buffer_real_r_def(1) = vreal
    buffer_character_str_def(1) = vstr

    call broadcast( buffer_character_str_def, 2*str_def, 0 )
    call broadcast( buffer_character_str_max_filename, 1*str_max_filename, 0 )
    call broadcast( buffer_integer_i_def, 2, 0 )
    call broadcast( buffer_integer_i_long, 1, 0 )
    call broadcast( buffer_integer_i_native, 1, 0 )
    call broadcast( buffer_integer_i_short, 1, 0 )
    call broadcast( buffer_logical_l_def, 1, 0 )
    call broadcast( buffer_real_r_def, 2, 0 )
    call broadcast( buffer_real_r_double, 1, 0 )
    call broadcast( buffer_real_r_single, 1, 0 )

    dint = buffer_integer_i_def(2)
    dlog = buffer_logical_l_def(1) /= 0
    dreal = buffer_real_r_def(2)
    dstr = buffer_character_str_def(2)
    dummy_enum = buffer_integer_i_native(1)
    fstr = buffer_character_str_max_filename(1)
    lint = buffer_integer_i_long(1)
    lreal = buffer_real_r_double(1)
    sint = buffer_integer_i_short(1)
    sreal = buffer_real_r_single(1)
    vint = buffer_integer_i_def(1)
    vreal = buffer_real_r_def(1)
    vstr = buffer_character_str_def(1)



    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_test_namelist()

    implicit none


  end subroutine postprocess_test_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function test_is_loadable()

    implicit none

    logical :: test_is_loadable

    test_is_loadable = .not. namelist_loaded

  end function test_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function test_is_loaded()

    implicit none

    logical :: test_is_loaded

    test_is_loaded = namelist_loaded

  end function test_is_loaded

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
    vint = imdi
    vreal = real(rmdi,r_def)
    vstr = cmdi

    return
  end subroutine test_final


end module test_config_mod
        '''.strip()

        random.seed(1)
        uut = description.NamelistDescription('test')
        uut.add_value('vint', 'integer')
        uut.add_value('dint', 'integer', 'default')
        uut.add_value('sint', 'integer', 'short')
        uut.add_value('lint', 'integer', 'long')
        uut.add_value('dlog', 'logical', 'default')
        uut.add_value('vreal', 'real')
        uut.add_value('dreal', 'real', 'default')
        uut.add_value('sreal', 'real', 'single')
        uut.add_value('lreal', 'real', 'double')
        uut.add_string('vstr')
        uut.add_string('dstr', configure_string_length='default')
        uut.add_string('fstr', configure_string_length='filename')
        uut.add_enumeration('enum', enumerators=['one', 'two', 'three'])
        output_file = io.StringIO()
        uut.write_module(output_file)

        assert expected_source + '\n' == output_file.getvalue()

    _FIRST_EXPECTED_GROWING_SOURCE = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the test namelist.
!>
module test_config_mod

  use constants_mod, only: i_def, &
                           i_native
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key

  implicit none

  private
  public :: read_test_namelist, postprocess_test_namelist, &
            test_is_loadable, test_is_loaded, test_final

  integer(i_def), public, protected :: foo = imdi

  logical :: namelist_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_test_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank )

  end subroutine read_test_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    integer(i_def) :: buffer_integer_i_def(1)

    namelist /test/ foo

    integer(i_native) :: condition

    foo = imdi

    if (local_rank == 0) then

      read( file_unit, nml=test, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_integer_i_def(1) = foo

    call broadcast( buffer_integer_i_def, 1, 0 )

    foo = buffer_integer_i_def(1)



    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_test_namelist()

    implicit none


  end subroutine postprocess_test_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function test_is_loadable()

    implicit none

    logical :: test_is_loadable

    test_is_loadable = .not. namelist_loaded

  end function test_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function test_is_loaded()

    implicit none

    logical :: test_is_loaded

    test_is_loaded = namelist_loaded

  end function test_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine test_final()

    implicit none

    foo = imdi

    return
  end subroutine test_final


end module test_config_mod
        '''.strip() + '\n'

    _SECOND_EXPECTED_GROWING_SOURCE = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the test namelist.
!>
module test_config_mod

  use constants_mod, only: i_def, &
                           i_native, &
                           r_def
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key

  implicit none

  private
  public :: read_test_namelist, postprocess_test_namelist, &
            test_is_loadable, test_is_loaded, test_final

  real(r_def), public, protected :: bar = rmdi
  integer(i_def), public, protected :: foo = imdi

  logical :: namelist_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_test_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank )

  end subroutine read_test_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    integer(i_def) :: buffer_integer_i_def(1)
    real(r_def) :: buffer_real_r_def(1)

    namelist /test/ bar, &
                    foo

    integer(i_native) :: condition

    bar = rmdi
    foo = imdi

    if (local_rank == 0) then

      read( file_unit, nml=test, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_real_r_def(1) = bar
    buffer_integer_i_def(1) = foo

    call broadcast( buffer_integer_i_def, 1, 0 )
    call broadcast( buffer_real_r_def, 1, 0 )

    bar = buffer_real_r_def(1)
    foo = buffer_integer_i_def(1)



    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_test_namelist()

    implicit none


  end subroutine postprocess_test_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function test_is_loadable()

    implicit none

    logical :: test_is_loadable

    test_is_loadable = .not. namelist_loaded

  end function test_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function test_is_loaded()

    implicit none

    logical :: test_is_loaded

    test_is_loaded = namelist_loaded

  end function test_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine test_final()

    implicit none

    bar = real(rmdi,r_def)
    foo = imdi

    return
  end subroutine test_final


end module test_config_mod
        '''.strip() + '\n'

    def test_module_write_growing(self):
        output_file = io.StringIO()

        uut = description.NamelistDescription('test')
        uut.add_value('foo', 'integer')
        uut.write_module(output_file)

        assert output_file.getvalue() == self._FIRST_EXPECTED_GROWING_SOURCE

        output_file = io.StringIO()
        uut.add_value('bar', 'real', 'default')
        uut.write_module(output_file)

        assert output_file.getvalue() == self._SECOND_EXPECTED_GROWING_SOURCE

    _EXPECTED_ENUM_ONLY_SOURCE = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the enum namelist.
!>
module enum_config_mod

  use constants_mod, only: i_native, &
                           str_def
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key

  implicit none

  private
  public :: value_from_key, key_from_value, &
            read_enum_namelist, postprocess_enum_namelist, &
            enum_is_loadable, enum_is_loaded, enum_final

  integer(i_native), public, parameter :: value_one = 138
  integer(i_native), public, parameter :: value_three = 870
  integer(i_native), public, parameter :: value_two = 584

  integer(i_native), public, protected :: value = emdi

  logical :: namelist_loaded = .false.

  character(str_def), parameter :: value_key(3) &
          = [character(len=str_def) :: 'one', &
                                       'three', &
                                       'two']

  integer(i_native), parameter :: value_value(3) &
          = [138_i_native, &
             870_i_native, &
             584_i_native]

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_native) function value_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    if (key == unset_key) then
      write( log_scratch_space, '(A)') &
          'Missing key for value enumeration in enum namelist.'
      value_from_key = emdi
      call log_event( log_scratch_space, LOG_LEVEL_WARNING )
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

    integer(i_native), intent(in) :: value

    integer(i_native) :: value_index

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
  !>
  subroutine read_enum_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank, &
                        value )

  end subroutine read_enum_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, &
                            dummy_value )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_native), intent(out) :: dummy_value

    integer(i_native) :: buffer_integer_i_native(1)

    character(str_def) :: value

    namelist /enum/ value

    integer(i_native) :: condition

    value = unset_key

    if (local_rank == 0) then

      read( file_unit, nml=enum, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      dummy_value = value_from_key( value )

    end if

    buffer_integer_i_native(1) = dummy_value

    call broadcast( buffer_integer_i_native, 1, 0 )

    dummy_value = buffer_integer_i_native(1)



    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_enum_namelist()

    implicit none


  end subroutine postprocess_enum_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function enum_is_loadable()

    implicit none

    logical :: enum_is_loadable

    enum_is_loadable = .not. namelist_loaded

  end function enum_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function enum_is_loaded()

    implicit none

    logical :: enum_is_loaded

    enum_is_loaded = namelist_loaded

  end function enum_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine enum_final()

    implicit none

    value = emdi

    return
  end subroutine enum_final


end module enum_config_mod
        '''.strip() + '\n'

    def test_enumeration_only(self):
        random.seed(1)
        uut = description.NamelistDescription('enum')
        uut.add_enumeration('value', enumerators=['one', 'two', 'three'])
        output_file = io.StringIO()
        uut.write_module(output_file)

        assert output_file.getvalue() == self._EXPECTED_ENUM_ONLY_SOURCE

    _EXPECTED_TWO_ENUM_SOURCE = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the twoenum namelist.
!>
module twoenum_config_mod

  use constants_mod, only: i_native, &
                           str_def
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key

  implicit none

  private
  public :: first_from_key, key_from_first, &
            second_from_key, key_from_second, &
            read_twoenum_namelist, postprocess_twoenum_namelist, &
            twoenum_is_loadable, twoenum_is_loaded, twoenum_final

  integer(i_native), public, parameter :: first_one = 138
  integer(i_native), public, parameter :: first_three = 870
  integer(i_native), public, parameter :: first_two = 584
  integer(i_native), public, parameter :: second_ay = 824
  integer(i_native), public, parameter :: second_bee = 785
  integer(i_native), public, parameter :: second_see = 65

  integer(i_native), public, protected :: first = emdi
  integer(i_native), public, protected :: second = emdi

  logical :: namelist_loaded = .false.

  character(str_def), parameter :: first_key(3) &
          = [character(len=str_def) :: 'one', &
                                       'three', &
                                       'two']
  character(str_def), parameter :: second_key(3) &
          = [character(len=str_def) :: 'ay', &
                                       'bee', &
                                       'see']

  integer(i_native), parameter :: first_value(3) &
          = [138_i_native, &
             870_i_native, &
             584_i_native]
  integer(i_native), parameter :: second_value(3) &
          = [824_i_native, &
             785_i_native, &
             65_i_native]

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_native) function first_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    if (key == unset_key) then
      write( log_scratch_space, '(A)') &
          'Missing key for first enumeration in twoenum namelist.'
      first_from_key = emdi
      call log_event( log_scratch_space, LOG_LEVEL_WARNING )
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

    integer(i_native), intent(in) :: value

    integer(i_native) :: value_index

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
  integer(i_native) function second_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    if (key == unset_key) then
      write( log_scratch_space, '(A)') &
          'Missing key for second enumeration in twoenum namelist.'
      second_from_key = emdi
      call log_event( log_scratch_space, LOG_LEVEL_WARNING )
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

    integer(i_native), intent(in) :: value

    integer(i_native) :: value_index

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
  !>
  subroutine read_twoenum_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank, &
                        first, &
                        second )

  end subroutine read_twoenum_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, &
                            dummy_first, &
                            dummy_second )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_native), intent(out) :: dummy_first
    integer(i_native), intent(out) :: dummy_second

    integer(i_native) :: buffer_integer_i_native(2)

    character(str_def) :: first
    character(str_def) :: second

    namelist /twoenum/ first, &
                       second

    integer(i_native) :: condition

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

    buffer_integer_i_native(1) = dummy_first
    buffer_integer_i_native(2) = dummy_second

    call broadcast( buffer_integer_i_native, 2, 0 )

    dummy_first = buffer_integer_i_native(1)
    dummy_second = buffer_integer_i_native(2)



    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_twoenum_namelist()

    implicit none


  end subroutine postprocess_twoenum_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function twoenum_is_loadable()

    implicit none

    logical :: twoenum_is_loadable

    twoenum_is_loadable = .not. namelist_loaded

  end function twoenum_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function twoenum_is_loaded()

    implicit none

    logical :: twoenum_is_loaded

    twoenum_is_loaded = namelist_loaded

  end function twoenum_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine twoenum_final()

    implicit none

    first = emdi
    second = emdi

    return
  end subroutine twoenum_final


end module twoenum_config_mod
                       '''.strip() + '\n'

    def test_more_than_one_enumeration(self):
        random.seed(1)
        uut = description.NamelistDescription('twoenum')
        uut.add_enumeration('first', enumerators=['one', 'two', 'three'])
        uut.add_enumeration('second', enumerators=['ay', 'bee', 'see'])
        output_file = io.StringIO()
        uut.write_module(output_file)

        assert output_file.getvalue() == self._EXPECTED_TWO_ENUM_SOURCE

    _EXPECTED_MODULE_WRITE_COMPUTED_SOURCE = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the teapot namelist.
!>
module teapot_config_mod

  use constants_mod, only: i_native, &
                           r_def
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key
  use fridge_config_mod, only: milk

  implicit none

  private
  public :: read_teapot_namelist, postprocess_teapot_namelist, &
            teapot_is_loadable, teapot_is_loaded, teapot_final

  real(r_def), public, protected :: bar = rmdi
  real(r_def), public, protected :: bif = rmdi
  real(r_def), public, protected :: foo = rmdi
  real(r_def), public, protected :: fum = rmdi

  logical :: namelist_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_teapot_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank )

  end subroutine read_teapot_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    real(r_def) :: buffer_real_r_def(2)

    namelist /teapot/ foo, &
                      fum

    integer(i_native) :: condition

    bar = rmdi
    bif = rmdi
    foo = rmdi
    fum = rmdi

    if (local_rank == 0) then

      read( file_unit, nml=teapot, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_real_r_def(1) = foo
    buffer_real_r_def(2) = fum

    call broadcast( buffer_real_r_def, 2, 0 )

    foo = buffer_real_r_def(1)
    fum = buffer_real_r_def(2)


    bar = foo ** 2
    if ( any([fum, foo, milk] == imdi) .or. &
         any([fum, foo, milk] == rmdi) ) then
      bif = rmdi
    else
      bif = fum + (foo ** 2) / milk
    end if

    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_teapot_namelist()

    implicit none


  end subroutine postprocess_teapot_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function teapot_is_loadable()

    implicit none

    logical :: teapot_is_loadable

    teapot_is_loadable = .not. namelist_loaded

  end function teapot_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function teapot_is_loaded()

    implicit none

    logical :: teapot_is_loaded

    teapot_is_loaded = namelist_loaded

  end function teapot_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine teapot_final()

    implicit none

    bar = real(rmdi,r_def)
    bif = real(rmdi,r_def)
    foo = real(rmdi,r_def)
    fum = real(rmdi,r_def)

    return
  end subroutine teapot_final


end module teapot_config_mod
        '''.strip() + '\n'

    def test_module_write_computed(self):
        output_file = io.StringIO()

        uut = description.NamelistDescription('teapot')
        uut.add_value('foo', 'real', 'default')
        uut.add_value('fum', 'real', 'default')
        uut.add_computed('bar', 'real', 'default', calculation=['foo ** 2'])
        uut.add_computed('bif', 'real', 'default',
                         calculation=[
                             'namelist:teapot=fum + ' +
                             '(namelist:teapot=foo ** 2) / ' +
                             'namelist:fridge=milk'])
        uut.write_module(output_file)

        assert output_file.getvalue() \
            == self._EXPECTED_MODULE_WRITE_COMPUTED_SOURCE

    _EXPECTED_MODULE_WRITE_CONSTANT_SOURCE = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the cheese namelist.
!>
module cheese_config_mod

  use constants_mod, only: i_native, &
                           r_def
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, FUDGE, imdi, rmdi, unset_key

  implicit none

  private
  public :: read_cheese_namelist, postprocess_cheese_namelist, &
            cheese_is_loadable, cheese_is_loaded, cheese_final

  real(r_def), public, protected :: fred = rmdi
  real(r_def), public, protected :: wilma = rmdi

  logical :: namelist_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_cheese_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank )

  end subroutine read_cheese_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    real(r_def) :: buffer_real_r_def(1)

    namelist /cheese/ fred

    integer(i_native) :: condition

    fred = rmdi
    wilma = rmdi

    if (local_rank == 0) then

      read( file_unit, nml=cheese, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_real_r_def(1) = fred

    call broadcast( buffer_real_r_def, 1, 0 )

    fred = buffer_real_r_def(1)


    wilma = fred * FUDGE

    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_cheese_namelist()

    implicit none


  end subroutine postprocess_cheese_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function cheese_is_loadable()

    implicit none

    logical :: cheese_is_loadable

    cheese_is_loadable = .not. namelist_loaded

  end function cheese_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function cheese_is_loaded()

    implicit none

    logical :: cheese_is_loaded

    cheese_is_loaded = namelist_loaded

  end function cheese_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine cheese_final()

    implicit none

    fred = real(rmdi,r_def)
    wilma = real(rmdi,r_def)

    return
  end subroutine cheese_final


end module cheese_config_mod
        '''.strip() + '\n'

    def test_module_write_constant(self):
        output_file = io.StringIO()

        uut = description.NamelistDescription('cheese')
        uut.add_usage('FUDGE', module='constants_mod')
        uut.add_value('fred', 'real', 'default')
        uut.add_computed('wilma', 'real', 'default',
                         calculation=['fred * FUDGE'])
        uut.write_module(output_file)

        assert output_file.getvalue() \
            == self._EXPECTED_MODULE_WRITE_CONSTANT_SOURCE

    _EXPECTED_MODULE_WRITE_ARRAY_SOURCE = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the aerial namelist.
!>
module aerial_config_mod

  use constants_mod, only: i_def, &
                           i_native, &
                           r_def, &
                           str_def
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key
  use wibble_mod, only: esize

  implicit none

  private
  public :: read_aerial_namelist, postprocess_aerial_namelist, &
            aerial_is_loadable, aerial_is_loaded, aerial_final

  integer(i_native), parameter, public :: max_array_size = 100

  character(str_def), public, protected :: absolute(5) = cmdi
  integer(i_def), public, protected, allocatable :: inlist(:)
  integer(i_native), public, protected :: lsize = imdi
  real(r_def), public, protected, allocatable :: outlist(:)
  integer(i_def), public, protected, allocatable :: unknown(:)

  logical :: namelist_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_aerial_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank )

  end subroutine read_aerial_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    integer(i_native) :: buffer_integer_i_native(1)

    namelist /aerial/ absolute, &
                      inlist, &
                      lsize, &
                      outlist, &
                      unknown

    integer(i_native) :: condition

    allocate( inlist(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "inlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    allocate( outlist(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "outlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    allocate( unknown(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "unknown"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    absolute = cmdi
    inlist = imdi
    lsize = imdi
    outlist = rmdi
    unknown = imdi

    if (local_rank == 0) then

      read( file_unit, nml=aerial, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_integer_i_native(1) = lsize

    call broadcast( buffer_integer_i_native, 1, 0 )

    lsize = buffer_integer_i_native(1)



    call broadcast( absolute, size(absolute, 1)*str_def, 0 )
    call broadcast( inlist, size(inlist, 1), 0 )
    call broadcast( outlist, size(outlist, 1), 0 )
    call broadcast( unknown, size(unknown, 1), 0 )

    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_aerial_namelist()

    implicit none

    integer(i_native) :: condition
    integer(i_native) :: size
    integer(i_def), allocatable :: new_inlist(:)
    real(r_def), allocatable :: new_outlist(:)
    integer(i_def), allocatable :: new_unknown(:)
    integer(i_native) :: index_unknown

    condition = 0
    size      = 0

    size = lsize
    allocate( new_inlist(size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "inlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_inlist(:size) = inlist(:size)
    call move_alloc( new_inlist, inlist )
    if (allocated(new_inlist)) deallocate( new_inlist)
    size = esize
    allocate( new_outlist(size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "outlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_outlist(:size) = outlist(:size)
    call move_alloc( new_outlist, outlist )
    if (allocated(new_outlist)) deallocate( new_outlist)
    do index_unknown=ubound(unknown, 1), 1, -1
      if (unknown(index_unknown) /= imdi) exit
    end do
    size = index_unknown
    allocate( new_unknown(size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "unknown"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_unknown(:size) = unknown(:size)
    call move_alloc( new_unknown, unknown )
    if (allocated(new_unknown)) deallocate( new_unknown)

  end subroutine postprocess_aerial_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function aerial_is_loadable()

    implicit none

    logical :: aerial_is_loadable

    aerial_is_loadable = .not. namelist_loaded

  end function aerial_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function aerial_is_loaded()

    implicit none

    logical :: aerial_is_loaded

    aerial_is_loaded = namelist_loaded

  end function aerial_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine aerial_final()

    implicit none

    absolute = cmdi
    lsize = imdi

    if ( allocated(inlist) ) deallocate(inlist)
    if ( allocated(outlist) ) deallocate(outlist)
    if ( allocated(unknown) ) deallocate(unknown)

    return
  end subroutine aerial_final


end module aerial_config_mod
        '''.strip() + '\n'

    def test_module_write_array(self):
        output_file = io.StringIO()

        uut = description.NamelistDescription('aerial')
        uut.add_usage('esize', module='wibble_mod')
        uut.add_value('lsize', 'integer', 'native')
        uut.add_string('absolute', bounds=['5'])
        uut.add_value('inlist', 'integer', bounds=['lsize'])
        uut.add_value('outlist', 'real', bounds=['esize'])
        uut.add_value('unknown', 'integer', bounds=[':'])
        uut.write_module(output_file)

        assert output_file.getvalue() \
            == self._EXPECTED_MODULE_WRITE_ARRAY_SOURCE


class TestNamelistConfigDescription():
    @staticmethod
    def _compile_dictionary(namelists):

        dictionary = {}

        for namelist in namelists:

            parameters = {}
            for parameter in namelist.get_parameters():

                parameters[parameter.name] = \
                    [parameter.fortran_type.intrinsic_type,
                     parameter.fortran_type.kind]
                if parameter.get_configure_type() == 'enumeration':
                    parameters[parameter.name].extend(
                        list(parameter.mapping.keys()))
                elif parameter.get_configure_type() == 'computed':
                    parameters[parameter.name].append(
                        parameter.computation)
                elif parameter.get_configure_type() == 'array':
                    if isinstance(parameter.bounds, list):
                        bounds = parameter.bounds[0]
                    else:
                        bounds = parameter.bounds
                    parameters[parameter.name].append('(' + bounds + ')')

            dictionary[namelist.get_namelist_name()] = parameters

        return dictionary

    def test_parser_good_file(self, tmp_path: Path):
        input_file = tmp_path / 'test.ini'
        input_file.write_text(dedent('''
            [namelist:fred]

            [namelist:fred=first_thing]
            type=character

            [namelist:fred=second]
            type=integer

            [namelist:fred=filename]
            type=character
            !string_length=filename

            [namelist:fred=choices]
            !enumeration=true
            values=foo, bar, baz, qux
            '''))

        picker_command = [PICKER_EXE, input_file]
        out = Popen(picker_command, cwd=input_file.parent)
        out.wait()

        json_file = input_file.with_suffix('.json')

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(uut.process_config(json_file))

        assert {'fred': {'choices':     ['integer', 'i_native',
                                         'foo', 'bar', 'baz', 'qux'],
                         'filename':    ['character', 'str_max_filename'],
                         'first_thing': ['character', 'str_def'],
                         'second':      ['integer', 'i_def']}} == result

    def test_only_enumeration(self, tmp_path: Path):
        input_file = tmp_path / 'input.ini'
        input_file.write_text(dedent('''
            [namelist:barney]

            [namelist:barney=stuff]
            !enumeration=true
            values=one, two, three
            '''))

        picker_command = [PICKER_EXE, str(input_file)]
        out = Popen(picker_command, cwd=input_file.parent)
        out.wait()

        json_file = input_file.with_suffix('.json')

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(uut.process_config(json_file))

        assert {'barney': {'stuff': ['integer', 'i_native',
                                     'one', 'two', 'three']}} == result

    _EXPECTED_NON_RANDOM_ENUMERATOR_SOURCE = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the telly namelist.
!>
module telly_config_mod

  use constants_mod, only: i_native, &
                           str_def
  use log_mod,       only: log_event, log_scratch_space &
                         , LOG_LEVEL_ERROR, LOG_LEVEL_WARNING, LOG_LEVEL_INFO
  use mpi_mod,       only: broadcast
  use mpi,           only: MPI_SUCCESS

  use constants_mod, only: cmdi, emdi, imdi, rmdi, unset_key

  implicit none

  private
  public :: tubbies_from_key, key_from_tubbies, &
            read_telly_namelist, postprocess_telly_namelist, &
            telly_is_loadable, telly_is_loaded, telly_final

  integer(i_native), public, parameter :: tubbies_inky = 3
  integer(i_native), public, parameter :: tubbies_lala = 1
  integer(i_native), public, parameter :: tubbies_po = 2

  integer(i_native), public, protected :: tubbies = emdi

  logical :: namelist_loaded = .false.

  character(str_def), parameter :: tubbies_key(3) &
          = [character(len=str_def) :: 'inky', &
                                       'lala', &
                                       'po']

  integer(i_native), parameter :: tubbies_value(3) &
          = [3_i_native, &
             1_i_native, &
             2_i_native]

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_native) function tubbies_from_key( key )

    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    if (key == unset_key) then
      write( log_scratch_space, '(A)') &
          'Missing key for tubbies enumeration in telly namelist.'
      tubbies_from_key = emdi
      call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      return
    end if

    key_index = 1
    do
      if (trim(tubbies_key(key_index)) == trim(key)) then
        tubbies_from_key = tubbies_value(key_index)
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(tubbies_key, 1)) then
          write( log_scratch_space, &
              '("Key ''", A, "'' not recognised for telly tubbies")' ) &
              trim(adjustl(key))
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function tubbies_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> @param[in] value Enumeration value.
  !>
  character(str_def) function key_from_tubbies( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: value_index

    value_index = 1
    do
      if (tubbies_value(value_index) == emdi) then
        key_from_tubbies = unset_key
        return
      else if (tubbies_value(value_index) == value) then
        key_from_tubbies = tubbies_key(value_index)
        return
      else
        value_index = value_index + 1
        if (value_index > ubound(tubbies_key, 1)) then
          write( log_scratch_space, &
                 '("Value ", I0, " is not in telly tubbies")' ) value
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function key_from_tubbies

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_telly_namelist( file_unit, local_rank )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    call read_namelist( file_unit, local_rank, &
                        tubbies )

  end subroutine read_telly_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, &
                            dummy_tubbies )

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_native), intent(out) :: dummy_tubbies

    integer(i_native) :: buffer_integer_i_native(1)

    character(str_def) :: tubbies

    namelist /telly/ tubbies

    integer(i_native) :: condition

    tubbies = unset_key

    if (local_rank == 0) then

      read( file_unit, nml=telly, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      dummy_tubbies = tubbies_from_key( tubbies )

    end if

    buffer_integer_i_native(1) = dummy_tubbies

    call broadcast( buffer_integer_i_native, 1, 0 )

    dummy_tubbies = buffer_integer_i_native(1)



    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_telly_namelist()

    implicit none


  end subroutine postprocess_telly_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function telly_is_loadable()

    implicit none

    logical :: telly_is_loadable

    telly_is_loadable = .not. namelist_loaded

  end function telly_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function telly_is_loaded()

    implicit none

    logical :: telly_is_loaded

    telly_is_loaded = namelist_loaded

  end function telly_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine telly_final()

    implicit none

    tubbies = emdi

    return
  end subroutine telly_final


end module telly_config_mod
    '''.strip() + '\n'

    def test_non_randomized_enumerations(self):
        output_file = io.StringIO()

        uut = description.NamelistDescription('telly', fix_enum=True)
        uut.add_enumeration('tubbies', enumerators=['lala', 'po', 'inky'])
        uut.write_module(output_file)

        assert output_file.getvalue() \
            == self._EXPECTED_NON_RANDOM_ENUMERATOR_SOURCE

    def test_non_enumeration_no_type(self, tmp_path: Path):
        input_file = tmp_path / 'test.init'
        input_file.write_text(dedent('''
            [namelist:barney]

            [namelist:barney=stuff]
            values=one, two, three
            '''))

        picker_command = [PICKER_EXE, input_file]
        out = Popen(picker_command, cwd=input_file.parent)
        out.wait()

        json_file = input_file.with_suffix('.json')

        uut = description.NamelistConfigDescription()
        with pytest.raises(description.NamelistDescriptionException):
            uut.process_config(json_file)

    def test_no_member_type(self, tmp_path: Path):
        input_file = tmp_path / 'input.init'
        input_file.write_text(dedent('''
            [namelist:santa]

            [namelist:santa=elf]
            length=:
            '''))

        picker_command = [PICKER_EXE, str(input_file)]
        out = Popen(picker_command, cwd=input_file.parent)
        out.wait()

        json_file = input_file.with_suffix('.json')

        uut = description.NamelistConfigDescription()
        with pytest.raises(description.NamelistDescriptionException):
            uut.process_config(json_file)

    def test_computed_fields(self, tmp_path: Path):
        input_file = tmp_path / 'input.ini'
        input_file.write_text(dedent('''
            [namelist:teapot]

            [namelist:teapot=foo]
            type=real

            [namelist:teapot=fum]
            type=real

            [!namelist:teapot=bar]
            type=real
            expression=namelist:teapot=foo ** 2

            [!namelist:teapot=baz]
            type=real
            expression=source:constants_mod=PI * foo

            [!namelist:teapot=dosh]
            type=real
            expression=namelist:fridge=milk + (namelist:teapot=foo ** 2) - (source:constants_mod=PI * namelist:teapot=fum)
            '''))

        picker_command = [PICKER_EXE, input_file]
        out = Popen(picker_command, cwd=input_file.parent)
        out.wait()

        json_file = input_file.with_suffix('.json')

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(uut.process_config(json_file))

        assert {'teapot': {'foo': ['real', 'r_def'],
                           'fum': ['real', 'r_def'],
                           'bar': ['real', 'r_def', 'foo ** 2'],
                           'baz': ['real', 'r_def', 'PI * foo'],
                           'dosh': ['real', 'r_def',
                                    'milk + (foo ** 2) - (PI * fum)']}} \
            == result


    def test_constant_in_computed(self, tmp_path: Path):
        input_file = tmp_path / 'input.ini'
        input_file.write_text(dedent('''
            [namelist:cheese]

            [namelist:cheese=fred]
            type=real

            [!namelist:cheese=wilma]
            type=real
            expression=fred * source:constants_mod=FUDGE
            '''))

        picker_command = [PICKER_EXE, str(input_file)]
        out = Popen(picker_command, cwd=input_file.parent)
        out.wait()

        json_file = input_file.with_suffix('.json')

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(uut.process_config(json_file))

        assert {'cheese': {'fred':  ['real', 'r_def'],
                'wilma': ['real', 'r_def', 'fred * FUDGE']}} == result

    def test_array_fields(self, tmp_path: Path):
        input_file = tmp_path / 'input.ini'
        input_file.write_text(dedent('''
            [namelist:aerial]

            [namelist:aerial=fred]
            type=real

            [namelist:aerial=wilma]
            type=real
            length=:
            !bounds=source:constants_mod=FUDGE

            [namelist:aerial=betty]
            type=logical
            length=:
            !bounds=fred

            [namelist:aerial=dino]
            type=integer
            length=:
            !bounds=namelist:sugar=TABLET

            [namelist:aerial=bambam]
            type=integer
            length=:
            '''))

        picker_command = [PICKER_EXE, str(input_file)]
        out = Popen(picker_command, cwd=input_file.parent)
        out.wait()

        json_file = input_file.with_suffix('.json')

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(
            uut.process_config(json_file))

        assert result == \
            {'aerial': {'bambam': ['integer', 'i_def', '(:)'],
                        'betty':  ['logical', 'l_def', '(fred)'],
                        'fred':   ['real', 'r_def'],
                        'wilma':  ['real', 'r_def', '(FUDGE)'],
                        'dino':   ['integer', 'i_def', '(TABLET)']}}
