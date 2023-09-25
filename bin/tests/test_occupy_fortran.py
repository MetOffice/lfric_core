from pathlib import Path
from textwrap import dedent

from ..modules.occupy_fortran import entry

program_text = dedent("""
    program test_program
      integer, public :: global_var
      real, public, save :: save_var
      type :: test_type
        integer :: type_var
      end type test_type
      type(test_type), public :: global_type
      type(test_type), public, save :: save_type
    end program test_program
    """)


def test_program(tmp_path: Path):
    """
    Ensures that globals may exist in programs. There can only ever be one
    program so globals are acceptable.
    """
    test_file = tmp_path / 'program.f90'
    test_file.write_text(program_text)

    dirty_list, clean_list, not_considered = entry([test_file])

    assert clean_list == [test_file]
    assert dirty_list == []
    assert not_considered == []


module_text = dedent("""
    module test_module
      integer, public :: global_var
      type :: test_type
        integer :: type_var
      end type test_type
      type(test_type), public :: global_type
    contains
      subroutine biscuits()
        integer, save :: save_local
      end subroutine biscuits
    end module test_module
    """)


def test_module(tmp_path: Path):
    """
    Ensures that globals are flagged in modules.
    """
    test_file = tmp_path / 'module.f90'
    test_file.write_text(module_text)

    dirty_list, clean_list, not_considered = entry([test_file])

    assert clean_list == []
    assert not_considered == []
    assert len(dirty_list) == 1
    assert dirty_list[0].filename == test_file
    assert len(dirty_list[0].dirt) == 3
    assert dirty_list[0].dirt[0].line_number == 3
    assert dirty_list[0].dirt[0].fortran_type == 'integer'
    assert dirty_list[0].dirt[0].variable_name == 'global_var'
    assert dirty_list[0].dirt[1].line_number == 7
    assert dirty_list[0].dirt[1].fortran_type == 'test_type'
    assert dirty_list[0].dirt[1].variable_name == 'global_type'
    assert dirty_list[0].dirt[2].line_number == 10
    assert dirty_list[0].dirt[2].fortran_type == 'integer'
    assert dirty_list[0].dirt[2].variable_name == 'save_local'
