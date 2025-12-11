#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
from pathlib import Path
from textwrap import dedent

import pytest

from dependerator.analyser import FortranAnalyser
from dependerator.database import FortranDependencies, SQLiteDatabase


class TestFortranAnalyser:
    @pytest.fixture
    def database(self, tmp_path_factory):
        filename = tmp_path_factory.mktemp("db-", True) / "test.db"
        database = SQLiteDatabase(filename)
        return FortranDependencies(database)

    def test_continuation_lines(self, database, tmp_path: Path):
        """
        Ensure continuation lines are handled correctly.
        """
        database.add_program("stock", Path("oxo.f90"))
        database.add_compile_dependency("stock", "widgit")
        database.add_link_dependency("stock", "widgit")
        database.add_compile_dependency("stock", "thingy")
        database.add_link_dependency("stock", "thingy")
        database.add_module("beef", Path("cow.f90"))
        database.add_module("pork", Path("pig.f90"))

        test_filename = tmp_path / "cont.f90"
        test_filename.write_text(
            dedent("""
            subroutine thingy(cheese, meat, &
                              teapot, fishslice, &
                             )

            end subroutine thingy

            subroutine widgit(grunk, &
                              ! This comment should be ignored
                              grank)
              use beef
              implicit none
            end subroutine widgit

            subroutine old_school(bibble, &
                                & bobble)
              use pork
              implicit none
            end subroutine old_school
            """)
        )

        uut = FortranAnalyser([], database)
        uut.analyse(test_filename)

        assert sorted(database.get_program_units()) == [
            ("beef", Path("cow.f90")),
            ("old_school", test_filename),
            ("pork", Path("pig.f90")),
            ("stock", Path("oxo.f90")),
            ("thingy", test_filename),
            ("widgit", test_filename),
        ]

        assert sorted(database.get_compile_dependencies()) == [
            (
                "stock",
                Path("oxo.f90"),
                "program",
                "thingy",
                test_filename,
                "procedure",
            ),
            (
                "stock",
                Path("oxo.f90"),
                "program",
                "widgit",
                test_filename,
                "procedure",
            ),
            (
                "widgit",
                test_filename,
                "procedure",
                "beef",
                Path("cow.f90"),
                "module",
            ),
        ]

        dependencies = list(database.get_link_dependencies("widgit"))
        assert [
            ("widgit", test_filename, "beef", Path("cow.f90"))
        ] == dependencies

    def test_exclamation_in_string(self, database, tmp_path: Path):
        """
        Ensure an exclamation mark inside a string doesn't end the line.
        """
        test_filename = tmp_path / "cont.f90"
        test_filename.write_text(
            dedent("""
            call log_event("Gungho: " // "stopping program! ", &
            &LOG_LEVEL_ERROR)
            """)
        )

        uut = FortranAnalyser([], database)
        uut.analyse(test_filename)

        programs = list(database.get_programs())
        assert [] == programs

        dependencies = list(database.get_compile_dependencies())
        assert [] == dependencies

    def test_procedure(self, database, tmp_path: Path):
        """
        Procedure as program unit.
        """
        test_filename = tmp_path / "test.f90"
        test_filename.write_text(
            dedent("""
           subroutine empty_sub()
           end subroutine empty_sub

           subroutine one_sub(argument)
           end subroutine one_sub

           real function cosd(degree)
           end function cosd

           function alternative() result(thing)
           end function alternative
           """)
        )

        uut = FortranAnalyser([], database)
        uut.analyse(test_filename)

        assert sorted(database.get_program_units()) == [
            ("alternative", test_filename),
            ("cosd", test_filename),
            ("empty_sub", test_filename),
            ("one_sub", test_filename),
        ]

        assert sorted(database.get_compile_dependencies()) == []

    def test_analyse_program(self, database, tmp_path: Path):
        """
        Includes disparate case to ensure case insensitivity.
        """
        database.add_module("constants_mod", Path("constants_mod.f90"))
        database.add_module("trumpton_mod", Path("trumpton_mod.f90"))

        test_filename = tmp_path / "test.f90"
        test_filename.write_text(
            dedent("""
            program fOo
              use constAnts_mod, only : i_def
              use trumpton_Mod, only : hew, pew, barney, mcgrey, &
                                       cuthbirt, dibble, grub
              implicit none
            end program fOo
            """)
        )

        uut = FortranAnalyser([], database)
        uut.analyse(test_filename)

        programs = list(database.get_programs())
        assert ["foo"] == programs

        dependencies = list(database.get_compile_dependencies())
        assert [
            (
                "foo",
                test_filename,
                "program",
                "constants_mod",
                Path("constants_mod.f90"),
                "module",
            ),
            (
                "foo",
                test_filename,
                "program",
                "trumpton_mod",
                Path("trumpton_mod.f90"),
                "module",
            ),
        ] == sorted(dependencies)

        dependencies = list(database.get_link_dependencies("foo"))
        assert [
            ("foo", test_filename, "constants_mod", Path("constants_mod.f90")),
            ("foo", test_filename, "trumpton_mod", Path("trumpton_mod.f90")),
        ] == sorted(dependencies)

    def test_analyse_module(self, database, tmp_path: Path):
        """
        Includes disparate case to ensure case insensitivity.
        """
        uut = FortranAnalyser([], database)

        test_filename = tmp_path / "test.f90"
        test_filename.write_text(
            dedent("""
            module foO
              ! Ignore this
              use consTants_mod, only : i_def
              use trumPton_mod, only : hew, pew, barney, mcgrey, &
                                       cuthbirt, dibble, grub
              implicit none
              private
            contains
          end module foO
          module truMpton_mod
          end module truMpton_mod
          """)
        )
        uut.analyse(test_filename)

        other_filename = tmp_path / "other.f90"
        other_filename.write_text(
            dedent("""
            module coNstants_mod
            end module coNstants_mod
            """)
        )
        uut.analyse(other_filename)

        assert sorted(database.get_program_units()) == [
            ("constants_mod", other_filename),
            ("foo", test_filename),
            ("trumpton_mod", test_filename),
        ]

        assert sorted(database.get_compile_dependencies("foo")) == [
            (
                "foo",
                test_filename,
                "module",
                "constants_mod",
                other_filename,
                "module",
            ),
            (
                "foo",
                test_filename,
                "module",
                "trumpton_mod",
                test_filename,
                "module",
            ),
        ]

        dependencies = list(database.get_link_dependencies("foo"))
        assert [
            ("foo", test_filename, "constants_mod", other_filename),
            ("foo", test_filename, "trumpton_mod", test_filename),
        ] == sorted(dependencies)

    def test_intra_module(self, database: FortranDependencies, tmp_path: Path):
        """
        Ensures Intra-module dependencies are handled correctly.
        """
        external_path = tmp_path / "external_mod.f90"
        database.add_module("external_mod", external_path)
        other_path = tmp_path / "other_mod.f90"
        database.add_module("other_mod", other_path)

        test_filename = tmp_path / "beef.f90"
        test_filename.write_text(
            dedent(
                """
                module my_mod
                  use external_mod, only : external_proc
                  use library_mod
                end module

                program beef
                  use my_mod, only: thing_proc
                  use other_mod, only: other_thing
                end program
                """
            )
        )

        test_unit = FortranAnalyser(["library_mod"], database)
        test_unit.analyse(test_filename)

        assert list(database.get_programs()) == ["beef"]
        assert sorted(list(database.get_compile_dependencies())) == [
            (
                "beef",
                test_filename,
                "program",
                "my_mod",
                test_filename,
                "module",
            ),
            (
                "beef",
                test_filename,
                "program",
                "other_mod",
                other_path,
                "module",
            ),
            (
                "my_mod",
                test_filename,
                "module",
                "external_mod",
                external_path,
                "module",
            ),
        ]
        assert sorted(list(database.get_link_dependencies("beef"))) == [
            ("beef", test_filename, "my_mod", test_filename),
            ("beef", test_filename, "other_mod", other_path),
            ("my_mod", test_filename, "external_mod", external_path),
        ]

    def test_analyse_submodule(self, database, tmp_path):
        """
        This test also includes disparate case to ensure case insensitivity is
        enforced.
        """
        uut = FortranAnalyser([], database)

        parent_filename = tmp_path / "parent.f90"
        parent_filename.write_text(
            dedent("""
            module Parent
              implicit none
              private
              type, public :: test_type
              contains
                procedure foo
                procedure bar
                procedure baz
              end type test_type
              interface
                module subroutine foo(this, cheese)
                  class(test_type), intent(inout) :: this
                  real,             intent(in)    :: cheese
                end subroutine foo
                module subroutine bar(this, teapot)
                  class(test_type), intent(inout) :: this
                  character(*),     intent(in)    :: teapot
                end subroutine bar
                type(something) module function baz(this)
                  class(test_type), intent(in) :: this
                end function baz
              end interface
            end module Parent
            submodule (pArent) chIld3
              implicit none
            contains
              type(something) module function baz(this)
                class(test_type), intent(in) :: this
              end function baz
            end submodule chIld3
            """)
        )
        uut.analyse(parent_filename)

        child1_filename = tmp_path / "child1.f90"
        child1_filename.write_text(
            dedent("""
            submodule(paRent) Child1
              implicit none
              type :: secondary_type
              contains
                procedure baz
              end type secondary_type
              interface
                module subroutine baz(this, wibble)
                  class(secondary_type), intent(inout) :: this
                  real,                  intent(in)    :: wibble
                end subroutine baz
              end interface
            contains
              module subroutine foo(this, cheese)
                implicit none
                class(test_type), intent(inout) :: this
                real,             intent(in)    :: cheese
                type(secondary_type) :: thang
                thang = secondary_type()
                write(6, *) cheese
                call thang%baz(12.7)
              end subroutine foo
            end submodule Child1
            """)
        )
        uut.analyse(child1_filename)

        child2_filename = tmp_path / "child2.f90"
        child2_filename.write_text(
            dedent("""
            submodule (parent) cHild2
              implicit none
            contains
              module procedure bar
                implicit none
                write(6, *) teapot
              end procedure bar
            end submodule cHild2
            """)
        )
        uut.analyse(child2_filename)

        child3_filename = tmp_path / "child3.f90"
        child3_filename.write_text(
            dedent("""
            submodule (parEnt:chilD1) grandChild
              implicit none
            contains
              module subroutine baz(this, wibble)
                implicit none
                class(secondary_type), intent(inout) :: this
                real,                  intent(in)    :: wibble
                write(6, *) wibble
              end subroutine baz
            end submodule grandChild
            """)
        )
        uut.analyse(child3_filename)

        assert sorted(database.get_program_units()) == [
            ("child1", child1_filename),
            ("child2", child2_filename),
            ("child3", parent_filename),
            ("grandchild", child3_filename),
            ("parent", parent_filename),
        ]

        assert sorted(database.get_compile_dependencies("parent")) == [
            (
                "child1",
                child1_filename,
                "submodule",
                "parent",
                parent_filename,
                "module",
            ),
            (
                "child2",
                child2_filename,
                "submodule",
                "parent",
                parent_filename,
                "module",
            ),
            (
                "child3",
                parent_filename,
                "submodule",
                "parent",
                parent_filename,
                "module",
            ),
            (
                "grandchild",
                child3_filename,
                "submodule",
                "child1",
                child1_filename,
                "submodule",
            ),
        ]

        assert sorted(list(database.get_link_dependencies("parent"))) == [
            ("child1", child1_filename, "grandchild", child3_filename),
            ("parent", parent_filename, "child1", child1_filename),
            ("parent", parent_filename, "child2", child2_filename),
            ("parent", parent_filename, "child3", parent_filename),
        ]

        assert sorted(list(database.get_link_dependencies("child1"))) == [
            ("child1", child1_filename, "grandchild", child3_filename)
        ]

    def test_function_in_module_name(self, database, tmp_path: Path):
        """
        Ensure the analyser isn't tripped up by naked global level procedures
        as program units.
        """
        uut = FortranAnalyser([], database)

        test_filename = tmp_path / "test.f90"
        test_filename.write_text(
            dedent("""
            module function_mod
              use constants_mod, only : i_def
              implicit none
              private
           contains
           end module function_mod
           """)
        )
        uut.analyse(test_filename)

        other_filename = tmp_path / "other.f90"
        other_filename.write_text(
            dedent("""
            module constants_mod
            end module constants_mod
            """)
        )
        uut.analyse(other_filename)

        assert sorted(database.get_program_units()) == [
            ("constants_mod", other_filename),
            ("function_mod", test_filename),
        ]

        assert sorted(database.get_compile_dependencies("function_mod")) == [
            (
                "function_mod",
                test_filename,
                "module",
                "constants_mod",
                other_filename,
                "module",
            )
        ]

        dependencies = list(database.get_link_dependencies("function_mod"))
        assert [
            ("function_mod", test_filename, "constants_mod", other_filename)
        ] == sorted(dependencies)

    def test_depends_on(self, database, tmp_path: Path):
        """
        The analyser has to be able to track dependencies using the deprecated
        "depends on:" comments of the UM.
        """
        database.add_procedure("flibble", "flibble.f90")
        test_filename = tmp_path / "test.f90"
        test_filename.write_text(
            dedent("""
            module some_mod

              use constants_mod, only : i_def

              implicit none

              ! Add in an interface block - this will test to make sure
              ! we don't pick up a spurious subroutine call
              interface
                subroutine wooble ()
                end subroutine
              end interface

              private

              ! Comments before the "depends on" shouldn't upset it.

              ! depends on: flibble.o

              ! depends on: wooble

            contains

            end module some_mod
            """)
        )

        other_filename = tmp_path / "other.f90"
        other_filename.write_text(
            dedent("""
            module constants_mod
            contains
              subroutine wooble
              end subroutine wooble
            end module constants_mod
            """)
        )

        depend_filename = tmp_path / "wooble.f90"
        depend_filename.write_text(
            dedent("""
            subroutine wooble
            end subroutine wooble
            """)
        )

        uut = FortranAnalyser([], database)
        uut.analyse(test_filename)
        uut.analyse(other_filename)
        uut.analyse(depend_filename)

        assert sorted(database.get_program_units()) == [
            ("constants_mod", other_filename),
            ("flibble", Path("flibble.f90")),
            ("some_mod", test_filename),
            ("wooble", depend_filename),
        ]
        assert database.get_compile_prerequisites("some_mod") == [
            "constants_mod",
            "flibble",
            "wooble",
        ]

        assert sorted(database.get_compile_dependencies("some_mod")) == [
            (
                "some_mod",
                test_filename,
                "module",
                "constants_mod",
                other_filename,
                "module",
            ),
            (
                "some_mod",
                test_filename,
                "module",
                "flibble",
                Path("flibble.f90"),
                "procedure",
            ),
            (
                "some_mod",
                test_filename,
                "module",
                "wooble",
                depend_filename,
                "procedure",
            ),
        ]

        dependencies = list(database.get_link_dependencies("some_mod"))
        assert [
            ("some_mod", test_filename, "constants_mod", other_filename),
            ("some_mod", test_filename, "flibble", Path("flibble.f90")),
            ("some_mod", test_filename, "wooble", depend_filename),
        ] == sorted(dependencies)

    def test_abstract_interface(self, database, tmp_path: Path):
        """
        The analyser must ignore abstract interface definitions. These are not
        program units.
        """
        first_filename = tmp_path / "test.f90"
        first_filename.write_text(
            dedent("""
            module first_mod

              implicit none

              private

              abstract interface
                subroutine thing_face()
                  implicit none
                end subroutine thing_face
              end interface

            contains

            end module first_mod
            """)
        )

        second_filename = tmp_path / "test2.f90"
        second_filename.write_text(
            dedent("""
            module second_mod

              implicit none

              private

              abstract interface
                subroutine thing_face()
                  implicit none
                end subroutine thing_face
              end interface

            contains

          end module second_mod
          """)
        )

        uut = FortranAnalyser([], database)
        uut.analyse(first_filename)
        uut.analyse(second_filename)

    def test_external(self, database, tmp_path):
        """
        Ensure "external" works as a dependency marker.
        """
        database.add_module("wibble", Path("wibble.f90"))
        database.add_module("bibble", Path("bibble.f90"))
        database.add_module("ibble", Path("ibble.f90"))
        database.add_module("gribble", Path("gribble.f90"))
        test_filename = tmp_path / "test.f90"
        test_filename.write_text(
            dedent("""
           program boo

             implicit none

             external ibble
             external wibble, bibble, gribble

             call wibble()
             call bibble()

           end program boo
           """)
        )

        uut = FortranAnalyser([], database)
        uut.analyse(test_filename)

        programs = list(database.get_programs())
        assert ["boo"] == programs

        dependencies = list(database.get_compile_dependencies())
        assert [
            (
                "boo",
                test_filename,
                "program",
                "bibble",
                Path("bibble.f90"),
                "module",
            ),
            (
                "boo",
                test_filename,
                "program",
                "gribble",
                Path("gribble.f90"),
                "module",
            ),
            (
                "boo",
                test_filename,
                "program",
                "ibble",
                Path("ibble.f90"),
                "module",
            ),
            (
                "boo",
                test_filename,
                "program",
                "wibble",
                Path("wibble.f90"),
                "module",
            ),
        ] == sorted(dependencies)

        dependencies = list(database.get_link_dependencies("boo"))
        assert [
            ("boo", test_filename, "bibble", Path("bibble.f90")),
            ("boo", test_filename, "gribble", Path("gribble.f90")),
            ("boo", test_filename, "ibble", Path("ibble.f90")),
            ("boo", test_filename, "wibble", Path("wibble.f90")),
        ] == sorted(dependencies)

    def test_openmp_sentinel(self, database, tmp_path):
        """
        Ensure OpenMP "sentinel" markup is handled.
        """
        database.add_module(
            "special_thread_sauce_mod", Path("special_thread_sauce_mod.f90")
        )
        test_filename = tmp_path / "sentinel.f90"
        test_filename.write_text(
            dedent(
                """
                module test_mod
                  !$ use special_thread_sauce_mod, only : thread_sauce
                  implicit none
                end module test_mod
                """
            )
        )
        test_unit = FortranAnalyser([], database)
        test_unit.analyse(test_filename)

        assert database.get_compile_prerequisites("test_mod") == [
            "special_thread_sauce_mod"
        ]

        dependencies = list(database.get_link_dependencies("test_mod"))
        assert sorted(dependencies) == [
            (
                "test_mod",
                test_filename,
                "special_thread_sauce_mod",
                Path("special_thread_sauce_mod.f90"),
            )
        ]
