#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################

from pathlib import Path

from pytest import fixture, raises

from dependerator.database import (
    DatabaseException,
    FileDependencies,
    FortranDependencies,
    SQLiteDatabase,
)


class TestDatabase:
    def test_all(self, tmp_path: Path):
        """
        Exercise basic database functions of the SQLite wrapper.
        """
        uut = SQLiteDatabase(tmp_path / "test.db")

        uut.ensure_table("test", ["'alpha' integer primary key"])

        result = uut.query("select * from test")
        assert 0 == len(result)

        result = uut.query("insert into 'test' values (1)")
        assert 0 == len(result)

        result = uut.query("select * from test")
        assert 1 == len(result)
        assert 1 == result[0][0]


class TestFileDependency:
    def test_all(self, tmp_path: Path):
        """
        Ensure file dependencies are correctly managed.
        """
        database = SQLiteDatabase(tmp_path / "file.db")
        uut = FileDependencies(database)

        result = uut.get_dependencies()
        assert [] == list(result)

        uut.add_file_dependency("foo.f90", "bar")
        result = uut.get_dependencies()
        assert [(Path("foo.f90"), [Path("bar")])] == list(result)

        uut.add_file_dependency(Path("foo.f90"), Path("baz"))
        uut.add_file_dependency(Path("qux.f90"), Path("bar"))
        result = uut.get_dependencies()
        assert [
            (Path("foo.f90"), [Path("bar"), Path("baz")]),
            (Path("qux.f90"), [Path("bar")]),
        ] == list(result)

        uut.remove_file("foo.f90")
        result = uut.get_dependencies()
        assert [(Path("qux.f90"), [Path("bar")])] == list(result)

        uut.remove_all_file_dependencies()
        result = uut.get_dependencies()
        assert [] == list(result)


class TestFortranDependency:
    @fixture
    def example_db(self, tmp_path: Path) -> FortranDependencies:
        database = FortranDependencies(SQLiteDatabase(tmp_path / "fortran.db"))
        database.add_program("foo", Path("foo.f90"))
        database.add_module("noonoo", Path("foo.f90"))
        database.add_module("bar", Path("bar.f90"))
        database.add_module("baz", Path("baz.f90"))
        database.add_module("qux", Path("qux.f90"))
        database.add_program("fred", Path("fred.f90"))
        database.add_module("wilma", Path("wilma.f90"))

        database.add_compile_dependency("foo", "bar")
        database.add_compile_dependency("foo", "baz")
        database.add_compile_dependency("foo", "noonoo")
        database.add_compile_dependency("bar", "qux")
        database.add_compile_dependency("noonoo", "wilma")
        database.add_compile_dependency("fred", "wilma")

        database.add_link_dependency("bar", "foo")
        database.add_link_dependency("baz", "foo")
        database.add_link_dependency("noonoo", "foo")
        database.add_link_dependency("qux", "bar")
        database.add_link_dependency("wilma", "fred")
        database.add_link_dependency("wilma", "noonoo")

        return database

    def test_add_program(self, tmp_path: Path):
        """
        Ensure that a program can be added and correctly retrieved.
        """
        database = SQLiteDatabase(tmp_path / "fortran.db")
        uut = FortranDependencies(database)
        uut.add_program("foo", Path("bar.f90"))

        result = uut.get_programs()
        assert ["foo"] == result

        # Should test that filename is correctly stored.

    def test_get_compile_prerequisites(self, example_db: FortranDependencies):
        assert sorted(example_db.get_compile_prerequisites("foo")) == [
            "bar",
            "baz",
            "noonoo",
        ]
        assert sorted(example_db.get_compile_prerequisites("fred")) == [
            "wilma"
        ]

    def test_get_compile_prerequisites_sub_tree(
        self, example_db: FortranDependencies
    ):
        assert sorted(example_db.get_compile_prerequisites("bar")) == ["qux"]

    def test_remove_source_file(self, example_db: FortranDependencies):
        """
        Ensure that file is removed correctly. References to it remain as they
        do in the source.
        """
        example_db.remove_file(Path("bar.f90"))

        assert list(example_db.get_programs()) == ["fred", "foo"]

        assert example_db.get_compile_prerequisites("foo") == [
            "bar",
            "baz",
            "noonoo",
        ]
        assert example_db.get_compile_prerequisites("fred") == ["wilma"]
        assert example_db.get_compile_prerequisites("noonoo") == ["wilma"]

        assert sorted(list(example_db.get_link_dependencies("baz"))) == [
            ("baz", Path("baz.f90"), "foo", Path("foo.f90"))
        ]
        assert sorted(list(example_db.get_link_dependencies("wilma"))) == [
            ("noonoo", Path("foo.f90"), "foo", Path("foo.f90")),
            ("wilma", Path("wilma.f90"), "fred", Path("fred.f90")),
            ("wilma", Path("wilma.f90"), "noonoo", Path("foo.f90")),
        ]

    def test_get_programs(self, example_db: FortranDependencies):
        """
        Ensure the list of programs is correctly retrieved.
        """
        programs = example_db.get_programs()
        assert ["fred", "foo"] == list(programs)

    def test_get_modules(self, example_db: FortranDependencies):
        """
        Ensure the list of modules is correctly retrieved.
        """
        modules = example_db.get_modules()
        assert [
            ("noonoo", Path("foo.f90")),
            ("bar", Path("bar.f90")),
            ("baz", Path("baz.f90")),
            ("qux", Path("qux.f90")),
            ("wilma", Path("wilma.f90")),
        ] == list(modules)

    def test_get_all_compile_dependencies(
        self, example_db: FortranDependencies
    ):
        """
        Ensure all file dependencies are correctly returned.
        """
        assert sorted(example_db.get_compile_dependencies()) == [
            (
                "bar",
                Path("bar.f90"),
                "module",
                "qux",
                Path("qux.f90"),
                "module",
            ),
            (
                "foo",
                Path("foo.f90"),
                "program",
                "bar",
                Path("bar.f90"),
                "module",
            ),
            (
                "foo",
                Path("foo.f90"),
                "program",
                "baz",
                Path("baz.f90"),
                "module",
            ),
            (
                "foo",
                Path("foo.f90"),
                "program",
                "noonoo",
                Path("foo.f90"),
                "module",
            ),
            (
                "fred",
                Path("fred.f90"),
                "program",
                "wilma",
                Path("wilma.f90"),
                "module",
            ),
            (
                "noonoo",
                Path("foo.f90"),
                "module",
                "wilma",
                Path("wilma.f90"),
                "module",
            ),
        ]

        assert [("bar", Path("bar.f90"), "foo", Path("foo.f90"))] == list(
            example_db.get_link_dependencies("bar")
        )
        assert [("baz", Path("baz.f90"), "foo", Path("foo.f90"))] == list(
            example_db.get_link_dependencies("baz")
        )
        assert sorted(list(example_db.get_link_dependencies("qux"))) == [
            ("bar", Path("bar.f90"), "foo", Path("foo.f90")),
            ("qux", Path("qux.f90"), "bar", Path("bar.f90")),
        ]
        assert sorted(list(example_db.get_link_dependencies("wilma"))) == [
            ("noonoo", Path("foo.f90"), "foo", Path("foo.f90")),
            ("wilma", Path("wilma.f90"), "fred", Path("fred.f90")),
            ("wilma", Path("wilma.f90"), "noonoo", Path("foo.f90")),
        ]

    def test_get_some_compile_dependencies(
        self, example_db: FortranDependencies
    ):
        """
        Ensures a sub-set of compile dependencies are correctly returned.
        """
        assert sorted(example_db.get_compile_dependencies("bar")) == [
            (
                "bar",
                Path("bar.f90"),
                "module",
                "qux",
                Path("qux.f90"),
                "module",
            )
        ]

    def test_get_submod_compile_dependencies(self, tmp_path: Path):
        """
        Ensures dependencies around a sub-module are correctly handled.
        """
        database = FortranDependencies(SQLiteDatabase(tmp_path / "fortran.db"))
        database.add_program("prog", Path("prog.f90"))
        database.add_module("parent", Path("parent.f90"))
        database.add_submodule("child1", Path("parent.f90"))
        database.add_submodule("child2", Path("child2.f90"))

        database.add_compile_dependency("prog", "parent")
        database.add_compile_dependency("child1", "parent")
        database.add_compile_dependency("child2", "parent")

        database.add_link_dependency("prog", "parent")
        database.add_link_dependency("parent", "child1")
        database.add_link_dependency("parent", "child2")

        assert sorted(database.get_compile_dependencies()) == [
            (
                "child1",
                Path("parent.f90"),
                "submodule",
                "parent",
                Path("parent.f90"),
                "module",
            ),
            (
                "child2",
                Path("child2.f90"),
                "submodule",
                "parent",
                Path("parent.f90"),
                "module",
            ),
            (
                "prog",
                Path("prog.f90"),
                "program",
                "parent",
                Path("parent.f90"),
                "module",
            ),
        ]

    def test_compile_dependencies_missing_prerequisite(self, tmp_path: Path):
        """
        Checks that the correct error is raised and message displayed if a
        prerequisite cannot be found for a module.
        """
        database = FortranDependencies(SQLiteDatabase(tmp_path / "fortran.db"))

        database.add_program("prog", Path("prog.f90"))
        database.add_module("parent", Path("parent.f90"))
        database.add_submodule("child2", Path("child2.f90"))

        database.add_compile_dependency("prog", "parent")
        database.add_compile_dependency("child2", "child1")

        with raises(DatabaseException) as err:
            sorted(database.get_compile_dependencies("child1"))
        # Assert prerequisite and units
        assert "Unable to find prerequisite 'child1' of 'child2'" in str(
            err.value
        )

    @staticmethod
    def test_duplicate_module(example_db: FortranDependencies):
        """
        Ensure that an exception is thrown if an attempt is made to add a
        module which is already in the database.
        """
        with raises(DatabaseException) as caught:
            example_db.add_module("qux", Path("cheese/qux.f90"))
        assert "qux" == caught.value.module
        assert Path("cheese/qux.f90") == caught.value.filename
