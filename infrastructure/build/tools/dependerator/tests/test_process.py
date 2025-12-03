#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################

from pathlib import Path
from typing import Tuple

import pytest

from dependerator.database import (
    FileDependencies,
    FortranDependencies,
    SQLiteDatabase,
)
from dependerator.process import FortranProcessor


class TestFortranProcessor:
    @pytest.fixture
    def databases(
        self, tmp_path: Path
    ) -> Tuple[FortranDependencies, FileDependencies]:
        """
        Creates an example dependencies database.
        """
        filename = tmp_path / "fortran.db"
        database = SQLiteDatabase(filename)
        fortran_db = FortranDependencies(database)

        fortran_db.add_program("foo", Path("foo.f90"))
        fortran_db.add_module("bar", Path("bits/bar.f90"))
        fortran_db.add_module("baz", Path("bits/baz.f90"))
        fortran_db.add_module("qux", Path("bobs/qux.f90"))
        fortran_db.add_module("corge", Path("bobs/grault.f90"))
        fortran_db.add_procedure("quux", Path("quux.f90"))
        fortran_db.add_program("fred", Path("fred.f90"))

        fortran_db.add_compile_dependency("foo", "bar")
        fortran_db.add_compile_dependency("bar", "baz")
        fortran_db.add_compile_dependency("qux", "baz")
        fortran_db.add_compile_dependency("corge", "bar")
        fortran_db.add_compile_dependency("foo", "quux")
        fortran_db.add_compile_dependency("fred", "quux")

        fortran_db.add_link_dependency("foo", "bar")
        fortran_db.add_link_dependency("bar", "baz")
        fortran_db.add_link_dependency("baz", "qux")
        fortran_db.add_link_dependency("corge", "bar")
        fortran_db.add_link_dependency("foo", "quux")
        fortran_db.add_link_dependency("fred", "quux")

        return fortran_db, FileDependencies(database)

    def test_compile_dependencies(self, databases):
        """
        Ensures a full list of compile-time dependencies can be generated.
        """
        uut = FortranProcessor(databases[0], Path("objects"), Path("modules"))
        uut.determine_compile_file_dependencies(databases[1])

        assert list(databases[1].get_dependencies()) == [
            (Path("modules/bits/bar.mod"), [Path("objects/bits/bar.o")]),
            (Path("modules/bits/baz.mod"), [Path("objects/bits/baz.o")]),
            (Path("modules/bobs/corge.mod"), [Path("objects/bobs/grault.o")]),
            (Path("modules/bobs/qux.mod"), [Path("objects/bobs/qux.o")]),
            (Path("objects/bits/bar.o"), [Path("modules/bits/baz.mod")]),
            (
                Path("objects/foo.o"),
                [Path("modules/bits/bar.mod"), Path("objects/quux.o")],
            ),
            (Path("objects/fred.o"), [Path("objects/quux.o")]),
        ]

    def test_compile_dependencies_module_objects(self, databases):
        """
        Ensures that module information in object files is supported.
        """
        uut = FortranProcessor(databases[0], Path("objects"), Path("modules"))
        uut.determine_compile_file_dependencies(
            databases[1], object_modules=True
        )

        assert list(databases[1].get_dependencies()) == [
            (Path("objects/bits/bar.o"), [Path("objects/bits/baz.o")]),
            (
                Path("objects/foo.o"),
                [Path("objects/bits/bar.o"), Path("objects/quux.o")],
            ),
            (Path("objects/fred.o"), [Path("objects/quux.o")]),
        ]

    def test_link_dependencies_programs(self, databases):
        """
        Ensures a list of all objects per program can be fetched.
        """
        uut = FortranProcessor(databases[0], Path("objects"), Path("modules"))
        result = list(uut.determine_link_dependencies())

        assert result == [
            (
                Path("objects/fred"),
                Path("objects/fred.o"),
                [Path("objects/quux.o")],
            ),
            (
                Path("objects/foo"),
                Path("objects/foo.o"),
                [
                    Path("objects/bits/bar.o"),
                    Path("objects/bits/baz.o"),
                    Path("objects/bobs/qux.o"),
                    Path("objects/quux.o"),
                ],
            ),
        ]

    def test_link_dependencies_sub_tree(self, databases):
        """
        Tests getting a link set for a subtree (or library).
        """
        test_unit = FortranProcessor(databases[0], Path("obj"), Path("mod"))
        assert list(test_unit.determine_link_dependencies("corge")) == [
            (
                Path("obj/corge"),
                Path("obj/bobs/grault.o"),
                [
                    Path("obj/bits/bar.o"),
                    Path("obj/bits/baz.o"),
                    Path("obj/bobs/qux.o"),
                ],
            )
        ]
