#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Unit tests for Extended Namelist Specific Object generator.
"""

from pathlib import Path

import configurator.extended_namelist_type as ExtendedNml

HERE = Path(__file__).resolve().parent


class TestExtendedNml:
    """
    Tests generation of extended namelist specific object.
    """

    def test_write_one_of_each(self, tmp_path: Path):
        # pylint: disable=no-self-use
        """
        Generating extended namelist object with one of each
        component member type.
        """
        output_file = tmp_path / "one_of_each_mod.f90"
        uut = ExtendedNml.NamelistDescription("one_of_each")

        uut.add_value("vint", "integer")
        uut.add_value("dint", "integer", "default")
        uut.add_value("sint", "integer", "short")
        uut.add_value("lint", "integer", "long")
        uut.add_value("dlog", "logical", "default")
        uut.add_value("vreal", "real")
        uut.add_value("dreal", "real", "default")
        uut.add_value("sreal", "real", "single")
        uut.add_value("lreal", "real", "double")
        uut.add_value("treal", "real", "second")
        uut.add_string("vstr")
        uut.add_string("dstr", configure_string_length="default")
        uut.add_string("fstr", configure_string_length="filename")
        uut.add_enumeration("enum", enumerators=["one", "two", "three"])

        uut.write_module(output_file)

        expected_file = HERE / "one_each_mod.f90"
        assert (
            expected_file.read_text(encoding="ascii")
            == output_file.read_text(encoding="ascii") + "\n"
        )
