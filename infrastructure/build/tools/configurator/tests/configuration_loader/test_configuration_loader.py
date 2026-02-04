#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Unit tests configuration loader generator.
"""

from pathlib import Path

import configurator.configurationloader as loader

HERE = Path(__file__).resolve().parent


class TestLoader:
    """
    Tests generation of configuration loader.
    """

    def test_with_content(self, tmp_path: Path):  # pylint: disable=no-self-use
        """
        Generating configuration loader.
        """
        uut = loader.ConfigurationLoader("content_mod")
        uut.add_namelist("foo")
        output_file = tmp_path / "content_mod.f90"
        uut.write_module(output_file)

        expected_file = HERE / "content_mod.f90"
        assert output_file.read_text(
            encoding="ascii"
        ) + "\n" == expected_file.read_text(encoding="ascii")
