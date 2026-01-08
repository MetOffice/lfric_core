#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Unit tests Application Configuration Object generator.
"""

from pathlib import Path

import configurator.config_type as AppConfig

HERE = Path(__file__).resolve().parent


class TestAppConfig:
    """
    Tests generation of application configuration object.
    """

    def test_with_content(self, tmp_path: Path):  # pylint: disable=no-self-use
        """
        Generating application configuration object.
        """
        uut = AppConfig.AppConfiguration("config_mod")
        uut.add_namelist("foo", duplicate=False)
        uut.add_namelist("bar", duplicate=True)
        uut.add_namelist("moo", duplicate=False)
        uut.add_namelist("pot", duplicate=True)
        output_file = tmp_path / "content_mod.f90"
        uut.write_module(output_file)

        expected_file = HERE / "content_mod.f90"
        assert output_file.read_text(
            encoding="ascii"
        ) + "\n" == expected_file.read_text(encoding="ascii")

        output_file = tmp_path / "bar_nml_iterator_mod.f90"
        expected_file = HERE / "bar_nml_iterator_mod.f90"
        assert output_file.read_text(
            encoding="ascii"
        ) + "\n" == expected_file.read_text(encoding="ascii")

        output_file = tmp_path / "pot_nml_iterator_mod.f90"
        expected_file = HERE / "pot_nml_iterator_mod.f90"
        assert output_file.read_text(
            encoding="ascii"
        ) + "\n" == expected_file.read_text(encoding="ascii")
