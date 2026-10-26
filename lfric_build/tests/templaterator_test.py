#!/usr/bin/env python3

##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Author: J. Henrichs, Bureau of Meteorology

"""
This module tests the Templaterator tool.
"""

from pathlib import Path
from unittest.mock import patch

import pytest

from templaterator import Templaterator  # adjust import as needed


@pytest.fixture(name="templaterator")
def templaterator_setup(tmp_path: Path) -> Templaterator:
    """
    :returns: A dummy Templaterator object (in tmp_path,
        which does not exist, but it's all we need for these tests since
        all tests are mocked.
    """
    return Templaterator(tmp_path / "Templaterator.py")


def test_init(templaterator: Templaterator,
              tmp_path: Path) -> None:
    """
    Test the constructor.
    """
    assert templaterator.name == "Templaterator"
    assert templaterator.exec_path == tmp_path / "Templaterator.py"
    assert templaterator.exec_name == "Templaterator.py"


def test_check_available(templaterator: Templaterator) -> None:
    """
    Test the check_available function.
    """
    with patch("fab.tools.tool.Tool.run", return_value=0) as mock_run:
        assert templaterator.check_available() is True
    mock_run.assert_called_once_with(additional_parameters="-h")

    with patch("fab.tools.tool.Tool.run",
               side_effect=RuntimeError()) as mock_run:
        assert templaterator.check_available() is False
    mock_run.assert_called_once_with(additional_parameters="-h")


def test_process_call(templaterator: Templaterator,
                      tmp_path: Path) -> None:
    """
    Test that execution passes on the right parameter to the
    Templaterator script.
    """
    input_template = tmp_path / "input.txt"
    output_file = tmp_path / "output.txt"
    key_values = {"A": "1", "B": "2"}

    with patch("fab.tools.tool.Tool.run", return_value=0) as mock_run:
        templaterator.process(input_template, output_file, key_values)

        expected_params = [
            input_template,
            "-o",
            output_file,
            "-s A=1",
            "-s B=2",
        ]

    mock_run.assert_called_once_with(additional_parameters=expected_params)
