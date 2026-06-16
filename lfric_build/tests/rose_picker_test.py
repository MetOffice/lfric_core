##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Author J. Henrichs, Bureau of Meteorology

"""
This module tests rose_picker_tool.
"""

import os
from pathlib import Path
from unittest.mock import patch

from rose_picker import RosePicker


def test_get_rose_picker_check_available() -> None:
    """
    Test RosePicker's check_available.
    """
    rose_picker = RosePicker()
    with patch("fab.tools.tool.Tool.run", return_value=True) as mock_run:
        assert rose_picker.check_available()
    mock_run.assert_called_once_with("-help")

    with patch.object(RosePicker, "run", side_effect=RuntimeError) as mock_run:
        assert not rose_picker.check_available()
    mock_run.assert_called_once_with("-help")


def test_get_rose_picker_execute() -> None:
    """
    Test RosePicker's check_available.
    """
    rose_picker = RosePicker()
    rose_meta_conf = Path("rose_meta_conf")
    directory = Path("/some/dir")
    p1 = Path("/path1")
    p2 = Path("/path12")
    with patch("fab.tools.tool.Tool.run", return_value=0) as mock_run, \
         patch.object(os, "environ", {}):
        rose_picker.execute(rose_meta_conf, directory, include_paths=[p1, p2])

    # Rose picker prepends the existing python path, separated by ":".
    # Since python path is not set, there will be a leading ":""
    mock_run.assert_called_once_with(
        additional_parameters=[rose_meta_conf, "-directory", directory,
                               "-include_dirs", p1, "-include_dirs", p2, ])
