##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Author J. Henrichs, Bureau of Meteorology
# Author J. Lyu, Bureau of Meteorology

"""
This module contains a function that returns a working version of a
rose_picker tool. It can either be a version installed in the system,
or otherwise a checked-out version in the fab-workspace will be used.
If required, a version of rose_picker will be checked out.
"""

import logging
from pathlib import Path

from fab.api import Tool

logger = logging.getLogger(__name__)


class RosePicker(Tool):
    '''This implements rose_picker as a Fab tool.

    :param path: the path to the rose picker binary.
    '''
    def __init__(self):
        super().__init__("rose_picker", exec_name="rose_picker",
                         availability_option="-help")

    def execute(self,
                rose_meta_conf: Path,
                directory: Path,
                include_paths: list[Path]) -> None:
        '''

        :param rose_meta_conf: Path to the metadata file to load.
        :param directory: Path to the output directory.
        :param include_paths: List of include directories which are
            searched for inherited metadata files.
        '''
        params = [rose_meta_conf, "-directory", directory]
        for inc_path in include_paths:
            params.extend(["-include_dirs", inc_path])

        super().run(additional_parameters=params)
