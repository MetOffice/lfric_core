#!/usr/bin/env python3

##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Author: J. Henrichs, Bureau of Meteorology
# Author: J. Lyu, Bureau of Meteorology

"""
A FAB build script for applications/skeleton. It relies on
the LFRicBase class contained in the infrastructure directory.
"""

import logging
from pathlib import Path
import sys

from fab.steps.grab.folder import grab_folder

# We need to import the base class:
sys.path.insert(0, str(Path(__file__).parents[2] / "lfric_build"))

from lfric_base import LFRicBase  # noqa: E402


class FabSkeleton(LFRicBase):
    """
    A Fab-based build script for skeleton. It relies on the LFRicBase class
    to implement the actual functionality, and only provides the required
    source files.

    :param name: The name of the application.
    """

    def __init__(self, name: str = "skeleton") -> None:
        super().__init__(name=name)
        # Store the root of this apps for later
        this_file = Path(__file__).resolve()
        self._this_root = this_file.parent

    def grab_files_step(self) -> None:
        """
        Grabs the required source files and optimisation scripts.
        """
        super().grab_files_step()
        dirs = ['applications/skeleton/source/']

        # pylint: disable=redefined-builtin
        for dir in dirs:
            grab_folder(self.config, src=self.lfric_core_root / dir,
                        dst_label='')

        # Copy the optimisation scripts into a separate directory
        grab_folder(self.config, src=self._this_root / "optimisation",
                    dst_label='optimisation')

    def get_rose_meta(self) -> Path:
        """
        :returns: the rose-meta.conf path.
        """
        return (self._this_root / 'rose-meta' / 'lfric-skeleton' / 'HEAD' /
                'rose-meta.conf')


# -----------------------------------------------------------------------------
if __name__ == '__main__':

    logger = logging.getLogger('fab')
    logger.setLevel(logging.DEBUG)
    fab_skeleton = FabSkeleton()
    fab_skeleton.build()
