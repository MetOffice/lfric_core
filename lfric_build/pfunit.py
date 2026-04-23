##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT
# which you should have received as part of this distribution
##############################################################################

"""This file contains the Rsync class for synchronising file trees.
"""

import logging
import os
from pathlib import Path

from fab.tools.tool import Tool
from fab.tools.category import Category


logger = logging.getLogger(__name__)


class PfUnit(Tool):
    """
    This is a class to encapsulate pFUnit. It relies on the environment
    variable $PFUNIT to indicate the location of the source code.
    This is required since besides .mod files and executable, it also
    contains the source code for a Fortran driver program .
    It assumes that pFUnit's preprocessor `funitproc` is in $PFUNIT/bin.
    """

    def __init__(self):
        pfunit_home = os.environ.get("PFUNIT", "")
        if not pfunit_home:
            logger.error("$PFUNIT not defined in an environment, testing will"
                         "likely not work.")
        self._pfunit_home = Path(pfunit_home)

        # Create a new category of pFUnit
        Category("PFUNIT")
        exec_name = self._pfunit_home / "bin" / "funitproc"
        super().__init__("funitproc", exec_name=exec_name,
                         category=Category.PFUNIT,
                         availability_option="-v")

    def get_root_path(self) -> Path:
        """
        :returns: the root path of pFUnit.
        """
        return self._pfunit_home

    def get_include_path(self) -> Path:
        """
        :returns: the include directory for PFUnit.
        """
        return self._pfunit_home / "include"

    def get_driver_f90(self) -> str:
        """
        :returns: the content of pFUnit's driver.F90 file.
        """
        driver_path = self._pfunit_home / "include" / "driver.F90"
        with driver_path.open("r", encoding='utf-8') as f:
            driver_f90 = f.read()
        return driver_f90

    def process(self, pf_path: Path,
                f90_out_path: Path):
        """
        Processes the .pf file to create an output f90 file.

        :param pf_path: the input path.
        :param f90_out_path: destination path.
        """
        return self.run(additional_parameters=[pf_path, f90_out_path])
