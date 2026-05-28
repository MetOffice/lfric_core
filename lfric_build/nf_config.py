##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT
# which you should have received as part of this distribution
##############################################################################

"""This file contains the class to interface with NetCDF's nf-config script.
"""

from typing import List

from fab.tools.category import Category
from fab.tools.tool import Tool


class NfConfig(Tool):
    '''This class interfaces with NetCDF's nf-config tool. It is not added
    to the ToolRepository, it is intended for site-specific configurations
    to make it easier to query for NetCDF settings.
    '''

    def __init__(self):
        super().__init__("nf-config", "nf-config", Category.MISC)

    def get_compiler_flags(self) -> List[str]:
        """
        :returns: the compilation flags to use for NetCDF.
        """
        flags = self.run(additional_parameters=["--fflags"])
        return flags.split()

    def get_linker_flags(self) -> List[str]:
        """
        :returns: the linker flags to use for NetCDF.
        """
        flags = self.run(additional_parameters=["--flibs"])
        return flags.split()
