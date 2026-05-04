#! /usr/bin/env python3

##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

'''
This module contains a setup NIWA's XC-50
'''

import os
from typing import cast

from fab.api import BuildConfig, Category, Linker, ToolRepository

from default.config import Config as DefaultConfig


class Config(DefaultConfig):
    '''
    This config class sets specific flags for NIWA's XC-50
    '''

    def __init__(self):
        super().__init__()
        tr = ToolRepository()
        tr.set_default_compiler_suite("intel-classic")

    def setup_cray(self, build_config: BuildConfig) -> None:
        '''
        First call the base class to get all default options.
        See the file ../default/setup_cray.py for the current
        default. Then the NIWA's XC-50 specific flags are added.
        The linker is also updated.

        :param build_config: the Fab build config instance from which
            required parameters can be taken.
        :type build_config: :py:class:`fab.BuildConfig`
        '''
        super().setup_cray(build_config)
        tr = ToolRepository()
        ftn = tr.get_tool(Category.FORTRAN_COMPILER, "crayftn-ifort")
        # Add any flags you want to have:
        ftn.add_flags([f"-I{os.environ['EBROOTXIOS']}/inc"])

        # Update the linker. This is what the default sets up
        # (except NetCDF, which is defined using nf-config, and
        # should likely work the way it is):
        linker = tr.get_tool(Category.LINKER, "linker-crayftn-ftn")
        linker = cast(Linker, linker)   # make mypy happy

        # The first parameter specifies the internal name for libraries,
        # followed by a list of linker options. If you should need additional
        # library paths, you could e.g. use:
        # linker.add_lib_flags("yaxt", ["-L", "/my/path/to/yaxt", "-lyaxt",
        #                               "-lyaxt_c"])
        # Make sure to not use a space as ONE parameter ("-L /my/lib"),
        # you have to specify them as two separate list elements

        linker.add_lib_flags("yaxt", ["-lyaxt", "-lyaxt_c"])
        linker.add_lib_flags("xios", ["-lxios"])
        linker.add_lib_flags("hdf5", ["-lhdf5"])
        linker.add_lib_flags("shumlib", ["-lshum"])
        linker.add_lib_flags("vernier", ["-lvernier_f", "-lvernier_c",
                                         "-lvernier"])
