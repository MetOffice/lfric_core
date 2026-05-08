#! /usr/bin/env python3

##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

'''
This module contains the default configuration for NCI. It will be invoked
by the Baf scripts. This script:
- sets intel-classic as the default compiler suite to use.
- Adds the tau compiler wrapper as (optional) compilers to the ToolRepository.
'''

from fab.api import Category, ToolRepository

from default.config import Config as DefaultConfig


class Config(DefaultConfig):
    '''
    For NCI, make intel the default, and add the Tau wrapper.
    '''

    def __init__(self):
        super().__init__()
        tr = ToolRepository()
        tr.set_default_compiler_suite("intel-classic")

        # ATM we don't use a shell when running a tool, and as such
        # we can't directly use "$()" as parameter. So query these values using
        # Fab's shell tool (doesn't really matter which shell we get, so just
        # ask for the default):
        shell = tr.get_default(Category.SHELL)
        # We must remove the trailing new line, and create a list:
        nc_flibs = shell.run(additional_parameters=["-c", "nf-config --flibs"],
                             capture_output=True).strip().split()
        linker = tr.get_tool(Category.LINKER, "linker-tau-ifort")

        # Setup all linker flags:
        linker.add_lib_flags("netcdf", nc_flibs)
        linker.add_lib_flags("yaxt", ["-lyaxt", "-lyaxt_c"])
        linker.add_lib_flags("xios", ["-lxios"])
        linker.add_lib_flags("hdf5", ["-lhdf5"])
        linker.add_lib_flags("shumlib", ["-lshum"])

        # Always link with C++ libs
        linker.add_post_lib_flags(["-lstdc++"])
