#! /usr/bin/env python3

##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

'''This module contains a setup for NCAS-EX (archer2)
'''

from typing import cast

from fab.api import BuildConfig, Category, Linker, ToolRepository

from default.config import Config as DefaultConfig


class Config(DefaultConfig):
    '''This config class sets specific flags for NCAS-EX (archer2)
    '''

    def __init__(self):
        super().__init__()
        tr = ToolRepository()
        tr.set_default_compiler_suite("gnu")

    def setup_cray(self, build_config: BuildConfig):
        '''First call the base class to get all default options.
        See the file ../default/setup_cray.py for the current
        default.
        Very likely, linker options need to be changed:
        '''
        super().setup_cray(build_config)
        tr = ToolRepository()
        ftn = tr.get_tool(Category.FORTRAN_COMPILER, "gfortran")
        # Any gfortran on Cray's EX need this flag in order to
        # compile mpi_mod:
        ftn.add_flags(["-fallow-argument-mismatch"])

        # Update the linker. This is what the default sets up
        # (except NetCDF, which is defined using nf-config, and
        # should likely work the way it is):
        linker = tr.get_tool(Category.LINKER, "linker-gfortran")
        linker = cast(Linker, linker)   # make mypy happy

        # Cray's don't have nf-config. Till we have figured out the
        # proper solution, hard-code some flags that might work
        # with a plain gfortran build in a spack environment:
        linker.add_lib_flags("netcdf", ["-lnetcdff", "-lnetcdf",
                                        "-lnetcdf", "-lm"])
        # That's pretty much the default:
        linker.add_lib_flags("yaxt", ["-lyaxt", "-lyaxt_c"])
        linker.add_lib_flags("xios", ["-lxios"])
        linker.add_lib_flags("hdf5", ["-lhdf5"])
        linker.add_lib_flags("shumlib", ["-lshum"])
        linker.add_lib_flags("vernier", ["-lvernier_f", "-lvernier_c",
                                         "-lvernier"])
