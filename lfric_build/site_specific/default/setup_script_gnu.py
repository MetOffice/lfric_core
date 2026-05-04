#!/usr/bin/env python3

##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

'''
This file contains a function that sets the default flags for all
GNU based compilers and linkers in the ToolRepository.

This function gets called from the default site-specific config file
'''

import argparse
from typing import cast

from fab.api import BuildConfig, Category, Compiler, Linker, ToolRepository


def setup_script_gnu(build_config: BuildConfig,
                     args: argparse.Namespace) -> None:
    # pylint: disable=unused-argument
    '''
    Defines the default flags for all GNU compilers and linkers.

    :param build_config: the Fab build config instance from which
        required parameters can be taken.
    :param args: all command line options
    '''

    tr = ToolRepository()
    gfortran = tr.get_tool(Category.FORTRAN_COMPILER, "gfortran")

    if not gfortran.is_available:
        gfortran = tr.get_tool(Category.FORTRAN_COMPILER, "mpif90-gfortran")
        if not gfortran.is_available:
            return
    gfortran = cast(Compiler, gfortran)

    if gfortran.get_version() < (4, 9):
        raise RuntimeError(f"GFortran is too old to build LFRic. "
                           f"Must be at least 4.9.0, it is "
                           f"'{gfortran.get_version_string()}'.")

    # The base flags
    # ==============

    # TODO: It should use -Werror=conversion, but:
    # Most lfric_atm dependencies contain code with implicit lossy
    # conversions.
    # This should be restricted to only the files/directories
    # that need it, but this needs Fab updates.

    gfortran.add_flags(
        ['-ffree-line-length-none', '-Wall', '-g',
         '-Werror=character-truncation',
         '-Werror=unused-value',
         '-Werror=tabs',
         '-std=f2008',
         '-fdefault-real-8',
         '-fdefault-double-8',
         ],
        "base")

    # TODO - Remove the -fallow-arguments-mismatch flag when MPICH no longer
    #        fails to build as a result of its mismatched arguments (see
    #        ticket summary for #2549 for reasoning).
    if gfortran.get_version() >= (10, 0):
        gfortran.add_flags("-fallow-argument-mismatch", "base")

    runtime = ["-fcheck=all", "-ffpe-trap=invalid,zero,overflow"]
    init = ["-finit-integer=31173",  "-finit-real=snan",
            "-finit-logical=true", "-finit-character=85"]
    # Full debug
    # ==========
    gfortran.add_flags(runtime + ["-O0"] + init, "full-debug")

    # Fast debug
    # ==========
    gfortran.add_flags(runtime + ["-Og"], "fast-debug")

    # Production
    # ==========
    gfortran.add_flags(["-Ofast"], "production")

    # Set up the linker
    # =================
    # This will implicitly affect all gfortran based linkers, e.g.
    # linker-mpif90-gfortran will use these flags as well.
    linker = tr.get_tool(Category.LINKER, f"linker-{gfortran.name}")
    linker = cast(Linker, linker)

    # ATM we don't use a shell when running a tool, and as such
    # we can't directly use "$()" as parameter. So query these values using
    # Fab's shell tool (doesn't really matter which shell we get, so just
    # ask for the default):
    shell = tr.get_default(Category.SHELL)

    try:
        # We must remove the trailing new line, and create a list:
        nc_flibs = shell.run(additional_parameters=["-c", "nf-config --flibs"],
                             capture_output=True).strip().split()
    except RuntimeError:
        nc_flibs = []

    linker.add_lib_flags("netcdf", nc_flibs)
    linker.add_lib_flags("yaxt", ["-lyaxt", "-lyaxt_c"])
    linker.add_lib_flags("xios", ["-lxios"])
    linker.add_lib_flags("hdf5", ["-lhdf5"])
    linker.add_lib_flags("shumlib", ["-lshum"])
    linker.add_lib_flags("vernier", ["-lvernier_f", "-lvernier_c",
                                     "-lvernier"])

    # Always link with C++ libs
    linker.add_post_lib_flags(["-lstdc++"], "base")
