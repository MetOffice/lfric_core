#!/usr/bin/env python3

##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

'''
This file contains a function that sets the default flags for all
Intel llvm based compilers and linkers in the ToolRepository (ifx, icx).

This function gets called from the default site-specific config file
'''

import argparse
from typing import cast

from fab.api import BuildConfig, Category, Compiler, Linker, ToolRepository


def setup_script_intel_llvm(build_config: BuildConfig,
                            args: argparse.Namespace) -> None:
    # pylint: disable=unused-argument, too-many-locals
    '''
    Defines the default flags for all Intel llvm compilers.

    :param build_config: the Fab build config instance from which
        required parameters can be taken.
    :param args: all command line options
    '''

    tr = ToolRepository()
    ifx = tr.get_tool(Category.FORTRAN_COMPILER, "ifx")
    ifx = cast(Compiler, ifx)

    if not ifx.is_available:
        ifx = tr.get_tool(Category.FORTRAN_COMPILER, "mpif90-ifx")
        ifx = cast(Compiler, ifx)
        if not ifx.is_available:
            return

    # The base flags
    # ==============
    # The following flags will be applied to all modes:
    ifx.add_flags(["-stand", "f08"],               "base")
    ifx.add_flags(["-g", "-traceback"],            "base")
    # With -warn errors we get externals that are too long. While this
    # is a (usually safe) warning, the long externals then causes the
    # build to abort. So for now we cannot use `-warn errors`
    ifx.add_flags(["-warn", "all"],                "base")

    # By default turning interface warnings on causes "genmod" files to be
    # created. This adds unnecessary files to the build so we disable that
    # behaviour.
    ifx.add_flags(["-gen-interfaces", "nosource"], "base")

    # Full debug
    # ==========
    ifx.add_flags(["-check", "all", "-fpe0"], "full-debug")
    ifx.add_flags(["-O0", "-ftrapuv"],        "full-debug")

    # Fast debug
    # ==========
    ifx.add_flags(["-O2", "-fp-model=strict"], "fast-debug")

    # Production
    # ==========
    ifx.add_flags(["-O3", "-xhost"], "production")

    # Set up the linker
    # =================
    # This will implicitly affect all ifx based linkers, e.g.
    # linker-mpif90-ifx will use these flags as well.
    linker = tr.get_tool(Category.LINKER, f"linker-{ifx.name}")
    linker = cast(Linker, linker)    # Make mypy happy
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
    linker.add_post_lib_flags(["-lstdc++"])
