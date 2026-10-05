#!/usr/bin/env python3

##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

'''
This file contains a function that sets the default flags for all
Intel classic based compilers in the ToolRepository (ifort, icc).

This function gets called from the default site-specific config file
'''

import argparse
from typing import cast

from fab.api import BuildConfig, Category, Compiler, Linker, ToolRepository

from nf_config import NfConfig


def setup_script_intel_classic(build_config: BuildConfig,
                               args: argparse.Namespace) -> None:
    # pylint: disable=unused-argument, too-many-locals
    '''
    Defines the default flags for all Intel classic compilers and linkers.

    :param build_config: the Fab build config instance from which
        required parameters can be taken.
    :param args: all command line options
    '''

    tr = ToolRepository()
    ifort = tr.get_tool(Category.FORTRAN_COMPILER, "ifort")
    ifort = cast(Compiler, ifort)

    if not ifort.is_available:
        # This can happen if ifort is not in path (in spack environments).
        # To support this common use case, see if mpif90-ifort is available,
        # and initialise this otherwise.
        ifort = tr.get_tool(Category.FORTRAN_COMPILER, "mpif90-ifort")
        ifort = cast(Compiler, ifort)
        if not ifort.is_available:
            # Since some flags depends on version, the code below requires
            # that the intel compiler actually works.
            return

    # The base flags
    # ==============
    # The following flags will be applied to all modes:
    ifort.add_flags(["-stand", "f08"],               "base")
    ifort.add_flags(["-g", "-traceback"],            "base")
    # With -warn errors we get externals that are too long. While this
    # is a (usually safe) warning, the long externals then causes the
    # build to abort. So for now we cannot use `-warn errors`
    ifort.add_flags(["-warn", "all"],                "base")

    # By default turning interface warnings on causes "genmod" files to be
    # created. This adds unnecessary files to the build so we disable that
    # behaviour.
    ifort.add_flags(["-gen-interfaces", "nosource"], "base")

    # The "-assume realloc-lhs" switch causes Intel Fortran prior to v17 to
    # actually implement the Fortran2003 standard. At version 17 it becomes the
    # default behaviour.
    if ifort.get_version() < (17, 0):
        ifort.add_flags(["-assume", "realloc-lhs"], "base")

    # Full debug
    # ==========
    # ifort.mk: bad interaction between array shape checking and
    # the matmul" intrinsic in at least some iterations of v19.
    if (19, 0, 0) <= ifort.get_version() < (19, 1, 0):
        runtime_flags = ["-check", "all,noshape", "-fpe0"]
    else:
        runtime_flags = ["-check", "all", "-fpe0"]
    ifort.add_flags(runtime_flags,        "full-debug")
    ifort.add_flags(["-O0", "-ftrapuv"],  "full-debug")

    # Fast debug
    # ==========
    ifort.add_flags(["-O2", "-fp-model=strict"], "fast-debug")

    # Production
    # ==========
    ifort.add_flags(["-O3", "-xhost"], "production")

    # Set up the linker
    # =================
    # This will implicitly affect all ifort based linkers, e.g.
    # linker-mpif90-ifort will use these flags as well.
    linker = tr.get_tool(Category.LINKER, f"linker-{ifort.name}")
    linker = cast(Linker, linker)

    nf_config = NfConfig()
    if nf_config.is_available:
        # If not available, the site-specific setup must define netcdf
        linker.add_lib_flags("netcdf", nf_config.get_linker_flags())

    linker.add_lib_flags("yaxt", ["-lyaxt", "-lyaxt_c"])
    linker.add_lib_flags("xios", ["-lxios"])
    linker.add_lib_flags("hdf5", ["-lhdf5"])
    linker.add_lib_flags("shumlib", ["-lshum"])
    linker.add_lib_flags("vernier", ["-lvernier_f", "-lvernier_c",
                                     "-lvernier"])

    # This likely needs adjusting, pfunit required fargparse and gftl
    linker.add_lib_flags("pfunit", ["-lfunit", "-lpfunit",
                                    "-lfargparse", "lgftl-shared-v2"])

    # Always link with C++ libs
    linker.add_post_lib_flags(["-lstdc++"])
