#!/usr/bin/env python3

##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

'''
This file contains a function that sets the default flags for the NVIDIA
compilers and linkers in the ToolRepository.

This function gets called from the default site-specific config file
'''

import argparse
from typing import cast

from fab.api import BuildConfig, Category, Compiler, Linker, ToolRepository

from nf_config import NfConfig


def setup_script_nvidia(build_config: BuildConfig,
                        args: argparse.Namespace) -> None:
    # pylint: disable=unused-argument
    '''
    Defines the default flags for nvfortran.

    :param build_config: the Fab build config instance from which
        required parameters can be taken.
    :param args: all command line options
    '''

    tr = ToolRepository()
    nvfortran = tr.get_tool(Category.FORTRAN_COMPILER, "nvfortran")
    nvfortran = cast(Compiler, nvfortran)

    if not nvfortran.is_available:
        nvfortran = tr.get_tool(Category.FORTRAN_COMPILER, "mpif90-nvfortran")
        nvfortran = cast(Compiler, nvfortran)
        if not nvfortran.is_available:
            return

    # The base flags
    # ==============
    flags = ["-Mextend",           # 132 characters line length
             "-g", "-traceback",
             "-r8",                # Default 8 bytes reals
             "-O0",                # No optimisations
             ]

    lib_flags = ["-c++libs"]

    # Handle accelerator options:
    if args.openacc or args.openmp:
        host = args.host.lower()
    else:
        # Neither openacc nor openmp specified
        host = ""

    if args.openacc:
        if host == "gpu":
            flags.extend(["-acc=gpu", "-gpu=managed"])
            lib_flags.extend(["-aclibs", "-cuda"])
        else:
            # CPU
            flags.extend(["-acc=cpu"])
    elif args.openmp:
        if host == "gpu":
            flags.extend(["-mp=gpu", "-gpu=managed"])
            lib_flags.append("-cuda")
        else:
            # OpenMP on CPU, that's already handled by Fab
            pass

    nvfortran.add_flags(flags, "base")

    # Full debug
    # ==========
    nvfortran.add_flags(["-O0", "-fp-model=strict"], "full-debug")

    # Fast debug
    # ==========
    nvfortran.add_flags(["-O2", "-fp-model=strict"], "fast-debug")

    # Production
    # ==========
    nvfortran.add_flags(["-O4"], "production")

    # Set up the linker
    # =================
    # This will implicitly affect all nvfortran based linkers, e.g.
    # linker-mpif90-nvfortran will use these flags as well.
    linker = tr.get_tool(Category.LINKER, f"linker-{nvfortran.name}")
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
    linker.add_post_lib_flags(lib_flags)
