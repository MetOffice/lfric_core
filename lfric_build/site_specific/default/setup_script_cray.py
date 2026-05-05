#!/usr/bin/env python3

##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

'''
This file contains a function that sets the default flags for the Cray
compilers and linkers in the ToolRepository.

This function gets called from the default site-specific config file
'''

import argparse
from typing import cast

from fab.api import BuildConfig, Category, Compiler, Linker, ToolRepository


def setup_script_cray(build_config: BuildConfig,
                      args: argparse.Namespace) -> None:
    # pylint: disable=unused-argument
    '''
    Defines the default flags for ftn.

    :param build_config: the Fab build config instance from which
        required parameters can be taken.
    :param args: all command line options
    '''

    tr = ToolRepository()
    ftn = tr.get_tool(Category.FORTRAN_COMPILER, "crayftn-ftn")
    ftn = cast(Compiler, ftn)

    if not ftn.is_available:
        return

    # The base flags
    # ==============
    flags = ["-g", "-G0", "-m", "0",    # ?
             "-M", "E664,E7208,E7212",  # ?
             "-en",                     # Fortran standard
             "-ef",                     # use lowercase module names!Important!
             "-hnocaf",                 # Required for linking with C++
             ]

    # Handle accelerator options:
    if args.openacc or args.openmp:
        host = args.host.lower()
    else:
        # Neither openacc nor openmp specified
        host = ""

    if args.openacc:
        if host == "gpu":
            flags.extend(["-h acc"])
        else:
            # CPU
            flags.extend(["-h acc"])
    elif args.openmp:
        if host == "gpu":
            flags.extend([])
        else:
            # OpenMP on CPU, that's already handled by Fab
            pass

    ftn.add_flags(flags, "base")

    # Full debug
    # ==========
    ftn.add_flags(["-Ktrap=fp",    # floating point checking
                   "-R", "bcdps",  # bounds, array shape, collapse,
                                   # pointer, string checking
                   "-O0"],         # No optimisation
                  "full-debug")
    if ftn.get_version() >= (15, 0):
        ftn.add_flags(["-G0"], "full-debug")
    else:
        ftn.add_flags(["-Gfast"], "full-debug")

    # Fast debug
    # ==========
    ftn.add_flags(["-O2", "-hflex_mp=strict"], "fast-debug")
    if ftn.get_version() >= (15, 0):
        ftn.add_flags(["-G2"], "fast-debug")
    else:
        ftn.add_flags(["-Gfast"], "fast-debug")

    # Production
    # ==========
    ftn.add_flags(["-O3", "-hipa3", "-m", "3"], "production")

    # Set up the linker
    # =================
    linker = tr.get_tool(Category.LINKER, f"linker-{ftn.name}")
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

    linker.add_post_lib_flags(["-lcraystdc++"])

    # Using the GNU compiler on Crays for now needs the additional
    # flag -fallow-argument-mismatch to compile mpi_mod.f90
    ftn = tr.get_tool(Category.FORTRAN_COMPILER, "crayftn-gfortran")
    ftn.add_flags("-fallow-argument-mismatch")
