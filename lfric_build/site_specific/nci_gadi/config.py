#! /usr/bin/env python3

##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

'''
This module contains the default configuration for NCI. It will be invoked
by the Fab scripts. This script sets intel-llvm as the default compiler
suite to use, and adds the required site-specific linker and include flag.
'''

import os

from fab.api import BuildConfig, Category, Linker, ToolRepository

from default.config import Config as DefaultConfig


class Config(DefaultConfig):
    '''
    For NCI, make intel the default, and setup link paths for gnu,
    intel-classic and intel-llvm.
    '''

    def __init__(self):
        super().__init__()
        tr = ToolRepository()
        tr.set_default_compiler_suite("intel-llvm")

    def setup_gnu(self, build_config: BuildConfig) -> None:
        '''
        This method sets up the Gnu compiler and linker flags.
        For now call an external function, since it is expected that
        this configuration can be very lengthy (once we support
        compiler modes).

        :param build_config: the Fab build configuration instance
        '''
        super().setup_gnu(build_config)
        tr = ToolRepository()
        linker = tr.get_tool(Category.LINKER, "linker-gfortran")
        # Add netcdf and pfunit flags
        self._setup_linker(linker)

    def setup_intel_classic(self, build_config: BuildConfig) -> None:
        '''
        This method sets up the Gnu compiler and linker flags.
        For now call an external function, since it is expected that
        this configuration can be very lengthy (once we support
        compiler modes).

        :param build_config: the Fab build configuration instance
        '''
        super().setup_intel_classic(build_config)
        tr = ToolRepository()
        linker = tr.get_tool(Category.LINKER, "linker-ifort")
        # Add netcdf and pfunit flags
        self._setup_linker(linker)

    def setup_intel_llvm(self, build_config: BuildConfig) -> None:
        '''
        This method sets up the Gnu compiler and linker flags.
        For now call an external function, since it is expected that
        this configuration can be very lengthy (once we support
        compiler modes).

        :param build_config: the Fab build configuration instance
        '''
        super().setup_intel_llvm(build_config)
        tr = ToolRepository()
        compiler = tr.get_tool(Category.FORTRAN_COMPILER, "ifx")
        if not self.args.no_test:
            # TODO: path-specific flags required here.
            # pfunit driver triggers an error " #7977: The type of the
            # function reference does not match the type of the function
            # definition"
            # for the call suite%addTest(skeleton_test_suite()) in the
            # driver. This flag should only be set for
            # unit-test/driver_<NAME>.f90:
            compiler.add_flags(["-warn", "nointerfaces"], "base")

        linker = tr.get_tool(Category.LINKER, "linker-ifx")
        # Add netcdf and pfunit flags
        self._setup_linker(linker)

    def _setup_linker(self, linker: Linker) -> None:
        """
        Generic setup of a linker (since in the NCI container environments
        the paths are actually the same, independent of the compiler used).
        This adds netcdf and pfunit definitions.

        :param linker: the linker instance to setup
        """
        nf_config = NfConfig()
        if nf_config.is_available:        
            # If not available, the site-specific setup must define netcdf
            linker.add_lib_flags("netcdf", nf_config.get_linker_flags(),
                                 silent_replace=True)

        tr = ToolRepository()
        pfunit = tr.get_tool(Category.PFUNIT, "funitproc")
        pfunit_root = pfunit.get_root_path()
        spack_view = os.environ.get("SPACK_ENV_VIEW", "")

        linker.add_lib_flags(
            "pfunit",
            [f"-L{pfunit_root}/lib", "-lfunit", "-lpfunit",
             f"-L{spack_view}/FARGPARSE-1.7/lib/", "-lfargparse",
             f"-L{spack_view}/GFTL_SHARED-1.8/lib", "-lgftl-shared-v2",
             ], silent_replace=True)
