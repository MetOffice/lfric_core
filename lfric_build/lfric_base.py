##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Author: J. Henrichs, Bureau of Meteorology
# Author: J. Lyu, Bureau of Meteorology

"""
This is an OO basic interface to FAB. It allows the typical LFRic
applications to only modify very few settings to have a working FAB build
script.
"""

import argparse
from pathlib import Path
import sys
from typing import List, Optional, Iterable, Union

from fab.api import (ArtefactSet, BuildConfig, Category, Exclude, grab_folder,
                     Include, input_to_output_fpath, preprocess_x90, psyclone,
                     step, SuffixFilter)
from fab.fab_base.fab_base import FabBase

from configurator import configurator
from templaterator import Templaterator


class LFRicBase(FabBase):
    '''
    This is the base class for all LFRic FAB scripts.

    :param name: the name to be used for the workspace. Note that
        the name of the compiler will be added to it.
    :param root_symbol: the symbol (or list of symbols) of the main
        programs. Defaults to the parameter `name` if not specified.

    '''
    # pylint: disable=too-many-instance-attributes
    def __init__(self, name: str,
                 root_symbol: Optional[Union[List[str], str]] = None
                 ):

        # List of all precision preprocessor symbols and their default.
        # Used to add corresponding command line options, and then to define
        # the preprocessor definitions. Note that precision_other
        # becomes RDEF.
        self._all_precisions = [("precision_other", "64"),
                                ("R_SOLVER_PRECISION", "32"),
                                ("R_TRAN_PRECISION", "64"),
                                ("R_BL_PRECISION", "64")]

        super().__init__(name)

        this_file = Path(__file__)
        # The root directory of the LFRic Core
        self._lfric_core_root = this_file.parents[1]

        # If the user wants to overwrite the default root symbol (which
        # is `name`):
        if root_symbol:
            self.set_root_symbol(root_symbol)

        # Many PSyclone scripts use module(s) from this directory. Additional
        # paths might need to be added later.
        self._add_python_paths = [str(self.lfric_core_root / "infrastructure" /
                                      "build" / "psyclone")]
        linker = self.config.tool_box.get_tool(Category.LINKER,
                                               mpi=self.config.mpi,
                                               openmp=self.config.openmp,
                                               enforce_fortran_linker=True)
        try:
            linker.get_lib_flags("netcdf")
        except RuntimeError as err:
            msg = (f"LFRic needs NetCDF, but the linker '{linker.name}' "
                   f"has no NetCDF library setting defined. Aborting.")
            raise RuntimeError(msg) from err

    def define_command_line_options(
            self,
            parser: Optional[argparse.ArgumentParser] = None
            ) -> argparse.ArgumentParser:
        '''
        This adds LFRic specific command line options to the base class
        define_command_line_option. Currently, precision-related options
        are added.

        :param parser: optional a pre-defined argument parser.

        :returns: the argument parser with the LFRic specific options added.
        '''
        parser = super().define_command_line_options()

        parser.add_argument(
            '--no-xios', action="store_true", default=False,
            help="Disable compilation with XIOS.")

        # Precision related command line arguments
        # ----------------------------------------
        group = parser.add_argument_group(
            title="Precisions",
            description="Arguments related to setting the floating "
                        "point precision.")

        for prec_name, default in self._all_precisions:
            lower_name = prec_name.lower()
            if prec_name == "precision_other":
                help_msg = "Precision for other floating point values."
            else:
                help_msg = f"Precision for '{prec_name}'."
            group.add_argument(
                f'--{lower_name}', type=str, choices=['32', '64'],
                default=default, help=help_msg)

        return parser

    @property
    def lfric_core_root(self) -> Path:
        '''
        :returns: the root directory of the LFRic core repository.
        '''
        return self._lfric_core_root

    def setup_site_specific_location(self):
        '''
        This method adds the required directories for site-specific
        configurations to the Python search path. We want to add the
        directory where this lfric_base class is located, and not the
        directory in which the application script is (which is what
        baf base would set up).
        '''
        this_dir = Path(__file__).parent
        # We need to add the 'site_specific' directory to the path, so
        # each config can import from 'default' (instead of having to
        # use 'site_specific.default', which would hard-code the name
        # `site_specific` in more scripts).
        sys.path.insert(0, str(this_dir / "site_specific"))

    def define_preprocessor_flags_step(self) -> None:
        '''
        This method overwrites the base class define_preprocessor_flags.
        It uses add_preprocessor_flags to set up preprocessing flags for LFRic
        applications. This includes:
        - various floating point precision related directives
        - Use of XIOS (if not disabled using --no-xios command line option)
        - Disabling MPI (if disabled using --no-mpi)
        '''
        preprocessor_flags: List[str] = []

        # Check all required precision defines
        for prec_name, _ in self._all_precisions:
            # Check if a value was specified on the command line:
            value = getattr(self.args, prec_name.lower())
            if prec_name == "precision_other":
                preprocessor_flags.append(f"-DRDEF_PRECISION={value}")
            else:
                preprocessor_flags.append(f"-D{prec_name}={value}")

        # core/components/lfric-xios/build/import.mk
        if not self.args.no_xios:
            preprocessor_flags.append('-DUSE_XIOS')

        if not self.config.mpi:
            preprocessor_flags.append("-DNO_MPI")

        self.add_preprocessor_flags(preprocessor_flags)

    def get_linker_flags(self) -> List[str]:
        '''
        This method overwrites the base class get_liner_flags. It passes the
        libraries that LFRic uses to the linker. Currently, these libraries
        include yaxt, xios, netcdf and hdf5.

        :returns: list of flags for the linker.
        '''
        libs = ['yaxt', 'xios', 'netcdf', 'hdf5']
        return libs + super().get_linker_flags()

    def grab_files_step(self) -> None:
        '''
        This method overwrites the base class grab_files_step. It includes all
        the LFRic core directories that are commonly required for building
        LFRic applications. It also grabs the psydata directory for profiling,
        if required.
        '''
        dirs = ['infrastructure/source/',
                'components/driver/source/',
                'components/inventory/source/',
                'components/science/source/',
                'components/lfric-xios/source/',
                ]

        # pylint: disable=redefined-builtin
        for dir in dirs:
            grab_folder(self.config, src=self.lfric_core_root / dir,
                        dst_label='')

        # Copy the PSyclone Config file into a separate directory
        dir = "etc"
        grab_folder(self.config, src=self.lfric_core_root / dir,
                    dst_label='psyclone_config')

    def find_source_files_step(
            self,
            path_filters: Optional[Iterable[Union[Exclude, Include]]] = None
            ) -> None:
        '''
        This method overwrites the base class find_source_files_step.
        It first calls the configurator_step to set up the configurator.
        Then it finds all the source files in the LFRic core directories,
        excluding the unit tests. Finally, it calls the templaterator_step.

        :param path_filters: optional list of path filters to be passed to
            Fab find_source_files, default is None.
        '''
        self.configurator_step()

        path_filter_list = list(path_filters) if path_filters else []
        path_filter_list.append(Exclude('unit-test', '/test/'))
        super().find_source_files_step(path_filters=path_filter_list)

        self.templaterator_step(self.config)

    @step
    def configurator_step(
            self,
            include_paths: Optional[list[Path]] = None) -> None:
        '''
        This method first gets the rose meta data information by calling
        get_rose_meta. If the rose meta data is available, it then get the
        rose picker tool by calling the get_rose_picker. Finally, it runs
        the LFRic configurator with the LFRic core and apps sources by calling
        configurator.

        :param include_paths: optional additional include paths
        '''
        rose_meta = self.get_rose_meta()
        if rose_meta:
            # Ideally we would want to get all source files created in
            # the build directory, but then we need to know the list of
            # files to add them to the list of files to process. Instead,
            # we create the files in the source directory, and find them
            # there later.
            include_paths = include_paths or []
            configurator(self.config, lfric_core_source=self.lfric_core_root,
                         rose_meta_conf=rose_meta,
                         include_paths=include_paths)

    @step
    def templaterator_step(self, config: BuildConfig) -> None:
        '''
        This method runs the LFRic templaterator Fab tool.

        :param config: the Fab build configuration
        :type config: :py:class:`fab.BuildConfig`
        '''
        base_dir = self.lfric_core_root / "infrastructure" / "build" / "tools"

        templaterator = Templaterator(base_dir/"Templaterator")
        config.artefact_store["template_files"] = set()
        t90_filter = SuffixFilter(ArtefactSet.INITIAL_SOURCE_FILES,
                                  [".t90", ".T90"])
        template_files = t90_filter(config.artefact_store)
        templ_r32 = {"kind": "real32", "type": "real"}
        templ_r64 = {"kind": "real64", "type": "real"}
        templ_i32 = {"kind": "int32", "type": "integer"}
        # Don't bother with parallelising this, it's fast
        for template_file in template_files:
            out_dir = input_to_output_fpath(config=config,
                                            input_path=template_file).parent
            out_dir.mkdir(parents=True, exist_ok=True)
            template_stem = template_file.stem.removesuffix("_mod")
            for key_values in [templ_r32, templ_r64, templ_i32]:
                out_file = (out_dir /
                            f"{template_stem}_{key_values['kind']}_mod.f90")
                templaterator.process(template_file, out_file,
                                      key_values=key_values)
                # Add the newly created file to the set of
                # Fortran files to compile
                config.artefact_store.add(ArtefactSet.FORTRAN_COMPILER_FILES,
                                          out_file)

    def get_rose_meta(self) -> Optional[Path]:
        '''
        This method returns the path to the rose meta data config file.
        Currently, it returns none. It's up to the LFRic applications to
        overwrite if required.
        '''
        return None

    def analyse_step(
            self,
            ignore_dependencies: Optional[Iterable[str]] = None,
            find_programs: bool = False
            ) -> None:
        '''
        The method overwrites the base class analyse_step.
        For LFRic, it first runs the preprocess_x90_step and then runs
        psyclone_step. Finally, it calls Fab's analyse for dependency
        analysis, ignoring the third party modules that are commonly
        used by LFRic.
        '''
        if ignore_dependencies is None:
            ignore_dependencies = []
        # core/infrastructure/build/import.mk
        ignore_dep_list = list(ignore_dependencies)
        ignore_dep_list += ['netcdf', 'mpi', 'mpi_f08', 'yaxt']
        # From core/components/lfric-xios/build/import.mk
        ignore_dep_list += ['xios', 'icontext', 'mod_wait']

        self.preprocess_x90_step()
        self.psyclone_step(ignore_dependencies=ignore_dep_list)
        super().analyse_step(
            ignore_dependencies=ignore_dep_list,
            find_programs=find_programs)

    def preprocess_x90_step(self) -> None:
        """
        Invokes the Fab preprocess step for all X90 files.
        """
        # TODO: Fab does not support path-specific flags for X90 files.
        preprocess_x90(self.config,
                       common_flags=self.preprocess_flags_common)

    def psyclone_step(
            self,
            ignore_dependencies: Optional[Iterable[str]] = None,
            additional_parameters: Optional[list[str]] = None
            ) -> None:
        '''
        This method runs Fab's psyclone. It first sets the psyclone
        command line arguments by calling get_psyclone_config to get the
        PSyclone configuration file. Additional flags can be set in the
        PSyclone tool. Finally, Fab's psyclone is called with the Fab build
        configuration, the kernel root directory, the transformation script
        got through calling `get_transformation_script`, the api, and the
        additional psyclone command line arguments.

        :param ignore_dependencies:
        :param additional_parameters: optional additional parameter for the
            PSyclone.
        '''
        psyclone_cli_args = ["--config", self.get_psyclone_config()]
        if additional_parameters:
            psyclone_cli_args.extend(additional_parameters)

        # To avoid impacting other code, store the original search path
        old_sys_path = sys.path[:]
        sys.path.extend(self._add_python_paths)
        psyclone(self.config, kernel_roots=[(self.config.build_output /
                                             "kernel")],
                 transformation_script=self.get_transformation_script,
                 api="lfric",
                 cli_args=psyclone_cli_args,
                 ignore_dependencies=ignore_dependencies)
        sys.path = old_sys_path

    def get_psyclone_config(self) -> str:
        '''
        This method can be overwritten if an application needs to provide
        a modified psyclone config file (e.g. to enable additional
        debug options).

        :returns: the PSyclone config file as string.
        '''
        return str(self.config.source_root / 'psyclone_config' /
                   'psyclone.cfg')

    def get_transformation_script(self, fpath: Path,
                                  config: BuildConfig) -> Optional[Path]:
        '''
        This method returns the path to the transformation script that PSyclone
        will use for each x90 file. It first checks if there is a specific
        transformation script for the x90 file. If not, it will see whether a
        global transformation script can be used.

        :param fpath: the path to the file being processed.
        :param config: the FAB BuildConfig instance.

        :returns: the transformation script to be used by PSyclone.
        '''
        # Newer LFRic versions have a psykal directory
        optimisation_path = (config.source_root / "optimisation" /
                             f"{self.site}-{self.platform}" / "psykal")
        relative_path = None
        # The soure file might be either in build_output (e.g. a preprocessed
        # .X90 file), or still in source (.x90 file). Check if the file
        # is in one of the two sub-trees, and use the relative path to
        # check if there is a file-specific optimisation script
        for base_path in [config.source_root, config.build_output]:
            try:
                relative_path = fpath.relative_to(base_path)
            except ValueError:
                # The file is not under the `base_path` - keep on checking
                pass

        if relative_path:
            # The file was under either source or build. Check if there
            # is a file-specific optimisation script:
            local_transformation_script = (optimisation_path /
                                           (relative_path.with_suffix('.py')))
            if local_transformation_script.exists():
                return local_transformation_script

        # No file-specific optimisation script found. Check for global.py:
        global_transformation_script = optimisation_path / 'global.py'
        if global_transformation_script.exists():
            return global_transformation_script
        return None
