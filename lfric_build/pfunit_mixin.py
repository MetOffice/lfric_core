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
from typing import List, Optional, Iterable, Union

from fab.api import (ArtefactSet, Category, Exclude, grab_folder, Include,
                     input_to_output_fpath, step)


class PfUnitMixin:
    '''
    This mixin adds support for pFUnit based testing. It also adds
    a command line option to disable testing. This class will also
    automatically detect if there is no unit-test directory and handle
    this case correctly.

    :param name: the name to be used for the workspace. Note that
        the name of the compiler will be added to it.
    :param app_dir: the base directory of the application.
    :param root_symbol: the symbol (or list of symbols) of the main
        programs. Defaults to the parameter `name` if not specified.

    '''

    # The new artefact set to use for pf files
    PF_SOURCE = "PF_SOURCE"

    # pylint: disable=too-many-instance-attributes
    def __init__(self, name: str,
                 app_dir: Path,
                 root_symbol: Optional[Union[List[str], str]] = None
                 ):

        self._has_test = False
        super().__init__(name, app_dir=app_dir, root_symbol=root_symbol)

    def define_command_line_options(
            self,
            parser: Optional[argparse.ArgumentParser] = None
            ) -> argparse.ArgumentParser:
        '''
        Adds an option to disable testing

        :param parser: optional a pre-defined argument parser.

        :returns: the argument parser with the LFRic specific options added.
        '''
        parser = super().define_command_line_options()

        parser.add_argument(
            '--no-test', action="store_true", default=False,
            help="Disable compilation of pFUnit tests.")

        return parser

    def get_linker_flags(self) -> List[str]:
        '''
        This method adds pFUnit as library if tests are available and enabled.

        :returns: list of flags for the linker.
        '''
        libs: list[Path] = []
        if self._has_test:
            # TODO: This implies that pfunit will be used when linking the
            #       actual app. We need improved support for path-specific
            #       flags in fab to specify a lib only to be used depending
            #       on output.
            libs.append('pfunit')
        return libs + super().get_linker_flags()

    def grab_files_step(self) -> None:
        '''
        This method adds files from APPS/unit-test if not disabled via the
        command line switch. If this directory exists (and testing is not
        disabled), it will set `self._has_test` to include compilation of
        pFUnit test files in the future build steps.
        '''

        super().grab_files_step()

        unit_test = "unit-test"
        # Check if there are unit tests
        if (not self.args.no_test) and (self.app_dir / unit_test).is_dir():
            grab_folder(self.config, src=self.app_dir / unit_test,
                        dst_label=unit_test)
            # Some tests also need the .f90 files from
            # components/science/unit-tests, but not the .pf files. So, only
            # pick the directories that contain f90 files (picking all files,
            # including .pf, would add these tests to each unit-test, and
            # besides being not intended, might not even compile in
            # lfric_apps).
            core_test_dir = (self.lfric_core_root / "components" / "science" /
                             "unit-test")
            dirs = set()
            for path in core_test_dir.rglob("*90"):
                dirs.add(path.parent)
            for path in dirs:
                # Store the files in the corresponding subdirectories (without
                # this when rsync-ing `a` and `a/b` you end up with duplicated
                # files).
                dst = path.relative_to(core_test_dir)
                grab_folder(self.config,
                            src=path,
                            dst_label=unit_test / dst)

            self._has_test = True

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
        super().find_source_files_step(path_filters=path_filters)
        if not self._has_test:
            # Don't do anything else if there are no test files (or testing
            # was explicitly disabled on the command line).
            return

        self.config.artefact_store[PfUnitMixin.PF_SOURCE] = set()
        self.config.artefact_store.copy_artefacts(
            ArtefactSet.INITIAL_SOURCE_FILES,
            PfUnitMixin.PF_SOURCE,
            suffixes=[".pf", ".PF"])

        pfunit = self.config.tool_box.get_tool(Category.PFUNIT)
        driver_f90 = pfunit.get_driver_f90()
        driver_f90 = driver_f90.replace("program main",
                                        f"program {self.name}_unit_test")

        out_driver = (self.config.build_output / "unit-test" /
                      f"driver_{self.name}.F90")
        out_driver.parent.mkdir(parents=True, exist_ok=True)
        with out_driver.open("w", encoding='utf-8') as f:
            f.write(driver_f90)

        self.config.artefact_store.add(ArtefactSet.FORTRAN_COMPILER_FILES,
                                       out_driver)

    @step
    def preprocess_pfunit_step(self) -> None:
        """
        Preprocess all .pf files with pfunit, and create test_list.inc
        to list all tests (which is required when preprocessing the
        pfunit driver program).
        """

        pf_files = self.config.artefact_store[PfUnitMixin.PF_SOURCE]
        pfunit = self.config.tool_box.get_tool(Category.PFUNIT)
        all_tests = []
        for pf_file in pf_files:
            all_tests.append(pf_file.stem)
            output_fpath = (input_to_output_fpath(config=self.config,
                                                  input_path=pf_file)
                            .with_suffix(".F90"))
            output_fpath.parent.mkdir(parents=True, exist_ok=True)
            pfunit.process(pf_file, output_fpath)
            self.config.artefact_store.add(ArtefactSet.FORTRAN_COMPILER_FILES,
                                           output_fpath)
        test_list = self.config.build_output / "unit-test" / "test_list.inc"
        with test_list.open("w", encoding="utf-8") as f:
            for test_name in all_tests:
                f.write(f"ADD_TEST_SUITE({test_name}_suite)\n")

        # TODO: That should be path-specific
        self.add_preprocessor_flags([f"-D_TEST_SUITES=\"{test_list.name}\"",
                                    "-I", str(pfunit.get_include_path())])
        self.root_symbols.append(f"{self.name}_unit_test")
        compiler = self.config.tool_box.get_tool(Category.FORTRAN_COMPILER)
        compiler.add_flags(["-I", str(pfunit.get_include_path())])

    def preprocess_fortran_step(self) -> None:
        """
        Calls Fab's preprocessing of all Fortran files. After preprocessing
        the sources, this implementation will then also pre-process the
        test files.
        """

        # We need to call preprocess_pfunit first, since it will create
        # the test_list.inc file
        if self._has_test:
            self.preprocess_pfunit_step()
        super().preprocess_fortran_step()

    def analyse_step(
            self,
            ignore_dependencies: Optional[Iterable[str]] = None,
            find_programs: bool = False
            ) -> None:
        '''
        The method overwrites the base class analyse_step and adds
        pfunit to be ignored in the analysis step.
        '''
        if ignore_dependencies is None:
            ignore_dependencies = []
        # core/infrastructure/build/import.mk
        ignore_dep_list = list(ignore_dependencies)
        if self._has_test:
            ignore_dep_list += ['pfunit']
        super().analyse_step(ignore_dependencies=ignore_dep_list,
                             find_programs=find_programs)
