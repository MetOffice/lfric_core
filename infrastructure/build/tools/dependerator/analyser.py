#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
# Generate a make file snippet holding build information about a source file.
#
# This snippet consists of a dependency list built up from which modules a
# file "use"s.
#
# Additionally a separate snippet file is created for each "program" unit found
# listing all the modules which go to make that program.
#
# These snippets may then be "include"ed into other make files.
#
"""
Examine Fortran source and build dependency information for use by "make".
"""

import logging
import os
import os.path
import re
import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from time import time
from typing import Dict, Generator, List, Optional, Tuple

from dependerator.database import FortranDependencies


class Analyser(ABC):
    """
    Interface for analysers.
    """

    @abstractmethod
    def analyse(self, source_filename: Path) -> None:
        """
        Examine a source file and store dependency information in the database.

        @param source_filename: File to analyse
        """
        pass


class FortranAnalyser(Analyser):
    ###########################################################################
    # Constructor
    #
    # Arguments:
    #   ignoreModules - Module names to ignore.
    #   database      - Backing store to hold details.
    #   preprocess_macros - Macro name is the key. Value may be None for
    #                            empty macros.
    #   preprocess_include_paths - Directories where inclusions will be saught.
    #
    def __init__(
        self,
        ignoreModules: List[str],
        database: FortranDependencies,
        preprocess_macros: Optional[Dict[str, Optional[str]]] = None,
        preprocess_include_paths: Optional[List[Path]] = None,
    ):
        self._ignoreModules = [str.lower(mod) for mod in ignoreModules]
        self._database = database
        self.__preprocess_macros = preprocess_macros or {}
        self.__preprocess_include_paths = preprocess_include_paths or []

        # The intrinsic Fortran modules
        self._ignoreModules.extend(
            [
                "iso_c_binding",
                "iso_fortran_env",
                "ieee_arithmetic",
                "ieee_exceptions",
                "ieee_features",
            ]
        )
        # The OpenMP libraries
        self._ignoreModules.extend(["omp_lib", "omp_lib_kinds"])

        fpp = os.getenv("FPP", None)
        if fpp is None:
            raise Exception("No Fortran preprocessor provided in $FPP")
        self._fpp = fpp.split()

        # Patterns to recognise scoping units
        #
        self._programPattern = re.compile(
            r"^\s*PROGRAM\s+(\S+)", flags=re.IGNORECASE
        )
        self._modulePattern = re.compile(
            r"^\s*MODULE\s+(?!(?:PROCEDURE|SUBROUTINE|FUNCTION)\s+)(\S+)",
            flags=re.IGNORECASE,
        )
        self._submodulePattern = re.compile(
            r"^\s*SUBMODULE\s*\((?:([^:]+):)?([^)]+)\)\s+(\S+)",
            flags=re.IGNORECASE,
        )
        self._subroutinePattern = re.compile(
            r"^\s*(MODULE\s+)?SUBROUTINE\s+([^(\s]+)", flags=re.IGNORECASE
        )
        self._functionPattern = re.compile(
            r"^\s*(?:TYPE\s*\(\s*\S+?\s*\)\s*|\S+\s*)?"
            r"(MODULE\s+)?FUNCTION\s+([^(\s]+)",
            flags=re.IGNORECASE,
        )
        self._endPattern = re.compile(
            r"^\s*END(?:\s+(\S+)(?:\s+(\S+))?)?", flags=re.IGNORECASE
        )

        # Patterns to recognise dependencies
        #
        self._usePattern = re.compile(
            r"^\s*USE\s+([^,\s]+)", flags=re.IGNORECASE
        )
        self.__openmp_use_pattern = re.compile(
            r"^\s*!\$\s+USE\s+([^,\s]+)", flags=re.IGNORECASE
        )
        self._externalPattern = re.compile(
            r"^\s*EXTERNAL\s+([^,\s]+(:?\s*,\s*[^,\s]+)*)", flags=re.IGNORECASE
        )
        self._external_attribute_pattern = re.compile(
            r"^\s*.+?external\s*?::\s*(.+)"
        )
        self._pFUnitPattern = re.compile(
            r'^\s*#\s+\d+\s+".*testSuites.inc"', flags=re.IGNORECASE
        )
        self._suitePattern = re.compile(r"^\s*ADD_TEST_SUITE\(\s*(\S+)\)")
        self._dependsPattern = re.compile(
            r"!\s*DEPENDS ON:\s*([^.\s]+)(.o)?", flags=re.IGNORECASE
        )

    ###########################################################################
    def analyse(self, source_filename: Path):
        """
        Scans a Fortran source file and harvest dependency information.

        @param source_filename: Fortran source file to be scanned.
        """
        logger = logging.getLogger(__name__)

        # Perform any necessary preprocessing
        #
        if source_filename.suffix in [".F90", ".X90"]:
            logging.getLogger(__name__).info(
                "  Preprocessing " + str(source_filename)
            )
            preprocess_command = self._fpp
            for path in self.__preprocess_include_paths:
                preprocess_command.append("-I" + str(path))
            for name, macro in self.__preprocess_macros.items():
                if macro:
                    preprocess_command.append(f"-D{name}={macro}")
                else:
                    preprocess_command.append("-D" + name)
            preprocess_command.append(str(source_filename))
            logging.getLogger(__name__).debug(preprocess_command)

            start_time = time()
            preprocessor = subprocess.Popen(
                preprocess_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                encoding="utf-8",
            )
            processed_source, errors = preprocessor.communicate()

            logger.debug(
                "Time to preprocess Fortran source: "
                + str(time() - start_time)
            )
            if preprocessor.returncode:
                logger.error(errors)
                raise subprocess.CalledProcessError(
                    preprocessor.returncode, " ".join(preprocess_command)
                )
        elif source_filename.suffix == ".f90":
            start_time = time()
            with source_filename.open("rt") as sourceFile:
                processed_source = sourceFile.read()
            logger.debug(
                "Time to read Fortran source: " + str(time() - start_time)
            )
        else:
            raise Exception(
                "File doesn't look like a Fortran file: "
                + str(source_filename)
            )

        def add_dependency(
            program_unit: str, prerequisite_unit: str, reverse_link=False
        ) -> None:
            """
            Hides the messiness of adding a dependency to the database.

            @param program_unit: Name of depender unit.
            @param prerequisite_unit: Name of dependee unit.
            @param reverse_link: Should the relationship be flipped?
                                 Used by submodules.
            """
            logger = logging.getLogger(__name__)
            prerequisite_unit = prerequisite_unit.lower()
            if prerequisite_unit in self._ignoreModules:
                logger.info("      - Ignored 3rd party prerequisite")
                return
            if prerequisite_unit in dependencies:
                logger.info("      - Ignore duplicate prerequisite")
                return

            dependencies.append(prerequisite_unit)
            self._database.add_compile_dependency(
                program_unit, prerequisite_unit
            )
            if reverse_link:
                self._database.add_link_dependency(
                    prerequisite_unit, program_unit
                )
            else:  # Normal link
                self._database.add_link_dependency(
                    program_unit, prerequisite_unit
                )

        def lines_of_code(
            source: str,
        ) -> Generator[Tuple[str, str, int], None, None]:
            """
            Reads lines from the file, concatenate at continuation markers and
            split comments off.

            TODO: This is complex. That complexity comes from the need
                  to preserve comments. This is needed to support "depends on"
                  comments. Ergo, once "depends on" is gone we can ignore
                  comments and this becomes a lot simpler.

            @param source: Fortran source code.
            @return:
            """
            line_number = 0
            code = ""
            comment = ""
            continuing = False
            for line in source.splitlines():  # Loop over every line in the
                # source
                line_number += 1
                state = "indent"  # Each line starts in the "indent" state
                index = -1
                code_start = 0
                code_end = len(line)
                comment_start = len(line)
                continuation = False
                for character in line:  # Scan every character in a line
                    index += 1
                    if state == "indent":
                        ##############################################
                        if character == "&":  # Start of continuation line
                            if not continuing:
                                message = (
                                    "Found continuation marker at "
                                    "start of line when there was none "
                                    "ending previous line"
                                )
                                raise Exception(message)
                        elif character == "!":  # Line contains only a comment
                            comment_start = index
                            state = "comment"
                        elif character != " ":  # Start of code located
                            code_start = index
                            state = "code"
                    elif state == "code":
                        ##############################################
                        if character == '"':  # String opened with double quote
                            state = "double"
                        elif character == "'":  # String opened with single
                            # quote
                            state = "single"
                        elif character == "&":  # Line continues on next line
                            code_end = index
                            continuing = True
                            continuation = True
                            state = "continue"
                        elif character == "!":  # The remainder of the line
                            # is a comment
                            comment_start = index
                            state = "comment"
                    elif state == "double":
                        ############################################
                        if character == '"':  # Quoted string has ended
                            state = "code"
                        elif character == "&":  # Line continues on next line
                            code_end = index
                            continuing = True
                            continuation = True
                            state = "continue"
                    elif state == "single":
                        ############################################
                        if character == "'":  # Quoted string has ended
                            state = "code"
                        elif character == "&":  # Line continues on next line
                            code_end = index
                            continuing = True
                            continuation = True
                            state = "continue"
                    elif state == "continue":
                        ##########################################
                        if character == "!":  # There is a comment after the
                            # continuation
                            state = "comment"
                    elif state == "comment":
                        ###########################################
                        if index - 1 < code_end:  # If we have not already
                            # ended the code
                            code_end = index - 1  # Mark it as ended.
                        break

                code += " " + line[code_start:code_end]
                comment += " " + line[comment_start:]

                if line[code_start:code_end].strip() and not continuation:
                    yield (code, comment, line_number)
                    code = ""
                    comment = ""
                    continuing = False

        # Scan file for dependencies.
        #
        logger.info("  Scanning " + str(source_filename))
        self._database.remove_file(source_filename)

        program_unit = None
        modules = []
        dependencies: List[str] = []
        scope_stack = []
        pfunit_driver = False
        for code, comment, line_number in lines_of_code(processed_source):
            match = self._programPattern.match(code)
            if match:
                program_unit = match.group(1).lower()
                logger.info("    Contains program: " + program_unit)
                self._database.add_program(program_unit, source_filename)
                scope_stack.append(("program", program_unit))
                continue

            match = self._modulePattern.match(code)
            if match:
                program_unit = match.group(1).lower()
                logger.info("    Contains module " + program_unit)
                modules.append(program_unit)
                self._database.add_module(program_unit, source_filename)
                scope_stack.append(("module", program_unit))
                continue

            match = self._submodulePattern.match(code)
            if match:
                ancestor_unit = match.group(1)
                if ancestor_unit:
                    ancestor_unit = ancestor_unit.lower()
                parent_unit = match.group(2).lower()
                program_unit = match.group(3).lower()

                message = "{}Contains submodule {} of {}".format(
                    " " * 4, program_unit, parent_unit
                )
                if ancestor_unit:
                    message = message + "({})".format(ancestor_unit)
                logger.info(message)

                self._database.add_submodule(program_unit, source_filename)
                add_dependency(program_unit, parent_unit, True)
                scope_stack.append(("submodule", program_unit))
                continue

            match = self._subroutinePattern.match(code)
            if match and len(scope_stack) == 0:
                # Only if this subroutine is a program unit.
                program_unit = match.group(2).lower()
                logger.info("    Contains subroutine " + program_unit)
                modules.append(program_unit)
                self._database.add_procedure(program_unit, source_filename)
                scope_stack.append(("subroutine", program_unit))
                continue

            match = self._functionPattern.match(code)
            if match and len(scope_stack) == 0:
                # Only if this function is a program unit.
                program_unit = match.group(2).lower()
                logger.info("    Contains function " + program_unit)
                modules.append(program_unit)
                self._database.add_procedure(program_unit, source_filename)
                scope_stack.append(("function", program_unit))
                continue

            match = self._endPattern.match(code)
            if match:
                end_scope = match.group(1)
                if end_scope is not None:
                    end_scope = end_scope.lower()
                end_unit = match.group(2)
                if end_unit is not None:
                    end_unit = end_unit.lower()

                try:
                    begin_scope, begin_unit = scope_stack[-1]

                    if end_scope == begin_scope and end_unit == begin_unit:
                        scope_stack.pop()
                        logger.debug(f"    End {end_scope} {end_unit}")
                    else:
                        logger.debug(
                            f"    Mismatched end {end_scope} {end_unit}"
                        )
                except IndexError:
                    message = (
                        'Mismatched begin/end. Found "{end}" but stack empty'
                    )
                    raise Exception(message.format(end=end_scope))

            match = self._usePattern.match(code)
            if match is not None:
                if program_unit is None:
                    raise Exception("Usage found before program unit.")
                module_name = match.group(1).lower()
                logger.info("    Depends on module " + module_name)
                add_dependency(program_unit, module_name)
                continue

            match = self.__openmp_use_pattern.match(comment)
            if match is not None:
                """
                There is no knowledge of whether OpenMP is actually turned on
                or not. This means that the dependency will always exist even
                if it shouldn't when OpenMP is off.
                """
                if program_unit is None:
                    raise Exception("OpenMP found before program unit.")
                module_name = match.group(1).lower()
                logger.info(f"    Depends on module {module_name} with OpenMP")
                add_dependency(program_unit, module_name)
                continue

            match = self._externalPattern.match(code)
            if match is not None:
                if program_unit is None:
                    raise Exception("External found before program unit.")
                module_names = [
                    moduleName.strip()
                    for moduleName in match.group(1).split(",")
                ]
                for module_name in module_names:
                    if module_name is not None:
                        logger.info("    Depends on external " + module_name)
                        add_dependency(program_unit, module_name)
                continue

            match = self._external_attribute_pattern.match(code)
            if match:
                module_names = [
                    module_name.strip()
                    for module_name in match.group(1).split(",")
                ]
                for module_name in module_names:
                    if module_name is not None:
                        logger.info("    Depends on external " + module_name)
                        add_dependency(program_unit, module_name)
                continue

            match = self._pFUnitPattern.match(code)
            if not pfunit_driver and match is not None:
                if program_unit is None:
                    raise Exception("pFUnit driver found before program unit.")
                logger.info("    Is driver")
                pfunit_driver = True

                start_time = time()
                include_filename = source_filename.parent / "testSuites.inc"
                with include_filename.open("rt") as includeFile:
                    for line in includeFile:
                        match = self._suitePattern.match(line)
                        if match is not None:
                            test_generator_function = match.group(1)
                            test_module = test_generator_function.replace(
                                "_suite", ""
                            )
                            logger.info(
                                "      Depends on module " + test_module
                            )
                            self._database.add_compile_dependency(
                                program_unit, test_module
                            )
                            self._database.add_link_dependency(
                                program_unit, test_module
                            )
                message = "Time to read pFUnit driver include file: {0}"
                logging.getLogger(__name__).debug(
                    message.format(time() - start_time)
                )
                continue

            for match in self._dependsPattern.finditer(comment):
                if program_unit is None:
                    raise Exception("Dependency found before program unit.")
                name = match.group(1).lower()
                logger.info(
                    "    %s depends on call to %s " % (program_unit, name)
                )
                add_dependency(program_unit, name)
