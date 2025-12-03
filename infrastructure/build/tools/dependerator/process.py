#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
# Process previously analysed dependency database. For fun and profit!

import logging
from pathlib import Path
from typing import Generator, List, Optional, Set, Tuple

from dependerator.database import FileDependencies, FortranDependencies


###############################################################################
# Process dependency database.
#
class FortranProcessor:
    ###########################################################################
    # Constructor.
    #
    def __init__(
        self,
        database: FortranDependencies,
        object_directory: Path,
        module_directory: Optional[Path],
    ):
        self.__database = database
        self.__object_directory = object_directory
        self.__module_directory = module_directory

    ###########################################################################
    # Examine the program unit dependecies and work out the file dependencies.
    #
    # :param file_store: FileDependencies object to accept computed
    #                    dependencies.
    # :param object_modules: Whether the compiler stores module information in
    #                        object files.

    def determine_compile_file_dependencies(
        self, file_store: FileDependencies, object_modules=False
    ):
        logging.getLogger(__name__).info(
            "Removing old file compile dependencies"
        )

        file_store.remove_all_file_dependencies()

        logging.getLogger(__name__).info(
            "Determining file compile dependencies..."
        )

        # we have 2 types of dependency:
        #
        # 1) module file dependencies: a module's .mod file depends on it's
        # source files .o file (ommitted of course if module information is
        # stored in object files)
        #
        if not object_modules:
            if self.__module_directory is None:
                raise Exception(
                    "Cannot determine compile dependencies when "
                    "no module directory is specified and modules "
                    "are not in object files."
                )

            for module, file_path in self.__database.get_modules():
                module_leafname = module + ".mod"
                module_path = (
                    self.__module_directory
                    / file_path.parent
                    / module_leafname
                )

                source_object_path = (
                    self.__object_directory / file_path.with_suffix(".o")
                )

                file_store.add_file_dependency(module_path, source_object_path)

        # 2) object file dependencies: %.o files depend on .mod files for
        # modules used within %.[fF]90 (or on object files if module
        # information stored in object files)
        #
        for (
            unit,
            unit_path,
            unit_type,
            prereq,
            prereq_path,
            prereq_type,
        ) in self.__database.get_compile_dependencies():
            unit_object_path = self.__object_directory / unit_path.with_suffix(
                ".o"
            )

            prereq_object_path = (
                self.__object_directory / prereq_path.with_suffix(".o")
            )

            if prereq_type == "module":
                message = f"{unit} depends on module {prereq}"
                logging.getLogger(__name__).info(message)

                if object_modules:
                    file_store.add_file_dependency(
                        unit_object_path, prereq_object_path
                    )
                else:  # not object_modules
                    if self.__module_directory is None:
                        raise Exception(
                            "Cannot determine compile dependencies when no "
                            "module directory is specified and modules are "
                            "not in object files."
                        )

                    prereq_leaf = prereq + ".mod"
                    prereq_module_path = (
                        self.__module_directory
                        / prereq_path.parent
                        / prereq_leaf
                    )

                    file_store.add_file_dependency(
                        unit_object_path, prereq_module_path
                    )

            if prereq_type == "procedure":
                message = f"{unit} depends on procedure {prereq}"
                logging.getLogger(__name__).info(message)

                file_store.add_file_dependency(
                    unit_object_path, prereq_object_path
                )

    ###########################################################################
    # Determine all program units needed to build each program.
    #
    # TODO: Once we have a more recent version of SQLite we could consider
    # doing this at the database level.
    #
    def determine_link_dependencies(
        self, root_unit: Optional[str] = None
    ) -> Generator[Tuple[Path, Path, List[Path]], None, None]:
        if root_unit is not None:
            roots = [root_unit]
        else:
            roots = self.__database.get_programs()

        for root in roots:
            logging.getLogger(__name__).info("Root {0}".format(root))

            root_object_file: Optional[Path] = None
            objects: Set[Path] = set()
            for (
                unit,
                unit_file,
                prereq,
                prereq_file,
            ) in self.__database.get_link_dependencies(root):
                if unit == root:
                    root_object_file = (
                        self.__object_directory / unit_file.with_suffix(".o")
                    )
                objects.add(
                    self.__object_directory / prereq_file.with_suffix(".o")
                )

            if root_object_file is None:
                raise Exception(f"Root object '{root}' not found")

            yield (
                self.__object_directory / root,
                root_object_file,
                sorted(objects),
            )
