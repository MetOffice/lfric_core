#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
# Manages a database of dependency information.

import logging
import sqlite3
from abc import ABC, abstractmethod
from collections import defaultdict
from pathlib import Path
from time import time
from typing import Dict, Generator, List, Optional, Tuple


##############################################################################
# Databases throw this exception.
#
class DatabaseException(Exception):
    def __init__(
        self,
        message: str,
        cause: Optional[Exception] = None,
        module: Optional[str] = None,
        filename: Optional[Path] = None,
    ):
        super().__init__(message, cause)
        self.module = module
        self.filename = filename


##############################################################################
# Basic backend database functionality.
#
class _Database(ABC):
    ##########################################################################
    # Makes sure a described table exists in the database.
    #
    # Arguments:
    #   name    - String by which the table will be known.
    #   columns - List of name/type/modifiers tuples.
    #
    @abstractmethod
    def ensure_table(self, name, columns):
        pass

    ##########################################################################
    # Passes an SQL query to the backend.
    #
    # \deprecated Clearly this breaks all encapsulation. It should be replaced
    #             with a series of objects from which to build an abstract
    #             query.
    #
    # Arguments:
    #   query - Potentially multi-line SQL query.
    #
    @abstractmethod
    def query(self, query):
        pass


##############################################################################
# Database backend.
#
class SQLiteDatabase(_Database):
    ##########################################################################
    # Default constructor.
    #
    # Arguments:
    #   filename - The filename of the database.
    #
    def __init__(self, filename: Path):
        super().__init__()

        start_time = time()
        self._database = sqlite3.connect(str(filename), timeout=5.0)
        self._database.row_factory = sqlite3.Row
        message = "Time to initialise database: {0}"
        logging.getLogger(__name__).debug(message.format(time() - start_time))

    ###########################################################################
    # Destructor.
    #
    def __del__(self):
        start_time = time()
        self._database.commit()
        self._database.close()
        message = "Time to finalise database: {0}"
        logging.getLogger(__name__).debug(message.format(time() - start_time))

    ###########################################################################
    # Creates a table if it does not already exist.
    #
    # Arguments:
    #   name    - String by which table will be identified.
    #   columns - List of lists of column definition clauses.
    #             e.g. [['foo', 'integer', 'primary key'],
    #                   ['bar', 'integer']]
    #
    def ensure_table(self, name, columns):
        column_definitions = []
        for columnDetails in columns:
            column_definitions.append(" ".join(columnDetails))
        query = "CREATE TABLE IF NOT EXISTS {} ( {} )"
        start_time = time()
        with self._database:
            column_list = ", ".join(column_definitions)
            self._database.execute(query.format(name, column_list))
        message = "Time to ensure database table: {0} [{1}]"
        logging.getLogger(__name__).debug(
            message.format(time() - start_time, name)
        )

    ###########################################################################
    # Execute an SQL query against the database.
    #
    # Arguments:
    #   query - Either a string containing a single SQL instruction or a list
    #           of strings, each containing a single SQL instruction.
    #
    def query(self, query):
        if isinstance(query, list):
            query = "; ".join(query) + ";"
        query = " ".join(query.split())  # This wheeze collapses whitespace
        start_time = time()
        try:
            cursor = self._database.cursor()
            if query.count(";") > 1:
                cursor.executescript(query)
            else:
                cursor.execute(query)
            return cursor.fetchall()
        except sqlite3.IntegrityError as ex:
            raise DatabaseException("Database error: ", ex)
        finally:
            message = "Time to query database: {0} [{1}]"
            logging.getLogger(__name__).debug(
                message.format(time() - start_time, query)
            )


###############################################################################
# Basic file dependencies.
#
class FileDependencies(object):
    ###########################################################################
    # Default constructor.
    #
    # Arguments:
    #   database - Database object to use as backend.
    #
    def __init__(self, database):
        self._database = database
        self._database.ensure_table(
            "file_dependency",
            (
                ("file", "TEXT", "NOT NULL"),
                ("prerequisite", "TEXT", "NOT NULL"),
            ),
        )

    ###########################################################################
    # Remove a file and all its dependencies from the database.
    #
    # Arguments:
    #   filename - The filename as it appears in the database.
    #
    def remove_file(self, filename):
        query = 'DELETE FROM file_dependency WHERE file="{}"'.format(filename)
        self._database.query(query)

    ###########################################################################
    # Remove all files and dependencies from the database.
    #
    def remove_all_file_dependencies(self):
        query = "DELETE FROM file_dependency"
        self._database.query(query)

    ###########################################################################
    # Add a dependency relationship to the database.
    #
    # Arguments:
    #   filename     - A filename string.
    #   prerequisite - The filename string of a file which the first depends
    #                  on.
    #
    def add_file_dependency(self, filename, prerequisite):
        query = "INSERT INTO file_dependency VALUES ('{}','{}')"
        self._database.query(query.format(filename, prerequisite))

    ###########################################################################
    # Get all the file dependency relationships.
    #
    # Arguments:
    # Return:
    #   A generator yielding (filename, filename) tuples.
    #
    def get_dependencies(self):
        query = "SELECT * FROM file_dependency ORDER BY file"
        result = self._database.query(query)

        last_file = None
        prerequisites = []
        for row in result:
            if row["file"] != last_file:
                if prerequisites:
                    yield Path(last_file), prerequisites
                last_file = row["file"]
                prerequisites = []

            prerequisites.append(Path(row["prerequisite"]))

        if prerequisites:
            yield Path(last_file), prerequisites


###############################################################################
# Fortran dependencies.
#
class FortranDependencies(object):
    ###########################################################################
    # Default constructor.
    #
    # Arguments:
    #   database - Database object to use as backend.
    #
    def __init__(self, database):
        self._database = database

        self._database.ensure_table(
            "fortran_unit_type", [("type", "TEXT", "PRIMARY KEY")]
        )
        self._database.query(
            [
                'INSERT OR IGNORE INTO fortran_unit_type VALUES("program")',
                'INSERT OR IGNORE INTO fortran_unit_type VALUES("module")',
                'INSERT OR IGNORE INTO fortran_unit_type VALUES("submodule")',
                'INSERT OR IGNORE INTO fortran_unit_type VALUES("procedure")',
            ]
        )

        self._database.ensure_table(
            "fortran_dependency_type", [("type", "TEXT", "PRIMARY KEY")]
        )
        self._database.query(
            [
                "INSERT OR IGNORE INTO fortran_dependency_type "
                'VALUES("compile")',
                'INSERT OR IGNORE INTO fortran_dependency_type VALUES("link")',
            ]
        )

        self._database.ensure_table(
            "fortran_program_unit",
            (
                ("unit", "TEXT", "PRIMARY KEY"),
                ("file", "TEXT", "NOT NULL"),
                ("type", "REFERENCES fortran_unit_type(type)"),
            ),
        )
        self._database.ensure_table(
            "fortran_unit_dependency",
            (
                ("unit", "TEXT", "NOT NULL"),
                ("prerequisite", "TEXT", "NOT NULL"),
                ("type", "REFERENCES fortran_dependency_type(type)"),
            ),
        )

    def remove_file(self, filename: Path) -> None:
        """
        Removes a Fortran source file and all its dependencies.

        @param filename: As it appears in the database.
        """
        query = [
            f'''
            CREATE TEMPORARY TABLE _units AS SELECT unit
            FROM fortran_program_unit
            WHERE file="{filename}"
            ''',
            """
            DELETE FROM fortran_unit_dependency
            WHERE unit IN _units
            """,
            "DROP TABLE _units",
            f'DELETE FROM fortran_program_unit WHERE file="{filename}"',
        ]
        self._database.query(query)

    def add_program(self, name: str, filename: Path) -> None:
        """
        Adds a program to the database.

        @param name: Program's program unit.
        @param filename: Source file in which program was found.
        @return:
        """
        # Changes are transacted to ensure other processes can't find the
        # database with half a program.
        query = (
            "INSERT INTO fortran_program_unit VALUES ( '{}', '{}', 'program' )"
        )
        self._database.query(query.format(name, filename))

    def add_module(self, name: str, filename: Path) -> None:
        """
        Adds a module to the database.

        :param name: module's program unit name.
        :param filename: source file in which the modules is found.
        """
        try:
            _ = self._database.query(
                "INSERT INTO fortran_program_unit "
                f"VALUES ( '{name}', '{str(filename)}', 'module' )"
            )
        except DatabaseException as ex:
            raise DatabaseException(
                f"Unable to add module '{name}' from '{filename}': {ex}",
                module=name,
                filename=filename,
            )

    def add_submodule(self, name: str, filename: Path) -> None:
        """
        Adds a sub-module to the database.

        :param name: Sub-modules program unit name.
        :param filename: Source file in which sub-module is found.
        """
        try:
            _ = self._database.query(
                "INSERT INTO fortran_program_unit "
                f"VALUES ( '{name}', '{str(filename)}', 'submodule' )"
            )
        except DatabaseException as ex:
            raise DatabaseException(
                f"Unable to add sub-module '{name}' from '{filename}': {ex}",
                module=name,
                filename=filename,
            )

    def add_procedure(self, name: str, filename: Path) -> None:
        """
        Adds a program-unit procedure to the database.

        @param name: Procedure's program unit name.
        @param filename: Source file in which the procedure was found.
        """
        try:
            query = (
                "INSERT INTO fortran_program_unit VALUES ( '{name}', "
                "'{filename}', 'procedure' )"
            )
            self._database.query(query.format(name=name, filename=filename))
        except DatabaseException as ex:
            new_exception = DatabaseException(
                f"Unable to add procedure '{name}' from '{filename}': {ex}"
            )
            new_exception.module = name
            new_exception.filename = filename
            raise new_exception

    def get_program_units(self) -> List[Tuple[str, Path]]:
        """
        Gets a list of all program units in the database.

        @return: Program unit name and containing file.
        """
        query = "SELECT * FROM fortran_program_unit"
        rows = self._database.query(query)
        return [(row["unit"], Path(row["file"])) for row in rows]

    def add_compile_dependency(self, unit: str, prerequisite: str) -> None:
        """
        Adds a compile dependency to the database.

        @param unit: Depender unit name.
        @param prerequisite: Dependee unit name.
        """
        query = (
            "INSERT INTO fortran_unit_dependency VALUES ( '{}', '{}', "
            "'compile' )"
        )
        self._database.query(query.format(unit, prerequisite))

    def get_compile_prerequisites(self, unit: str) -> List[str]:
        """
        Gets a list of prerequisites for a program unit.
        """
        query = (
            "SELECT prerequisite "
            "FROM fortran_unit_dependency "
            f"WHERE unit='{unit}' AND type='compile'"
        )
        rows = self._database.query(query)
        return [row["prerequisite"] for row in rows]

    def add_link_dependency(self, unit: str, prerequisite: str) -> None:
        """
        Adds a link dependency to the database.

        @param unit: Depender unit name.
        @param prerequisite: Dependee unit name.
        """
        query = (
            "INSERT INTO fortran_unit_dependency VALUES ( '{}', '{}', 'link' )"
        )
        self._database.query(query.format(unit, prerequisite))

    def get_programs(self) -> List[str]:
        """
        Gets all the programs from the database.

        @return: program names
        """
        query = (
            "SELECT unit FROM fortran_program_unit "
            + " WHERE type='program' ORDER BY unit DESC"
        )
        rows = self._database.query(query)
        return [row["unit"] for row in rows]

    def get_modules(self) -> List[Tuple[str, Path]]:
        """
        Gets all the modules from the database.

        @return: Module name and containing file.
        """
        query = (
            "SELECT unit, file FROM fortran_program_unit "
            + " WHERE type='module' OR type='submodule'"
        )
        rows = self._database.query(query)
        return [(row["unit"], Path(row["file"])) for row in rows]

    def get_link_dependencies(
        self, program_unit: str
    ) -> Generator[Tuple[str, Path, str, Path], None, None]:
        """
        Gets program unit dependencies for a program unit.

        This could be done (and was for a while) as a single SELECT statement.
        The problem with that is it becomes hard to know when nothing is
        returned because there is nothing to return and when it happens because
        a dependency is missing from the database.

        @param program_unit: Desired link target name.
        @return: program unit name, containing file,
                 prerequisite unit name, containing file
        """
        # We cache program unit details for performance.
        #
        unit_cache = UnitCache(self._database)

        processed: List[str] = []
        candidates = [program_unit]
        while candidates:
            unit = candidates.pop(0)

            processed.append(unit)
            rows = self._database.query(
                f'''
                SELECT unit, prerequisite FROM fortran_unit_dependency
                WHERE unit="{unit}" and type="link"
                ORDER BY unit
                '''
            )
            for unit, prerequisite in rows:
                if prerequisite in processed:
                    continue
                candidates.append(prerequisite)

                unit_filename, _ = unit_cache.details(unit)
                try:
                    prerequiste_filename, _ = unit_cache.details(prerequisite)
                except DatabaseException:
                    raise DatabaseException(
                        f"Program unit '{unit}' requires '{prerequisite}' "
                        "but it was not found in the database"
                    )
                yield (
                    unit,
                    Path(unit_filename),
                    prerequisite,
                    Path(prerequiste_filename),
                )

    def get_compile_dependencies(
        self, root: Optional[str] = None
    ) -> Generator[Tuple[str, Path, str, str, Path, str], None, None]:
        """
        Gets all program unit dependencies from database for compilation.

        This could be done (and was for a while) as a single SELECT statement.
        The problem with that is it becomes hard to know when nothing is
        returned because there is nothing to return and when it happens because
        a dependency is missing from the database.

        @param: root  Program unit to start descent with. If not specified,
                      all programs are used.

        @return: unit name, unit filename, unit type,
                 prerequisite name, prerequisite filename, prerequisite type
        """
        if root is None:
            units = self.get_programs()
        else:
            units = [root]

        # Sub-modules are rare but need to be checked every time round the
        # dependency loop. As such it makes sense to cache them and save a
        # database access per iteration.
        #
        submodule_cache = SubmoduleCache(self._database)

        # To further save database lookups we cache unit name to details
        # mapping.
        #
        unit_cache = UnitCache(self._database)

        units_seen: List[str] = []
        while units:
            unit = units.pop()

            if unit in units_seen:
                continue
            units_seen.append(unit)

            units.extend(submodule_cache.submodules(unit))

            dep_rows = self._database.query(
                "SELECT unit, prerequisite, type "
                "FROM fortran_unit_dependency "
                f"WHERE (type = 'compile' AND unit = '{unit}')"
            )
            for dep_row in dep_rows:
                unit_file, unit_type = unit_cache.details(dep_row["unit"])

                try:
                    prerequisite_file, prerequisite_type = unit_cache.details(
                        dep_row["prerequisite"]
                    )
                except DatabaseException:
                    raise DatabaseException(
                        "Unable to find prerequisite "
                        f"'{dep_row['prerequisite']}' of '{dep_row['unit']}'"
                    )

                units.append(dep_row["prerequisite"])
                yield (
                    unit,
                    Path(unit_file),
                    unit_type,
                    dep_row["prerequisite"],
                    Path(prerequisite_file),
                    prerequisite_type,
                )


# I'm not sure either of the caches are achieving much additional performance
# but I don't have time to investigate futher.
#
class SubmoduleCache:
    def __init__(self, database: SQLiteDatabase):
        self.__cache: Dict[str, List[str]] = defaultdict(list)
        sub_rows = database.query(
            "SELECT u.prerequisite AS module, u.unit AS submodule "
            "FROM fortran_unit_dependency AS u, fortran_program_unit AS s "
            "WHERE u.unit=s.unit "
            "AND u.type='compile' "
            "AND s.type='submodule'"
        )
        for sub_row in sub_rows:
            self.__cache[sub_row["module"]].append(sub_row["submodule"])

    def submodules(self, unit: str) -> List[str]:
        if unit in self.__cache:
            return self.__cache[unit]
        else:
            return []


class UnitCache:
    def __init__(self, database: SQLiteDatabase):
        self.__database = database
        self.__cache: Dict[str, Tuple[str, str]] = {}

    def details(self, unit: str):
        if unit not in self.__cache:
            unit_rows = self.__database.query(
                "SELECT file, type FROM fortran_program_unit "
                f"WHERE unit='{unit}'"
            )
            if not unit_rows:
                raise DatabaseException(f"Unable to find unit '{unit}'")
            self.__cache[unit] = unit_rows[0]
        return self.__cache[unit]
