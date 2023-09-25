##############################################################################
# (c) Crown copyright 2023 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Crufty tool to detect global variables in Fortran source. This is a stop-gap
to aid the eradication of globals until Stylist can do it properly.
"""
from errno import ENOENT
from logging import getLogger
from os import strerror
from pathlib import Path
from re import compile as re_compile
from typing import Callable, List, Optional, Sequence, Tuple

from fparser.common.readfortran import FortranFileReader  # type: ignore
from fparser.two.Fortran2003 import (Attr_Spec,  # type: ignore
                                     Declaration_Type_Spec,
                                     Entity_Decl,
                                     Intrinsic_Type_Spec,
                                     Main_Program,
                                     Module,
                                     Name,
                                     Program,
                                     Type_Declaration_Stmt)
from fparser.two.parser import ParserFactory  # type: ignore
from fparser.two.utils import get_child, walk  # type: ignore


# pylint: disable=too-few-public-methods
class Dirt:
    """
    Describes a dirty file.
    """
    def __init__(self,
                 line_number: int, fortran_type: str, variable_name: str):
        self.line_number = line_number
        self.fortran_type = fortran_type.lower()
        self.variable_name = variable_name


# pylint: disable=too-few-public-methods
class DirtyFile:
    """
    List of dirty files.
    """
# pylint: disable = redefined-outer-name
    def __init__(self, filename: Path) -> None:
        self.filename = filename
        self.dirt: List[Dirt] = []

    def add_dirt(self, line_number: int, fortran_type: str,
                 variable_name: str):
        """
        Extends the list of dirty files.
        """
        self.dirt.append(Dirt(line_number, fortran_type, variable_name))


def __find_declarations(root: Program,
                        dirty_file: DirtyFile,
                        handlers: Sequence[Callable[[DirtyFile,
                                                     Program,
                                                     int,
                                                     Sequence[str],
                                                     str,
                                                     Sequence[str]], None]]) \
        -> None:
    for declaration in walk(root, Type_Declaration_Stmt):
        if isinstance(declaration.parent.parent, Main_Program):
            continue

        line_number = declaration.item.span[0]
        attributes = [str(attribute).lower()
                      for attribute in walk(declaration, Attr_Spec)]
        intrinsic_type = get_child(declaration, Intrinsic_Type_Spec)
        user_type = get_child(declaration, Declaration_Type_Spec)
        actual_type = intrinsic_type or get_child(user_type, Name)
        names = [str(entity) for entity in walk(declaration, Entity_Decl)]

        for handler in handlers:
            handler(dirty_file, declaration.parent.parent,
                    line_number, attributes, str(actual_type), names)


def __process_file(filename: Path,
                   tree_handler: Sequence[Callable[[DirtyFile,
                                                    Program,
                                                    int,
                                                    Sequence[str],
                                                    str,
                                                    Sequence[str]], None]]) \
        -> Optional[DirtyFile]:
    file_tally = DirtyFile(filename)

    reader = FortranFileReader(str(filename))
    # pylint: disable = redefined-outer-name
    parser = ParserFactory().create(std='f2008')
    tree = parser(reader)
    __find_declarations(tree, file_tally, tree_handler)

    # pylint: disable = no-else-return
    if len(file_tally.dirt) == 0:
        return None
    else:
        return file_tally


def __find_globals(dirty_file: DirtyFile,
                   parent: Program,
                   line_number: int,
                   attributes: Sequence[str],
                   fortran_type: str,
                   names: Sequence[str]) -> None:
    if not isinstance(parent, Module):
        return

    if 'parameter' in attributes:
        return

    for variable_name in names:
        dirty_file.dirt.append(Dirt(line_number, fortran_type, variable_name))


def __find_saved(dirty_file: DirtyFile,
                 parent: Program,
                 line_number: int,
                 attributes: Sequence[str],
                 fortran_type: str,
                 names: Sequence[str]) -> None:
    if isinstance(parent, (Module, Main_Program)):
        return

    if 'save' not in attributes:
        return

    for variable_name in names:
        dirty_file.dirt.append(Dirt(line_number, fortran_type, variable_name))


__LOG_MESSAGE = "{filename}: {fortran_type}: {names}"
__FORTRAN_EXTENSION_PATTERN = re_compile(r'\.[FfXx]90')


# pylint: disable=redefined-outer-name
def entry(file_objects: List[Path]) -> Tuple[List[DirtyFile],
                                             List[Path],
                                             List[Path]]:
    """
    Descend file tree processing files.
    """
    dirty_list: List[DirtyFile] = []  # pylint: disable=redefined-outer-name
    clean_list: List[Path] = []  # pylint: disable=redefined-outer-name
    not_considered: List[Path] = []  # pylint: disable=redefined-outer-name

    while len(file_objects) > 0:
        file_object = file_objects.pop()
        if not file_object.exists():
            raise FileNotFoundError(ENOENT, strerror(ENOENT), file_object)
        if file_object.is_dir():
            getLogger('occupyfortran').debug("Descending into %s", file_object)
            file_objects.extend(file_object.iterdir())
        else:  # Object is a file
            if __FORTRAN_EXTENSION_PATTERN.match(file_object.suffix):
                getLogger('occupyfortran').debug("Processing %s", file_object)
                report = __process_file(file_object,
                                        [__find_globals, __find_saved])
                if report is None:
                    clean_list.append(file_object)
                else:  # File has dirt
                    if report is not None:
                        dirty_list.append(report)
            else:
                getLogger('occupyfortran').debug("Ignoring %s", file_object)
                not_considered.append(file_object)

    return dirty_list, clean_list, not_considered
