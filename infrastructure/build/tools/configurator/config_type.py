#!/usr/bin/env python3
##############################################################################
# (C) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Generates Fortran source for application specifc configuration object.
"""

from pathlib import Path
from typing import List

import jinja2


##############################################################################
class AppConfiguration:
    """
    Fortran source to load configuration namelists.
    """

    def __init__(self, module_name: str):
        self._engine = jinja2.Environment(
            loader=jinja2.PackageLoader("configurator", "templates")
        )
        self._module_name = module_name
        self._namelists: List[str] = []
        self._duplicates: List[bool] = []

    def add_namelist(self, name:str, duplicate:bool) -> None:
        """
        Registers a namelist name with the loader.

        :param name: Name to register.
        :param duplicate: Is this namelist allowed multiple instances.
        """
        self._namelists.append(name)
        self._duplicates.append(duplicate)

    def write_module(self, module_file: Path) -> None:
        """
        Stamps out the Fortran source.

        :param module_file: Filename to use.
        """
        inserts = {
            "moduleName": self._module_name,
            "namelists": self._namelists,
            "duplicates": self._duplicates,
        }

        template = self._engine.get_template("config_type.f90.jinja")
        module_file.write_text(template.render(inserts))

        for i, duplicate in enumerate(self._duplicates):
             if duplicate:
                 iterator_file = self._namelists[i] + '_nml_iterator_mod.f90'
                 iterator_filepath = module_file.parent.joinpath(iterator_file)
                 template = self._engine.get_template("namelist_iterator_type.f90.jinja")
                 iterator_filepath.write_text(template.render({"listname":self._namelists[i]}))
