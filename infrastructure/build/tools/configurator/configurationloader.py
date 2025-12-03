#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
"""
Generates Fortran source for loading all the configuration namelists.
"""

from pathlib import Path
from typing import List

import jinja2


##############################################################################
class ConfigurationLoader:
    """
    Fortran source to load configuration namelists.
    """

    def __init__(self, module_name: str):
        self._engine = jinja2.Environment(
            loader=jinja2.PackageLoader("configurator", "templates")
        )
        self._module_name = module_name
        self._namelists: List[str] = []

    def add_namelist(self, name: str) -> None:
        """
        Registers a namelist name with the loader.

        :param name: Name to register.
        """
        self._namelists.append(name)

    def write_module(self, module_file: Path) -> None:
        """
        Stamps out the Fortran source.

        :param module_file: Filename to use.
        """
        inserts = {
            "moduleName": self._module_name,
            "namelists": self._namelists,
        }

        template = self._engine.get_template("loader.f90.jinja")
        module_file.write_text(template.render(inserts))
