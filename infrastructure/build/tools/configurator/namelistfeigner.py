#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
"""
Manage the feign functions for namelists used in testing.
"""

import collections
from pathlib import Path
from typing import Dict, List, Sequence

import jinja2

from configurator import jinjamacros
from configurator.namelistdescription import NamelistDescription, _Property


##############################################################################
class NamelistFeigner:
    """
    Generate all the feigner functions.
    """

    def __init__(self, module_name: str):
        """
        :param module_name: Name the generated Fortran module will have.
        """
        self._module_name = module_name

        self._engine = jinja2.Environment(
            loader=jinja2.PackageLoader("configurator", "templates"),
            extensions=["jinja2.ext.do"],
        )
        self._engine.filters["decorate"] = jinjamacros.decorate_macro

        self._namelists: Dict[str, NamelistDescription] = (
            collections.OrderedDict()
        )

    def add_namelist(self, namelists: Sequence[NamelistDescription]):
        """
        Registers namelists with this feigner generator.

        :param namelist: Namelists to add.
        """
        for item in namelists:
            self._namelists[item.get_namelist_name()] = item

    def write_module(self, module_file: Path) -> None:
        """
        Writes Fortran source file containing feign functions.

        :param module_file: Filename to create.
        """
        enumerations = collections.defaultdict(list)
        kinds = set(["i_def"])
        namelists: List[str] = []
        parameters: Dict[str, List[_Property]] = {}
        character_arrays = None
        non_character_arrays = None

        for namelist in self._namelists.values():
            namelists.append(namelist.get_namelist_name())
            parameters[namelist.get_namelist_name()] = []
            for param in namelist.get_parameters():
                if param.get_configure_type() == "enumeration":
                    enumerations[namelist.get_namelist_name()].append(
                        param.name
                    )
                if param.get_configure_type() != "computed":
                    parameters[namelist.get_namelist_name()].append(param)
                    kinds.add(param.fortran_type.kind)
                if param.get_configure_type() == "array":
                    if param.fortran_type.intrinsic_type == "character":
                        character_arrays = True
                    else:
                        non_character_arrays = True

        inserts = {
            "enumerations": enumerations,
            "kinds": kinds,
            "modulename": self._module_name,
            "namelists": namelists,
            "parameters": parameters,
            "string_arrays": character_arrays,
            "non_string_arrays": non_character_arrays,
        }

        template = self._engine.get_template("feign_config.f90.jinja")
        module_file.write_text(template.render(inserts))
