#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
"""
Turns namelist descriptions into namelist modules.
"""

import collections
import json
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
from zlib import crc32

import jinja2

from configurator import jinjamacros


##############################################################################
class NamelistDescriptionException(Exception):
    """
    Thrown for problems in the namelist.
    """

    pass  # pylint: disable=unnecessary-pass


##############################################################################
class FortranType:
    """
    Represents a Fortran type.

    Implements the singleton pattern such that there is only one object per
    type.
    """

    _singletonMap: Dict[str, Dict[str, Dict[str, "FortranType"]]] = {}

    def __init__(self, intrinsic_type: str, kind: str, write_format: str):
        """
        :param intrinsic_type: One of "integer", "real", etc.
        :param kind: Name of data type kind.
        :param write_format: Formatting string for this type.
        """
        self.intrinsic_type = intrinsic_type
        self.kind = kind
        self.write_format = write_format

    def declaration(self) -> str:
        """
        Gets the type designator used by declarations in source files.
        """
        return f"{self.intrinsic_type}({self.kind})"

    def label(self) -> str:
        """
        Gets a label for this type.
        """
        return f"{self.intrinsic_type}_{self.kind}"

    def __lt__(self, other):
        return self.declaration() < other.declaration()

    def __eq__(self, other):
        return self.declaration() == other.declaration()

    def __key(self):
        return (self.intrinsic_type, self.kind, self.write_format)

    def __hash__(self):
        return hash(self.__key())

    @classmethod
    def instance(cls, intrinsic_type, kind, write_format) -> "FortranType":
        """
        Gets the singleton object for a given type.
        """
        if intrinsic_type not in cls._singletonMap:
            cls._singletonMap[intrinsic_type] = {}

        if kind not in cls._singletonMap[intrinsic_type]:
            cls._singletonMap[intrinsic_type][kind] = {}

        if write_format not in cls._singletonMap[intrinsic_type][kind]:
            cls._singletonMap[intrinsic_type][kind][write_format] = cls(
                intrinsic_type, kind, write_format
            )

        return cls._singletonMap[intrinsic_type][kind][write_format]


##############################################################################
class _Property(ABC):
    """
    Root of all namelist fields.

    .. todo:: This interface is used externally so shouldn't be "private."
    """

    def __init__(self, name: str, fortran_type: FortranType):
        """
        :param name: Identifying name.
        :param fortran_type: field's Fortran type.
        """
        self.name = name
        self.fortran_type = fortran_type

    def required_kinds(self) -> List[str]:
        """
        Gets Fortran kind of this field.
        """
        return [self.fortran_type.kind]

    @abstractmethod
    def get_configure_type(self) -> str:
        """
        Gets the configuration meta-data type of this field.
        """
        raise NotImplementedError()

    @property
    @abstractmethod
    def missing_data_indicator(self) -> str:
        """
        Gets the value used to indicate an unset field.
        """
        raise NotImplementedError()


##############################################################################
class _String(_Property):
    """
    Namelist string field.
    """

    _fortranStringMap = {"default": "str_def", "filename": "str_max_filename"}

    def __init__(self, name: str, length: Optional[str] = None):
        """
        :param name: Identifying name.
        :param length: String length is a name which resolves to a length.
        """
        if not length:
            length = "default"

        super().__init__(
            name,
            FortranType.instance(
                "character", self._fortranStringMap[length], "A"
            ),
        )

    def get_configure_type(self) -> str:
        return "string"

    @property
    def missing_data_indicator(self) -> str:
        return "cmdi"


##############################################################################
class _Enumeration(_Property):
    """
    Namelist enumeration field.
    """

    def __init__(self, name: str, keyDictionary: Dict[str, int]):
        """
        :param name: Identifying name.
        :param keyDictionary: Mapping of enumerator to representation.
        """
        super().__init__(name, FortranType.instance("integer", "i_def", "I0"))

        self.mapping = keyDictionary
        self.inverse_mapping = {
            value: key for key, value in self.mapping.items()
        }
        self.first_key = self.inverse_mapping[min(self.inverse_mapping.keys())]

    def required_kinds(self):
        return [self.fortran_type.kind, "str_def"]

    def get_configure_type(self):
        return "enumeration"

    @property
    def missing_data_indicator(self):
        return "emdi"


##############################################################################
class _Scalar(_Property):
    """
    Namelist scalar value field.
    """

    _fortranKindMap = {
        "character": {"default": "str_def", "filename": "str_max_filename"},
        "logical": {"default": "l_def", "native": "l_native"},
        "integer": {
            "default": "i_def",
            "short": "i_short",
            "medium": "i_medium",
            "long": "i_long",
        },
        "real": {
            "default": "r_def",
            "native": "r_native",
            "single": "r_single",
            "double": "r_double",
            "second": "r_second",
        },
    }

    _fortranFormatMap = {
        "character": "A",
        "logical": "L2",
        "integer": "I0",
        "real": "E14.7",
    }

    _fortranMissingDataIndicator = {
        "character": "cmdi",
        "logical": ".false.",
        "integer": "imdi",
        "real": "rmdi",
    }

    def __init__(
        self,
        name: str,
        configure_type: str,
        configure_kind: Optional[str] = None,
    ):
        """
        :param name: Identifying name.
        :param configure_type: Configuration type identifier.
        :param configure_kind: Configuration kind identifier.
        """
        if not configure_kind:
            configure_kind = "default"

        if configure_type == "string":
            configure_type = "character"

        super().__init__(
            name,
            FortranType.instance(
                configure_type,
                self._fortranKindMap[configure_type][configure_kind],
                self._fortranFormatMap[configure_type],
            ),
        )
        self._mdi = self._fortranMissingDataIndicator[configure_type]

    def get_configure_type(self):
        return "scalar"

    @property
    def missing_data_indicator(self):
        return self._mdi


##############################################################################
class _Computed(_Scalar):
    """
    Namelist computed value field.
    """

    def __init__(
        self,
        name: str,
        configure_type: str,
        computation: str,
        configure_kind: Optional[str],
        dereferenced_list_vars: Optional[Sequence[str]] = None,
    ):
        # pylint: disable=too-many-arguments
        """
        :param name: Identifying name.
        :param configure_type: Configuration type identifier.
        :param configure_kind: Configuration kind identifier.
        :param computation: Fortran expression.
        :param derefernced_list_vars: Fields needed from other namelists.
        """
        super().__init__(name, configure_type, configure_kind)
        self.computation = computation
        self.dereferenced_list_vars = dereferenced_list_vars

    def get_configure_type(self):
        return "computed"


##############################################################################
class _Array(_Property):
    """
    Namelist array field.
    """

    def __init__(self, name: str, contentProperty: _Property, bounds: str):
        """
        :param name: Identifying name.
        :param contentProperty: Description of array elements.
        :param bounds: Description of array size.
        """
        super().__init__(name, contentProperty.fortran_type)
        self.content = contentProperty

        if "," in bounds:
            message = "Only 1D arrays allowed in configuration: {}"
            raise NamelistDescriptionException(message.format(bounds))

        if ":" in bounds and bounds.strip() != ":":
            lower, upper = bounds.split(":")

            if lower.strip() not in ["1", ""]:
                message = (
                    "Only lower bound of 1 is allowed in configuration: {}"
                )
                raise NamelistDescriptionException(message.format(bounds))

            self.bounds = upper
        else:
            self.bounds = bounds

    def get_configure_type(self):
        return "array"

    @property
    def missing_data_indicator(self):
        return self.content.missing_data_indicator

    def is_immediate_size(self) -> bool:
        """
        :return: True if array size is a fixed number.
        """
        if self.bounds.isdigit():
            return True

        return False

    def is_deferred_size(self):
        """
        :return: True if array size is dependent on another field.
        """
        if not self.bounds[0].isdigit() and self.bounds[0] != ":":
            return True

        return False

    def is_arbitrary_size(self):
        """
        :return: True if array size is unspecified.
        """
        if self.bounds[0] == ":":
            return True

        return False


##############################################################################
class NamelistDescription:
    """
    Describes a namelist and its contained fields.
    """

    def __init__(
        self,
        listname: str,
        multiple_instances_allowed: bool = False,
        instance_key_member: Optional[str] = None,
    ):
        """
        :param listname: Identifying name.
        """
        self._listname = listname
        self._multiple_instances_allowed = multiple_instances_allowed
        self._instance_key_member = instance_key_member

        self._engine = jinja2.Environment(
            loader=jinja2.PackageLoader("configurator", "templates"),
            extensions=["jinja2.ext.do"],
        )
        self._engine.filters["decorate"] = jinjamacros.decorate_macro

        self._parameters: Dict[str, _Property] = collections.OrderedDict()
        self._module_usage = collections.defaultdict(set)
        self._module_usage["constants_mod"] = set(
            ["cmdi", "emdi", "unset_key", "imdi", "rmdi", "str_def"]
        )

    def get_namelist_name(self) -> str:
        """
        :return: Namelist identifier.
        """
        return self._listname

    def get_module_name(self) -> str:
        """
        :return: Namelist loader Fortran module name.
        """
        return self._listname + "_config_mod"

    def add_enumeration(self, name: str, enumerators: Sequence[str]) -> None:
        """
        Adds an enumerated field to the namelist.

        .. warning::
            This routine will becomes stuck in an infinite loop if asked
            to handle an enumeration with 2^31 enumerators.

        :param name: Identifying name.
        :param enumerators:
        """
        if not isinstance(enumerators, list):
            message = "Expected list of enumerators"
            raise NamelistDescriptionException(message)

        key_dict: Dict[str, int] = collections.OrderedDict()
        for key in enumerators:
            # Hash collisions are always possible and uniqueness is essential
            # for our enumerators. This is a simple way of ensuring that
            # uniqueness. Obviously it will get in an infinite loop if there
            # are more than 2^32 things to deal with but that seems unlikely.
            #
            # Furthermore everything is limited to 2^31 as Fortran integers are
            # always signed.
            #
            value = crc32(bytes(name + key, encoding="ascii")) & 0x7FFFFFFF
            while value in key_dict.values():
                value = (value + 1) & 0x7FFFFFFF
            key_dict[key] = value

        self._parameters[name] = _Enumeration(name, key_dict)

    def add_usage(self, name: str, module: str) -> None:
        """
        Makes this namelist loading module depend on another Fortran module
        for values used in computed fields.

        :param name: Variable name.
        :param module: Module name.
        """
        self._module_usage[module].add(name)

    def add_string(
        self,
        name: str,
        configure_string_length: Optional[str] = None,
        bounds: Optional[str] = None,
    ) -> None:
        """
        Adds a scalar or array string field to the namelist.

        :param name: Field name.
        :param configure_string_length: Length of string is a label which
            resolves to a length.
        :param bounds: Either a length, slice or naked colon.
        """
        new_parameter = _String(name, configure_string_length)

        if bounds:
            dereffed_bounds, _ = self._dereference_expression(bounds)
            self._parameters[name] = _Array(
                name, new_parameter, dereffed_bounds
            )
        else:
            self._parameters[name] = new_parameter

    def add_value(
        self,
        name: str,
        configure_type: str,
        configure_kind: Optional[str] = None,
        bounds: Optional[str] = None,
    ) -> None:
        """
        Adds a scalar or array field of type logical, integer or real to the
        namelist.

        :param name: Field name.
        :param configure_type: type identifier.
        :param configure_kind: kind identifier.
        :param bounds: Either a length, slice or naked colon.
        """
        new_parameter = _Scalar(name, configure_type, configure_kind)
        if bounds:
            dereffed_bounds, _ = self._dereference_expression(bounds)
            self._parameters[name] = _Array(
                name, new_parameter, dereffed_bounds
            )
        else:
            self._parameters[name] = new_parameter

    def add_computed(
        self,
        name: str,
        configure_type: str,
        calculation: str,
        configure_kind: Optional[str] = None,
    ) -> None:
        """
        Adds a computed field to the namelist.

        :param name: Field name.
        :param configure_type: type identifier.
        :param configure_kind: kind identifier.
        :param colculation: Fortran expression.
        """
        calculation, dereferenced_list_vars = self._dereference_expression(
            calculation
        )
        self._parameters[name] = _Computed(
            name,
            configure_type,
            calculation,
            configure_kind,
            dereferenced_list_vars=dereferenced_list_vars,
        )

    def get_parameters(self) -> List[_Property]:
        """
        Gets all the properties associated with this namelist.
        """
        return list(self._parameters.values())

    def write_module(self, file_object: Path) -> None:
        """
        Generates Fortran module source and writes it to a file.

        :param file_object: Filename to write to.
        """
        if not self._parameters:
            message = (
                "Cannot write a module to load an empty namelist ("
                + self._listname
                + ")"
            )
            raise NamelistDescriptionException(message)

        all_kinds = set(["i_def"])
        lone_kind_index = {}
        lone_kind_tally: Dict[FortranType, int] = collections.defaultdict(int)
        namelist = []

        for name, parameter in self._parameters.items():
            all_kinds.update(parameter.required_kinds())

            if not isinstance(parameter, _Computed) and not isinstance(
                parameter, _Array
            ):
                lone_kind_tally[parameter.fortran_type] += 1
                lone_kind_index[name] = lone_kind_tally[parameter.fortran_type]

            if not isinstance(parameter, _Computed):
                namelist.append(parameter.name)

        inserts = {
            "all_kinds": all_kinds,
            "arrays": [
                parameter.name
                for parameter in self._parameters.values()
                if isinstance(parameter, _Array)
            ],
            "allocatables": [
                parameter.name
                for parameter in self._parameters.values()
                if (
                    isinstance(parameter, _Array)
                    and not parameter.is_immediate_size()
                )
            ],
            "enumerations": [
                parameter.name
                for parameter in self._parameters.values()
                if isinstance(parameter, _Enumeration)
            ],
            "listname": self._listname,
            "multiple_instances_allowed": self._multiple_instances_allowed,
            "instance_key_member": self._instance_key_member,
            "lonekindindex": lone_kind_index,
            "lonekindtally": lone_kind_tally,
            "namelist": namelist,
            "parameters": self._parameters,
            "use_from": self._module_usage,
        }

        template = self._engine.get_template("namelist.f90.jinja")
        file_object.write_text(template.render(inserts))

    def _dereference_expression(
        self, expression: str
    ) -> Tuple[str, List[str]]:
        """
        Resolve field references in an expression.

        :param expression: Fortran expression containing field references.
        :result: Expression with references resolved and a list of namelist
                 fields involved.
        """
        str_dict = {
            "namelist": {
                "regexString": r"namelist:(\w*)=(\w*)",
                "removalString": r"namelist:\w*=",
                "moduleSuffix": "_config_mod",
            },
            "source": {
                "regexString": r"source:(\w*)=(\w*)",
                "removalString": r"source:\w*=",
                "moduleSuffix": "",
            },
        }
        result = expression

        dereferenced_list_vars: List[str] = []

        for key, value in str_dict.items():
            use_variables = re.findall(value["regexString"], result)
            if use_variables is not None:
                n_vars = len(use_variables)

                for i_var in range(0, n_vars):
                    list_name = use_variables[i_var][0]
                    var_name = use_variables[i_var][1]

                    if use_variables[i_var][0] != self._listname:
                        module_name = f"{list_name}{value['moduleSuffix']}"
                        self.add_usage(var_name, module=module_name)

                    if key == "namelist":
                        dereferenced_list_vars.append(var_name)

            result = re.sub(value["removalString"], "", result)

        if len(dereferenced_list_vars) == 0:
            dereferenced_list_vars = []

        return result, dereferenced_list_vars

    def add_member(self, member_name: str, meta_dict: Dict[str, str]) -> None:
        # pylint: disable=too-many-branches
        """
        Processes one field entry from the metadata and adds the appropriate
        property to this namelist.

        :param member_name: Identifying name.
        :param meta_dict: Field description.
        """
        meta_keys = list(meta_dict.keys())
        string_length: Optional[str] = None
        xtype: str = ""
        xkind: Optional[str] = None
        xbounds: Optional[str] = None

        if "string_length" in meta_keys:
            string_length = meta_dict["string_length"]

        if "kind" in meta_keys:
            xkind = meta_dict["kind"]

        if "type" in meta_keys:
            xtype = meta_dict["type"]
            if isinstance(xtype, str):
                xtype = xtype.replace("character", "string")

        elif (
            "enumeration" not in meta_keys
            or meta_dict["enumeration"] == "false"
        ):
            message = (
                "namelist:"
                + self._listname
                + "="
                + member_name
                + ": Non-enumeration metadata requires "
                + "a type definition"
            )
            raise NamelistDescriptionException(message)

        # Determining array bounds if any.
        if "length" in meta_keys:
            xlength = meta_dict["length"]

            if xlength == ":":
                if "bounds" in meta_keys:
                    xbounds = meta_dict["bounds"]
                else:
                    xbounds = ":"

            elif isinstance(int(xlength), int):
                xbounds = xlength

        # Generating Enumerators from metadata
        # These are not dependant on xtype being specified
        if "enumeration" in meta_keys and meta_dict["enumeration"] == "true":
            key_values = meta_dict["values"]
            if all(isinstance(item, str) for item in key_values):
                key_values = key_values.replace("\n", "")
                key_values = key_values.replace(" ", "")
                key_values = key_values.replace("'", "")
                keys = key_values.split(",")

                enumeration_keys = [
                    re.sub(r"namelist:", "", member) for member in keys
                ]

                self.add_enumeration(member_name, enumerators=enumeration_keys)

        # Check to see if member is a derived variable
        elif "expression" in meta_keys:
            expression_string = meta_dict["expression"]
            self.add_computed(
                member_name,
                xtype,
                configure_kind=xkind,
                calculation=expression_string,
            )

        elif xtype == "string":
            self.add_string(
                member_name,
                configure_string_length=string_length,
                bounds=xbounds,
            )
        else:
            self.add_value(
                member_name, xtype, configure_kind=xkind, bounds=xbounds
            )


###############################################################################
class NamelistConfigDescription:  # pylint: disable=too-few-public-methods
    """
    Manages the JSON representation of the configuration metadata.
    """

    @staticmethod
    def process_config(nml_config_file: Path) -> List[NamelistDescription]:
        """
        Loads the file and dissects it.
        :param nml_config_file: Input JSON file.
        """
        with open(nml_config_file, encoding="utf8") as config_file:
            namelist_config = json.load(config_file)

        result = []

        for listname in namelist_config.keys():
            multiple_instances_allowed = False
            instance_key_member = None
            if (
                "multiple_instances_allowed"
                in namelist_config[listname].keys()
            ):
                multiple_instances_allowed = namelist_config[listname][
                    "multiple_instances_allowed"
                ]
                instance_key_member = namelist_config[listname][
                    "instance_key_member"
                ]

            description = NamelistDescription(
                listname, multiple_instances_allowed, instance_key_member
            )
            members_dict = namelist_config[listname]["members"]

            for member in sorted(members_dict.keys()):
                meta_dict = members_dict[member]
                description.add_member(member, meta_dict)

            result.append(description)

        return result
