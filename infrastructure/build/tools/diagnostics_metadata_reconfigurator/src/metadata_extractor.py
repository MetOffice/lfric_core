##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Populates a Metadata object from data in a rose suite and a JSON file
containing immutable metadata
"""

import copy
import hashlib
import json
import re
import logging
from pathlib import Path
from typing import TextIO, Union

from entities import Field, FieldGroup, OutputStream, OutputStreamField
from diagnostics_metadata_collection import Metadata

KEY_VALUE_RE = re.compile(r'^([a-zA-Z_]+)=\'?(\w+)\'?')
FIELD_CONFIG_RE = re.compile(r'^\[field_config:([a-zA-Z_]+):([a-zA-Z_]+)\]$')
FIELD_RE = re.compile(
    r'^([a-zA-Z]+(?:_[a-zA-Z]+)*__[a-zA-Z]+(?:_[a-zA-Z]+)*)=(true|false)$')
ADDITIONAL_INPUT_RE = re.compile(
        r'^[a-zA-Z_]+__[a-zA-Z_]+__([a-zA-Z]+)=(true|false)$'
)
BASE_OUTPUT_RE = re.compile(r'^\[output_stream\(([0-9]+)\)\]$')
OUTPUT_FIELD_RE = re.compile(
    r'^\[output_stream\(([0-9]+)\):field\(([0-9]+)\)\]$')

REQUIRED_METADATA = ['unique_id', 'long_name', 'standard_name',
                     'units', 'grid_ref', 'function_space', 'mesh_id', 'order',
                     'io_driver', 'data_type']
LOGGER = logging.getLogger("reconfigurator.extractor")


def create_md5_checksum(obj) -> str:
    """
    :param obj: Object to hash
    :return: String containing checksum
    """
    obj_str = json.dumps(obj, sort_keys=True)
    checksum_hex = hashlib.md5(obj_str.encode('utf-8')).hexdigest()
    return f'md5: {checksum_hex}'


class MetadataExtractor:
    """
    Get fields and output streams which are configured in a Rose suite and
    extract their immutable metadata from a JSON file
        Usage:
        >>> extractor = MetadataExtractor('path/to/suite/', 'immutable.json')
        >>> extractor.extract_metadata()
    """

    def __init__(self, diagnostic_config_path: str,
                 immutable_data_file: Union[str, 'Path']):
        LOGGER.debug("Initialising metadata object")
        self._metadata = Metadata()
        immutable_data_path = Path(immutable_data_file).expanduser()
        self._immutable_metadata = self._get_immutable_data(
                immutable_data_path)

        self._rose_app_path = Path(diagnostic_config_path).expanduser()

        if not self._rose_app_path.exists():
            raise IOError("Could not find rose-app.conf")

    @staticmethod
    def _get_immutable_data(file_path) -> dict:
        """
        Returns immutable metadata from the JSON file given.
        Validates the MD5 checksum to ensure the file has not been
        modified by hand after it was generated.

        :param file_path: Path to the immutable data JSON file
        :return: Multidimensional dictionary containing immutable metadata
        """
        with open(file_path, 'r') as file_handle:
            LOGGER.debug("Opening immutable metadata")
            file_data = json.load(file_handle)
            if 'checksum' not in file_data:
                raise KeyError("Can't find checksum "
                               "to validate immutable data")

            # generate a checksum from a copy of the data without the checksum
            # and check it against the original
            immutable_metadata = copy.deepcopy(file_data)
            del immutable_metadata['checksum']
            checksum = create_md5_checksum(immutable_metadata)
            if file_data['checksum'] != checksum:
                raise RuntimeError("Immutable data has been modified by hand\n"
                                   f"Expected checksum: "
                                   f"{file_data['checksum']}")
            LOGGER.info("Immutable metadata checksum valid")
            return immutable_metadata

    def _add_immutable_field_metadata(self, field: Field) -> Field:
        """
        Populates the given metadata field with required metadata attributes
        taken from the immutable metadata

        :param field: Object containing field's data
        """
        # get the current field's immutable metadata
        section_name, group_name = field.field_group_id.split('__')
        metadata_section = self._immutable_metadata["meta_data"][section_name][
            "groups"][group_name]["fields"][field.unique_id]

        # add the immutable metadata to the field
        for attr in metadata_section:
            if attr in REQUIRED_METADATA:
                setattr(field, attr, metadata_section[attr])

        # hard-code following attributes for now
        field.grid_ref = 'half_level_face_grid'
        field.order = 0
        field.mesh_id = 1
        field.io_driver = 'WRITE_FIELD_FACE'
        return field

    def _parse_field_config(self, file_pointer: TextIO, section_name: str,
                            group_name: str) -> str:
        """
        Parse the [field_config:...] section of the rose-app the file pointer
        is pointing to and set the active status and checksum status of each
        metadata field according to the values found

        :param file_pointer: IO object for the rose-app.conf file
        :param section_name: String containing the name the metadata
                               section currently being processed
        :param group_name: String containing the name of the metadata group
                            currently being processed
        :return: The line the file pointer is currently pointing to so
                 calling method can continue parsing where this method left off
        """
        field_group_id = section_name + "__" + group_name
        field_group = FieldGroup(field_group_id)
        self._metadata.add_field_group(field_group)

        # process fields within group
        line = file_pointer.readline()

        # stop if file pointer reaches the start of the next section
        while line and line[0] != '[':
            f_match = FIELD_RE.search(line)
            if f_match:
                field_id, active = f_match.groups()
                field = Field(field_id, field_group_id,
                              active=active.lower() == 'true')

                # process additional configuration for field
                # currently only checksum is supported
                line = file_pointer.readline()
                match = ADDITIONAL_INPUT_RE.search(line)

                while line and line != '[' and match:
                    input_name, value = match.groups()

                    if input_name == 'checksum':
                        field.checksum = value.lower() == 'true'
                    else:
                        raise ValueError("Unrecognised additional field input "
                                         f"'{input_name}' found for field "
                                         f"'{field_id}'")
                    line = file_pointer.readline()
                    match = ADDITIONAL_INPUT_RE.search(line)

                field = self._add_immutable_field_metadata(field)
                self._metadata.add_field(field)
            else:
                line = file_pointer.readline()
        return line

    def _parse_output_stream(self, file_pointer: TextIO, stream_id: int) \
            -> str:
        """
        Read through an [output_stream(x)] section of the rose-app and store
        information about the output_stream being declared

        :param file_pointer: IO object pointing to the output_stream(x) section
        :param stream_id: Identifier for current output_stream being processed
        :return: Current line of the file pointer to allow calling method to
                  continue parsing where this method left off
        """
        line = file_pointer.readline()
        stream = OutputStream(stream_id)

        # stop if file pointer reaches the start of the next section
        while line and line[0] != '[':
            # check if line has key-value pair
            match = KEY_VALUE_RE.search(line)

            if match:
                key, value = match.groups()
                recognised_keys = {'name': 'name', 'timestep': 'timestep'}

                if key in recognised_keys:
                    setattr(stream, recognised_keys[key], value)
                else:
                    LOGGER.warning("Key '%s' unrecognised in "
                                   "[output_stream(%s)]", key, stream_id)
                    LOGGER.debug("Recognised keys are: %s",
                                 ', '.join(recognised_keys.keys()))
            line = file_pointer.readline()
        self._metadata.add_output_stream(stream)
        return line

    def _parse_output_stream_field(self, file_pointer: TextIO, stream_id: int,
                                   output_stream_field_id: int) -> str:
        """
        Read an output_stream(x):field(y) section (in the rose-app) and store
        info about the output stream field in the Metadata object

        :param file_pointer: IO object pointing to the start of an
                              output_stream(x):field(y) section
        :param stream_id: Identifier for current output_stream being processed
        :param output_stream_field_id: Identifier for output_stream field
                                being processed
        :return: Current line of the file pointer to allow calling method to
                  continue parsing where this method left off
        """
        output_stream_field = OutputStreamField(output_stream_field_id)
        line = file_pointer.readline()

        # stop if file pointer reaches the start of the next section
        while line and line[0] != '[':
            # check if line has key-value pair
            match = KEY_VALUE_RE.search(line)

            if match:
                key, value = match.groups()
                recognised_keys = {'id': 'field_ref', 'temporal': 'temporal'}

                if key in recognised_keys:
                    setattr(output_stream_field, recognised_keys[key], value)
                else:
                    LOGGER.warning("Key '%s' unrecognised in "
                                   "[output_stream(%s):field(%s)]",
                                   key, stream_id, output_stream_field_id)
                    LOGGER.debug("Recognised keys are: %s",
                                 ', '.join(recognised_keys.keys()))
            line = file_pointer.readline()

        # Check whether field has been configured
        try:
            self._metadata.get_field(output_stream_field.field_ref)
        except KeyError:
            LOGGER.warning("Output stream field '%s' has no field config",
                           output_stream_field.field_ref)

        self._metadata.add_output_stream_field(output_stream_field, stream_id)
        return line

    def _parse_rose_app(self):
        """
        Read through the rose-app.conf to find which fields are present in the
        current Rose suite configuration, store their data, and store
        information about which output streams they are destined for
        """
        with open(self._rose_app_path, 'r') as file_pointer:
            LOGGER.debug("Opening rose-app.conf")
            line = file_pointer.readline()

            # File pointer is advanced by each individual method as it
            # processes its own section (updated line is returned by each
            # method pointing to start of next section)
            while line:
                # if line is start of a new section
                if line[0] == '[':
                    LOGGER.debug("Parsing section: %s", line[:-1])

                    # matches [field_config:(section):(field_group)]
                    match = FIELD_CONFIG_RE.search(line)
                    if match:
                        line = self._parse_field_config(
                                file_pointer,
                                match.group(1),
                                match.group(2)
                        )
                        continue

                    # matches [output_stream(x)]
                    match = BASE_OUTPUT_RE.search(line)
                    if match:
                        line = self._parse_output_stream(
                                file_pointer,
                                int(match.group(1))
                        )
                        continue

                    # matches [output_stream(x):field(y)]
                    match = OUTPUT_FIELD_RE.search(line)
                    if match:
                        line = self._parse_output_stream_field(
                                file_pointer,
                                int(match.group(1)),
                                int(match.group(2))
                        )
                        continue

                line = file_pointer.readline()

            # Log number of fields in configs
            LOGGER.info("Found %s fields in %s field groups",
                        len(self._metadata.get_fields()),
                        len(self._metadata.get_field_groups()))

            # Log number of fields in output streams
            stream_fields = 0
            for stream in self._metadata.get_output_streams():
                stream_fields += len(stream._stream_fields)
            LOGGER.info("Found %s fields in %s output streams", stream_fields,
                        len(self._metadata.get_output_streams()))

    def extract_metadata(self):
        """
        Populate a Metadata object with data extracted from the rose-app file
        and the immutable data

        :return: A newly populated Metadata object
        """
        self._parse_rose_app()

        return self._metadata
