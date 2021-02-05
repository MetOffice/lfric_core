##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""Class for storing and operating on metadata fields and their output streams.
"""
from typing import Dict

from entities import Field, FieldGroup, OutputStream, OutputStreamField


class Metadata:
    """Store data for fields and output streams"""

    def __init__(self):
        self._field_groups: Dict[str, FieldGroup] = {}
        self._output_streams: Dict[int, OutputStream] = {}
        self._fields: Dict[str, Field] = {}

    def get_field(self, field_id: str) -> Field:
        """
        :param field_id: Identifier for field
        :return: A Field object with the ID given
        """
        return self._fields[field_id]

    def get_fields(self) -> [Field]:
        """:return: A sorted list of fields"""
        return sorted(self._fields.values(),
                      key=lambda field: field.unique_id)

    def get_field_groups(self) -> [FieldGroup]:
        """:return: A sorted list of fields contained in the field group"""
        return sorted(self._field_groups.values(),
                      key=lambda field_group: field_group.name)

    def get_output_streams(self) -> [OutputStream]:
        """:return: A sorted list of all output streams"""
        return sorted(self._output_streams.values(),
                      key=lambda streams: streams.unique_id)

    def add_output_stream(self, stream: OutputStream):
        """Add an OutputStream object to the collection of output streams"""
        self._output_streams.update({stream.unique_id: stream})

    def add_field_group(self, field_group: FieldGroup):
        """Add a FieldGroup to the collection of field groups"""
        self._field_groups.update({field_group.name: field_group})

    def add_field(self, field: Field):
        """
        Add a Field object to the collection of fields and the ID to the
        given field group
        """
        self._fields.update({field.unique_id: field})
        self._field_groups[field.field_group_id].add_field(field)

    def add_output_stream_field(self, output_stream_field: OutputStreamField,
                                stream_unique_id: int):
        """
        Add an output stream field to the output stream corresponding to the
        given output stream ID
        """
        self._output_streams[stream_unique_id].add_field(output_stream_field)
