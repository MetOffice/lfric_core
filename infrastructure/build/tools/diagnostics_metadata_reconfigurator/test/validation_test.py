##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""Test metadata_validator"""

import logging
import unittest

from entities import Field, FieldGroup
from diagnostics_metadata_collection import Metadata
from metadata_validator import validate_metadata, InvalidMetadataError

logging.disable(logging.CRITICAL)  # disable logging to neaten test output


class ValidationTest(unittest.TestCase):
    def test_validate_metadata_valid(self):
        """
        Check that function validate_metadata does not return an error when
        passed valid diagnostic field metadata
        """

        metadata = Metadata()

        field_group_name = "TestFieldGroup"
        metadata.add_field_group(FieldGroup(name=field_group_name))

        # Fields with valid naming metadata
        valid_field1 = Field(unique_id="valid_field1",
                             field_group_id=field_group_name,
                             long_name='valid_long_name1')
        valid_field2 = Field(unique_id="valid_field2",
                             field_group_id=field_group_name,
                             standard_name="valid_standard_name2")
        valid_field3 = Field(unique_id="valid_field3",
                             field_group_id=field_group_name,
                             long_name='valid_long_name3',
                             standard_name="valid_standard_name")

        metadata.add_field(valid_field1)
        metadata.add_field(valid_field2)
        metadata.add_field(valid_field3)

        try:
            validate_metadata(metadata)
        except InvalidMetadataError:
            self.fail('An unexpected InvalidMetadataError occurred')

    def test_validate_metadata_invalid_names(self):
        """
        Check that metadata is invalid when long name and standard name are
        not specified, or set to either None or an empty string
        """

        metadata = Metadata()

        field_group_name = "TestFieldGroup"
        metadata.add_field_group(FieldGroup(name=field_group_name))

        invalid_field1 = Field(unique_id="invalid_field1",
                               field_group_id=field_group_name,)
        invalid_field2 = Field(unique_id="invalid_field2",
                               field_group_id=field_group_name,
                               long_name="",
                               standard_name="")
        invalid_field3 = Field(unique_id="invalid_field3",
                               field_group_id=field_group_name,
                               long_name=None,
                               standard_name=None)

        valid_field1 = Field(unique_id="valid_field",
                             field_group_id=field_group_name,
                             long_name='valid_long_name')
        valid_field2 = Field(unique_id="valid_field2",
                             field_group_id=field_group_name,
                             standard_name="valid_standard_name")

        metadata.add_field(invalid_field1)
        metadata.add_field(invalid_field2)
        metadata.add_field(invalid_field3)
        metadata.add_field(valid_field1)
        metadata.add_field(valid_field2)

        # Test wording of error message, beware if it is reformatted
        expected = r"\['invalid_field1', 'invalid_field2', 'invalid_field3'\]"
        self.assertRaisesRegex(InvalidMetadataError, expected,
                               validate_metadata, metadata)

    def test_validate_metadata_invalid_names_force(self):
        """
        Check that metadata is invalid when long name and standard name are
        not specified or set to either None or an empty string, but that no
        error when is raised when force is set to True
        """

        metadata = Metadata()

        field_group_name = "TestFieldGroup"
        metadata.add_field_group(FieldGroup(name=field_group_name))

        invalid_field1 = Field(unique_id="invalid_field1",
                               field_group_id=field_group_name,)
        invalid_field2 = Field(unique_id="invalid_field2",
                               field_group_id=field_group_name,
                               long_name="",
                               standard_name="")
        invalid_field3 = Field(unique_id="invalid_field3",
                               field_group_id=field_group_name,
                               long_name=None,
                               standard_name=None)

        valid_field1 = Field(unique_id="valid_field",
                             field_group_id=field_group_name,
                             long_name='valid_long_name')
        valid_field2 = Field(unique_id="valid_field2",
                             field_group_id=field_group_name,
                             standard_name="valid_standard_name")

        metadata.add_field(invalid_field1)
        metadata.add_field(invalid_field2)
        metadata.add_field(invalid_field3)
        metadata.add_field(valid_field1)
        metadata.add_field(valid_field2)

        # Test wording of error message, beware if it is reformatted
        expected = r"\['invalid_field1', 'invalid_field2', 'invalid_field3'\]"
        self.assertRaisesRegex(InvalidMetadataError, expected,
                               validate_metadata, metadata)

        try:
            validate_metadata(metadata, True)
        except InvalidMetadataError:
            self.fail('An unexpected InvalidMetadataError occurred')


if __name__ == '__main__':
    unittest.main()
