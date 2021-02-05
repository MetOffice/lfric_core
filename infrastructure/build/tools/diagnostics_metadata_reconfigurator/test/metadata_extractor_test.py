##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
import unittest
from pathlib import Path

from metadata_extractor import MetadataExtractor


TEST_DIR = Path(__file__).parent
IMMUTABLE_DATA_PATH = TEST_DIR / Path('input/LFRic_meta_data_test.JSON')
IMMUTABLE_DATA_NO_CHECKSUM_PATH = TEST_DIR / Path(
    'input/LFRic_meta_data_no_checksum.JSON')
IMMUTABLE_DATA_BAD_CHECKSUM_PATH = TEST_DIR / Path(
    'input/LFRic_meta_data_bad_checksum.JSON')
ROSE_SUITE_PATH = TEST_DIR / Path('input/rose-suite/rose-app.conf')
BAD_ADDITIONAL_INPUT_ROSE_SUITE_PATH = TEST_DIR / Path(
    'input/rose-suite-bad-additional-input/rose-app.conf')
BAD_ROSE_SUITE_PATH = TEST_DIR / Path('input/rose-app.conf')


class TestExtractor(unittest.TestCase):
    def test_extractor(self):
        immutable_metadata = {
            "meta_data": {
                "section_name": {
                    "groups": {
                        "field_group_1": {
                            "fields": {
                                "section_name__field_1": {
                                    "_unique_id": "section_name__field_1",
                                    "units": "units_1"},
                                "section_name__field_2": {
                                    "_unique_id": "section_name__field_2",
                                    "units": "units_2"
                                }
                            }
                        },
                        "field_group_2": {
                            "fields": {
                                "section_name__field_3": {
                                    "_unique_id": "section_name__field_3",
                                    "units": "units_3"
                                }
                            }
                        }
                    }
                }
            }
        }
        extractor = MetadataExtractor(ROSE_SUITE_PATH, IMMUTABLE_DATA_PATH)
        assert (extractor._immutable_metadata == immutable_metadata)

    def test_extractor_no_checksum(self):
        with self.assertRaises(KeyError):
            MetadataExtractor(ROSE_SUITE_PATH, IMMUTABLE_DATA_NO_CHECKSUM_PATH)

    def test_extractor_incorrect_checksum(self):
        with self.assertRaises(RuntimeError):
            MetadataExtractor(ROSE_SUITE_PATH,
                              IMMUTABLE_DATA_BAD_CHECKSUM_PATH)

    def test_extractor_no_rose_app_conf(self):
        with self.assertRaises(IOError):
            MetadataExtractor(BAD_ROSE_SUITE_PATH, IMMUTABLE_DATA_PATH)

    def test_extractor_bad_additional_input(self):
        with self.assertRaises(ValueError):
            extractor = MetadataExtractor(BAD_ADDITIONAL_INPUT_ROSE_SUITE_PATH,
                                          IMMUTABLE_DATA_PATH)
            extractor.extract_metadata()
