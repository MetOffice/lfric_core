#!/usr/bin/env python3
###############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
###############################################################################
# Some of the content of this file has been produced with the assistance of
# Met Office Github Copilot Enterprise.
"""Test harness for the field_cache_test integration test.
Validates that a constant-valued test field is correctly created and
printed by the Fortran test executable. Parses the output to extract
field values for each MPI rank and compares against expected values.

"""
from sys import argv as sys_argv
from testframework import MpiTest, TestEngine, TestFailed
from field_output_utils import uniform_field_data, extract_field_and_compare


class Field_CacheTest(MpiTest):
    """Test class for validating field_cache test field creation and output."""

    def __init__(self):
        """Initialize the test case.
        Configures the test to run the Fortran executable with 6 MPI processes.
        """
        super().__init__([sys_argv[1], '64'], processes=6)

    def test(self, returncode: int, out: str, err: str) -> str:
        """Validate test output.
        Parses the program output to extract field values for each MPI rank,
        then compares against expected constant values (all 1.0). Verifies
        that the field was correctly created and printed across all processes.
        Args:
            returncode: Exit code from the Fortran executable
            out: Standard output from the executable
            err: Standard error from the executable
        Returns:
            Success message if test passes
        Raises:
            TestFailed: If returncode != 0 or field values don't match expected
        """
        # Check that the executable exited successfully
        if returncode != 0:
            message = f"Test program failed with exit code: {returncode}"
            raise TestFailed(message,
                             stdout=out,
                             stderr=err)

        # Setup expected field values for comparison
        before_cache_expected = uniform_field_data(6, [1.1] * 27)
        after_change_expected = uniform_field_data(6, [2.2] * 27)

        # Extract and compare field output
        # This parses the output, checks for duplicates, and compares values
        extract_field_and_compare('field_before_caching', before_cache_expected, out, err)
        extract_field_and_compare('field_after_caching', after_change_expected, out, err)
        extract_field_and_compare('field_after_reset', before_cache_expected, out, err)
        return "Field_Cache field test passed"


if __name__ == '__main__':
    TestEngine.run(Field_CacheTest())
