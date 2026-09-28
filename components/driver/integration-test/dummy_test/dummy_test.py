#!/usr/bin/env python3
###############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
###############################################################################
# Some of the content of this file has been produced with the assistance of
# Met Office Github Copilot Enterprise.

"""Test harness for the dummy_test integration test.

Validates that a constant-valued test field is correctly created and
printed by the Fortran test executable. Parses the output to extract
field values for each MPI rank and compares against expected values.

TODO: Replace this with a proper test of the field caching system.
"""

from re import compile as re_compile
from sys import argv as sys_argv
from testframework import MpiTest, TestEngine, TestFailed
from typing import Dict

class DummyFieldTest(MpiTest):
    """Test class for validating dummy test field creation and output."""

    def __init__(self):
        """Initialize the test case.

        Configures the test to run the Fortran executable with 6 MPI processes
        and sets up a regular expression pattern to parse field values from
        the program output. The pattern extracts rank number and field values.
        """
        super().__init__([sys_argv[1], '64'], processes=6)
        field_to_test = 'dummy_field'

        # Pattern to match output like: "dummy_field (0) = 1.0 1.0 1.0 ..."
        self.__pattern = re_compile(rf'{field_to_test} \((\d+)\) = (.*)')

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
        # Expected field values: constant 1.0 across 27 degrees of freedom per rank
        expected = {
            0: [1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000],
            1: [1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000],
            2: [1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000],
            3: [1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000],
            4: [1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000],
            5: [1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000]
        }

        # Check that the executable exited successfully
        if returncode != 0:
            message = f"Test program failed with exit code: {returncode}"
            raise TestFailed(message,
                             stdout=out,
                             stderr=err)

        # Parse output to extract field values for each rank
        result: Dict[float] = {}
        for line in out.splitlines():
            match = self.__pattern.match(line)
            if match:
                # Convert matched values to list of floats
                result[int(match.group(1))] = [float(value)
                                               for value
                                               in match.group(2).split()]

        # Verify field values match expected values for each rank
        for rank in range(0, 5):
            if result[rank] != expected[rank]:
                message = ("Dummy field does not match expectation for rank "
                           + str(rank))
                raise TestFailed(message, stdout=out, stderr=err)

        return "Dummy field test passed"


if __name__ == '__main__':
    TestEngine.run(DummyFieldTest())