#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
"""
Test the Fortran code for handling command line arguments. This can't be done
using the unit test framework due to the requirement to interfact with the
command line.
"""


import sys

from testframework import Test, TestEngine, TestFailed


##############################################################################
class cli_mod_normal_test(Test):
    """
    Tests the case where everything is normal and a filename is passed on the
    command line.
    """

    def __init__(self):
        self._INJECT = "onwards/waffles.nml"
        super().__init__([sys.argv[1], self._INJECT])

    def test(self, returncode, out, err):
        if returncode != 0:
            raise TestFailed(
                "Unexpected failure of test executable: {code}".format(
                    code=returncode
                )
            )

        if out.strip() != self._INJECT:
            raise TestFailed(
                'Expected filename "{expected}" but found "{found}"'.format(
                    expected=self._INJECT, found=out.strip()
                )
            )

        return "Filename extracted from command line"


##############################################################################
class cli_mod_too_few_test(Test):
    """
    Tests the case where nothing is passed on the command line. This should be
    an error as a filename is always needed.
    """

    def test(self, returncode, out, err):
        if returncode == 0:
            raise TestFailed("Unexpected success with no arguments")

        return "Command line with no arguments returned with error"


##############################################################################
class cli_mod_too_many_test(Test):
    """
    Tests the case where more than one argument is passed on the command line.
    This should be an error as only a filename is needed.
    """

    def __init__(self):
        super().__init__(
            [sys.argv[1], "onwards/waffles.nml", "2"]
        )

    def test(self, returncode, out, err):
        if returncode == 0:
            raise TestFailed("Unexpected success with 2 arguments")

        return "Command line with 2 arguments returned with error"


##############################################################################
class cli_mod_help_test(Test):
    """
    Tests the case where the "-help" switch is specified.
    """

    def __init__(self):
        super().__init__([sys.argv[1], "-help"])

    def test(self, returncode, out, err):
        if returncode == 0:
            raise TestFailed('Unexpected success with "-help"')

        return 'Command line with "-help" returned an error'


##############################################################################
class cli_mod_h_test(Test):
    """
    Tests the case where the "-h" switch is specified.
    """

    def __init__(self):
        super().__init__([sys.argv[1], "-h"])

    def test(self, returncode, out, err):
        if returncode == 0:
            raise TestFailed('Unexpected success with "-h"')

        return 'Command line with "-h" returned an error'


##############################################################################
class cli_mod_other_help_test(Test):
    """
    Tests the case where the "--help" switch is specified.
    """

    def __init__(self):
        super().__init__([sys.argv[1], "--help"])

    def test(self, returncode, out, err):
        if returncode == 0:
            raise TestFailed('Unexpected success with "--help"')

        return 'Command line with "--help" returned an error'


##############################################################################
if __name__ == "__main__":
    TestEngine.run(cli_mod_normal_test())
    TestEngine.run(cli_mod_too_few_test())
    TestEngine.run(cli_mod_too_many_test())
    TestEngine.run(cli_mod_help_test())
    TestEngine.run(cli_mod_h_test())
    TestEngine.run(cli_mod_other_help_test())
