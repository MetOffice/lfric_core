#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################


from sys import exit

from .exception import TestFailed


class TestEngine:
    """
    Handles the running of test cases.

    TODO: This is just a skeleton which provides standardised result reporting.
          In the future it could deal with registering and dispatching tests
          as well.
    """

    @staticmethod
    def run(testcase):
        """
        Runs the test case and reports the result in a standardised form.
        """
        try:
            success = testcase.performTest()
            message = "[PASS] {test}: {message}"
            print(
                message.format(test=type(testcase).__name__, message=success)
            )
        except TestFailed as ex:
            message = "[FAIL] {test}: {message}"
            exit(message.format(test=type(testcase).__name__, message=str(ex)))
