#!/usr/bin/env python3
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################


class TestFailed(Exception):
    def __init__(
        self, message, return_code=None, stdout=None, stderr=None, log=None
    ):
        super().__init__(message)
        self._return_code = return_code
        self._stdout = stdout
        self._stderr = stderr
        self._log = log

    def __str__(self):
        string = super().__str__()

        if self._return_code:
            string += "\n** Return code: {code}".format(code=self._return_code)

        if self._stdout:
            string += "\n** Standard out:\n{out}".format(out=self._stdout)

        if self._stderr:
            string += "\n** Standard error:\n{err}".format(err=self._stderr)

        if self._log:
            string += "\n** Log file:\n{log}".format(log=self._log)

        return string
