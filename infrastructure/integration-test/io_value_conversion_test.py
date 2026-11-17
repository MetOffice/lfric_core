#!/usr/bin/env python3
##############################################################################
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""Exercise numeric conversion with the real I/O types and logging module."""

from testframework import Test, TestEngine, TestFailed


class IoValueConversionTest(Test):
    """Check that conversions, overflow warnings, and direct I/O writes work."""

    def test(self, returncode, out, err):
        if returncode != 0:
            raise TestFailed(f"Conversion test failed ({returncode}):\n{out}\n{err}")
        if out.strip() != "Numeric I/O conversion tests passed" or err.strip():
            raise TestFailed(f"Unexpected conversion test output:\n{out}\n{err}")
        return "Numeric I/O conversion tests passed"


if __name__ == "__main__":
    TestEngine.run(IoValueConversionTest())
