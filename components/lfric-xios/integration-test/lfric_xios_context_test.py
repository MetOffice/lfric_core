#!/usr/bin/env python3
##############################################################################
# (C) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Simple test which initialises a minimal `lfric_xios_context_type` object and
then destroys it. This will also create an attached XIOS context.
"""

from testframework import TestEngine, TestFailed
from xiostest import LFRicXiosTest
import sys

###############################################################################
class LfricXiosContextTest(LFRicXiosTest):
    """
    Tests the lfric_xios_context_type by creating and destroying it
    """

    def __init__(self):
        super().__init__(command=[sys.argv[1], "context.nml"], processes=1)
        self.gen_config( "context.nml", "context.nml", {} )

    def test(self, returncode: int, out: str, err: str):
        """
        Test the output of the context test
        """

        if returncode != 0:
            raise TestFailed(
                f"Unexpected failure of test executable: {returncode}\n"
                "stderr:\n"
                f"{err}"
            )

        for xios_out in self.xios_out:
            if not xios_out.exists():
                raise TestFailed("XIOS context log output not found")
            if (
                "-> info : Client side context is finalized\n"
                not in xios_out.contents
            ):
                raise TestFailed("XIOS context not finalised")

        for xios_err in self.xios_err:
            if not xios_err.exists():
                raise TestFailed("XIOS context log err not found")

        return "XIOS context initialised and destroyed successfully"


##############################################################################
if __name__ == "__main__":
    TestEngine.run(LfricXiosContextTest())
