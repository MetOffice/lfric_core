#!/usr/bin/env python3
##############################################################################
# (C) Crown copyright 2025 Met Office. All rights reserved.
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
from pathlib import Path


###############################################################################
class LfricXiosTimeReadTest(LFRicXiosTest):  # pylint: disable=too-few-public-methods
    """
    Tests the lfric_xios_context_type by creating and destroying it
    """

    def __init__(self, nprocs: int):
        super().__init__(command=[sys.argv[1], "context.nml"], processes=nprocs)
        self.gen_data('temporal_data.cdl', 'lfric_xios_time_read_data.nc')
        self.gen_config( "context.nml", "context.nml", {} )
        self.nprocs = nprocs

    def test(self, returncode: int, out: str, err: str):
        """
        Test the output of the context test
        """

        if returncode == 1:
            raise TestFailed("Unexpected failure of test executable")
        else:
            s_or_no_s = ""
            if self.nprocs > 1:
                s_or_no_s = "s"
            return f"Successful read of time data on {self.nprocs} MPI rank{s_or_no_s}"


##############################################################################
if __name__ == "__main__":
    TestEngine.run(LfricXiosTimeReadTest(1))
    TestEngine.run(LfricXiosTimeReadTest(2))
