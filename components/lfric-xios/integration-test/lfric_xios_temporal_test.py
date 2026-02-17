#!/usr/bin/env python3
##############################################################################
# (C) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
A set of tests which exercise the temporal reading functionality provided by
the LFRic-XIOS component.
The tests cover the reading of a piece of non-cyclic temporal data with data
points ranging from 15:01 to 15:10 in 10 1-minute intervals. The model start
time is changed to change how the model interacts with the data.
"""
from testframework import TestEngine, TestFailed
from xiostest import LFRicXiosTest
from pathlib import Path
import sys

###############################################################################
class LfricXiosFullNonCyclicTest(LFRicXiosTest):  # pylint: disable=too-few-public-methods
    """
    Tests the LFRic-XIOS temporal reading functionality for a full set of non-cyclic data
    """

    def __init__(self):
        super().__init__(command=[sys.argv[1], "non_cyclic_full.nml"], processes=1)
        self.gen_data('temporal_data.cdl', 'lfric_xios_temporal_input.nc')
        self.gen_config( "non_cyclic_base.nml", "non_cyclic_full.nml", {} )

    def test(self, returncode: int, out: str, err: str):
        """
        Test the output of the context test
        """

        if returncode != 0:
            print(out)
            raise TestFailed(f"Unexpected failure of test executable: {returncode}\n" +
                             f"stderr:\n" +
                             f"{err}")
        if not self.nc_data_match(Path(self.test_working_dir, 'lfric_xios_temporal_input.nc'),
                                  Path(self.test_working_dir, 'lfric_xios_temporal_output.nc'),
                                  'temporal_field'):
            raise TestFailed("Output data does not match input data for same time values")

        return "Reading full set of non-cylic data okay..."


class LfricXiosNonCyclicHighFreqTest(LFRicXiosTest):  # pylint: disable=too-few-public-methods
    """
    Tests the LFRic-XIOS temporal reading functionality for a full set of
    non-cyclic data at higher frequency than the input data
    """

    def __init__(self):
        super().__init__(command=[sys.argv[1], "non_cyclic_high_freq.nml"], processes=1)
        self.gen_data('temporal_data.cdl', 'lfric_xios_temporal_input.nc')
        self.gen_config( "non_cyclic_base.nml", "non_cyclic_high_freq.nml", {"dt":10.0} )

    def test(self, returncode: int, out: str, err: str):
        """
        Test the output of the context test
        """

        if returncode != 0:
            print(out)
            raise TestFailed(f"Unexpected failure of test executable: {returncode}\n" +
                             f"stderr:\n" +
                             f"{err}")
        if not self.nc_data_match(Path(self.test_working_dir, 'lfric_xios_temporal_input.nc'),
                                  Path(self.test_working_dir, 'lfric_xios_temporal_output.nc'),
                                  'temporal_field'):
            raise TestFailed("Output data does not match input data for same time values")

        return "Reading full set of non-cylic data at higher model frequency okay..."


class LfricXiosPartialNonCyclicTest(LFRicXiosTest):  # pylint: disable=too-few-public-methods
    """
    Tests the LFRic-XIOS temporal reading functionality for a partial set of non-cyclic data
    (starting half-way through)
    """

    def __init__(self):
        super().__init__(command=[sys.argv[1], "non_cyclic_mid.nml"], processes=1)
        self.gen_data('temporal_data.cdl', 'lfric_xios_temporal_input.nc')
        self.gen_config( "non_cyclic_base.nml", "non_cyclic_mid.nml", {'calendar_start':'2024-01-01 15:01:00'} )

    def test(self, returncode: int, out: str, err: str):
        """
        Test the output of the context test
        """

        if returncode != 0:
            raise TestFailed(f"Unexpected failure of test executable: {returncode}\n" +
                             f"stderr:\n" +
                             f"{err}")

        if not self.nc_data_match(Path(self.test_working_dir, 'lfric_xios_temporal_input.nc'),
                                  Path(self.test_working_dir, 'lfric_xios_temporal_output.nc'),
                                  'temporal_field'):
            raise TestFailed("Output data does not match input data for same time values")

        return "Reading partial set of non-cylic data okay..."


class LfricXiosNonCyclicFutureTest(LFRicXiosTest):  # pylint: disable=too-few-public-methods
    """
    Tests the LFRic-XIOS reading for non-cyclic data in the future (expected failure)
    """

    def __init__(self):
        super().__init__(command=[sys.argv[1], "non_cyclic_future.nml"], processes=1)
        self.gen_data('temporal_data.cdl', 'lfric_xios_temporal_input.nc')
        self.gen_config( "non_cyclic_base.nml", "non_cyclic_future.nml", {'calendar_start':'2024-01-01 10:00:00',
                          'calendar_origin':'2024-01-01 10:00:00'} )

    def test(self, returncode: int, out: str, err: str):
        """
        Test the output of the context test
        """

        expected_error_code = "ERROR: Context must start within data time window for non-cyclic temporal data"

        if returncode == 1:
            errorcode = err.split("\n")[0].split("0:")[1]
            if not errorcode == expected_error_code:
                raise TestFailed("Incorrect error handling of non-cyclic future data")
        else:
            raise TestFailed("Unexpected non-failure of test executable")


        return "Expected error for future non-cyclic data reading..."


class LfricXiosNonCyclicPastTest(LFRicXiosTest):  # pylint: disable=too-few-public-methods
    """
    Tests the LFRic-XIOS reading for non-cyclic data in the future (expected failure)
    """

    def __init__(self):
        super().__init__(command=[sys.argv[1], "non_cyclic_past.nml"], processes=1)
        self.gen_data('temporal_data.cdl', 'lfric_xios_temporal_input.nc')
        self.gen_config( "non_cyclic_base.nml", "non_cyclic_past.nml", {'calendar_start':'2024-02-01 10:00:00',
                          'calendar_origin':'2024-02-01 10:00:00'} )

    def test(self, returncode: int, out: str, err: str):
        """
        Test the output of the context test
        """

        expected_error_code = "ERROR: Context must start within data time window for non-cyclic temporal data"

        if returncode == 1:
            errorcode = err.split("\n")[0].split("0:")[1]
            if not errorcode == expected_error_code:
                raise TestFailed("Incorrect error handling of non-cyclic past data")
        else:
            raise TestFailed("Unexpected non-failure of test executable")


        return "Expected error for past non-cyclic data reading..."


##############################################################################
if __name__ == "__main__":
    TestEngine.run(LfricXiosFullNonCyclicTest())
    TestEngine.run(LfricXiosNonCyclicHighFreqTest())
    TestEngine.run(LfricXiosPartialNonCyclicTest())
    TestEngine.run(LfricXiosNonCyclicFutureTest())
    TestEngine.run(LfricXiosNonCyclicPastTest())