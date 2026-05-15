#!/usr/bin/env python3
##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
A set of tests which exercise the temporal reading functionality provided by
the LFRic-XIOS component.
"""
from testframework import TestEngine, TestFailed
from xiostest import LFRicXiosTest
from pathlib import Path
import sys

###############################################################################
class LfricXiosFullCyclicTest(LFRicXiosTest):  # pylint: disable=too-few-public-methods
    """
    Tests the LFRic-XIOS temporal reading functionality for a full set of cyclic data
    """

    def __init__(self):
        super().__init__(command=[sys.argv[1], "cyclic_full.nml"], processes=1)
        self.gen_data('temporal_data.cdl', 'lfric_xios_cyclic_input.nc')
        self.gen_config( 'cyclic_base.nml', 'cyclic_full.nml', {} )

    def test(self, returncode: int, out: str, err: str):
        """
        Test the output of the full cyclic test
        """

        if returncode != 0:
            print(out)
            raise TestFailed(f"Unexpected failure of test executable: {returncode}\n" +
                             f"stderr:\n" +
                             f"{err}")

        self.plot_output(Path(self.test_working_dir, 'lfric_xios_cyclic_input.nc'),
                         Path(self.test_working_dir, 'lfric_xios_cyclic_output.nc'),
                         'temporal_field')

        if not self.nc_data_match(Path(self.test_working_dir, 'lfric_xios_cyclic_input.nc'),
                                  Path(self.test_working_dir, 'lfric_xios_cyclic_output.nc'),
                                  'temporal_field'):
            raise TestFailed("Output data does not match input data for same time values")

        return "Reading full set of cyclic data okay..."


class LfricXiosFutureCyclicTest(LFRicXiosTest):  # pylint: disable=too-few-public-methods
    """
    Tests the LFRic-XIOS temporal reading functionality when data is in the future
    """

    def __init__(self):
        super().__init__(command=[sys.argv[1], "cyclic_future.nml"], processes=1)
        self.gen_data('temporal_data.cdl', 'lfric_xios_cyclic_input.nc')
        self.gen_config( 'cyclic_base.nml', 'cyclic_future.nml',
                         {"calendar_start":'2024-01-01 14:55:00'} )

    def test(self, returncode: int, out: str, err: str):
        """
        Test the output of the future cyclic test
        """

        expected_error_code = "ERROR: I/O context must start after data time " \
                              "window when reading cyclic temporal data"

        if returncode == 1:
            errorcode = err.split("\n")[0].split("0:")[1]
            if not errorcode == expected_error_code:
                raise TestFailed("Incorrect error handling of cyclic future data")
        else:
            raise TestFailed("Unexpected non-failure of test executable")

        return "Expected error for future cyclic data reading..."


class LfricXiosPastCyclicTest(LFRicXiosTest):  # pylint: disable=too-few-public-methods
    """
    Tests the LFRic-XIOS temporal reading functionality when data is in the past
    """

    def __init__(self):
        super().__init__(command=[sys.argv[1], "cyclic_future.nml"], processes=1)
        self.gen_data('temporal_data.cdl', 'lfric_xios_cyclic_input.nc')
        self.gen_data('cyclic_past_kgo.cdl', 'cyclic_past_kgo.nc')
        self.gen_config( 'cyclic_base.nml', 'cyclic_future.nml',
                         {"calendar_start":'2025-01-01 14:55:00'} )

    def test(self, returncode: int, out: str, err: str):
        """
        Test the output of the past cyclic test
        """

        if returncode != 0:
            print(out)
            raise TestFailed(f"Unexpected failure of test executable: {returncode}\n" +
                             f"stderr:\n" +
                             f"{err}")
        if not self.nc_data_match(Path(self.test_working_dir, 'cyclic_past_kgo.nc'),
                                  Path(self.test_working_dir, 'lfric_xios_cyclic_output.nc'),
                                  'temporal_field'):
            raise TestFailed("Output data does not match expected values")

        return "Reading full set of cyclic data from the past okay..."


class LfricXiosCyclicHighFreqTest(LFRicXiosTest):  # pylint: disable=too-few-public-methods
    """
    Tests the LFRic-XIOS temporal reading functionality for a full set of
    cyclic data at higher frequency than the input data
    """

    def __init__(self):
        super().__init__(command=[sys.argv[1], "cyclic_high_freq.nml"], processes=1)
        self.gen_data('temporal_data.cdl', 'lfric_xios_cyclic_input.nc')
        self.gen_data('cyclic_high_freq_kgo.cdl', 'cyclic_high_freq_kgo.nc')
        self.gen_config( 'cyclic_base.nml', 'cyclic_high_freq.nml',
                         {"dt":10.0,
                          "timestep_end":'150'} )

    def test(self, returncode: int, out: str, err: str):
        """
        Test the output of the high frequency cyclic test
        """

        if returncode != 0:
            print(out)
            raise TestFailed(f"Unexpected failure of test executable: {returncode}\n" +
                             f"stderr:\n" +
                             f"{err}")

        self.plot_output(Path(self.test_working_dir, 'lfric_xios_cyclic_input.nc'),
                         Path(self.test_working_dir, 'lfric_xios_cyclic_output.nc'),
                         'temporal_field')

        if not self.nc_data_match(Path(self.test_working_dir, 'cyclic_high_freq_kgo.nc'),
                                  Path(self.test_working_dir, 'lfric_xios_cyclic_output.nc'),
                                  'temporal_field'):
            raise TestFailed("Output data does not match expected values")

        return "Reading full set of cyclic data from the past okay..."


class LfricXiosCyclicNonSyncTest(LFRicXiosTest):  # pylint: disable=too-few-public-methods
    """
    Tests the LFRic-XIOS temporal reading functionality when model timesteps do not match data timesteps
    """

    def __init__(self):
        super().__init__(command=[sys.argv[1], "cyclic_non_sync.nml"], processes=1)
        self.gen_data('temporal_data.cdl', 'lfric_xios_cyclic_input.nc')
        self.gen_data('non_sync_kgo.cdl', 'non_sync_kgo.nc')
        self.gen_config( 'cyclic_base.nml', 'cyclic_non_sync.nml',
                         {"dt":10.0,
                          "calendar_start":"2024-01-01 15:03:20",
                          "timestep_end":"30"} )

    def test(self, returncode: int, out: str, err: str):
        """
        Test the output of the non-synchronised cyclic test
        """

        if returncode != 0:
            print(out)
            raise TestFailed(f"Unexpected failure of test executable: {returncode}\n" +
                             f"stderr:\n" +
                             f"{err}")

        self.plot_output(Path(self.test_working_dir, 'lfric_xios_cyclic_input.nc'),
                         Path(self.test_working_dir, 'lfric_xios_cyclic_output.nc'),
                         'temporal_field')

        if not self.nc_data_match(Path(self.test_working_dir, 'non_sync_kgo.nc'),
                                  Path(self.test_working_dir, 'lfric_xios_cyclic_output.nc'),
                                  'temporal_field'):
            raise TestFailed("Output data does not match expected values")

        return "Reading non-synchronised cyclic data okay..."



##############################################################################
if __name__ == "__main__":
    TestEngine.run(LfricXiosFullCyclicTest())
    TestEngine.run(LfricXiosFutureCyclicTest())
    TestEngine.run(LfricXiosPastCyclicTest())
    TestEngine.run(LfricXiosCyclicHighFreqTest())
    TestEngine.run(LfricXiosCyclicNonSyncTest())