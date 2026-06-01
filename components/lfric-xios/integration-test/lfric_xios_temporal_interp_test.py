#!/usr/bin/env python3
##############################################################################
# (C) Crown copyright Met Office. All rights reserved.
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


class LfricXiosFullInterpTest(LFRicXiosTest):
    """
    Tests the LFRic-XIOS temporal reading functionality for a full set of
    non-cyclic data
    """

    def __init__(self):
        super().__init__(command=[sys.argv[1], "cyclic_high_freq.nml"],
                         processes=1)
        self.gen_data('temporal_data.cdl', 'lfric_xios_interp_input.nc')
        self.gen_data('cyclic_high_freq_kgo.cdl', 'cyclic_high_freq_kgo.nc')
        self.gen_config("cyclic_base.nml", "cyclic_high_freq.nml",
                        {"dt": 10.0, "timestep_end": '150'})

    def test(self, returncode: int, out: str, err: str):
        """
        Test the output of the interpolation test
        """

        if returncode != 0:
            print(out)
            raise TestFailed("Unexpected failure of test executable: "
                             f"{returncode}\n stderr:\n {err}")

        self.plot_output(
                Path(self.test_working_dir, 'lfric_xios_interp_input.nc'),
                Path(self.test_working_dir, 'lfric_xios_interp_output.nc'),
                'temporal_field')

        if not self.nc_data_match(
                Path(self.test_working_dir, 'lfric_xios_interp_input.nc'),
                Path(self.test_working_dir, 'lfric_xios_interp_output.nc'),
                'temporal_field'):
            raise TestFailed("Output data does not match input data for same "
                             "time values")

        return "Reading and interpolating data okay..."


class LfricXiosNonSyncInterpTest(LFRicXiosTest):
    """
    Tests the LFRic-XIOS temporal reading functionality for a full set of
    non-cyclic data
    """

    def __init__(self):
        super().__init__(command=[sys.argv[1], "cyclic_high_freq.nml"],
                         processes=1)
        self.gen_data('temporal_data.cdl', 'lfric_xios_interp_input.nc')
        self.gen_data('cyclic_high_freq_kgo.cdl', 'cyclic_high_freq_kgo.nc')
        self.gen_config("cyclic_base.nml", "cyclic_high_freq.nml",
                        {"dt": 10.0,
                         "calendar_start": "2024-01-01 15:03:20",
                         "timestep_end": '30'})

    def test(self, returncode: int, out: str, err: str):
        """
        Test the output of the interpolation test
        """

        if returncode != 0:
            print(out)
            raise TestFailed("Unexpected failure of test executable: "
                             f"{returncode}\n stderr:\n {err}")

        self.plot_output(
                Path(self.test_working_dir, 'lfric_xios_interp_input.nc'),
                Path(self.test_working_dir, 'lfric_xios_interp_output.nc'),
                'temporal_field')

        if not self.nc_data_match(
                Path(self.test_working_dir, 'lfric_xios_interp_input.nc'),
                Path(self.test_working_dir, 'lfric_xios_interp_output.nc'),
                'temporal_field'):
            raise TestFailed("Output data does not match input data for same "
                             "time values")

        return "Reading and interpolating data okay..."


if __name__ == "__main__":

    TestEngine.run(LfricXiosFullInterpTest())
    TestEngine.run(LfricXiosNonSyncInterpTest())
