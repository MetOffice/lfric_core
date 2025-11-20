#!/usr/bin/env python3
##############################################################################
# (C) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
import os
import subprocess
from pathlib import Path
import sys
from typing import List

from testframework import MpiTest
import xarray as xr


##############################################################################
class LFRicXiosTest(MpiTest):
    """
    Base for LFRic-XIOS integration tests.
    """

    def __init__(self, command=sys.argv[1], processes=1):
        super().__init__(command, processes)
        self.xios_out: List[XiosOutput] = []
        self.xios_err: List[XiosOutput] = []

    def gen_data(self, source: Path, dest: Path):
        """
        Create input data files from CDL formatted text.
        """
        proc = subprocess.Popen(
            ['ncgen', '-k', 'nc4', '-o', f'{dest}', f'{source}'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            )
        _, err = proc.communicate()
        if proc.returncode != 0:
            raise Exception("Test data generation failed:\n" + f"{err}")
        
    def gen_config(self, config_source: Path, config_out: Path, new_config: dict):
        """
        Create an LFRic configuration namelist.
        """
        config_in = open(config_source, 'r')
        config = config_in.readlines()
        for key in new_config.keys():
            for i in range(len(config)):
                if key in config[i]:
                    config[i] = f"  {key}={new_config[key]}\n"
        config_in.close()

        f = open(config_out, "w")
        for line in config:
            f.write(line)
        f.close()            

    def nc_kgo_check(self, output: Path, kgo: Path):
        """
        Compare output files with nccmp.
        """
        proc = subprocess.Popen(
            ['nccmp', '-Fdm', '--exclude=Mesh2d', '--tolerance=0.000001', f'{output}', f'{kgo}'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            )
        _, err = proc.communicate()

        return proc.returncode, err
    
    def nc_data_match(self, in_file: Path, out_file: Path, varname: str):
        """
        Contextually compare output data.
        """
        ds_in = xr.open_dataset(in_file, engine='netcdf4', decode_timedelta=False)
        ds_out = xr.open_dataset(out_file, engine='netcdf4', decode_timedelta=False)

        comparison_window = [max(min(ds_out['time'].values), min(ds_in['time'].values)),
                            min(max(ds_out['time'].values), max(ds_in['time'].values))]

        ds_in_comp = ds_in.sel(time=slice(comparison_window[0], comparison_window[1]))
        ds_out_comp = ds_out.sel(time=slice(comparison_window[0], comparison_window[1]))

        result = [(ds_in_comp['time'] == ds_out_comp['time']).values.all(),
                (ds_in_comp[varname] == ds_out_comp[varname]).values.all()]

        return all(result)

    def post_execution(self, return_code):
        """
        Cache XIOS logging output for analysis.
        """

        for proc in range(self._processes):
            self.xios_out.append(XiosOutput(f"xios_client_{proc}.out"))
            self.xios_err.append(XiosOutput(f"xios_client_{proc}.err"))


class XiosOutput:
    """
    Simple class to hold XIOS output log information
    """

    def __init__(self, filename):
        self.path: Path = Path(os.getcwd()) / Path(filename)

        with open(self.path, "rt") as handle:
            self.contents = handle.read()

    def exists(self):
        """
        Checks if log output file exists
        """
        return os.path.exists(self.path)
