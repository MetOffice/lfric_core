#!/usr/bin/env python3
##############################################################################
# (C) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
from pathlib import Path
import os
import subprocess
import sys
import shutil
from typing import List, Optional

from testframework import MpiTest
import xarray as xr


##############################################################################
class LFRicXiosTest(MpiTest):
    """
    Base for LFRic-XIOS integration tests.
    """

    def __init__(self, command=sys.argv[1], processes:int=1, iodef_file: Optional[Path]="iodef.xml"):

        self.iodef_file = Path(iodef_file)

        super().__init__(command, processes)

        self.xios_out: List[XiosOutput] = []
        self.xios_err: List[XiosOutput] = []

        # Setup test working directory
        self.test_top_level = Path(os.getcwd())
        self.resources_dir = self.test_top_level / "resources"
        self.test_working_dir = self.test_top_level / "working" / type(self).__name__
        self.test_working_dir.mkdir(parents=True, exist_ok=True)

        # Create symlink to test executable in working directory
        executable = self.test_working_dir / command[0].split('/')[-1]
        if not executable.exists():
            command_path = Path(command[0])
            executable.symlink_to(command_path)

        # Change to test working directory
        os.chdir(self.test_working_dir)


    def gen_data(self, source: Path, dest: Path):
        """
        Create input data files from CDL formatted text. Looks for source file
        in resources/data directory and generates dest file in test working directory.
        """
        dest_path = Path(self.test_working_dir) / dest
        source_path = Path(self.resources_dir, 'data') / source
        dest_path.unlink(missing_ok=True)

        proc = subprocess.run(
            ['ncgen', '-k', 'nc4', '-o', f'{dest_path}', f'{source_path }'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            )
        if proc.returncode != 0:
            raise Exception("Test data generation failed:\n" + f"{proc.stderr}")


    def gen_config(self, config_source: Path, config_out: Path, new_config: dict):
        """
        Create an LFRic configuration namelist. Looks for source file
        in resources/configs directory and generates dest file in test working directory.
        """
        filename = Path(self.resources_dir, 'configs', config_source)
        config = filename.read_text().splitlines()
        for key in new_config.keys():
            for i in range(len(config)):
                if key in config[i]:
                    if type(new_config[key]) == str:
                        config[i] = f"  {key}='{new_config[key]}'\n"
                    else:
                        config[i] = f"  {key}={new_config[key]}\n"

        Path(self.test_working_dir, config_out).write_text('\n'.join(config) + '\n')


    def performTest(self):
        """
        Removes any old log files and runs the executable.
        """

        # Handle iodef file
        self.iodef_file.unlink(missing_ok=True)
        shutil.copy(self.resources_dir / self.iodef_file, self.test_working_dir / "iodef.xml")

        return super().performTest()


    def nc_kgo_check(self, output: Path, kgo: Path):
        """
        Compare output files with nccmp.
        """
        proc = subprocess.run(
            ['nccmp', '-Fdm', '--exclude=Mesh2d', '--tolerance=0.000001', f'{output}', f'{kgo}'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            )

        return proc.returncode, proc.stderr
    
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
            self.xios_out.append(XiosOutput(self.test_working_dir / f"xios_client_{proc}.out"))
            self.xios_err.append(XiosOutput(self.test_working_dir / f"xios_client_{proc}.err"))

        # Return to top level directory
        os.chdir(self.test_top_level)


class XiosOutput:
    """
    Simple class to hold XIOS output log information
    """

    def __init__(self, filename):
        self.path: Path = Path(filename)

        with open(self.path, "rt") as handle:
            self.contents = handle.read()

    def exists(self):
        """
        Checks if log output file exists
        """
        return os.path.exists(self.path)
