#!/usr/bin/env python3
##############################################################################
# (C) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Test framework for LFRic-XIOS integration tests. Provides a base class which
sets up the test environment and provides utility functions for generating
input data, generating configuration files, checking output data against KGO
files, and plotting input and output data for visual comparison.
"""
from pathlib import Path
import os
import subprocess
import sys
import shutil
from typing import List, Optional

from testframework import MpiTest # pylint: disable=no-name-in-module


##############################################################################
class LFRicXiosTest(MpiTest):
    """
    Base for LFRic-XIOS integration tests.
    """

    def __init__( self, command=sys.argv[1],
                        processes:int=1,
                        iodef_file: Optional[Path]="iodef.xml" ):

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
            check=False,
            )
        if proc.returncode != 0:
            raise Exception("Test data generation failed:\n" + f"{proc.stderr}\n")

    def gen_config(self, config_source: Path, config_out: Path, new_config: dict):
        """
        Create an LFRic configuration namelist. Looks for source file
        in resources/configs directory and generates dest file in test working directory.
        """
        filename = Path(self.resources_dir, 'configs', config_source)
        config = filename.read_text(encoding="utf-8").splitlines()
        for key in new_config.keys():
            for i, line in enumerate(config):
                if key in line:
                    if isinstance(new_config[key], str):
                        config[i] = f"  {key}='{new_config[key]}'"
                    else:
                        config[i] = f"  {key}={new_config[key]}"

        Path(self.test_working_dir, config_out).write_text('\n'.join(config) + '\n',
                                                           encoding="utf-8")

    def performTest(self): # pylint: disable=invalid-name ; This needs to be fixed in the base class
        """
        Removes any old log files and runs the executable.
        """

        # Handle iodef file
        self.iodef_file.unlink(missing_ok=True)
        shutil.copy(self.resources_dir / self.iodef_file, self.test_working_dir / "iodef.xml")

        return super().performTest()

    @classmethod
    def nc_kgo_check(cls, output: Path, kgo: Path):
        """
        Compare output files with nccmp.
        """
        proc = subprocess.run(
            ['nccmp',
             '-Fdm',
             '--exclude=Mesh2d,Mesh2d_face_edges,Mesh2d_face_links', # We use a unit test mesh in the integration tests, so the values for these connectivity fields are incorrect. pylint: disable=line-too-long
             '--tolerance=0.000001',
             f'{output}',
             f'{kgo}'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            )

        kgo_check_okay = (proc.returncode == 0)
        if not kgo_check_okay:
            print(f"{proc.stderr}\n")

        return kgo_check_okay

    def plot_output(self, in_file: Path, out_file: Path, varname: str):
        """
        Visually compare input and output data. If the environment variable
        PLOT_TEST_OUTPUT is not set as True, No plot will be generated. This
        routine depends on the matplotlib package.
        """

        if os.environ.get('PLOT_TEST_OUTPUT', False):
            sys.path.append(str((Path(__file__).parent.parent.parent /
                                 "integration-test" / "tools")))
            from plot_output import plot_test_output # pylint: disable=import-outside-toplevel, import-error

            plot_path = self.test_working_dir / f"{type(self).__name__}.png"
            plot_test_output(in_file, out_file, varname, plot_path)

    def post_execution(self, return_code: int): # pylint: disable=unused-argument
        """
        Cache XIOS logging output for analysis.
        """

        for proc in range(self._processes):
            self.xios_out.append(XiosOutput(self.test_working_dir / f"xios_client_{proc}.out"))
            self.xios_err.append(XiosOutput(self.test_working_dir / f"xios_client_{proc}.err"))

        # Return to top level directory
        os.chdir(self.test_top_level)


class XiosOutput: # pylint: disable=too-few-public-methods
    """
    Simple class to hold XIOS output log information
    """

    def __init__(self, filename):
        self.path: Path = Path(filename)

        self.contents = self.path.read_text(encoding="utf-8")

    def exists(self):
        """
        Checks if log output file exists
        """
        return os.path.exists(self.path)
