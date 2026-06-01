##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Author: J. Henrichs, Bureau of Meteorology
# Author: J. Lyu, Bureau of Meteorology

'''
This module contains the Fab Templaterator class.
'''

import logging
from pathlib import Path
from typing import Dict, List, Union

from fab.api import Tool

logger = logging.getLogger('fab')


class Templaterator(Tool):
    '''This implements the LFRic templaterator as a Fab tool.
    It can check whether the templaterator is available and
    creates command line options for it to run.

    :param Path exec_name: the path to the templaterator binary.
    '''
    def __init__(self, exec_name: Path):
        # Remove suffix as a name
        super().__init__(exec_name.stem,
                         exec_name=exec_name)

    def check_available(self) -> bool:
        '''
        :returns bool: whether templaterator works by running
            `Templaterator -help`.
        '''
        try:
            super().run(additional_parameters="-h")
        except RuntimeError:
            return False

        return True

    def process(self, input_template: Path,
                output_file: Path,
                key_values: Dict[str, str]) -> None:
        """
        This wrapper runs the Templaterator, which replaces the
        give keys in the input template with the value in the
        `key_values` dictionary. The new file is written to the
        specified output file.

        :param input_template: the path to the input template.
        :param output_file: the output file path.
        :param key_values: the keys and values for the keys to
            define as a dictionary.
        """
        replace_list = [f"-s {key}={value}"
                        for key, value in key_values.items()]
        params: List[Union[str, Path]]
        params = [input_template, "-o", output_file]
        params.extend(replace_list)

        super().run(additional_parameters=params)
