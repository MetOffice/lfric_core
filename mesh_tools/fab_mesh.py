#!/usr/bin/env python3
# ##############################################################################
#  (c) Crown copyright Met Office. All rights reserved.
#  For further details please refer to the file COPYRIGHT
#  which you should have received as part of this distribution
# ##############################################################################

'''A FAB build script for mesh_tools. It relies on the LFRicBase class
contained in the infrastructure directory.
'''

import argparse
import logging
from pathlib import Path
import sys
from typing import Optional

from fab.steps.grab.folder import grab_folder

# We need to import the base class:
sys.path.insert(0, str(Path(__file__).parents[1] / "lfric_build"))

from lfric_base import LFRicBase  # noqa: E402


class FabMeshTool(LFRicBase):

    def __init__(self,
                 name: str = "mesh_tools",
                 root_symbols: list[str] = []):

        app_dir = Path(__file__).parent
        super().__init__(name=name, app_dir=app_dir)
        self.set_root_symbols(root_symbols)

    def define_command_line_options(
            self,
            parser: Optional[argparse.ArgumentParser] = None
            ) -> argparse.ArgumentParser:
        '''
        Overwrite to change the default of the psyclone-control
        option to be the dummy one here, which returns will indicate
        to run PSyclone without a script.

        :param parser: optional a pre-defined argument parser. If not
            specified, a new instance will be created.
        '''
        parser = super().define_command_line_options(parser)

        control = self.app_dir / "psyclone_control.yaml"
        parser.set_defaults(psyclone_control=[str(control)])
        return parser

    def grab_files_step(self):
        super().grab_files_step()
        dirs = ['mesh_tools/source/']

        # pylint: disable=redefined-builtin
        for dir in dirs:
            grab_folder(self.config, src=self.lfric_core_root / dir,
                        dst_label='')

    def get_rose_meta(self):
        return (self.lfric_core_root / 'mesh_tools' / 'rose-meta' /
                'lfric-mesh_tools' / 'HEAD' / 'rose-meta.conf')


# -----------------------------------------------------------------------------
if __name__ == '__main__':

    logger = logging.getLogger('fab')
    logger.setLevel(logging.DEBUG)
    fab_mesh_tool = FabMeshTool(name="mesh_tools",
                                root_symbols=[
                                    'cubedsphere_mesh_generator',
                                    'planar_mesh_generator',
                                    'summarise_ugrid'])
    fab_mesh_tool.build()
