import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version31_32 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


"""
Copy this template and complete to add your macro

class vnXX_txxx(MacroUpgrade):
    # Upgrade macro for <TICKET> by <Author>

    BEFORE_TAG = "vnX.X"
    AFTER_TAG = "vnX.X_txxx"

    def upgrade(self, config, meta_config=None):
        # Add settings
        return config, self.reports
"""

class vn32_t386(MacroUpgrade):
    # Upgrade macro for #386 by Christine Johnson

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t386"

    def upgrade(self, config, meta_config=None):
        # Add settings
        topology = self.get_setting_value(config, ["namelist:base_mesh", "topology"])
        partitioner = self.get_setting_value(
            config, ["namelist:partitioning", "partitioner"]
        )
        if topology == "'non_periodic'" or partitioner == "'planar'":
            self.change_setting_value(
                 config,
                ["namelist:base_mesh", "prime_mesh_name"],
                "'l0_planar'",
            )
        else:
           self.change_setting_value(
                 config,
                ["namelist:base_mesh", "prime_mesh_name"],
                "'l0_cubedsphere'",
            ) 
        return config, self.reports
