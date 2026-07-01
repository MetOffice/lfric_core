import re
import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version30_31 import *


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


class vn31_t238(MacroUpgrade):
    """Upgrade macro for ticket #238 by Thomas Bendall."""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t238"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-driver
        self.add_setting(
            config, ["namelist:finite_element", "coord_space"], "'Wchi'"
        )
        coord_order = self.get_setting_value(
            config, ["namelist:finite_element", "coord_order"]
        )
        self.add_setting(
            config,
            ["namelist:finite_element", "coord_order_nonprime"],
            coord_order,
        )
        return config, self.reports


class vn31_t324(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn3.1_t238"
    AFTER_TAG = "vn3.1_t324"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-driver
        # Only add in new configuration settings if the namelists
        # are already present
        #
        if config.get(["namelist:partitioning"]) is not None:
            self.add_setting(
                config, ["namelist:partitioning", "inner_halo_tiles"], ".false."
            )
            self.add_setting(
                config, ["namelist:partitioning", "tile_size_x"], "1"
            )
            self.add_setting(
                config, ["namelist:partitioning", "tile_size_y"], "1"
            )
        if config.get(["namelist:multigrid"]) is not None:
            self.add_setting(
                config,
                ["namelist:multigrid", "coarsen_multigrid_tiles"],
                ".false.",
            )
            self.add_setting(
                config, ["namelist:multigrid", "max_tiled_multigrid_level"], "1"
            )

        return config, self.reports
