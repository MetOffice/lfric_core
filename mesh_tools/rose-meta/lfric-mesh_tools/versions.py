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


class vn31_t270(MacroUpgrade):
    """Upgrade macro for ticket #270 by cjohnson-pi."""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t270"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-mesh_tools
        nml = "namelist:planar_mesh"
        apply_stretch_transform = self.get_setting_value(
            config, [nml, "apply_stretch_transform"]
        )
        if apply_stretch_transform == ".true.":
            self.add_setting(
                config,
                [nml, "stretch_function"],
                "'inflation'",
            )
        else:
            self.add_setting(
                config,
                [nml, "stretch_function"],
                "'uniform'",
            )
        self.remove_setting(config, [nml, "apply_stretch_transform"])

        return config, self.reports
