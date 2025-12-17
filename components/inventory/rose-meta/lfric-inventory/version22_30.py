import sys

from metomi.rose.upgrade import MacroUpgrade

from .version21_22 import *


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


class vn22_t4661(MacroUpgrade):
    """Upgrade macro for ticket #4661 by Denis Sergeev."""

    BEFORE_TAG = "vn2.2"
    AFTER_TAG = "vn2.2_t4661"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-driver
        self.add_setting(config, ["namelist:extrusion", "eta_values"], "''")
        return config, self.reports


class vn22_t4020(MacroUpgrade):
    """Upgrade macro for ticket #4020 by Andrew Coughtrie."""

    BEFORE_TAG = "vn2.2_t4661"
    AFTER_TAG = "vn2.2_t4020"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-driver
        self.add_setting(
            config, ["namelist:io", "end_of_run_checkpoint"], ".true."
        )
        self.add_setting(config, ["namelist:io", "checkpoint_times"], "")
        return config, self.reports


class vn22_t34(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn2.2_t4020"
    AFTER_TAG = "vn3.0"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-driver
        # Blank Upgrade Macro
        # Commands From: rose-meta/lfric-inventory
        # Blank Upgrade Macro
        return config, self.reports
