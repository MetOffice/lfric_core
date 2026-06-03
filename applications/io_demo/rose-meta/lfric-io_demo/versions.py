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


class vn31_t232(MacroUpgrade):
    # Upgrade macro for 232 by Ed Hone

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t232"

    def upgrade(self, config, meta_config=None):
        """Add new io_demo namelist"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"namelist:io",
            r"namelist:io" + "\n" + " namelist:io_demo",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        self.add_setting(config, ["namelist:io_demo"])

        """Move multifile_io setting from io namelist to io_demo"""
        self.remove_setting(config, ["namelist:io", "multifile_io"])
        self.add_setting(
            config, ["namelist:io_demo", "multifile_io"], ".false."
        )
        self.add_setting(
            config, ["namelist:io_demo", "benchmark_sleep_time"], 0
        )
        self.add_setting(
            config, ["namelist:io_demo", "io_benchmark"], ".false."
        )
        self.add_setting(config, ["namelist:io_demo", "n_benchmark_fields"], 0)
        return config, self.reports


class vn31_t330(MacroUpgrade):
    # Upgrade macro for 330 by Ed Hone

    BEFORE_TAG = "vn3.1_t232"
    AFTER_TAG = "vn3.1_t330"

    def upgrade(self, config, meta_config=None):
        """Add new io_demo namelist"""
        source = self.get_setting_value(
            config, ["file:configuration.nml", "source"]
        )
        source = re.sub(
            r"namelist:extrusion",
            r"namelist:extrusion" + "\n" + " namelist:files",
            source,
        )
        self.change_setting_value(
            config, ["file:configuration.nml", "source"], source
        )
        self.add_setting(config, ["namelist:files"])
        self.add_setting(config, ["namelist:files", "temporal_file_path"], "")

        self.add_setting(
            config, ["namelist:io_demo", "temporal_reading"], ".false."
        )
        return config, self.reports


class vn31_t238(MacroUpgrade):
    """Upgrade macro for ticket #238 by Thomas Bendall."""

    BEFORE_TAG = "vn3.1_t330"
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
