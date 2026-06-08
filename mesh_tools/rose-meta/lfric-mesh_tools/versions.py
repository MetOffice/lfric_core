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


class vn31_t11(MacroUpgrade):
    """Upgrade macro for ticket TTTT by Unknown."""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t11"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-mesh_tools
        geometry = self.get_setting_value(config, ["namelist:mesh", "geometry"])
        n_meshes = self.get_setting_value(config, ["namelist:mesh", "n_meshes"])
        if geometry == "planar":
            if n_meshes == "4":
                self.add_setting(
                    config,
                    ["mesh_maps"],
                    [
                        "'regional_primary':'regional_coarse_l1',
                         'regional_coarse_l1':'regional_coarse_l2',
                         'regional_coarse_l2':'regional_coarse_l3'",
                    ],
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["mesh_names"],
                    [
                        "'regional_primary',
                         'regional_coarse_l1',
                         'regional_coarse_l2',
                         'regional_coarse_l3'",
                    ],
                    forced=True,
                )
            elif n_meshes == "3":
                self.add_setting(
                    config,
                    ["mesh_maps"],
                    [
                        "'regional_primary':'regional_coarse_l1',
                         'regional_coarse_l1':'regional_coarse_l2'",
                    ],
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["mesh_names"],
                    [
                        "'regional_primary',
                         'regional_coarse_l1',
                         'regional_coarse_l2'",
                    ],
                    forced=True,
                )
            elif n_meshes == "2":
                self.add_setting(
                    config,
                    ["mesh_maps"],
                    ["'regional_primary':'regional_coarse_l1'"],
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["mesh_names"],
                    ["'regional_primary','regional_coarse_l1'"],
                    forced=True,
                )
            else:
                self.add_setting(
                    config, ["mesh_names"], ["'regional_primary'"], forced=True
                )
        else:
            if n_meshes == "4":
                self.add_setting(
                    config,
                    ["mesh_maps"],
                    [
                        "'global_primary:global_coarse_l1',
                         'global_coarse_l1:global_coarse_l2',
                         'global_coarse_l2:global_coarse_l3'",
                    ],
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["mesh_names"],
                    [
                        "'global_primary',
                         'global_coarse_l1',
                         'global_coarse_l2',
                         'global_coarse_l3'",
                    ],
                    forced=True,
                )
            elif n_meshes == "3":
                self.add_setting(
                    config,
                    ["mesh_maps"],
                    [
                        "'global_primary:global_coarse_l1',
                         'global_coarse_l1:global_coarse_l2'",
                    ],
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["mesh_names"],
                    ["'global_primary,global_coarse_l1,global_coarse_l2'"],
                    forced=True,
                )
            elif n_meshes == "2":
                self.add_setting(
                    config,
                    ["mesh_maps"],
                    ["'global_primary:global_coarse_l1'"],
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["mesh_names"],
                    ["'global_primary','global_coarse_l1'"],
                    forced=True,
                )
            else:
                self.add_setting(
                    config, ["mesh_names"], ["'global_primary'"], forced=True
                )

        return config, self.reports
