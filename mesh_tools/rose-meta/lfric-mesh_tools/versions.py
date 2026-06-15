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
        n_meshes = self.get_setting_value(config, ["namelist:mesh", "n_meshes"])
        if config.get_value(["namelist:planar_mesh"]) is not None:
            if n_meshes == "4":
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_maps"],
                    "'planar_l0:planar_l1','planar_l1:planar_l2','planar_l2:planar_l3'",
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'planar_l0','planar_l1','planar_l2','planar_l3'",
                    forced=True,
                )
            elif n_meshes == "3":
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_maps"],
                    "'planar_l0:planar_l1','planar_l1:planar_l2'",
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'planar_l0','planar_l1','planar_l2'",
                    forced=True,
                )
            elif n_meshes == "2":
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_maps"],
                    "'planar_l0:planar_l1'",
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'planar_l0','planar_l1'",
                    forced=True,
                )
            else:
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'planar_l0'",
                    forced=True,
                )
        else:
            if n_meshes == "4":
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_maps"],
                    "'cubedsphere_l0:cubedsphere_l1','cubedsphere_l1:cubedsphere_l2','cubedsphere_l2:cubedsphere_l3'",
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'cubedsphere_l0','cubedsphere_l1','cubedsphere_l2','cubedsphere_l3'",
                    forced=True,
                )
            elif n_meshes == "3":
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_maps"],
                    "'cubedsphere_l0:cubedsphere_l1','cubedsphere_l1:cubedsphere_l2'",
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'cubedsphere_l0','cubedsphere_l1','cubedsphere_l2'",
                    forced=True,
                )
            elif n_meshes == "2":
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_maps"],
                    "'cubedsphere_l0:cubedsphere_l1'",
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'cubedsphere_l0','cubedsphere_l1'",
                    forced=True,
                )
            else:
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'cubedsphere_l0'",
                    forced=True,
                )

        topology = self.get_setting_value(
            config, ["namelist:base_mesh", "topology"] )
        if topology == "'non_periodic'":
            self.update_setting(
                config, ["namelist:base_mesh", "prime_mesh_name" ],
                "'planar_l0'")

        return config, self.reports
