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
                    "'l0_planar:l1_planar','l1_planar:l2_planar','l2_planar:l3_planar','l0_planar:l3_planar",
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'l0_planar','l1_planar','l2_planar','l3_planar'",
                    forced=True,
                )
            elif n_meshes == "3":
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_maps"],
                    "'l0_planar:l1_planar','l1_planar:l2_planar'",
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'l0_planar','l1_planar','l2_planar'",
                    forced=True,
                )
            elif n_meshes == "2":
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_maps"],
                    "'l0_planar:l1_planar'",
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'l0_planar','l1_planar'",
                    forced=True,
                )
            else:
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'l0_planar'",
                    forced=True,
                )
        else:
            if n_meshes == "4":
                edge_cells = self.get_setting_value(
                    config, ["namelist:cubedsphere_mesh", "edge_cells"]
                )
                if edge_cells == "48,24,12,6":
                    self.add_setting(
                        config,
                        ["namelist:mesh", "mesh_maps"],
                        "'l0_cubedsphere:l1_cubedsphere','l1_cubedsphere:l2_cubedsphere','l2_cubedsphere:l3_cubedsphere','l0_cubedsphere:l2_cubedsphere'",
                        forced=True,
                    )
                else:
                    self.add_setting(
                        config,
                        ["namelist:mesh", "mesh_maps"],
                        "'l0_cubedsphere:l1_cubedsphere','l1_cubedsphere:l2_cubedsphere','l2_cubedsphere:l3_cubedsphere'",
                        forced=True,
                    )
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'l0_cubedsphere','l1_cubedsphere','l2_cubedsphere','l3_cubedsphere'",
                    forced=True,
                )
            elif n_meshes == "3":
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_maps"],
                    "'l0_cubedsphere:l1_cubedsphere','l1_cubedsphere:l2_cubedsphere'",
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'l0_cubedsphere','l1_cubedsphere','l2_cubedsphere'",
                    forced=True,
                )
            elif n_meshes == "2":
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_maps"],
                    "'l0_cubedsphere:l1_cubedsphere'",
                    forced=True,
                )
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'l0_cubedsphere','l1_cubedsphere'",
                    forced=True,
                )
            else:
                self.add_setting(
                    config,
                    ["namelist:mesh", "mesh_names"],
                    "'l0_cubedsphere'",
                    forced=True,
                )
        self.change_setting_value(
            config, ["namelist:planar_mesh", "lbc_parent_mesh"], "'l0_planar'"
        )
        self.change_setting_value(
            config,
            ["namelist:stretch_transform", "transform_mesh"],
            "'l0_planar'",
        )

        return config, self.reports
