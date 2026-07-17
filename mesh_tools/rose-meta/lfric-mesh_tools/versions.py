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
    """Upgrade macro for ticket #386 by Christine Johnson."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t386"

    def upgrade(self, config, meta_config=None):

        # Renames meshes (mesh names and mesh maps) as
        # planar geometry: l0_planar, l1_planar, l2_planar and l3_planar
        # cubedsphere geometry: l0_cubedsphere, l1_cubedsphere, l2_cubedsphere, l3_cubedsphere
        # and depending on the number of meshes
        
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
                
        # Point to the correct parent mesh for LAM meshes
        
        self.change_setting_value(
            config, ["namelist:planar_mesh", "lbc_parent_mesh"], "'l0_planar'"
        )
        self.change_setting_value(
            config,
            ["namelist:stretch_transform", "transform_mesh"],
            "'l0_planar'",
        )

        return config, self.reports
