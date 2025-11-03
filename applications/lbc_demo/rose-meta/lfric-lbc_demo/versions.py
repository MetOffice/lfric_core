import sys
from metomi.rose.upgrade import MacroUpgrade

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

class vn22_t4231(MacroUpgrade):
    # Upgrade macro for #4231 by Ricky Wong
    # Ticket #4231 adds a new application, "lbc_demo"
    # This ugrade macro does nothing but set the
    # appropriate "After_tag"

    BEFORE_TAG = "vn2.2"
    AFTER_TAG = "vn2.2_t4231"

    def upgrade(self, config, meta_config=None):
        # Add settings
        return config, self.reports
"""
