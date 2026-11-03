##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
export PROJECT_SOURCE = $(CORE_ROOT_DIR)/components/coupling/source

export IGNORE_DEPENDENCIES += mod_oasis xios
# For static linking the xios library must be linked before the Oasis libaries.
# To guarantee this, "xios" is listed here, even though it is not used by
# coupling
export EXTERNAL_STATIC_LIBRARIES += xios psmile.MPI1 mct mpeu scrip

.PHONY: import-coupling
import-coupling:
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk SOURCE_DIR=$(PROJECT_SOURCE)
