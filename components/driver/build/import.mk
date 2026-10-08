##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
export PROJECT_SOURCE = $(CORE_ROOT_DIR)/components/driver/source

.PHONY: import-driver
import-driver:
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk SOURCE_DIR=$(PROJECT_SOURCE)
	# Generate PSyclone kernel code and PSy-layer code
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
            SOURCE_DIR=$(PROJECT_SOURCE) \
            OPTIMISATION_PATH=$(OPTIMISATION_PATH)
