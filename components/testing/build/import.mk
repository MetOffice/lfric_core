##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Build configuration for importing the testing component
#
# This makefile defines the build process for the testing component,
# including source extraction, PSyclone code generation, and installation
# of support files.

export PROJECT_SOURCE = $(CORE_ROOT_DIR)/components/testing/source
PROJECT_SUPPORT_DIR = $(CORE_ROOT_DIR)/components/testing/support

.PHONY: import-testing
import-testing:
	# Extract Fortran source files from project directories
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk SOURCE_DIR=$(PROJECT_SOURCE)
	# Generate PSyclone kernel code and PSy-layer code
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
            SOURCE_DIR=$(PROJECT_SOURCE) \
            OPTIMISATION_PATH=$(OPTIMISATION_PATH)
	# Copy support files and utilities to binary directory
	$Qmkdir -p $(BIN_DIR)
	$Qrsync -a $(PROJECT_SUPPORT_DIR)/ $(BIN_DIR)