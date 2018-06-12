##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Master make file for LFRic mesh tools project.
# Targets provided our detailed below...
#
# all: (default) Complete build and test the mesh tools.
#      Intended for interactive use.
# build: Build all executables and any supporting files
# test: Run test battery including unit tests and others
# clean: Delete all final products and working files
#
# The following variables may be specified to modify the build process...
#
# WORKING_DIR: Path to scratch space in which intermediate files will be
#              placed. This should be somewhere with good "many small
#              files" performance, i.e. probably not Lustre.
#              Default: ./working
# VERBOSE: Set in order to see actual commands issued by the build system.
# PROFILE: Set to a string representing a package of compiler options.
#          Potential profiles are 'full-debug', 'fast-debug' and 'production'.
#          Default: fast-debug
# LINK_TYPE: Either 'static' or 'dynamic'.
#            Default: dynamic
#
##############################################################################

PROJECT_NAME = mesh_tools

PROFILE ?= fast-debug

export IGNORE_DEPENDENCIES = netcdf MPI yaxt pfunit_mod planet_constants_mod xios mod_wait
export EXTERNAL_DYNAMIC_LIBRARIES = netcdff netcdf hdf5 yaxt yaxt_c

PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
ROOT_DIR    ?= ../

.PHONY: default
default: build unit-tests
	$(Q)echo > /dev/null

include $(ROOT_DIR)/infrastructure/build/lfric.mk

export FFLAGS += $(OPENMP_ARG) $(LFRIC_INTEL_FIX_ARG)
export LDFLAGS += $(OPENMP_ARG)

##############################################################################

.PHONY: test-suite
test-suite: SUITE_BASE_NAME = $(notdir $(realpath $(shell pwd)/..))-$(PROJECT_NAME)
test-suite: SUITE_CONFIG = rose-stem
test-suite: launch-test-suite
	$(Q)echo > /dev/null

##############################################################################
# Build
#
.PHONY: build
ifeq "$(PROFILE)" "full-debug"
build: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) $(FFLAGS_INIT) \
        $(FFLAGS_RUNTIME) $(FFLAGS_NO_OPTIMISATION) $(STANDARDS_CHECK_ARG)
else ifeq "$(PROFILE)" "fast-debug"
build: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) \
        $(FFLAGS_SAFE_OPTIMISATION) $(STANDARDS_CHECK_ARG)
else ifeq "$(PROFILE)" "production"
build: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) $(FFLAGS_RISKY_OPTIMISATION)
else
  $(error Unrecognised profile "$(PROFILE)". Must be one of full-debug, fast-debug or production)
endif
build: PROJECT      = mesh_tools
build: SOURCE_DIR   = source
build: WORKING_DIR := $(WORKING_DIR)/mesh_tools
build: BIN_DIR     ?= $(PROJECT_DIR)/bin
build: PROGRAMS    := $(basename $(notdir $(shell find $(SOURCE_DIR) -maxdepth 1 -name '*.[Ff]90' -print)))
build: generate-configuration extract
	$(Q)$(MAKE) $(QUIET_ARG) -C $(LFRIC_INFRASTRUCTURE) extract WORKING_DIR=$(realpath $(WORKING_DIR))
	$(Q)$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$(Q)$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk \
	                             BIN_DIR=$(BIN_DIR) PROGRAMS="$(PROGRAMS)"

.PHONY: generate-configuration
generate-configuration: SOURCE_DIR = source
generate-configuration: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk PROJECT=$(PROJECT) \
	                                           SOURCE_DIR=$(SOURCE_DIR) \
	                                           WORKING_DIR=$(WORKING_DIR)

.PHONY: extract
extract: export SOURCE_DIR = source
extract: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk

##############################################################################
# Unit tests
#
.PHONY: unit-tests
unit-tests: export PROJECT      = mesh_tools
unit-tests: export SOURCE_DIR   = unit-test
unit-tests: export WORKING_DIR := $(WORKING_DIR)
unit-tests: export BIN_DIR     ?= $(PROJECT_DIR)/test
unit-tests: export PROGRAMS     = mesh_tools_unit_tests
unit-tests: export EXTERNAL_STATIC_LIBRARIES += pfunit
unit-tests: export FFLAGS      += $(FFLAGS_DEBUG) $(FFLAGS_NO_OPTIMISATION) $(FFLAGS_INIT)
unit-tests: run-unit-tests
	$(Q)echo > /dev/null

.PHONY: generate-unit-tests
generate-unit-tests: extract-unit-tests
	$(Q)$(MAKE) $(QUIET_ARG) -C $(LFRIC_INFRASTRUCTURE) extract WORKING_DIR=$(realpath $(WORKING_DIR))
	$(Q)$(MAKE) $(QUIET_ARG) extract
	$(Q)$(MAKE) $(QUIET_ARG) generate-configuration
	$(Q)$(MAKE) $(QUIET_ARG) pfunit

.PHONY: extract-unit-tests
extract-unit-tests:
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk

##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"mesh tools work space")
	$(Q)-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"mesh tools binaries")
	$(Q)-if [ -d bin ] ; then rm -r bin ; fi
	$(call MESSAGE,Removing,"mesh tools tests")
	$(Q)-if [ -d test ] ; then rm -r test ; fi
