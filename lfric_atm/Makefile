##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

PROJECT_NAME = lfric_atm

export PROFILE ?= fast-debug
export LFRIC_TARGET_PLATFORM ?= generic-x86

# This top level makefile is very order sensitive. Source code extraction and
# generation must happen in a certain order. Due to this we turn off
# multithreading for this file only. Any called recursively (i.e. with $(MAKE))
# run in parallel. Unless they specify NOTPARALLEL as well.
#
.NOTPARALLEL:

export IGNORE_DEPENDENCIES = netcdf MPI yaxt pfunit_mod xios mod_wait
export EXTERNAL_DYNAMIC_LIBRARIES = yaxt yaxt_c netcdff netcdf hdf5
export EXTERNAL_STATIC_LIBRARIES = xios

export UMPHYSICS_TARGETS = extract-um-physics
export UMPHYSICS_CLEAN   = clean-um-physics
export UM_FCM=$(UM_PHYS_DIR)/
export OPTIMISATION_PATH ?= $(wildcard optimisation/$(LFRIC_TARGET_PLATFORM))

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
export ROOT_DIR    ?= ..
export SUITE_GROUP ?= developer
export SUITE_GROUP_NAME ?= $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)-.*

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

.PHONY: default
default: build
	$(Q)echo > /dev/null

include $(ROOT_DIR)/infrastructure/build/lfric.mk
include build/project.mk

##############################################################################
# Launch gscan to monitor suites from this suite-group
#
.PHONY: launch-suite-gscan
launch-suite-gscan: gscan_processes := $(shell ps --no-headers -o command -C cylc-gscan)
launch-suite-gscan: ALWAYS
	$(Q)-if [[ "$(gscan_processes)" != *"--name=$(SUITE_GROUP_NAME)"* ]]; then \
          cylc gscan  $(DOUBLE_VERBOSE_ARG) --name=$(SUITE_GROUP_NAME) &           \
          usleep 1                                                                ;\
        fi

##############################################################################
# Test
#
.PHONY: test-suite
test-suite: SUITE_BASE_NAME = $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)
test-suite: SUITE_CONFIG = rose-stem
test-suite: launch-suite-gscan launch-test-suite
	$(Q)echo > /dev/null

##############################################################################
# Build
#
.PHONY: build
build: export BIN_DIR     ?= $(PROJECT_DIR)/bin
build: export COMPILE_OPTIONS = $(PWD)/build/compiler_options.mk
build: export CXX_LINK     = YES
build: export PROGRAMS    := $(basename $(notdir $(shell find source -maxdepth 1 -name '*.[Ff]90' -exec egrep -l "^\s*program" {} \;)))
build: export PROJECT      = lfric_atm
build: export SCRATCH_DIR := $(WORKING_DIR)/um-scratch
build: export WORKING_DIR := $(WORKING_DIR)/lfric_atm
ifeq "$(PROFILE)" "full-debug"
build: export FFLAG_GROUPS = DEBUG WARNINGS INIT RUNTIME NO_OPTIMISATION
else ifeq "$(PROFILE)" "fast-debug"
build: export FFLAG_GROUPS = DEBUG WARNINGS SAFE_OPTIMISATION
else ifeq "$(PROFILE)" "production"
build: export FFLAG_GROUPS = DEBUG WARNINGS RISKY_OPTIMISATION
else
  $(error Unrecognised profile "$(PROFILE)". Must be one of full-debug, fast-debug or production)
endif
build: ALWAYS
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
                  SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	          SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting Gungho dynamical core)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	         SOURCE_DIR=$(ROOT_DIR)/gungho/source \
	         WITHOUT_PROGRAMS=1
	$(call MESSAGE,====================================)
	$(call MESSAGE,Psycloning Gungho dynamical core)
	$(call MESSAGE,====================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	        SOURCE_DIR=$(ROOT_DIR)/gungho/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
                  SOURCE_DIR=$(ROOT_DIR)/gungho/source \
                  META_FILE_DIR=$(META_FILE_DIR) \
                  FIX_ENUMS_OPT=$(FIX_ENUMS_OPT)
	$(call MESSAGE,Extracting,"UM physics" )
	$Q$(MAKE) $(QUIET_ARG) -f $(PROJECT_DIR)/build/extract_physics.mk
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(ROOT_DIR)/um_physics/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	          SOURCE_DIR=$(ROOT_DIR)/um_physics/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(ROOT_DIR)/socrates/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	          SOURCE_DIR=$(ROOT_DIR)/socrates/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(ROOT_DIR)/jules/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	          SOURCE_DIR=$(ROOT_DIR)/jules/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
                  SOURCE_DIR=source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
                  SOURCE_DIR=source
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk

##############################################################################
# Unit tests
#
.PHONY: unit-tests
unit-tests: export ADDITIONAL_EXTRACTION = $(ROOT_DIR)/infrastructure/source
unit-tests: export BIN_DIR ?= $(PROJECT_DIR)/test
unit-tests: export EXTERNAL_STATIC_LIBRARIES += pfunit
unit-tests: export FFLAG_GROUPS += DEBUG NO_OPTIMISATION INIT
unit-tests: export PROGRAMS = lfric_atm_unit_tests
unit-tests: export PROJECT = lfric_atm 
unit-tests: export SOURCE_DIR = source
unit-tests: export TEST_DIR = unit-test
unit-tests: export WORKING_DIR := $(WORKING_DIR)/unit-tests
unit-tests: do-unit-tests
	$(Q)echo > /dev/null


##############################################################################
# Integration tests
#
.PHONY: integration-tests
integration-tests: export ADDITIONAL_EXTRACTION = $(ROOT_DIR)/infrastructure/source
integration-tests: export BIN_DIR ?= $(PROJECT_DIR)/test
integration-tests: export PROJECT =lfric_atm 
integration-tests: export SOURCE_DIR = source
integration-tests: export TEST_DIR = unit-test
integration-tests: export WORKING_DIR := $(WORKING_DIR)/integration-tests
integration-tests: do-integration-tests
	$(Q)echo > /dev/null


##############################################################################
# Clean
#
.PHONY: clean
clean: partial-clean
	$(call MESSAGE,Removing,"UM work space")
	$(Q)-rm -r $(WORKING_DIR)

.PHONY: partial-clean pclean
partial-clean pclean: ALWAYS
	$(call MESSAGE,Removing,"LFRic atmosphere  workspace")
	$(Q)-rm -r $(WORKING_DIR)/lfric_atm
	$(call MESSAGE,Removing,"LFRic atmosphere binaries")
	$(Q)-if [ -d bin ] ; then rm -r bin ; fi

