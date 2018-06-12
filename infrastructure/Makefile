##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Master make file for LFRic infrastructure project.
# Targets provided our detailed below...
#
# all: (default) Complete build and test the infrastructure.
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

PROJECT_NAME = infrastructure

PROFILE ?= fast-debug

export IGNORE_DEPENDENCIES = netcdf MPI yaxt pfunit_mod planet_constants_mod xios mod_wait
export EXTERNAL_DYNAMIC_LIBRARIES = netcdff netcdf hdf5 yaxt yaxt_c

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

.PHONY: default
default: unit-tests integration-tests
	$(Q)echo > /dev/null

.PHONY: documentation doc docs
documentation doc docs: document-uml document-api
	$(Q)echo > /dev/null

.PHONY: test-suite
test-suite: SUITE_BASE_NAME = $(notdir $(realpath $(shell pwd)/..))-$(PROJECT_NAME)
test-suite: SUITE_CONFIG = rose-stem
test-suite: launch-test-suite
	$(Q)echo > /dev/null

include build/lfric.mk

export FFLAGS  := $(FFLAGS) $(OPENMP_ARG) $(LFRIC_INTEL_FIX_ARG)
export LDFLAGS := $(LDFLAGS) $(OPENMP_ARG)

##############################################################################
# Extract infrastructure source
#
.PHONY: extract
extract: export SOURCE_DIR = source
extract: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk

##############################################################################
# Documentation
#
.PHONY: document-uml
document-uml: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/uml
document-uml: SOURCE_DIR    = documentation/uml
document-uml: WORKING_DIR  := $(WORKING_DIR)/uml
document-uml: uml-documentation
	$(Q)echo > /dev/null


.PHONY: document-api
document-api: PROJECT       = infrastructure
document-api: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/api
document-api: CONFIG_DIR    = documentation
document-api: SOURCE_DIR    = source
document-api: WORKING_DIR  := $(WORKING_DIR)/api
document-api: api-documentation
	$(Q)echo > /dev/null

##############################################################################
# Unit tests
#
.PHONY: unit-tests
unit-tests: export PROJECT      = infrastructure
unit-tests: export SOURCE_DIR   = unit-test
unit-tests: export WORKING_DIR := $(WORKING_DIR)/unit-tests
unit-tests: export BIN_DIR     ?= $(PROJECT_DIR)/test
unit-tests: export PROGRAMS     = infrastructure_unit_tests
unit-tests: export EXTERNAL_STATIC_LIBRARIES += pfunit
unit-tests: export FFLAGS      += $(FFLAGS_DEBUG) $(FFLAGS_NO_OPTIMISATION) $(FFLAGS_INIT)
unit-tests: run-unit-tests
	$(Q)echo > /dev/null

.PHONY: generate-unit-tests
generate-unit-tests: extract-unit-tests
	$(Q)$(MAKE) $(QUIET_ARG) extract
	$(Q)$(MAKE) $(QUIET_ARG) pfunit

.PHONY: extract-unit-tests
extract-unit-tests: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk
	$(Q)cp -r $(SOURCE_DIR)/data $(WORKING_DIR)/

##############################################################################
# Integration tests
#
.PHONY: integration-tests
integration-tests: export PROJECT      = infrastructure
integration-tests: export WORKING_DIR := $(WORKING_DIR)/integration-tests
integration-tests: export BIN_DIR     ?= $(PROJECT_DIR)/test
integration-tests: export SOURCE_DIR   = integration-test
integration-tests: $(addprefix integration-test-run/,$(basename $(shell find integration-test -name '*.[Ff]90' -exec egrep -l "^\s*program" {} \;)))
	$(Q)echo > /dev/null

.PHONY: generate-integration-tests
generate-integration-tests: extract-integration-test \
	                   configuration-integration-test
	$(Q)$(MAKE) $(QUIET_ARG) extract

.PHONY: configuration-integration-test
configuration-integration-test: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/lfric.mk configuration

.PHONY: extract-integration-test
extract-integration-test: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk

##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"infrastructure work space")
	$(Q)-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"infrastructure documents")
	$(Q)if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"infrastructure binaries")
	$(Q)-if [ -d test ] ; then rm -r test ; fi
