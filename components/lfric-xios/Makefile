##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Make file for LFRic-XIOS model component
# Targets provided our detailed below...
#
# all: (default) Complete build and test the LFRic-XIOS model component.

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

PROJECT_NAME = lfric-xios

PROFILE ?= fast-debug

export IGNORE_DEPENDENCIES = netcdf MPI yaxt pfunit_mod xios
export EXTERNAL_DYNAMIC_LIBRARIES = netcdff netcdf yaxt yaxt_c

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
export ROOT_DIR    ?= ../..

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/$(PROJECT_NAME)/$(META_VN)

.PHONY: default
default: unit-tests
	$(Q)echo > /dev/null

.PHONY: documentation doc docs
documentation doc docs: document-uml document-api
	$(Q)echo > /dev/null

include $(ROOT_DIR)/infrastructure/build/lfric.mk

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
document-api: PROJECT       = lfric-xios
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
unit-tests: export ADDITIONAL_EXTRACTION = "$(ROOT_DIR)/infrastructure/source"
unit-tests: export BIN_DIR ?= $(PROJECT_DIR)/test
unit-tests: export CXX_LINK = TRUE
unit-tests: export EXTERNAL_STATIC_LIBRARIES += pfunit
unit-tests: export FFLAG_GROUPS = DEBUG NO_OPTIMISATION INIT UNIT_WARNINGS
unit-tests: export FFLAGS += -D UNIT_TEST
unit-tests: export META_FILE_DIR = rose-meta/lfric-xios/HEAD
unit-tests: export PROGRAMS = lfric_xios_unit_tests
unit-tests: export PROJECT = lfric_xios
unit-tests: export SOURCE_DIR = source
unit-tests: export TEST_DIR = unit-test
unit-tests: export WORKING_DIR := $(WORKING_DIR)/unit-tests
unit-tests: do-unit-tests
	$(Q)echo > /dev/null

##############################################################################
# Integration tests
#
#.PHONY: integration-tests
#integration-tests: export BIN_DIR ?= $(PROJECT_DIR)/test
#integration-tests: export META_FILE_DIR = integration-test/configuration_test
#integration-tests: export PROJECT = lfric-xios
#integration-tests: export SOURCE_DIR = source
#integration-tests: export TEST_DIR = integration-test
#integration-tests: export WORKING_DIR := $(WORKING_DIR)/integration-tests
#integration-tests: do-integration-tests
#	$(Q)echo > /dev/null

##############################################################################
# Test suite rules
#
.PHONY: launch-suite-gscan
launch-suite-gscan: gscan_processes := $(shell ps --no-headers -o command -C cylc-gscan)
launch-suite-gscan: ALWAYS
	$(Q)-if [[ "$(gscan_processes)" != *"--name=$(SUITE_GROUP_NAME)"* ]]; then \
          cylc gscan  $(DOUBLE_VERBOSE_ARG) --name=$(SUITE_GROUP_NAME) &           \
          usleep 1                                                                ;\
        fi

.PHONY: test-suite
test-suite: SUITE_BASE_NAME = $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)
test-suite: SUITE_CONFIG = rose-stem
test-suite: launch-suite-gscan launch-test-suite
	$(Q)echo > /dev/null

print-%  : ; @echo $* = $($*)

##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"lfric-xios work space")
	$(Q)-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"lfric-xios documents")
	$(Q)if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"lfric-xios binaries")
	$(Q)-if [ -d test ] ; then rm -r test ; fi
