##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

PROJECT_NAME = gungho

PROFILE ?= fast-debug

export IGNORE_DEPENDENCIES = netcdf MPI yaxt pfunit_mod xios mod_wait
export EXTERNAL_DYNAMIC_LIBRARIES = yaxt yaxt_c netcdff netcdf hdf5 \
                                    $(CXX_RUNTIME_LIBRARY)
export EXTERNAL_STATIC_LIBRARIES = xios 


LFRIC_TARGET_PLATFORM ?= generic-x86
export OPTIMISATION_PATH ?= $(wildcard optimisation/$(LFRIC_TARGET_PLATFORM))

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
ROOT_DIR ?= ../

.PHONY: default
default: build unit-tests integration-tests
	$(Q)echo > /dev/null

.PHONY: documentation doc docs
documentation doc docs: document-uml document-latex document-api
	$(Q)echo > /dev/null

include $(ROOT_DIR)/infrastructure/build/lfric.mk

export FFLAGS  := $(FFLAGS) $(OPENMP_ARG) $(LFRIC_INTEL_FIX_ARG)
export LDFLAGS := $(LDFLAGS) $(OPENMP_ARG)

##############################################################################

.PHONY: test-suite
test-suite: SUITE_BASE_NAME = $(notdir $(realpath $(shell pwd)/..))-$(PROJECT_NAME)
test-suite: SUITE_CONFIG = rose-stem
test-suite: launch-test-suite
	$(Q)echo > /dev/null

##############################################################################
# Documentation
#
.PHONY: document-uml
document-uml: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/uml
document-uml: SOURCE_DIR    = documentation/uml
document-uml: WORKING_DIR  := $(WORKING_DIR)/uml
document-uml: uml-documentation
	$(Q)echo > /dev/null

.PHONY: document-latex
document-latex: DOCUMENT_DIR   = $(PROJECT_DIR)/documents
document-latex: SOURCE_DIR     = documentation
document-latex: DOCUMENTS      = $(shell find $(SOURCE_DIR) -name '*.latex' -print)
document-latex: WORKING_DIR   := $(WORKING_DIR)/latex
document-latex: TEX_STUFF      = documentation/tex
document-latex: COMMON_FIGURES = documentation/common-figures
document-latex: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/latex.mk \
                                 SOURCE_DIR=$(SOURCE_DIR) \
	                         WORKING_DIR=$(WORKING_DIR) \
	                         DOCUMENT_DIR=$(DOCUMENT_DIR) \
	                         TEX_STUFF=$(TEX_STUFF) \
	                         COMMON_FIGURES=$(COMMON_FIGURES) \
	                         DOCUMENTS="$(DOCUMENTS)"

.PHONY: document-api
document-api: PROJECT       = gungho
document-api: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/api
document-api: CONFIG_DIR    = documentation
document-api: SOURCE_DIR    = source
document-api: WORKING_DIR  := $(WORKING_DIR)/api
document-api: api-documentation
	$(Q)echo > /dev/null

##############################################################################
# Extract GungHo source
#
.PHONY: extract
extract: SOURCE_DIR = source
extract: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk SOURCE_DIR=$(SOURCE_DIR) \
	                                     WORKING_DIR=$(WORKING_DIR)

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
build: PROJECT      = gungho
build: SOURCE_DIR   = source
build: WORKING_DIR := $(WORKING_DIR)/gungho
build: BIN_DIR     ?= $(PROJECT_DIR)/bin
build: PROGRAMS    := $(basename $(notdir $(shell find $(SOURCE_DIR) -maxdepth 1 -name '*.[Ff]90' -print)))
build: generate-configuration generate-psykal
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

.PHONY: generate-psykal
generate-psykal: export SOURCE_DIR = source
generate-psykal: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) extract
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk WORKING_DIR=$(WORKING_DIR)

##############################################################################
# Unit tests
#
.PHONY: unit-tests
unit-tests: export PROJECT      = gungho
unit-tests: export SOURCE_DIR   = unit-test
unit-tests: export WORKING_DIR := $(WORKING_DIR)/unit-tests
unit-tests: export BIN_DIR     ?= $(PROJECT_DIR)/test
unit-tests: export PROGRAMS     = gungho_unit_tests
unit-tests: export EXTERNAL_STATIC_LIBRARIES += pfunit
unit-tests: export FFLAGS      += $(FFLAGS_DEBUG) $(FFLAGS_NO_OPTIMISATION) $(FFLAGS_INIT)
unit-tests:run-unit-tests
	$(Q)echo > /dev/null

.PHONY: generate-unit-tests
generate-unit-tests: extract-unit-tests
	$(Q)$(MAKE) $(QUIET_ARG) -C $(LFRIC_INFRASTRUCTURE) extract WORKING_DIR=$(realpath $(WORKING_DIR))
	$(Q)$(MAKE) $(QUIET_ARG) extract
	$(Q)$(MAKE) $(QUIET_ARG) generate-configuration
	$(Q)$(MAKE) $(QUIET_ARG) generate-psykal
	$(Q)$(MAKE) $(QUIET_ARG) pfunit

.PHONY: extract-unit-tests
extract-unit-tests:
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk
	$(Q)cp -r $(SOURCE_DIR)/data $(WORKING_DIR)

##############################################################################
# Integration tests
#
.PHONY: integration-tests
integration-tests: export PROJECT      = gungho
integration-tests: export SOURCE_DIR   = integration-test
integration-tests: export WORKING_DIR := $(WORKING_DIR)/integration-tests
integration-tests: export BIN_DIR     ?= $(PROJECT_DIR)/test
integration-tests: $(addprefix integration-test-run/,$(basename $(shell find integration-test -name '*.[Ff]90' -exec egrep -l "^\s*program" {} \;)))
	$(Q)echo > /dev/null

.PHONY: generate-integration-tests
generate-integration-tests: extract-integration-tests
	$(Q)$(MAKE) $(QUIET_ARG) -C $(LFRIC_INFRASTRUCTURE) extract WORKING_DIR=$(realpath $(WORKING_DIR))
	$(Q)$(MAKE) $(QUIET_ARG) extract
	$(Q)$(MAKE) $(QUIET_ARG) generate-configuration
	$(Q)$(MAKE) $(QUIET_ARG) generate-psykal

.PHONY: extract-integration-tests
extract-integration-tests: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk

##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"gungho work space")
	$(Q)-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"gungho documents")
	$(Q)if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"gungho binaries")
	$(Q)-rm -r bin test
