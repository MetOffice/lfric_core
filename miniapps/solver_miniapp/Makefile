##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# This top level makefile is very order sensitive. Source code extraction and
# generation must happen in a certain order. Due to this we turn off
# multithreading for this file only. Any called recursively (i.e. with $(MAKE))
# run in parallel. Unless they specify NOTPARALLEL as well.
#
.NOTPARALLEL:

export PROJECT_NAME = solver_miniapp
export PROJECT_DIR_NAME = solver_miniapp
PROFILE ?= fast-debug

export IGNORE_DEPENDENCIES = netcdf MPI yaxt pfunit_mod xios mod_wait
export EXTERNAL_DYNAMIC_LIBRARIES = yaxt yaxt_c netcdff netcdf hdf5
export EXTERNAL_STATIC_LIBRARIES = xios

LFRIC_TARGET_PLATFORM ?= generic-x86
export OPTIMISATION_PATH ?= $(wildcard optimisation/$(LFRIC_TARGET_PLATFORM))

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
export ROOT_DIR    ?= ../..

export SUITE_GROUP ?= developer
export SUITE_GROUP_NAME ?= $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)-.*

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

.PHONY: default
default: build unit-tests integration-tests
	$(Q)echo > /dev/null

.PHONY: documentation doc docs
documentation doc docs: document-uml document-latex document-api
	$(Q)echo > /dev/null

include $(ROOT_DIR)/infrastructure/build/lfric.mk

ANY_UNIT_TESTS := $(shell if [ -d unit-test ] ; then \
                            find unit-test -name *test.pf ; \
                          fi )
ANY_INTEGRATION_TESTS := $(shell if [ -d integration-test ] ; then \
                                   find integration-test -name '*.[Ff]90' -exec egrep -l "^\s*program" {} \; ; \
                                 fi )

##############################################################################
# Launch gscan to monitor suites from this suite-group
#
.PHONY: launch-suite-gscan
launch-suite-gscan: gscan_processes := $(shell ps --no-headers -o command -C cylc-gscan)
launch-suite-gscan: ALWAYS
	$(Q)-if [[ "$(gscan_processes)" != *"--name=$(SUITE_GROUP_NAME)"* ]]; then \
          cylc gscan $(DOUBLE_VERBOSE_ARG) --name=$(SUITE_GROUP_NAME) &            \
          usleep 1                                                                ;\
        fi

##############################################################################

.PHONY: test-suite
test-suite: SUITE_BASE_NAME = $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)
test-suite: SUITE_CONFIG = rose-stem
test-suite: launch-suite-gscan launch-test-suite
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
document-latex: DOCUMENT_DIR  ?= $(PROJECT_DIR)/documents
document-latex: SOURCE_DIR     = documentation
document-latex: DOCUMENTS      = $(shell find $(SOURCE_DIR) -name '*.latex' -print)
document-latex: WORKING_DIR   := $(WORKING_DIR)/latex
document-latex: TEX_STUFF      = documentation/tex
document-latex: COMMON_FIGURES = documentation/common-figures
document-latex: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/latex.mk       \
                                 SOURCE_DIR=$(SOURCE_DIR)         \
	                         WORKING_DIR=$(WORKING_DIR)       \
	                         DOCUMENT_DIR=$(DOCUMENT_DIR)     \
	                         TEX_STUFF=$(TEX_STUFF)           \
	                         COMMON_FIGURES=$(COMMON_FIGURES) \
	                         DOCUMENTS="$(DOCUMENTS)"

.PHONY: document-api
document-api: PROJECT       = $(PROJECT_NAME)
document-api: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/api
document-api: CONFIG_DIR    = documentation
document-api: SOURCE_DIR    = source
document-api: WORKING_DIR  := $(WORKING_DIR)/api
document-api: api-documentation
	$(Q)echo > /dev/null


##############################################################################
# Build
#
.PHONY: build
build: export BIN_DIR     ?= $(PROJECT_DIR)/bin
build: export CXX_LINK     = TRUE
build: export PROGRAMS    := $(basename $(notdir $(shell find source -maxdepth 1 -name '*.[Ff]90' -print)))
build: export PROJECT      = $(PROJECT_NAME)
build: export WORKING_DIR := $(WORKING_DIR)
ifeq "$(PROFILE)" "full-debug"
build: export FFLAG_GROUPS = DEBUG WARNINGS INIT RUNTIME NO_OPTIMISATION FORTRAN_STANDARD
else ifeq "$(PROFILE)" "fast-debug"
build: export FFLAG_GROUPS = DEBUG WARNINGS SAFE_OPTIMISATION FORTRAN_STANDARD
else ifeq "$(PROFILE)" "production"
build: export FFLAG_GROUPS = DEBUG WARNINGS RISKY_OPTIMISATION
else
  $(error Unrecognised profile "$(PROFILE)". Must be one of full-debug, fast-debug or production)
endif
build: ALWAYS
	$(call MESSAGE,===============================)
	$(call MESSAGE,Extracting LFRic infrastructure)
	$(call MESSAGE,===============================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	          SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source \
	          OPTIMISATION_PATH=$(OPTIMISATION_PATH)
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
	          SOURCE_DIR=$(ROOT_DIR)/gungho/source \
	          OPTIMISATION_PATH=$(OPTIMISATION_PATH)
	$(call MESSAGE,==============================)
	$(call MESSAGE,Extracting $(PROJECT_NAME))
	$(call MESSAGE,==============================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=source
	$(call MESSAGE,==============================)
	$(call MESSAGE,Psycloning $(PROJECT_NAME))
	$(call MESSAGE,==============================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	          SOURCE_DIR=source \
	          WORKING_DIR=$(WORKING_DIR) \
	          OPTIMISATION_PATH=$(OPTIMISATION_PATH)
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Processing namelist loaders for $(PROJECT))
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
	          SOURCE_DIR=source \
	          META_FILE_DIR=$(META_FILE_DIR) \
	          FIX_ENUMS_OPT=$(FIX_ENUMS_OPT)
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Analysing $(PROJECT) build dependencies)
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Compiling $(PROJECT))
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk \
	         BIN_DIR=$(BIN_DIR) PROGRAMS="$(PROGRAMS)"


##############################################################################
# Unit tests
#
.PHONY: unit-tests
unit-tests: export BIN_DIR ?= $(PROJECT_DIR)/test
unit-tests: export EXTERNAL_STATIC_LIBRARIES += pfunit
unit-tests: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_NO_OPTIMISATION) $(FFLAGS_INIT)
unit-tests: export IGNORE_DEPENDENCIES += pfunit_mod
unit-tests: export PROGRAMS = $(PROJECT_NAME)_unit_tests
unit-tests: export PROJECT = $(PROJECT_NAME)
unit-tests: export SOURCE_DIR = $(PROJECT_DIR)/source
unit-tests: export TEST_DIR = $(PROJECT_DIR)/unit-test
unit-tests: export WORKING_DIR := $(WORKING_DIR)
unit-tests: do-unit-tests
	$(Q)echo > /dev/null


##############################################################################
# Integration tests
#
.PHONY: integration-tests
integration-tests: export BIN_DIR ?= $(PROJECT_DIR)/test
integration-tests: export PROJECT = $(PROJECT_NAME)
integration-tests: export SOURCE_DIR = $(PROJECT_DIR)/source
integration-tests: export TEST_DIR = $(PROJECT_DIR)/integration-test
integration-tests: export WORKING_DIR := $(WORKING_DIR)
integration-tests: do-integration-tests
	$(Q)echo > /dev/null


##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"$(PROJECT_NAME) work space")
	$(Q)-if [ $(WORKING_DIR) != *[\*]* ] && [ -d $(WORKING_DIR) ] ; then rm -r $(WORKING_DIR) ; fi
	$(call MESSAGE,Removing,"$(PROJECT_NAME) documents")
	$(Q)-if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"$(PROJECT_NAME) binaries")
	$(Q)-if [ -d bin ] ; then rm -r bin ; fi
	$(Q)-if [ -d test ] ; then rm -r test ; fi
