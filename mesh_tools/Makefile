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

# Select option for PSyclone transformation. Select "none" for no transformation
export MESH_PSYCLONE_TRANSFORMATION ?= none

# This top level makefile is very order sensitive. Source code extraction and
# generation must happen in a certain order. Due to this we turn off
# multithreading for this file only. Any called recursively (i.e. with $(MAKE))
# run in parallel. Unless they specify NOTPARALLEL as well.
#
.NOTPARALLEL:

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
export CORE_ROOT_DIR    ?= ..

export INTERNAL_DEPENDENCIES = $(CORE_ROOT_DIR)/infrastructure \
                               $(CORE_ROOT_DIR)/components/science

export SUITE_GROUP ?= developer
export SUITE_GROUP_NAME ?= $(notdir $(realpath $(shell pwd)/$(CORE_ROOT_DIR)))-$(PROJECT_NAME)-.*

export PRE_PROCESS_MACROS := MESH_TOOLS

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

# Check PSyclone transformation option is valid
ifneq ("$(MESH_PSYCLONE_TRANSFORMATION)","none")
  export OPTIMISATION_PATH ?= $(wildcard optimisation/$(MESH_PSYCLONE_TRANSFORMATION))
  ifeq ("$(wildcard $(OPTIMISATION_PATH)/global.py)","")
    $(error PSyclone transformation dir "$(MESH_PSYCLONE_TRANSFORMATION)" does not have a global.py script)
  endif
endif

.PHONY: default
default: build unit-tests
	$(Q)echo > /dev/null

include $(CORE_ROOT_DIR)/infrastructure/build/lfric.mk
include $(INTERNAL_DEPENDENCIES:=/build/import.mk)

##############################################################################
# Build
#
.PHONY: build
build: BIN_DIR     ?= $(PROJECT_DIR)/bin
build: PROGRAMS    := $(basename $(notdir $(shell find source -maxdepth 1 -name '*.[Ff]90' -print)))
build: PROJECT      = mesh_tools
build: WORKING_DIR := $(WORKING_DIR)/mesh_tools

build: export LDFLAGS_GROUPS = OPENMP
ifeq "$(PROFILE)" "full-debug"
build: export FFLAG_GROUPS = OPENMP DEBUG WARNINGS INIT RUNTIME NO_OPTIMISATION FORTRAN_STANDARD
else ifeq "$(PROFILE)" "fast-debug"
build: export FFLAG_GROUPS = OPENMP DEBUG WARNINGS SAFE_OPTIMISATION FORTRAN_STANDARD
else ifeq "$(PROFILE)" "production"
build: export FFLAG_GROUPS = OPENMP DEBUG WARNINGS RISKY_OPTIMISATION
else
  $(error Unrecognised profile "$(PROFILE)". Must be one of full-debug, fast-debug or production)
endif
build: ALWAYS
	$(call MESSAGE,========================================)
	$(call MESSAGE,Importing internal dependencies...)
	$(call MESSAGE,========================================)
	$(Q)for SUBPROJECT in $(INTERNAL_DEPENDENCIES) ; do \
		$(MAKE) $(QUIET_ARG) -f $$SUBPROJECT/build/import.mk ; done
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting $(PROJECT_NAME))
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=source \
	          WORKING_DIR=$(WORKING_DIR)
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Generating $(PROJECT) namelist loaders)
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
	          PROJECT=$(PROJECT) \
	          SOURCE_DIR=source \
	          WORKING_DIR=$(WORKING_DIR) \
	          META_FILE_DIR=$(META_FILE_DIR)
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
unit-tests/%: export BIN_DIR       ?= $(PROJECT_DIR)/test
unit-tests/%: export IMPORT_PARTS   = $(CORE_ROOT_DIR)/infrastructure \
                                      $(CORE_ROOT_DIR)/components/science
unit-tests/%: export META_FILE_DIR  = rose-meta/lfric-mesh_tools/HEAD
unit-tests/%: export PROGRAMS       = mesh_tools_unit_tests
unit-tests/%: export PROJECT        = mesh_tools
unit-tests/%: export SOURCE_DIR     = source
unit-tests/%: export TEST_DIR       = unit-test
unit-tests/%: export WORKING_DIR   := $(WORKING_DIR)/unit-test
unit-tests: unit-tests/run


##############################################################################
# Integration tests
#
integration-tests:
	$(call MESSAGE,Testing,"There are no integration tests.")


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
