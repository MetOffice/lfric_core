##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

PROJECT_NAME = adjoint
PROFILE ?= fast-debug

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
export ROOT_DIR    ?= ..

export INTERNAL_DEPENDENCIES = $(ROOT_DIR)/infrastructure \
                               $(ROOT_DIR)/components/driver \
                               $(ROOT_DIR)/components/science \
                               $(ROOT_DIR)/components/inventory \
                               $(ROOT_DIR)/components/lfric-xios \
                               $(ROOT_DIR)/gungho \
                               $(ROOT_DIR)/linear

export SUITE_GROUP ?= developer
export SUITE_GROUP_NAME ?= $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)-.*

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

.PHONY: default
default: unit-tests integration-tests
	$(Q)echo > /dev/null

.PHONY: documentation doc docs
documentation doc docs: document-api
	$(Q)echo > /dev/null

include $(ROOT_DIR)/infrastructure/build/lfric.mk
include $(INTERNAL_DEPENDENCIES:=/build/import.mk)

##############################################################################
# Documentation
#
.PHONY: document-api
document-api: PROJECT       = adjoint
document-api: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/api
document-api: CONFIG_DIR    = documentation
document-api: SOURCE_DIR    = source
document-api: WORKING_DIR  := $(WORKING_DIR)/api
document-api: api-documentation
	$(Q)echo > /dev/null

.PHONY: document-uml
document-uml: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/uml
document-uml: SOURCE_DIR    = documentation/uml
document-uml: WORKING_DIR  := $(WORKING_DIR)/uml
document-uml: uml-documentation

.PHONY: document-latex
document-latex: DOCUMENT_DIR  ?= $(PROJECT_DIR)/documents
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

##############################################################################
# Unit tests
#
unit-tests:
	$(call MESSAGE,Testing,"There are no unit tests.")


##############################################################################
# Integration tests
#
integration-tests:
	$(call MESSAGE,Testing,"There are no integration tests.")


##############################################################################
# System tests
#
.PHONY: test-suite
test-suite: SUITE_BASE_NAME = $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)
test-suite: SUITE_CONFIG = rose-stem
test-suite: launch-suite-gscan launch-test-suite
	$(Q)echo > /dev/null

.PHONY: launch-suite-gscan
launch-suite-gscan: gscan_processes := $(shell ps --no-headers -o command -C cylc-gscan)
launch-suite-gscan: ALWAYS
	$(Q)-if [[ "$(gscan_processes)" != *"--name=$(SUITE_GROUP_NAME)"* ]]; \
        then \
          cylc gscan $(DOUBLE_VERBOSE_ARG) --name=$(SUITE_GROUP_NAME) & \
          usleep 1;\
        fi

##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"$(PROJECT_NAME) work space")
	$(Q)-if [ $(WORKING_DIR) != *[\*]* ] && [ -d $(WORKING_DIR) ] ; then rm -r $(WORKING_DIR) ; fi
	$(call MESSAGE,Removing,"$(PROJECT_NAME) documents")
	$(Q)if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"$(PROJECT_NAME) binaries")
	$(Q)-if [ -d test ] ; then rm -r test ; fi
