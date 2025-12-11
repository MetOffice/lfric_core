##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
#
# Include this file from your model make file in order to gain access to the
# LFRic build system. Include it at the end of the make file as it contains
# targets which you do not want to become the default target.
#
# Variables provided by including this file...
#
# LFRIC_BUILD: Path to the build system
# COMPILE_OPTIONS: File of target-specific compile options used in compile.mk
#
# Macros expected by the build system are as follows...
#
# APPS_ROOT_DIR: Absolute path to an lfric_apps directory (not expected for
#                core only tasks)
# CORE_ROOT_DIR: Absolute path to the project directory for lfric_core
# OPTIMISATION_PATH: Where PSyclone optimisation scripts may be found.
# WORKING_DIR: Path to scratch space in which intermediate files will be
#              placed. This should be somewhere with good "many small
#              files" performance, i.e. probably not Lustre.
#              Default: ./working
# VERBOSE: Set in order to see actual commands issued by the build system.
# PURGE_SUITES: Set to non-zero value to clean out exisiting rose suites of
#               same name. (this is also the default action)
#               Set to zero to not clean out existing rose suite.
# TEST_SUITE_TARGETS: Space separated list of target identifiers to be used
#                     when launching the test suite. Default is "meto-azspice
#                     meto-ex1a"
# SUITE_GROUP_ABRV: Set to a non-zero value to cause the names of the rose
#                   stem groups to always be abbreviated in the suite name.
#                   Set to zero to cause the names to always be unabbreviated.
#                   The default is abbreviated for multi-group runs, but
#                   unabbreviated for single-group runs.
#
##############################################################################

.SECONDEXPANSION:

# Ensure make offers the features we need...
#
$(info ** Make version $(MAKE_VERSION))
REQUIRED_FEATURES = else-if order-only second-expansion target-specific
ifneq ("x$(filter-out $(value .FEATURES), $(REQUIRED_FEATURES))", "x")
  $(error The build system requires support for the following from GMake: $(REQUIRED_FEATURES))
endif
#
# But of course one of the features is not covered by the .FEATURES list...
#
INT_MAKE_VERSION := $(shell echo $(MAKE_VERSION) \
                    | awk -F . '{ for (i=1; i<NF+1; ++i) if (i == 1) printf "%i", $$i; else printf "%02i", $$i }')
MAKE_AGE := $(shell if [[ "$(INT_MAKE_VERSION)" -lt "382" ]]; then echo old; fi )
ifneq ("x$(MAKE_AGE)", "x")
  $(error The build system requires at least GMake 3.82)
endif


# Default variables...
#
export WORKING_DIR ?= working
export PWD ?= $(shell pwd)

TEST_SUITE_TARGETS ?= meto-azspice meto-ex1a

# Make the build system available...
#
export LFRIC_BUILD := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

# Make the infrastructure available...
#
export LFRIC_INFRASTRUCTURE := $(realpath $(LFRIC_BUILD)/..)


# Attempt to identify Cray systems...
#
ifdef PE_ENV
  CRAY_ENVIRONMENT = true
  export CRAY_ENVIRONMENT
endif

ifdef NO_MPI
  export PRE_PROCESS_MACROS += NO_MPI=no_mpi
endif

# Currently default to legacy MPI, rather than use the not completely well-
# supported mpi_f08 (which would be the default in the absence of any setting)
#
ifndef USE_MPI_F08
  export PRE_PROCESS_MACROS += LEGACY_MPI
endif

ifdef USE_VERNIER
  export PRE_PROCESS_MACROS += VERNIER
endif

ifdef USE_TIMING_WRAPPER
  export PRE_PROCESS_MACROS += TIMING_ON
endif


# Set the default precision for reals
RDEF_PRECISION ?= 64
export PRE_PROCESS_MACROS += RDEF_PRECISION=$(RDEF_PRECISION)

# Set the r_solver precision for reals
R_SOLVER_PRECISION ?= 32
export PRE_PROCESS_MACROS += R_SOLVER_PRECISION=$(R_SOLVER_PRECISION)

# Set the r_tran precision for reals
R_TRAN_PRECISION ?= 64
export PRE_PROCESS_MACROS += R_TRAN_PRECISION=$(R_TRAN_PRECISION)

# Set the r_bl precision for reals
R_BL_PRECISION ?= 64
export PRE_PROCESS_MACROS += R_BL_PRECISION=$(R_BL_PRECISION)

# The compile options file overrides compile options based on file-name pattern matching.
# Use the miniapp-specific file if it exists. Otherwise use the infrastructure file.
ifeq (,$(wildcard $(PROJECT_DIR)/build/compile_options.mk))
  export COMPILE_OPTIONS := $(abspath $(CORE_ROOT_DIR)/infrastructure/build/compile_options.mk)
else
  export COMPILE_OPTIONS := $(abspath $(PROJECT_DIR)/build/compile_options.mk)
endif

# Set up verbose logging...
#
ifdef VERBOSE
  Q :=
  VERBOSE_ARG = -verbose
  SHORT_VERBOSE_ARG = -v
  DOUBLE_VERBOSE_ARG = --verbose
else
  Q := @
  QUIET_ARG = --quiet
  QUIET_ARG_SINGLE = -quiet
  VERBOSE_REDIRECT = >/dev/null
endif

export Q QUIET_ARG VERBOSE_REDIRECT

# Set flag to perform a fresh rose stem suite

CLEAN_OPT ?= '--new'
ifeq '$(PURGE_SUITES)' '0'
  CLEAN_OPT =
endif

# We only want to send terminal control characters if there is a terminal to
# interpret them...
#
_NOW = `date +%H:%M:%S`
ifneq 'x$(TERM)' 'x'
  MESSAGE = $(Q)printf "%s \\033[1m$(1)\\033[0m %s\n" $(_NOW) $(2)
else
  MESSAGE = $(Q)echo $(_NOW) *$(1)* $(2)
endif

# Set up some special macros for hard to escape characters
#
EMPTY :=
SPACE := $(EMPTY) # This comment highlights space character.
PERCENT := %
OPEN_PAREN := (
COMMA := ,

# Prerequisite for targets which should always be run.
#
.PHONY: ALWAYS
ALWAYS:

# The directory containing the target. Useful for order-only prerequisites to
# create that directory.
#
TARGET_DIR = $(patsubst $(PERCENT)/,$(PERCENT),$(dir $@))

# Default tool executables.
#
export INKSCAPE ?= inkscape

define ANNOUNCE_BODY

*******************************************************************************
Error

The Cylc 7 test suite is no longer in use, please use the Cylc 8 suite.

New commands:

export CYLC_VERSION=8
rose stem --group=developer
cylc play <working_copy_name>

*******************************************************************************

endef

##############################################################################
# Build UML documentation
#
# DOCUMENT_DIR - Directory in which documentation is placed.
# SOURCE_DIR   - Root directory of source tree.
# WORKING_DIR  - Location for temporary files.
#
.PHONY: uml-documentation
uml-documentation:
	$(Q)$(MAKE) $(QUIET_ARG) uml-pdfs DOCUMENT_DIR=$(DOCUMENT_DIR) SOURCE_DIR=$(SOURCE_DIR) WORKING_DIR=$(WORKING_DIR)

.PHONY: uml-pdfs
uml-pdfs: $$(patsubst $$(SOURCE_DIR)/$$(PERCENT).puml, \
                      $$(DOCUMENT_DIR)/$$(PERCENT).pdf, \
                      $$(wildcard $$(SOURCE_DIR)/*.puml))
	$(Q)echo >/dev/null

.PRECIOUS: $(DOCUMENT_DIR)/%.pdf
$(DOCUMENT_DIR)/%.pdf: $(DOCUMENT_DIR)/%.svg
	$(call MESSAGE,Translating,$(notdir $<))
	$Q$(INKSCAPE) $< --export-pdf=$@

.PRECIOUS: $(DOCUMENT_DIR)/%.svg
$(DOCUMENT_DIR)/%.svg: $(SOURCE_DIR)/%.puml \
                       $$(addprefix $$(SOURCE_DIR)/, \
                                    $$(shell sed -n -e 's/!include[ ]*\([^ \n]*\)/\1/p' $$(SOURCE_DIR)/$$*.puml))
	$(call MESSAGE,Generating,$(notdir $@))
	$(Q)mkdir -p $(DOCUMENT_DIR)
	$(Q)plantuml $(SHORT_VERBOSE_ARG) -tsvg -o $(abspath $(dir $@)) $(abspath $<)

##############################################################################
# Build API documentation
#
# PROJECT      - Name of the project for logging purposes.
# DOCUMENT_DIR - Directory in which documentation is placed.
# CONFIG_DIR   - Directory where Doxyfile can be found.
# SOURCE_DIR   - Root directory of source tree.
# WORKING_DIR  - Location for temporary files.
#
.PHONY: api-documentation
api-documentation: ALWAYS
	$(call MESSAGE,API,$(PROJECT))
	$(Q)mkdir -p $(DOCUMENT_DIR)
	$(Q)( cat $(CONFIG_DIR)/Doxyfile; \
	      echo INPUT=$(SOURCE_DIR) $(CONFIG_DIR)/$$(sed -n -e 's/\s*INPUT\s*=\s*//p' $(CONFIG_DIR)/Doxyfile); \
	      echo USE_MDFILE_AS_MAINPAGE=$(CONFIG_DIR)/$$(sed -n -e 's/\s*USE_MDFILE_AS_MAINPAGE\s*=\s*//p' $(CONFIG_DIR)/Doxyfile); \
	      echo OUTPUT_DIRECTORY=$(DOCUMENT_DIR) ) \
	    | doxygen - $(VERBOSE_REDIRECT)

##############################################################################
# Launch test suite
#
# SUITE_CONFIG    - Path to rose-stem directory.
# SUITE_BASE_NAME - Name for suites.
# SUITE_GROUP_NAME_ABRV - Name(s) of the rose stem group(s) with abbreviations applied.
#
.PHONY: launch-test-suite
launch-test-suite:
	$(error $(ANNOUNCE_BODY) )


##############################################################################
# Generate configuration source.
#
.PHONY: configuration
configuration:
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk


##############################################################################
# Extract parts of other projects. No source generation.
#
.PHONY: %/extract
%/extract:
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	            SOURCE_DIR=$* WORKING_DIR=$(WORKING_DIR)

##############################################################################
# Import another project including generated source.
#
.PHONY: %/import
%/import: $$*/build/import.mk  # If there are special import instructions.
	$Q$(MAKE) $(QUIET_ARG) -f $<

%/import: $$**/source  # In the absense of special instructions.
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$<
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
	          SOURCE_DIR=$<

##############################################################################
# Invoke PSyclone to generate PSy layer.
#
# Psyclone is called on the original source but that source may use other
# modules so extraction must be complete first.
#
.PHONY: %/psyclone
%/psyclone: $$(addsuffix /extract, $$*)
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone_psykal.mk \
	            SOURCE_DIR=$* \
	            WORKING_DIR=$(WORKING_DIR)


##############################################################################
# Run unit tests.
#
# We recurse into the make file here in order to reify all the target
# specific variables. They only appear in recipes and we need them to
# make decisions.
#
.PHONY: unit-tests/%
unit-tests/%: export LDFLAGS_GROUPS = OPENMP
unit-tests/%: export FFLAG_GROUPS = OPENMP DEBUG RUNTIME NO_OPTIMISATION INIT UNIT_WARNINGS
unit-tests/%:
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/tests.mk do-unit-test/$*


##############################################################################
# Run integration tests.
#
# We recurse into the make file here in order to reify all the target
# specific variables. They only appear in recipes and we need them to
# make decisions.
#
.PHONY: integration-tests/%
integration-tests/%: export LDFLAGS_GROUPS = OPENMP
integration-tests/%: export FFLAG_GROUPS = OPENMP DEBUG RUNTIME NO_OPTIMISATION INIT UNIT_WARNINGS
integration-tests/%:
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/tests.mk do-integration-tests/$*

###############################################################################
# End
