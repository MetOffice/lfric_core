##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
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
# ROOT_DIR: Absolute path to the project directory
# OPTIMISATION_PATH: Where PSyclone optimisation scripts may be found.
# WORKING_DIR: Path to scratch space in which intermediate files will be
#              placed. This should be somewhere with good "many small
#              files" performance, i.e. probably not Lustre.
#              Default: ./working
# VERBOSE: Set in order to see actual commands issued by the build system.
# PURGE_SUITES: Set in order clean out exisiting rose suites of same name.
#
# Plus the normal compiler macros...
#
# FC: Fortran compiler command
# FFLAGS: Fortran compiler flags including "-I" module file locations
# LD: Linker command, probably the same as FC
# MPILD: MPI linker command
# LDFLAGS: Linker flags including "-L" library locations
#
##############################################################################

.SECONDEXPANSION:

# Ensure make offers the features we need...
#
$(info ** Make version $(MAKE_VERSION))
ifeq ($(filter else-if,$(value .FEATURES)),)
  $(error The build system requires else-if support from GMake)
endif

# The default value of FC is almost always "f77" which is of no use to us.
# An empty FC is also of no use.
ifneq "$(or $(filter default, $(origin FC)), $(filter x, x$(FC)))" ""
  $(error The FC environment variable must be set to a Fortran compiler command)
endif

# Default variables...
#
export WORKING_DIR ?= working
export PWD ?= $(shell pwd)

# Make the build system available...
#
export LFRIC_BUILD := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

# Make the infrastructure available...
#
export LFRIC_INFRASTRUCTURE := $(realpath $(LFRIC_BUILD)/..)

# Make compiler macros available...
#
# Sometimes FC holds a full path which needs to be stripped off. It may also
# include a version number which also needs to go.
#
FORTRAN_COMPILER := $(firstword $(subst -, ,$(notdir $(FC))))

# Attempt to identify Cray systems...
#
ifdef PE_ENV
  CRAY_ENVIRONMENT = true
  ifeq '$(PE_ENV)' 'CRAY'
    FORTRAN_COMPILER = crayftn
  else ifeq '$(PE_ENV)' 'INTEL'
    FORTRAN_COMPILER = ifort
  else ifeq '$(PE_ENV)' 'GNU'
    FORTRAN_COMPILER = gfortran
  else ifeq '$(PE_ENV)' 'PGI'
    FORTRAN_COMPILER = pgfortran
  else
    $(error Unrecognised Cray programming environment)
  endif
endif

include $(LFRIC_BUILD)/fortran/$(FORTRAN_COMPILER).mk
export F_MOD_DESTINATION_ARG OPENMP_ARG

FFLAGS += $(FFLAGS_COMPILER)
export FFLAGS

# As XIOS is written in C++ we have to concern ourselves with
# compilers for that as well.
#
# We assume GCC is in use unless the CXX environment variable is set.
#
CXX ?= g++

# Try to work out the compiler from CXX which may contain a full path and a
# version number.
#
CXX_COMPILER := $(firstword $(subst -, ,$(notdir $(CXX))))

# Attempt to identify Cray systems...
ifdef CRAY_ENVIRONMENT
  ifeq '$(PE_ENV)' 'CRAY'
    CXX_COMPILER = craycc
  else ifeq '$(PE_ENV)' 'INTEL'
    CXX_COMPILER = icc
  else ifeq '$(PE_ENV)' 'GNU'
    CXX_COMPILER = g++
  else ifeq '$(PE_ENV)' 'PGI'
    CXX_COMPILER = pgc++
  else
    $(error Unrecognised Cray programming environment)
  endif
endif

include $(LFRIC_BUILD)/cxx/$(CXX_COMPILER).mk
export CXX_RUNTIME_LIBRARY

# The compile options file overrides compile options based on file-name pattern matching.
# Use the miniapp-specific file if it exists. Otherwise use the infrastructure file.
ifeq (,$(wildcard $(PROJECT_DIR)/build/compile_options.mk))
  export COMPILE_OPTIONS := $(abspath $(ROOT_DIR)/infrastructure/build/compile_options.mk)
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
  VERBOSE_REDIRECT = >/dev/null
endif

export Q QUIET_ARG VERBOSE_REDIRECT

# Set flag to perform a fresh rose stem suite
ifdef PURGE_SUITES
  CLEAN_OPT :='--new'
endif

# We only want to send terminal control characters if there is a terminal to
# interpret them...
#
ifneq 'x$(TERM)' 'x'
  MESSAGE = $(Q)printf "\\033[1m$(1)\\033[0m %s\n" $(2)
else
  MESSAGE = $(Q)echo *$(1)* $(2)
endif

# Set up some special macros for hard to escape characters
#
EMPTY :=
SPACE := $(EMPTY) # This comment highlights space character.
PERCENT := %
OPEN_PAREN := (

# Prerequisite for targets which should always be run.
#
.PHONY: ALWAYS
ALWAYS:

# The directory containing the target. Useful for order-only prerequisites to
# create that directory.
#
TARGET_DIR = $(patsubst $(PERCENT)/,$(PERCENT),$(dir $@))

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
	$(Q)inkscape $< --export-pdf=$@

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
#
.PHONY: launch-test-suite
launch-test-suite: SUITE_GROUP ?= developer
launch-test-suite: TEST_SUITE_TARGETS ?= $(error Please set the TEST_SUITE_TARGETS environment variable.)
launch-test-suite: SUITE_NAME = $(SUITE_BASE_NAME)-$$target-$(SUITE_GROUP)
launch-test-suite:
	$(Q)umask 022; for target in $(TEST_SUITE_TARGETS) ; do \
	echo Launching $(PROJECT_NAME) test suite against $$target with $(SUITE_GROUP) group ; \
	rose stem --name=$(SUITE_NAME) \
	          --config=$(SUITE_CONFIG) \
	          --opt-conf-key=$$target \
	          $(CLEAN_OPT) $(QUIET_ARG) \
	          --group=$(SUITE_GROUP) ; \
	done

##############################################################################
# Run unit tests.
#
.PHONY: run-unit-tests
run-unit-tests: generate-unit-tests
	$(Q)$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$(Q)$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk
	$(call MESSAGE,Running,$(PROGRAMS))
	$(Q)cd $(WORKING_DIR); mpiexec -n 1 $(BIN_DIR)/$(PROGRAMS) $(DOUBLE_VERBOSE_ARG)


##############################################################################
# Run integration tests.
#
.PHONY: integration-test-run/%
integration-test-run/%: export PYTHONPATH := $(PYTHONPATH):$(LFRIC_BUILD)
integration-test-run/%: generate-integration-tests
	$(Q)$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk PROGRAMS=$(notdir $*)
	$(Q)$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk \
	        PROGRAMS=$(notdir $*) FFLAGS="$(FFLAGS) $(FFLAGS_DEBUG) $(FFLAGS_RUNTIME)"
	$(call MESSAGE,Running,$*)
	$(Q)cd $(dir $*); \
	    ./$(notdir $(addsuffix .py,$*)) $(addprefix $(BIN_DIR)/,$(notdir $*))

##############################################################################
# Generate configuration source.
#
.PHONY: configuration
configuration:
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk

##############################################################################
# Generate pFUnit unit tests.
#
.PHONY: pfunit
pfunit:
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/pfunit.mk
