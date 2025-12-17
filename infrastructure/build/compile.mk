##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
#
# Locate and compile all Fortran source in the current directory.
#
# This is done as a standalone make file so that your base makefile does not
# require the dependency analysis. This is annoying if you issue a
# "make clean" on an already clean working copy. In that case the analysis
# will be performed in order that it then be deleted.
#
# The following variables may be specified to modify the compile process...
#
# ROOT: Project directory
# BIN_DIR: Path to directory for resulting executables.
#          Default: $(ROOT)/bin
# CXX_LINK: Set this macro to have the C++ runtime library linked to the
#           executable.
# EXTERNAL_LIBRARIES: Libraries required by link stage. These should be
#                     specified in static link order, even if you are linking
#                     dynamically.
# FFLAG_GROUPS: Space separated list of FFLAG_<group name> variables to use in
#               building up the FFLAGS variable. Only the group name is
#               specified.
# LINK_TYPE: 'static' or 'dynamic' linking.
#            Default: dynamic, except on Crays where it's static
# PROGRAMS: Names of programs to compile.
#           Default: Everything listed in programs.mk
# PRE_PROCESS_INCLUDE_DIRS: Space separated list of directories to search for
#                           inclusions.
# PRE_PROCESS_MACROS: Space separated list of macro definitions in the form
#                     NAME[=MACRO] to be passed to the compiler.
# PROJECT_MAKE_DIR: Used to locate project specific targets modifiers and
#                   such.
# COMPILE_OPTIONS: Name of an optional file that can be included to list
#                  project specific compile options
#
##############################################################################

.SECONDEXPANSION:

# Build a set of "-I" arguments to seach the whole object tree:
INCLUDE_ARGS := $(subst ./,-I,$(shell find . -mindepth 1 -type d -print)) \
                $(addprefix -I, $(PRE_PROCESS_INCLUDE_DIRS))

# Build a set of "-D" argument for any pre-processor macros
#
MACRO_ARGS := $(addprefix -D,$(PRE_PROCESS_MACROS))

include programs.mk

PROGRAMS ?= $(basename $(notdir $(PROG_OBJS)))

# Convert the program names to all caps, append "_OBJS" to the end and
# dereference that to get a list of all objects needed by all programs:
#
ALL_OBJECTS = $(foreach proj, $(shell echo $(PROGRAMS) | tr a-z A-Z), $($(proj)_OBJS))

.PHONY: applications
applications: FFLAGS_BASE = $(FFLAGS) $(foreach group, $(FFLAG_GROUPS), $(FFLAGS_$(group)))
applications: LDFLAGS_BASE = $(foreach group, $(LDFLAGS_GROUPS), $(LDFLAGS_$(group)))
applications: $(addprefix $(BIN_DIR)/, $(PROGRAMS))

.PHONY: libraries
libraries: FFLAGS_BASE = $(FFLAGS) $(foreach group, $(FFLAG_GROUPS), $(FFLAGS_$(group)))
libraries: LDFLAGS_BASE = $(foreach group, $(LDFLAGS_GROUPS), $(LDFLAGS_$(group)))
libraries: $(addsuffix .a, $(addprefix $(LIB_DIR)/lib, $(PROGRAMS)))

##############################################################################

include $(LFRIC_BUILD)/lfric.mk
include $(LFRIC_BUILD)/fortran.mk
include $(LFRIC_BUILD)/cxx.mk
-include $(COMPILE_OPTIONS)

BIN_DIR ?= $(ROOT)/bin
LIB_DIR ?= lib
MOD_DIR ?= mod

# Determine what performance tool to use.
#
# Technically the "time" tool we want to use is part of the operating system
# and not the kernel. It's not clear what "Linux" means in the context of the
# result from "uname" but it should be good enough for the moment.
#
ifeq ($(shell uname),Linux)
  $(warning linux)
  TIME_TOOL = /usr/bin/time -f "Compiled $<: Wallclock=%E, Highwater=%MKiB"
else
  $(warning not linux)
  TIME_TOOL = time
endif

# If the compiler produces module files, tell it where to put them
#
ifdef F_MOD_DESTINATION_ARG
  MODULE_DESTINATION_ARGUMENT = $(F_MOD_DESTINATION_ARG)$(MOD_DIR)
  MODULE_SOURCE_ARGUMENT = $(F_MOD_SOURCE_ARG)$(MOD_DIR)
else
$(error Compilers which do not produce module files are not currently supported)
endif

ifdef CXX_LINK
  EXTERNAL_DYNAMIC_LIBRARIES += $(CXX_RUNTIME_LIBRARY)
endif

ifdef CRAY_ENVIRONMENT
  $(warning Running on a Cray, selecting static linking)
  LINK_TYPE ?= static
else
  LINK_TYPE ?= dynamic
endif

# Work out what to do with external libraries.
#
ifeq "$(LINK_TYPE)" "static"
  override EXTERNAL_STATIC_LIBRARIES  := $(EXTERNAL_STATIC_LIBRARIES) $(EXTERNAL_DYNAMIC_LIBRARIES)
  override EXTERNAL_DYNAMIC_LIBRARIES :=
else ifeq "$(LINK_TYPE)" "dynamic"
  # Nothing further needs to be done.
else
  $(error Unrecognised LINK_TYPE. Must be either "static" or "dynamic")
endif

ifdef NO_MPI
  LINKER = $(FC)
else
  LINKER = $(LDMPI)
endif

##############################################################################

.PRECIOUS: $(BIN_DIR)/%
$(BIN_DIR)/%: %.o $$(LIB_DIR)/lib$$(*F).a
	$(call MESSAGE,Linking,$*)
	$(Q)mkdir -p $(@D)
	$(Q)$(LINKER) $(LDFLAGS) $(LDFLAGS_BASE) $(LDFLAGS_COMPILER) -o $@ $^ \
	            $(patsubst %,-l%,$(EXTERNAL_STATIC_LIBRARIES)) \
	            $(patsubst %,-l%,$(EXTERNAL_DYNAMIC_LIBRARIES))

.PRECIOUS: $(LIB_DIR)/lib%.a
$(LIB_DIR)/lib%.a: $$($$(shell basename $$* | tr a-z A-Z)_OBJS) | $(LIB_DIR)
	$(call MESSAGE,Archiving,$(@F))
	$(Q) ar -rcs $@ $^

.PRECIOUS: %.o
%.o: %.f90 | $(MOD_DIR)
	$(call MESSAGE,Compile,$<)
	$(Q)$(TIME_TOOL) $(FC) $(FFLAGS_BASE) $(FFLAGS_EXTRA)\
	          $(MODULE_DESTINATION_ARGUMENT) \
	          $(MODULE_SOURCE_ARGUMENT) \
	          $(INCLUDE_ARGS) -c -o $(basename $@).o $<
	$(call MESSAGE,Compiled,$<)

%.o: %.F90 | $(MOD_DIR)
	$(call MESSAGE,Pre-process and compile,$<)
	$(Q)$(TIME_TOOL) $(FC) $(FFLAGS_BASE) $(FFLAGS_EXTRA) \
	          $(MODULE_DESTINATION_ARGUMENT) \
	          $(MODULE_SOURCE_ARGUMENT) \
	          $(INCLUDE_ARGS) $(MACRO_ARGS) -c -o $(basename $@).o $<
	$(call MESSAGE,Compiled,$<)


#############################################################################
# Directories

# Sort removes duplicates in this list
$(sort $(LIB_DIR) $(MOD_DIR)):
	$(call MESSAGE,Creating,$@)
	$(Q)mkdir -p $@

#############################################################################
# Dependencies
#
include dependencies.mk
