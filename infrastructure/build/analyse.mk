##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
# Some of the content of this file has been produced with the assistance of
# Met Office Github Copilot Enterprise."
#
# Scan all Fortran source files in the current directory and build up
# dependency information.
#
# This is done as a stand alone make file as it appears the Python bindings
# to SQLite screw up the latter's multi-thread support. Thus this make file
# needs to operate in a single thread regime. We don't want to impose that
# restriction on the rest of the build system.
#
# All out of date source files are handed to a single invocation of the
# analyser rather than one invocation per file. Starting a Python interpreter,
# importing its modules and opening the database for every source file
# dominated the cost of this phase. The analyser parallelises the reading and
# preprocessing of source internally while keeping database access on a single
# thread.
#
# The following variables may be specified to modify behaviour:
#
# PRE_PROCESS_INCLUDE_DIRS: Space separated list of directories to search for
#                           inclusions.
# PRE_PROCESS_MACROS: Space separated list of macro definitions in the form
#                     NAME[=MACRO] to be passed to the preprocessor.
#
##############################################################################

.NOTPARALLEL:

DATABASE ?= dependencies.db

SOURCE_FILES := $(subst ./,,$(shell find . -name '*.[Ff]90' -print))

# Record of which sources have been analysed. Make compares the modification
# time of this file against each source so only those which have changed are
# presented to the analyser, through $?.
#
ANALYSED_STAMP = analysed.stamp

programs.mk: dependencies.mk
	$(call MESSAGE,Collating,$@)
	$(Q)$(LFRIC_BUILD)/tools/ProgramObjects $(VERBOSE_ARG) \
                                                -database $(DATABASE) \
	                                        -objectdir . $@

dependencies.mk: $(ANALYSED_STAMP)
	$(call MESSAGE,Building,$@)
	$(Q)$(LFRIC_BUILD)/tools/DependencyRules $(VERBOSE_ARG) \
                                                 -database $(DATABASE) \
	                                         -objectdir . \
	                                         -moduledir . \
	                                         $(DEPRULE_FLAGS) $@

IGNORE_ARGUMENTS = $(addprefix -ignore ,$(IGNORE_DEPENDENCIES))
INCLUDE_ARGUMENTS = $(addprefix -include , $(PRE_PROCESS_INCLUDE_DIRS))
MACRO_ARGUMENTS = $(addprefix -macro , $(PRE_PROCESS_MACROS))

# All out of date sources are passed to a single invocation of the analyser.
# Launching one Python interpreter per source file dominated the cost of this
# phase of the build.
#
# The analyser removes a file's existing entries before rescanning it so
# presenting only the changed subset leaves the database consistent.
#
$(ANALYSED_STAMP): $(SOURCE_FILES)
	$(call MESSAGE,Analysing,$(words $?) source files)
	$(Q)$(LFRIC_BUILD)/tools/DependencyAnalyser \
	    $(IGNORE_ARGUMENTS) $(INCLUDE_ARGUMENTS) $(MACRO_ARGUMENTS) \
	    $(VERBOSE_ARG) $(DATABASE) $?
	$(Q)touch $@

include $(LFRIC_BUILD)/lfric.mk
include $(LFRIC_BUILD)/fortran.mk
