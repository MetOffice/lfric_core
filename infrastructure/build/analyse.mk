##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
#
# Scan all Fortran source files in the current directory and build up
# dependency information.
#
# This is done as a stand alone make file as it appears the Python bindings
# to SQLite screw up the latter's multi-thread support. Thus this make file
# needs to operate in a single thread regime. We don't want to impose that
# restriction on the rest of the build system.
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
TOUCH_FILES = $(subst .F90,.t,$(subst .f90,.t,$(SOURCE_FILES)))

programs.mk: dependencies.mk
	$(call MESSAGE,Collating,$@)
	$(Q)$(LFRIC_BUILD)/tools/ProgramObjects $(VERBOSE_ARG) \
                                                -database $(DATABASE) \
	                                        -objectdir . $@

dependencies.mk: $(TOUCH_FILES)
	$(call MESSAGE,Building,$@)
	$(Q)$(LFRIC_BUILD)/tools/DependencyRules $(VERBOSE_ARG) \
                                                 -database $(DATABASE) \
	                                         -objectdir . \
	                                         -moduledir . \
	                                         $(DEPRULE_FLAGS) $@

IGNORE_ARGUMENTS = $(addprefix -ignore ,$(IGNORE_DEPENDENCIES))
INCLUDE_ARGUMENTS = $(addprefix -include , $(PRE_PROCESS_INCLUDE_DIRS))
MACRO_ARGUMENTS = $(addprefix -macro , $(PRE_PROCESS_MACROS))

%.t: %.f90
	$(call MESSAGE,Analysing,$<)
	$(Q)$(LFRIC_BUILD)/tools/DependencyAnalyser \
	    $(IGNORE_ARGUMENTS) $(VERBOSE_ARG) $(DATABASE) $<
	$(Q)touch $@

%.t: %.F90
	$(call MESSAGE,Analysing,$<)
	$(Q)$(LFRIC_BUILD)/tools/DependencyAnalyser \
	    $(IGNORE_ARGUMENTS) $(INCLUDE_ARGUMENTS) $(MACRO_ARGUMENTS) \
	    $(VERBOSE_ARG) $(DATABASE) $<
	$(Q)touch $@

include $(LFRIC_BUILD)/lfric.mk
include $(LFRIC_BUILD)/fortran.mk
