##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Preprocess Fortran source. Either from a source directory into the working
# directory or within the working directory.
#
# This file is intended to be "include"ed into other make files.
#
# It makes use of the following variables:
#
# FPP: Fortran preprocessor tool.
# FPPFLAGS: arguments passed to all FPP invocations.
# PRE_PROCESS_INCLUDE_DIRS: Space separated list of directories to search for
#                           inclusions.
# PRE_PROCESS_MACROS: Space separated list of macro definitions in the form
#                     NAME[=MACRO] to be passed to the compiler.
# SOURCE_DIR: Where the source code is coming from.
# WORKING_DIR: Scratch space for intermediate files.
#
###############################################################################

# Build a set of "-I" arguments to seach the whole object tree:
include_args := $(subst ./,-I,$(shell find . -mindepth 1 -type d -print)) \
                $(addprefix -I, $(PRE_PROCESS_INCLUDE_DIRS))

# Build a set of "-D" argument for any pre-processor macros
#
macro_args := $(addprefix -D,$(PRE_PROCESS_MACROS))

.PRECIOUS: $(WORKING_DIR)/%.f90
$(WORKING_DIR)/%.f90: $(SOURCE_DIR)/%.F90 | $(WORKING_DIR)
	$(call MESSAGE,Preprocessing source, $<)
	$Qmkdir -p $(dir $@)
	$Q$(FPP) $(FPPFLAGS) $(include_args) $(macro_args) $< $@

$(WORKING_DIR)/%.f90: $(WORKING_DIR)/%.F90
	$(call MESSAGE,Preprocessing intermediate, $<)
	$Q$(FPP) $(FPPFLAGS) $(include_args) $(macro_args) $< $@
