##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
#
# Run this make file to copy a source tree from SOURCE_DIR to WORKING_DIR
#
# Environment variables which affect this script:
#
# SOURCE_DIR: Where to look for files.
# WORKING_DIR: Where intermediate files end up.
# PRE_PROCESS_INCLUDE_DIRS: Space separated list of directories to search for
#                           inclusions.
# PRE_PROCESS_MACROS: Space separated list of macro definitions in the form
#                     NAME[=MACRO] to be passed to the compiler.
#
candidate_files := $(patsubst $(SOURCE_DIR)/%,%,$(shell find $(SOURCE_DIR) \( -name '*.[Ff]90' -o -name '*.h' \)))
fortran_files = $(filter %.f90,$(candidate_files)) $(patsubst %.F90,%.f90,$(filter %.F90, $(candidate_files)))
header_files = $(addprefix $(WORKING)/,$(filter %.h, $(candidate_files)))
.PHONY: files-to-extract
files-to-extract: $(addprefix $(WORKING_DIR)/,$(fortran_files)) \
                  $(addprefix $(WORKING_DIR)/,$(header_files))  \
                  | $(WORKING_DIR)

.PRECIOUS: $(WORKING_DIR)/%.f90
$(WORKING_DIR)/%.f90: $(SOURCE_DIR)/%.f90 | $(WORKING_DIR)
	$(call MESSAGE,Copying source,$<)
	$(Q)mkdir -p $(dir $@)
	$(Q)cp $< $@

.PRECIOUS: $(WORKING_DIR)/%.h
$(WORKING_DIR)/%.h: $(SOURCE_DIR)/%.h | $(WORKING_DIR)
	$(call MESSAGE,Copying source,$<)
	$(Q)mkdir -p $(dir $@)
	$(Q)cp $< $@

$(WORKING_DIR):
	$(call MESSAGE,Creating,$@)
	$(Q)mkdir -p $@

include $(LFRIC_BUILD)/preprocess.mk
include $(LFRIC_BUILD)/lfric.mk
include $(LFRIC_BUILD)/fortran.mk
-include $(COMPILE_OPTIONS)
