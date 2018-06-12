##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
#
# Run this make file to copy a source tree from SOURCE_DIR to WORKING_DIR
#
# Set WITHOUT_PROGRAMS in order to prevent files containing program program
# units from being extracted.
#
.PHONY: files-to-extract
NLD_FILES = $(shell find $(SOURCE_DIR) -path "$(SOURCE_DIR)/*/*" -name '*.nld' -printf '%P\n')
ifdef WITHOUT_PROGRAMS
MODULE_FILES = $(shell find $(SOURCE_DIR) -path "$(SOURCE_DIR)/*/*" -name '*.[Ff]90' -printf '%P\n')
files-to-extract: $(addprefix $(WORKING_DIR)/,$(MODULE_FILES)) \
                  $(addprefix $(WORKING_DIR)/,$(NLD_FILES)) | $(WORKING_DIR)
	$(Q)echo >/dev/null
else
files-to-extract: $(addprefix $(WORKING_DIR)/,$(shell find $(SOURCE_DIR) -name '*.[Ff]90' -printf '%P\n')) \
                  $(addprefix $(WORKING_DIR)/,$(NLD_FILES)) | $(WORKING_DIR)
	$(Q)echo >/dev/null
endif

.PRECIOUS: $(WORKING_DIR)/%.F90
$(WORKING_DIR)/%.F90: $(SOURCE_DIR)/%.F90 | $(WORKING_DIR)
	$(call MESSAGE,Copying source,$<)
	$(Q)mkdir -p $(dir $@)
	$(Q)cp $< $@

.PRECIOUS: $(WORKING_DIR)/%.f90
$(WORKING_DIR)/%.f90: $(SOURCE_DIR)/%.f90 | $(WORKING_DIR)
	$(call MESSAGE,Copying source,$<)
	$(Q)mkdir -p $(dir $@)
	$(Q)cp $< $@

.PRECIOUS: $(WORKING_DIR)/%.nld
$(WORKING_DIR)/%.nld: $(SOURCE_DIR)/%.nld | $(WORKING_DIR)
	$(call MESSAGE,Copying source,$<)
	$(Q)mkdir -p $(dir $@)
	$(Q)cp $< $@

$(WORKING_DIR):
	$(call MESSAGE,Creating,$@)
	$(Q)mkdir -p $@

include $(LFRIC_BUILD)/lfric.mk
