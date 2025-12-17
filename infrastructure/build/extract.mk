##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
#
# Run this make file to copy a source tree from SOURCE_DIR to WORKING_DIR
#
.PHONY: files-to-extract
files-to-extract: $(addprefix $(WORKING_DIR)/,$(shell find $(SOURCE_DIR) \( -name '*.[Ff]90' -o -name '*.h' \) -print | sed "s|$(SOURCE_DIR)/||")) \
                  | $(WORKING_DIR)

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

.PRECIOUS: $(WORKING_DIR)/%.h
$(WORKING_DIR)/%.h: $(SOURCE_DIR)/%.h | $(WORKING_DIR)
	$(call MESSAGE,Copying source,$<)
	$(Q)mkdir -p $(dir $@)
	$(Q)cp $< $@

$(WORKING_DIR):
	$(call MESSAGE,Creating,$@)
	$(Q)mkdir -p $@

include $(LFRIC_BUILD)/lfric.mk
