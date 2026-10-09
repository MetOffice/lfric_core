##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
#
# Run this make file to generate pFUnit source in WORKING_DIR from test
# descriptions in SOURCE_DIR.
#
# Since unit test source is preprocessed, will also need the preprocessor
# variables PRE_PROCESS_INCLUDE_DIRS and PRE_PROCESS_MACROS ang with FPP and
# FPPFLAGS.
#
PF_FILES = $(shell find $(SOURCE_DIR) -name '*.pf' -print | sed "s|$(SOURCE_DIR)||")
PPF_FILES = $(shell find $(SOURCE_DIR) -name '*.PF' -print | sed "s|$(SOURCE_DIR)||")

scratch_dir = $(abspath $(WORKING_DIR)/../pfunit)

.PHONY: prepare-pfunit
prepare-pfunit: $(patsubst %.pf,$(scratch_dir)/%.F90,$(PF_FILES)) \
                $(patsubst %.PF,$(scratch_dir)/%.F90,$(PPF_FILES)) \
                $(scratch_dir)/$(PROJECT)_unit_tests.F90
	$Q$(MAKE) -f $(LFRIC_BUILD)/lfric.mk \
	          SOURCE_DIR=$(scratch_dir) \
	          $(addsuffix /extract, $(scratch_dir))

include $(LFRIC_BUILD)/lfric.mk

.PRECIOUS: $(scratch_dir)/$(PROJECT)_unit_tests.F90
$(scratch_dir)/$(PROJECT)_unit_tests.F90: $(PFUNIT)/include/driver.F90 \
                                          $(scratch_dir)/$(TEST_LIST_FILE)
	$(call MESSAGE,Processing, "pFUnit driver source")
	$(Q)sed -e "s/program main/program $(basename $(notdir $@))/" \
	        <$< >$@

.PRECIOUS: $(scratch_dir)/$(TEST_LIST_FILE)
$(scratch_dir)/$(TEST_LIST_FILE): | $(scratch_dir)
	$(call MESSAGE,Collating, $@)
	$(Q)mkdir -p $(dir $@)
	$(Q)echo ! Tests to run >$@
	$Qfor test in $(basename $(notdir $(shell find $(SOURCE_DIR) -iname '*.pf'))); \
	      do echo ADD_TEST_SUITE\($${test}_suite\) >> $@; done

$(scratch_dir)/%.F90: $(SOURCE_DIR)/%.pf
	$(call MESSAGE,Generating unit test,$@)
	$(Q)mkdir -p $(dir $@)
	$(Q)$(PFUNIT)/bin/funitproc $(QUIET_ARG_SINGLE) $< $@

$(scratch_dir)/%.F90: $(WORKING_DIR)/%.pf
	$(call MESSAGE,Generating unit test, $@)
	$Qmkdir -p $(dir $@)
	$Q$(PFUNIT)/bin/funitproc $(QUIET_ARG_SINGLE) $< $@

$(scratch_dir)/%.pf: $(SOURCE_DIR)/%.PF
	$(call MESSAGE,Preprocessing unit test, $<)
	$Qmkdir -p $(dir $@)
	$Q$(FPP) $(addprefix -I, $(PRE_PROCESS_INCLUDE_DIRS)) \
	         $(addprefix -D, $(PRE_PROCESS_MACROS)) \
	         <$< >$@

$(scratch_dir):
	$(call MESSAGE,Creating,$@)
	$Qmkdir -p $@
