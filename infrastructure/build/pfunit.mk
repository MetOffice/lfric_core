##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
#
# Run this make file to generate pFUnit source in WORKING_DIR from test
# descriptions in SOURCE_DIR.
#
PF_FILES = $(shell find $(SOURCE_DIR) -name '*.pf' -print | sed "s|$(SOURCE_DIR)||")
PPF_FILES = $(shell find $(SOURCE_DIR) -name '*.PF' -print | sed "s|$(SOURCE_DIR)||")

.PHONY: prepare-pfunit
prepare-pfunit: $(patsubst %.pf,$(WORKING_DIR)/%.F90,$(PF_FILES)) \
                $(patsubst %.PF,$(WORKING_DIR)/%.F90,$(PPF_FILES)) \
                $(WORKING_DIR)/$(PROJECT)_unit_tests.F90

include $(LFRIC_BUILD)/lfric.mk

.PRECIOUS: $(WORKING_DIR)/$(PROJECT)_unit_tests.F90
$(WORKING_DIR)/$(PROJECT)_unit_tests.F90: $(PFUNIT)/include/driver.F90 \
                                          $(WORKING_DIR)/$(TEST_LIST_FILE)
	$(call MESSAGE,Processing, "pFUnit driver source")
	$(Q)sed -e "s/program main/program $(basename $(notdir $@))/" \
	        <$< >$@

.PRECIOUS: $(_TEST_SUITES)
$(WORKING_DIR)/$(TEST_LIST_FILE):
	$(call MESSAGE,Collating, $@)
	$(Q)mkdir -p $(dir $@)
	$(Q)echo ! Tests to run >$@
	$Qfor test in $(basename $(notdir $(shell find $(SOURCE_DIR) -iname '*.pf'))); \
	      do echo ADD_TEST_SUITE\($${test}_suite\) >> $@; done

$(WORKING_DIR)/%.F90: $(SOURCE_DIR)/%.pf
	$(call MESSAGE,Generating unit test,$@)
	$(Q)mkdir -p $(dir $@)
	$(Q)$(PFUNIT)/bin/funitproc $(QUIET_ARG_SINGLE) $< $@

$(WORKING_DIR)/%.F90: $(WORKING_DIR)/%.pf
	$(call MESSAGE, Generating unit test, $@)
	$Qmkdir -p $(dir $@)
	$Q$(PFUNIT)/bin/funitproc $(QUIET_ARG_SINGLE) $< $@

$(WORKING_DIR)/%.pf: $(SOURCE_DIR)/%.PF
	$(call MESSAGE, Preprocessing unit test, $<)
	$Qmkdir -p $(dir $@)
	$Q$(FPP) $(addprefix -I, $(PRE_PROCESS_INCLUDE_DIRS)) \
	         $(addprefix -D, $(PRE_PROCESS_MACROS)) \
	         <$< >$@
