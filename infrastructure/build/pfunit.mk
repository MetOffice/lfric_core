##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
#
# Run this make file to generate pFUnit source in WORKING_DIR from test
# descriptions in SOURCE_DIR.
#
PF_FILES = $(shell find $(SOURCE_DIR) -name '*.pf' -print | sed "s|$(SOURCE_DIR)||")

.PHONY: prepare-pfunit
prepare-pfunit: $(patsubst %.pf,$(WORKING_DIR)/%.F90,$(PF_FILES)) \
        $(WORKING_DIR)/$(PROJECT)_unit_tests.F90
	$(Q)echo >/dev/null

include $(LFRIC_BUILD)/lfric.mk

.PRECIOUS: $(WORKING_DIR)/$(PROJECT)_unit_tests.F90
$(WORKING_DIR)/$(PROJECT)_unit_tests.F90: $(PFUNIT)/include/driver.F90 \
                                         $(WORKING_DIR)/testSuites.inc
	$(call MESSAGE,Processing, "pFUnit driver source")
	$(Q)sed "s/program main/program $(basename $(notdir $@))/" <$< >$@

.PRECIOUS: $(WORKING_DIR)/testSuites.inc
$(WORKING_DIR)/testSuites.inc:
	$(call MESSAGE,Collating, $@)
	$(Q)mkdir -p $(dir $@)
	$(Q)echo ! Tests to run >$@
	$(Q)for test in $(basename $(notdir $(shell find $(SOURCE_DIR) -name '*.pf'))); do echo ADD_TEST_SUITE\($${test}_suite\) >> $@; done

$(WORKING_DIR)/%.F90: $(SOURCE_DIR)/%.pf
	$(call MESSAGE,Generating unit test,$@)
	$(Q)mkdir -p $(dir $@)
	$(Q)$(PFUNIT)/bin/pFUnitParser.py $< $@ $(VERBOSE_REDIRECT)
