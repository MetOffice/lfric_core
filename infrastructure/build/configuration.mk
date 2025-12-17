##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
#
# Run this make file to generate configuration found in SOURCE_DIR
# to WORKING_DIR. Uses PROJECT to know what to call master files.
#

export CONFIG_DIR=$(WORKING_DIR)/configuration

.PHONY: configuration_files
configuration_files: $(WORKING_DIR)/configuration_mod.f90 \
                     $(WORKING_DIR)/feign_config_mod.f90

.INTERMEDIATE: $(CONFIG_DIR)/rose-meta.json $(CONFIG_DIR)/config_namelists.txt
$(CONFIG_DIR)/rose-meta.json $(CONFIG_DIR)/config_namelists.txt: $(META_FILE_DIR)/rose-meta.conf
	$(call MESSAGE,Generating namelist configuration file.)
	$(Q)mkdir -p $(dir $@)
ifdef APPS_ROOT_DIR
	$(Q)rose_picker $(META_FILE_DIR)/rose-meta.conf    \
	                -directory $(CONFIG_DIR)           \
	                -include_dirs $(APPS_ROOT_DIR)/rose-meta     \
	                -include_dirs $(CORE_ROOT_DIR)/rose-meta
else
	$(Q)rose_picker $(META_FILE_DIR)/rose-meta.conf    \
	                -directory $(CONFIG_DIR)           \
	                -include_dirs $(CORE_ROOT_DIR)/rose-meta
endif
	# It's not clear why this is needed but as of 5/2/20 the diagnostic
	# application test suite fails without it.
	$(Q)sleep 20

.INTERMEDIATE: $(CONFIG_DIR)/build_config_loaders
$(CONFIG_DIR)/build_config_loaders: $(CONFIG_DIR)/rose-meta.json
	$(call MESSAGE,Generating namelist loading modules.)
	$(Q)$(LFRIC_BUILD)/tools/GenerateNamelist $(VERBOSE_ARG) \
                           $(CONFIG_DIR)/rose-meta.json          \
                           -directory $(CONFIG_DIR)
	$(Q)touch $(CONFIG_DIR)/build_config_loaders

# This recipe requires config_namelists.txt, although adding it to the dependencies
# causes a race condition when calling Make in parallel. The generation
# of config_namelists.txt is done at the same time as rose-meta.json, so the
# presense of config_namelists.txt is implied as true if rose-meta.json is present
.PRECIOUS: $(WORKING_DIR)/configuration_mod.f90 $(CONFIG_DIR)/%_config_mod.f90
$(WORKING_DIR)/configuration_mod.f90: $(CONFIG_DIR)/build_config_loaders
	$(call MESSAGE,Generating configuration loader module,$(notdir $@))
	$(Q)mkdir -p $(dir $@)
	$(Q)$(LFRIC_BUILD)/tools/GenerateLoader $(VERBOSE_ARG) $@ $(shell cat $(CONFIG_DIR)/config_namelists.txt)


.PRECIOUS: $(WORKING_DIR)/feign_config_mod.f90
$(WORKING_DIR)/feign_config_mod.f90: $(CONFIG_DIR)/rose-meta.json
	$(call MESSAGE,Generating namelist feigning module.)
	$(Q)mkdir -p $(dir $@)
	$(Q)$(LFRIC_BUILD)/tools/GenerateFeigns        \
                           $(CONFIG_DIR)/rose-meta.json \
                           -output $@


include $(LFRIC_BUILD)/lfric.mk
