##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
.SECONDEXPANSION:

# Unit tests
##############################################################################
UNIT_TEST_EXE = $(BIN_DIR)/$(firstword $(PROGRAMS))

UNIT_TEST_FILTER ?=
UNIT_TEST_FILTER_ARG =
ifdef UNIT_TEST_FILTER
  UNIT_TEST_FILTER_ARG = -f $(UNIT_TEST_FILTER)
endif

ifdef MPI_TESTS
  UNIT_TEST_PRE_PROCESS_MACROS = USE_MPI=YES
  LAUNCHER = mpiexec -n 6
else
  UNIT_TEST_PRE_PROCESS_MACROS = NO_MPI=no_mpi
  # It seems that the Cray 'ftn' wrapper always builds for MPI...
  #
  ifdef CRAY_ENVIRONMENT
    LAUNCHER = mpiexec -n 1
  endif
endif


UNIT_TEST_DATA_DIR = $(if $(wildcard $(TEST_DIR)/data), $(WORKING_DIR)/data)

.PHONY: do-unit-test/%
do-unit-test/run: $(UNIT_TEST_EXE) $(UNIT_TEST_DATA_DIR)
	$(call MESSAGE,Running,$(PROGRAMS))
ifdef UNIT_TEST_FILTER
	$(call MESSAGE,Filter,$(UNIT_TEST_FILTER))
endif
	$Qcd $(WORKING_DIR); \
	    $(LAUNCHER) $(UNIT_TEST_EXE) $(DOUBLE_VERBOSE_ARG) $(UNIT_TEST_FILTER_ARG)

# The addition of this target is a bit messy but it allows us to guarantee that
# no build will happen when running from a test suite.
#
do-unit-test/rerun: $(UNIT_TEST_DATA_DIR)
	$(call MESSAGE,Running,$(PROGRAMS))
	$Qcd $(WORKING_DIR); \
	    $(LAUNCHER) $(UNIT_TEST_EXE) $(DOUBLE_VERBOSE_ARG)

$(WORKING_DIR)/data: $(TEST_DIR)/data | $(WORKING_DIR)
	$(call MESSAGE,Copying test data,$<)
	$Qrsync -avz $(TEST_DIR)/data $(WORKING_DIR)/

$(WORKING_DIR):
	$(call MESSAGE,Creating $@)
	$Qmkdir -p $(WORKING_DIR)  # Ensure the target directory exists.

do-unit-test/build: $(UNIT_TEST_EXE)

$(UNIT_TEST_EXE): export EXTERNAL_STATIC_LIBRARIES += pfunit funit fargparse gftl-shared-v2
$(UNIT_TEST_EXE): export IGNORE_DEPENDENCIES += funit pfunit
$(UNIT_TEST_EXE): export TEST_LIST_FILE = test_list.inc
$(UNIT_TEST_EXE): export PRE_PROCESS_MACROS += $(UNIT_TEST_PRE_PROCESS_MACROS)
$(UNIT_TEST_EXE): export PRE_PROCESS_MACROS += _TEST_SUITES=\"$(TEST_LIST_FILE)\"
$(UNIT_TEST_EXE): do-unit-test/generate $(addsuffix /extract, $(TEST_DIR))
	$Qmkdir -p $(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/pfunit.mk \
	            SOURCE_DIR=$(TEST_DIR)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk

do-unit-test/generate: do-unit-test/get-source \
                       $(if $(META_FILE_DIR), configuration)

do-unit-test/get-source: $(addsuffix /import, $(IMPORT_PARTS)) \
                         $(addsuffix /extract, $(ADDITIONAL_EXTRACTION))
	$Q$(MAKE) -f $(LFRIC_BUILD)/lfric.mk \
	          $(addsuffix /import, $(PROJECT_DIR))

###############################################################################
# Integration tests
###############################################################################

ALL_INTEGRATION_TESTS = $(patsubst $(TEST_DIR)/%,%,$(basename                 \
                            $(shell find $(TEST_DIR) -name '*.[Ff]90'         \
                                         -exec egrep -l "^\s*program\s" {} \; \
                                         2>/dev/null)))
.PHONY: do-integration-tests/%
do-integration-tests/%: export PYTHONPATH    := $(PYTHONPATH):$(LFRIC_BUILD)
do-integration-tests/%: export PROGRAMS       = $(ALL_INTEGRATION_TESTS)

do-integration-tests/run: $(foreach test,$(ALL_INTEGRATION_TESTS),do-integration-tests/run/$(test))

do-integration-tests/run/%: $$(BIN_DIR)/$$*.py \
$$(if $$(realpath $$(TEST_DIR)/$$(dir $$*)/iodef.xml), $$(BIN_DIR)/$$(dir $$*)/iodef.xml) \
| do-integration-tests/build do-integration-tests/resources/$$*
	$(call MESSAGE,Running,$*)
	$Qcd $(dir $<); python3 ./$(notdir $<) $(BIN_DIR)/$*

# The addition of this target is a bit messy but it allows us to guarantee that
# no build will happen when running from a test suite.
do-integration-tests/rerun: $(foreach test,$(ALL_INTEGRATION_TESTS),do-integration-tests/rerun/$(test))

do-integration-tests/rerun/%: $$(BIN_DIR)/$$*.py \
$$(if $$(realpath $$(TEST_DIR)/$$(dir $$*)/iodef.xml), $$(BIN_DIR)/$$(dir $$*)/iodef.xml) \
| do-integration-tests/resources/$$*
	$(call MESSAGE,Rerunning,$*)
	$Qcd $(dir $<); python3 ./$(notdir $<) $(BIN_DIR)/$*

do-integration-tests/resources/%: \
| $$(if $$(realpath $$(TEST_DIR)/support/resources), do-integration-tests/resources/support/$$*)
	$(call MESSAGE,Harvesting, $*)
	$Qif [ -e $(TEST_DIR)/$(dir $*)resources ]; then rsync -a $(TEST_DIR)/$(dir $*)resources $(BIN_DIR)/$(dir $*); fi

do-integration-tests/resources/support/%: $(BIN_DIR)/resources
	$(call MESSAGE,Symlinking to, support/resources)
	$Qmkdir -p $(BIN_DIR)/$(dir $*)
	$Qln -sf $(BIN_DIR)/resources $(BIN_DIR)/$(dir $*)shared-resources

$(BIN_DIR)/resources:
	$(call MESSAGE,Harvesting, support/resources)
	$Qrsync -a $(TEST_DIR)/support/resources $(BIN_DIR)/

$(BIN_DIR)/%.py: $(TEST_DIR)/%.py
	$(call MESSAGE,Copying,$<)
	$Qmkdir -p $(dir $@)
	$Qcp $< $@

$(BIN_DIR)%iodef.xml: $(TEST_DIR)%iodef.xml
	$(call MESSAGE,Copying, $<)
	$Qmkdir -p $(dir $@)
	$Qcp $< $@

do-integration-tests/build: do-integration-tests/generate \
                            $(addsuffix /extract, $(TEST_DIR))
	$Qmkdir -p $(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/lfric.mk \
	          $(addsuffix /psyclone, $(TEST_DIR))
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk

do-integration-tests/generate: do-integration-test/get-source \
                               $(if $(META_FILE_DIR), configuration)

do-integration-test/get-source: $(addsuffix /import, $(IMPORT_PARTS)) \
                                $(addsuffix /extract, $(ADDITIONAL_EXTRACTION))
	$Q$(MAKE) -f $(LFRIC_BUILD)/lfric.mk \
	          $(addsuffix /import, $(PROJECT_DIR))

###############################################################################
# Utilities
###############################################################################

include $(LFRIC_BUILD)/lfric.mk
