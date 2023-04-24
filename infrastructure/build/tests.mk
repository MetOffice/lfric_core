##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Unit tests
##############################################################################
UNIT_TEST_EXE = $(BIN_DIR)/$(firstword $(PROGRAMS))

ifdef MPI_TESTS
  UNIT_TEST_PRE_PROCESS_MACROS = USE_MPI=YES
  LAUNCHER = mpiexec -n 4
else
  UNIT_TEST_PRE_PROCESS_MACROS = NO_MPI=no_mpi
  # It seems that the Cray 'ftn' wrapper always builds for MPI...
  #
  ifdef CRAY_ENVIRONMENT
    LAUNCHER = mpiexec -n 1
  endif
endif


.PHONY: do-unit-test/%
do-unit-test/run: $(UNIT_TEST_EXE)
	$(call MESSAGE,Running,$(PROGRAMS))
	$Qcd $(TEST_DIR); \
	    $(LAUNCHER) $(UNIT_TEST_EXE) $(DOUBLE_VERBOSE_ARG)

# The addition of this target is a bit messy but it allows us to guarantee that
# no build will happen when running from a test suite.
#
do-unit-test/rerun:
	$(call MESSAGE,Running,$(PROGRAMS))
	$Qcd $(TEST_DIR); \
	    $(LAUNCHER) $(UNIT_TEST_EXE) $(DOUBLE_VERBOSE_ARG)


do-unit-test/build: $(UNIT_TEST_EXE)

$(UNIT_TEST_EXE): export EXTERNAL_STATIC_LIBRARIES += pfunit
$(UNIT_TEST_EXE): do-unit-test/generate \
                  $(addsuffix /extract, $(TEST_DIR))
	$Qmkdir -p $(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/pfunit.mk \
	            SOURCE_DIR=$(TEST_DIR) WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$Q$(MAKE) $(QUIET_ARG) \
                    -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk \
                    PRE_PROCESS_MACROS="$(PRE_PROCESS_MACROS) $(UNIT_TEST_PRE_PROCESS_MACROS)"

# Ensure all extraction is performed before PSyclone otherwise kernel files may
# not have arrived when they are needed.
#
do-unit-test/generate: do-unit-test/get-source \
                       $(if $(META_FILE_DIR), configuration)
	$Q$(MAKE) -f $(LFRIC_BUILD)/lfric.mk           \
	          $(addsuffix /psyclone, $(SOURCE_DIR) \
	                                 $(ADDITIONAL_EXTRACTION))

do-unit-test/get-source: $(addsuffix /extract, $(SOURCE_DIR) \
                                               $(ADDITIONAL_EXTRACTION))

###############################################################################
# Integration tests
###############################################################################

ALL_INTEGRATION_TESTS = $(patsubst $(TEST_DIR)/%,%,$(basename                 \
                            $(shell find $(TEST_DIR) -name '*.[Ff]90'         \
                                         -exec egrep -l "^\s*program\s" {} \; \
                                         2>/dev/null)))
.PHONY: do-integration-tests/%
do-integration-tests/%: export PYTHONPATH  := $(PYTHONPATH):$(LFRIC_BUILD)
do-integration-tests/%: export PROGRAMS     = $(ALL_INTEGRATION_TESTS)
do-integration-tests/%: export TEST_RUN_DIR = $(BIN_DIR)/test_files

do-integration-tests/run: $(foreach test,$(ALL_INTEGRATION_TESTS),do-integration-tests/run/$(test))

do-integration-tests/run/%: do-integration-tests/build
	$(call MESSAGE,Running,$*)
	$Qcd $(TEST_RUN_DIR)/$(dir $*); \
	    ./$(notdir $(addsuffix .py,$*)) $(addprefix $(BIN_DIR)/,$*)

# The addition of this target is a bit messy but it allows us to guarantee that
# no build will happen when running from a test suite.
do-integration-tests/rerun: $(foreach test,$(ALL_INTEGRATION_TESTS),do-integration-tests/rerun/$(test))

do-integration-tests/rerun/%:
	$(call MESSAGE,Rerunning,$*)
	$Qcd $(TEST_RUN_DIR)/$(dir $*); \
	    ./$(notdir $(addsuffix .py,$*)) $(addprefix $(BIN_DIR)/,$*)

do-integration-tests/build: do-integration-tests/generate \
                           $(addsuffix /extract, $(TEST_DIR))
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk
	$Qrsync -a $(TEST_DIR)/ $(TEST_RUN_DIR)

# Extraction is performed before psyclone to ensure kernel files have arrived
# before they are needed.
#
do-integration-tests/generate: do-integration-test/get-source \
                               $(if $(META_FILE_DIR), configuration)
	$Q$(MAKE) -f $(LFRIC_BUILD)/lfric.mk           \
	          $(addsuffix /psyclone, $(SOURCE_DIR) \
	                                 $(TEST_DIR)   \
	                                 $(ADDITIONAL_EXTRACTION))

do-integration-test/get-source: $(addsuffix /extract, $(SOURCE_DIR) \
                                                      $(ADDITIONAL_EXTRACTION))

###############################################################################
# Utilities
###############################################################################

include $(LFRIC_BUILD)/lfric.mk
