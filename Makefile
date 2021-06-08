##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# LFRic project make file. This simply calls down into sub-project make files.
#
# Variables:
#   OPERATE_ON - Sub-projects which will be affected by operations.
#                default: infrastructure, mesh_tools and gungho
#   TEST_SUITE_TARGETS - Platforms to target with test suite.
#
#############################################################################

# Operate only on this list of sub-projects. May be overridden from the
# terminal.
#
OPERATE_ON ?= infrastructure                             \
              components/lfric-xios                      \
              mesh_tools gungho lfric_atm                \
              miniapps/skeleton                          \
              miniapps/diagnostics                       \
              miniapps/gravity_wave                      \
              miniapps/solver_miniapp                    \
              miniapps/io_dev                            \
              miniapps/transport                         \
              lfricinputs

export SUITE_GROUP ?= developer
export SUITE_GROUP_NAME ?= $(notdir $(realpath $(shell pwd)))-.*

##############################################################################
# Perform default action on each sub-project in OPERATE_ON list.
#
.PHONY: default
default: $(addprefix default/,$(OPERATE_ON))
	$(Q)echo >/dev/null

.PHONY: default/%
default/%: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -C $*

##############################################################################
# Perform the clean action on each sub-project in OPERATE_ON list.
#
.PHONY: clean
clean: $(addprefix clean/,$(OPERATE_ON))
	$(Q)echo >/dev/null

.PHONY: clean/%
clean/%: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -C $* clean

##############################################################################
# Lauanch gscan to monitor suites from this suite-group
#
.PHONY: launch-suite-gscan
launch-suite-gscan: gscan_processes := $(shell ps --no-headers -o command -C cylc-gscan)
launch-suite-gscan: ALWAYS
	$(Q)-if [[ "$(gscan_processes)" != *"--name=$(SUITE_GROUP_NAME)"* ]]; then \
          cylc gscan $(DOUBLE_VERBOSE_ARG) --name=$(SUITE_GROUP_NAME) &            \
          usleep 1                                                                ;\
        fi


##############################################################################
# Launch test suite for each sub-project in OPERATE_ON list.
#
.PHONY: test-suite
test-suite: launch-suite-gscan $(addprefix test-suite/,$(OPERATE_ON))
	$(Q)echo >/dev/null

.PHONY: test-suite/%
test-suite/%: ALWAYS
	$(Q)-$(MAKE) $(QUIET_ARG) -C $* test-suite TEST_SUITE_TARGETS="$(TEST_SUITE_TARGETS)"




##############################################################################

include infrastructure/build/lfric.mk
