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
OPERATE_ON ?= lfric_atm                  \
              miniapps/lfric_coupled     \
              gungho                     \
              infrastructure             \
              components/science         \
              components/driver          \
              components/lfric-xios      \
              components/jedi-lfric      \
              components/coupler-oasis   \
              components/inventory       \
              mesh_tools                 \
              linear                     \
              miniapps/gungho_model      \
              miniapps/linear_model      \
              miniapps/skeleton          \
              miniapps/jedi_lfric_tests  \
              miniapps/gravity_wave      \
              miniapps/simple_diffusion  \
              miniapps/solver_miniapp    \
              miniapps/shallow_water     \
              miniapps/io_dev            \
              miniapps/transport         \
              miniapps/multires_coupling \
              um_physics                 \
              socrates                   \
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
