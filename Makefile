##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# LFRic project make file. This simply calls down into sub-project make files.
#
# Variables:
#   OPERATE_ON - Sub-projects which will be affected by operations.
#
#############################################################################

# Operate only on this list of sub-projects. May be overridden from the
# terminal.
#
OPERATE_ON ?= infrastructure                \
              components/science            \
              components/driver             \
              components/lfric-xios         \
              components/inventory          \
              components/coupling           \
              mesh_tools                    \
              applications/skeleton         \
              applications/simple_diffusion \
              applications/lbc_demo         \
              applications/io_demo

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

include infrastructure/build/lfric.mk
