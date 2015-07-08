##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
#
# Python's support for SQLite seems to hobble its concurrency support so this
# stuff has to be done serially.

.NOTPARALLEL:

ROOT = ../..
include $(ROOT)/make/include.mk

DYNAMO_DEPENDENCY_DB ?= $(OBJ_DIR)/dependencies.db

OBJ_SUBDIRS = $(patsubst %,$(OBJ_DIR)/%,$(shell find * -type d))
ALL_SRC = $(shell find . -name "*.[Ff]90")
EXISTING_TOUCH = $(shell find $(OBJ_DIR) -name "*.t")
TOUCH_FILES = $(patsubst ./%.f90,$(OBJ_DIR)/%.t,$(patsubst ./%.F90,$(OBJ_DIR)/%.t,$(ALL_SRC)))
STALE_TOUCH = $(filter-out $(TOUCH_FILES), $(EXISTING_TOUCH) )
IGNORE_ARGUMENTS := $(patsubst %,-ignore %,$(IGNORE_DEPENDENCIES))

$(OBJ_DIR)/programs.mk: $(OBJ_DIR)/dependencies.mk | $(OBJ_DIR)
	@echo -e $(VT_BOLD)Building$(VT_RESET) $@
	$(Q)$(TOOL_DIR)/ProgramObjects -database $(DYNAMO_DEPENDENCY_DB) $@

$(OBJ_DIR)/dependencies.mk: $(TOUCH_FILES) $(STALE_TOUCH) | $(OBJ_DIR) $(OBJ_SUBDIRS)
	@echo -e $(VT_BOLD)Building$(VT_RESET) $@
	$(Q)$(TOOL_DIR)/DependencyRules -database $(DYNAMO_DEPENDENCY_DB) $@

.PHONY: unused-check
unused-check: $(OBJ_DIR)/dependencies.mk
	@echo -e $(VT_BOLD)Checking for unused source$(VT_RESET)
	$(Q)$(TOOL_DIR)/CheckForUnused -database $(DYNAMO_DEPENDENCY_DB) $(ALL_SRC)

$(OBJ_DIR)/%.t: %.F90 | $(OBJ_DIR) $(OBJ_SUBDIRS)
	@echo -e $(VT_BOLD)Analysing$(VT_RESET) $<
	$(Q)$(TOOL_DIR)/DependencyAnalyser $(IGNORE_ARGUMENTS) \
	                                   $(DYNAMO_DEPENDENCY_DB) $< && touch $@

$(OBJ_DIR)/%.t: %.f90 | $(OBJ_DIR) $(OBJ_SUBDIRS)
	@echo -e $(VT_BOLD)Analysing$(VT_RESET) $<
	$(Q)$(TOOL_DIR)/DependencyAnalyser $(IGNORE_ARGUMENTS) \
	                                   $(DYNAMO_DEPENDENCY_DB) $< && touch $@

$(OBJ_DIR) $(OBJ_SUBDIRS):
	@echo -e $(VT_BOLD)Creating$(VT_RESET) $@
	$(Q)mkdir -p $@

$(STALE_TOUCH):
	@echo -e $(VT_BOLD)Removing$(VT_RESET) $<
	$(Q)rm $<
