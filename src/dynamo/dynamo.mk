##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

AR ?= ar

ROOT     = ../..
BIN_DIR  = $(ROOT)/bin
DATABASE = $(OBJ_DIR)/dependencies.db

include $(ROOT)/make/include.mk
include $(OBJ_DIR)/programs.mk $(OBJ_DIR)/dependencies.mk

# Extract the program names from the program objects:
PROGRAMS = $(patsubst %.o,%,$(notdir $(PROG_OBJS)))

# Convert the program names to all caps, append "_OBJS" to the end and
# dereference that to get a list of all objects needed by all programs:
ALL_OBJS = $(foreach proj, $(shell echo $(PROGRAMS) | tr a-z A-Z), $($(proj)_OBJS))

# Remove the program objects from the list of all objects, leaving modules
ALL_MODULES = $(filter-out $(PROG_OBJ),$(ALL_OBJS))

.SECONDEXPANSION:

.PHONY: applications
applications: $(patsubst %,$(BIN_DIR)/%,$(PROGRAMS))

.PHONY: modules
modules: $(OBJ_DIR)/modules.a | $(OBJ_DIR)

$(BIN_DIR)/%: $(OBJ_DIR)/%.x | $(BIN_DIR)
	@echo -e $(VT_BOLD)Installing\$(VT_RESET) $@
	$(Q)cp $(OBJ_DIR)/$(notdir $@).x $@

# Directories

# Find all the subdirectories within the source directory:
SUBDIRS = $(shell find * -type d -prune)

# Mirror the source directory tree in the object directory:
OBJ_SUBDIRS = $(patsubst %,$(OBJ_DIR)/%/,$(SUBDIRS))

$(BIN_DIR) $(OBJ_DIR):
	@echo -e $(VT_BOLD)Creating$(VT_RESET) $@
	$(Q)mkdir -p $@

# Build Rules

# Build a set of "-I" arguments to seach the whole object tree:
INCLUDE_ARGS = -I$(OBJ_DIR) $(patsubst %, -I$(OBJ_DIR)/%, $(SUBDIRS))

$(OBJ_DIR)/%.o $(OBJ_DIR)/%.mod: %.F90 | $(OBJ_DIR)
	@echo -e \$(VT_BOLD)Compile$(VT_RESET) $<
	$(Q)$(FC) $(FPPFLAGS) $(FFLAGS) \
	          $(F_MOD_DESTINATION_ARG)$(OBJ_DIR)/$(dir $@) \
	          $(INCLUDE_ARGS) -c -o $(basename $@).o $<

$(OBJ_DIR)/%.o $(OBJ_DIR)/%.mod: %.f90 | $(OBJ_DIR)
	@echo -e \$(VT_BOLD)Compile$(VT_RESET) $<
	$(Q)$(FC) $(FPPFLAGS) $(FFLAGS) \
	          $(F_MOD_DESTINATION_ARG)$(OBJ_DIR)/$(dir $@) \
	          $(INCLUDE_ARGS) -c -o $(basename $@).o $<

$(OBJ_DIR)/modules.a: $(ALL_MODULES)
	$(Q)$(AR) -r $@ $^

$(OBJ_DIR)/%.x: $$($$(shell echo $$* | tr a-z A-Z)_OBJS)
	@echo -e \$(VT_BOLD)Linking$(VT_RESET) $@
	$(Q)$(MPI_LD) $(LDFLAGS) -o $@ \
	              $(patsubst %,-l%,$(EXTERNAL_DYNAMIC_LIBRARIES)) \
	              $^ \
	              $(patsubst %,-l%,$(EXTERNAL_STATIC_LIBRARIES))

# Dependencies

# It is important that the two dependency files are build sequentially
# otherwise they will fight over the dependency database.
#
$(OBJ_DIR)/programs.mk: $(OBJ_DIR)/dependencies.mk | $(OBJ_DIR)
	$(MAKE) -f examine.mk $(OBJ_DIR)/programs.mk OBJ_SUBDIRS="$(OBJ_SUBDIRS)"

$(OBJ_DIR)/dependencies.mk: ALWAYS | $(OBJ_DIR)
	$(MAKE) -f examine.mk $(OBJ_DIR)/dependencies.mk OBJ_SUBDIRS="$(OBJ_SUBDIRS)"

.PHONY: ALWAYS
ALWAYS:

# Special Rules

.PHONY: clean
clean:
	-rm -rf $(OBJ_DIR)
	-rm -f $(BIN_DIR)/*
