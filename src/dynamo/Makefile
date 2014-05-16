##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

IGNORE_DEPENDENCIES = netcdf
EXTERNAL_LIBRARIES = netcdff

ROOT = ../..
include $(ROOT)/Makefile.inc

SUBDIRS = $(shell find * -maxdepth 1 -type d)
IGNORE_DEPENDENCIES := $(patsubst %,-ignore %,$(IGNORE_DEPENDENCIES))

export SRC        = $(wildcard *.F90)
export OBJ_DIR    = $(BUILD_DIR)/dynamo
export DEP_FILE   = dependencies.mk

export LDFLAGS += $(patsubst %, -l%,$(EXTERNAL_LIBRARIES))

OBJ_SUBDIRS = $(patsubst %,$(OBJ_DIR)/%,$(SUBDIRS))

.PHONY: all
all: $(OBJ_DIR)/$(DEP_FILE)
	$(MAKE) -f dynamo.mk all

$(OBJ_DIR)/$(DEP_FILE): $(SRC) $(OBJ_DIR) $(OBJ_SUBDIRS)
	@$(TOOL_DIR)/DependencyAnalyser $(IGNORE_DEPENDENCIES) \
	                                $(OBJ_DIR) $(OBJ_DIR) $(OBJ_DIR) .

$(OBJ_DIR) $(OBJ_SUBDIRS):
	@echo "Creating $@"
	@mkdir -p $@

.PHONY: clean
clean:
	$(MAKE) -f dynamo.mk clean
