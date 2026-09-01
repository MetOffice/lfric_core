##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
#
# Run this make file to generate PSyKAl source in WORKING_DIR from algorithms
# and kernels in SOURCE_DIR. Transformation scripts are sought in
# OPTIMISATION_PATH.
#
# Set the DSL Method in use to collect the correct transformation files.
DSL = psykal
#

# Set default psyclone command additional options
PSYCLONE_PSYKAL_EXTRAS ?= -l all
#
# The command used to invoke PSyclone. By default this is the persistent
# server client (psyclone_client.py) which keeps a single PSyclone instance
# resident and dispatches jobs to a pool of pre-forked workers, avoiding the
# repeated cost of loading the Python libraries from disk. It is a drop-in
# replacement for the "psyclone" binary and falls back to it automatically if
# the server is unavailable. Override PSYCLONE=psyclone to bypass the server.
#
# The server's lifetime is pinned to the owning make process: lfric.mk exports
# PSYCLONE_OWNER_PID (the pid of the top-level make) and the server exits as
# soon as that process does, so no server survives the build that started it.
# Should that variable be unset the client determines the owner itself, by
# finding the outermost make process in its own ancestry.
PSYCLONE ?= $(LFRIC_BUILD)/psyclone/psyclone_client.py
#
# Number of pre-forked PSyclone worker processes. Sized to the build
# parallelism where known (MAKE_THREADS), otherwise to the number of
# available processors.
export PSYCLONE_WORKERS ?= $(if $(MAKE_THREADS),$(MAKE_THREADS),$(shell nproc))
#

ALGORITHM_F_FILES := $(patsubst $(SOURCE_DIR)/%.X90, \
                                $(WORKING_DIR)/%.f90, \
                                $(shell find $(SOURCE_DIR) -name '*.X90' -print))

ALGORITHM_f_FILES := $(patsubst $(SOURCE_DIR)/%.x90, \
                                $(WORKING_DIR)/%.f90, \
                                $(shell find $(SOURCE_DIR) -name '*.x90' -print))

DIRECTORIES := $(patsubst $(SOURCE_DIR)%,$(WORKING_DIR)%, \
                          $(shell find $(SOURCE_DIR) -type d -printf '%p/\n'))
PSYCLONE_CONFIG_FILE ?= $(CORE_ROOT_DIR)/etc/psyclone.cfg

.PHONY: psyclone
psyclone: $(ALGORITHM_F_FILES) $(ALGORITHM_f_FILES)

include $(LFRIC_BUILD)/lfric.mk
include $(LFRIC_BUILD)/fortran.mk

MACRO_ARGS := $(addprefix -D,$(PRE_PROCESS_MACROS))

# Where an override file exists in the "psy" directory we invoke PSyclone, then
# delete the resulting PSy source. The override has been copied as part of the
# rest of the source.
#
$(WORKING_DIR)/%.f90: \
$$(SOURCE_DIR)/psy/$$(notdir $$*)_psy.f90 $(WORKING_DIR)/%_psy.f90
	$(call MESSAGE,Removing,$*_psy.f90)
	$Qrm $(WORKING_DIR)/$*_psy.f90

# Where an optimisation script exists for a specific file, use it.
#
$(WORKING_DIR)/%.f90 $(WORKING_DIR)/%_psy.f90: \
$(WORKING_DIR)/%.x90 $$(OPTIMISATION_PATH)/$(DSL)/$$*.py | $$(dir $$@)
	$(call MESSAGE,PSyclone - local optimisation,$(subst $(SOURCE_DIR)/,,$<))
	$QPYTHONPATH=$(LFRIC_BUILD)/psyclone:$$PYTHONPATH $(PSYCLONE) -api lfric \
	           -d $(WORKING_DIR) \
	           --config $(PSYCLONE_CONFIG_FILE) \
	           -s $(OPTIMISATION_PATH)/$(DSL)/$*.py \
	           -okern $(WORKING_DIR)/kernel \
	           -oalg $(WORKING_DIR)/$*.f90 \
	           -opsy $(WORKING_DIR)/$*_psy.f90 \
	           $(PSYCLONE_PSYKAL_EXTRAS) \
	           $<

# Where a global optimisation script exists, use it.
#
$(WORKING_DIR)/%.f90 $(WORKING_DIR)/%_psy.f90: \
$(WORKING_DIR)/%.x90 $(OPTIMISATION_PATH)/$(DSL)/global.py | $$(dir $$@)
	$(call MESSAGE,PSyclone - global optimisation,$(subst $(SOURCE_DIR)/,,$<))
	$QPYTHONPATH=$(LFRIC_BUILD)/psyclone:$$PYTHONPATH $(PSYCLONE) -api lfric \
	           -d $(WORKING_DIR) \
	           --config $(PSYCLONE_CONFIG_FILE) \
	           -s $(OPTIMISATION_PATH)/$(DSL)/global.py \
	           -okern $(WORKING_DIR)/kernel \
	           -oalg  $(WORKING_DIR)/$*.f90 \
	           -opsy $(WORKING_DIR)/$*_psy.f90 \
	           $(PSYCLONE_PSYKAL_EXTRAS) \
	           $<

# Where no optimisation script exists, don't use it.
#
$(WORKING_DIR)/%.f90 $(WORKING_DIR)/%_psy.f90: \
$(WORKING_DIR)/%.x90 | $$(dir $$@)
	$(call MESSAGE,PSyclone,$(subst $(SOURCE_DIR)/,,$<))
	$QPYTHONPATH=$(LFRIC_BUILD)/psyclone:$$PYTHONPATH $(PSYCLONE) -api lfric \
	           -l all -d $(WORKING_DIR) \
	           --config $(PSYCLONE_CONFIG_FILE) \
	           -okern $(WORKING_DIR)/kernel \
	           -oalg  $(WORKING_DIR)/$*.f90 \
	           -opsy $(WORKING_DIR)/$*_psy.f90 \
	           $(PSYCLONE_PSYKAL_EXTRAS) \
	           $<

.PRECIOUS: $(WORKING_DIR)/%.x90
# Perform preprocessing for big X90 files.
#
ifeq ("$(FORTRAN_COMPILER)", "nvfortran")
$(WORKING_DIR)/%.x90: $(SOURCE_DIR)/%.X90 | $$(dir $$@)
	$(call MESSAGE,Preprocessing, $(subst $(SOURCE_DIR)/,,$<))
	$Q$(FPP) $(FPPFLAGS) $(MACRO_ARGS) -o $@ $<
else
$(WORKING_DIR)/%.x90: $(SOURCE_DIR)/%.X90 | $$(dir $$@)
	$(call MESSAGE,Preprocessing, $(subst $(SOURCE_DIR)/,,$<))
	$Q$(FPP) $(FPPFLAGS) $(MACRO_ARGS) $< $@
endif

# Little x90 files are just copied to the workspace.
#
$(WORKING_DIR)/%.x90: $(SOURCE_DIR)/%.x90 | $$(dir $$@)
	$(call MESSAGE,Copying, $(subst $(SOURCE_DIR)/,,$<))
	$Qcp $< $@

# Create directories in the workspace as needed.
#
$(DIRECTORIES):
	$(call MESSAGE,Creating,$@)
	$Qmkdir -p $@
