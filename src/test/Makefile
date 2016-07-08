##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

IGNORE_DEPENDENCIES = netcdf ESMF
EXTERNAL_LIBRARIES = esmf netcdff netcdf hdf5_fortran

export ROOT = ../..

include $(ROOT)/make/include.mk
include $(MAKE_DIR)/compiler.mk
include $(MAKE_DIR)/fortran/$(FORTRAN_COMPILER).mk
include $(MAKE_DIR)/mpi_link.mk

FFLAGS += -O0 -g

UNIT_TEST_LAUNCH_TIME = 2
UNIT_TEST_RUN_TIME    = 2

ifdef CRAY_ENVIRONMENT
LD_COMPILER_LIBRARIES = stdc++
endif

ifeq ($(FORTRAN_COMPILER), ifort)
# This will be required until Intel promote it to default behaviour. Possibly
# in ifort 17
FFLAGS += -assume realloc_lhs
endif

PFUNIT_INSTALL_DIR      = $(DYNAMO_BUILD_ROOT)/pfunit-install
PFUNIT_EXTRA_USAGE      = suite_fixture_mod
PFUNIT_EXTRA_INITIALIZE = set_up_suite
PFUNIT_EXTRA_FINALIZE   = tear_down_suite
FFLAGS += -I$(OBJ_DIR)/support

OBJ_DIR = $(DYNAMO_BUILD_ROOT)/tests
OBJ_SUBDIRS = $(patsubst %,$(OBJ_DIR)/%,$(filter-out data,$(shell find . -type d -not -path '*/.*' )))

DRIVER_OBJ=$(OBJ_DIR)/driver.o

TEST_SOURCES = $(shell find -name "*.pf")
TEST_OBJECTS = $(patsubst %.pf,$(OBJ_DIR)/%.o,$(TEST_SOURCES))

SUPPORT_SOURCE  = $(shell find -name "*.[Ff]90")
SUPPORT_OBJECTS = $(OBJ_DIR)/support/feign_config_mod.o \
                  $(addprefix $(OBJ_DIR)/,$(addsuffix .o,$(basename $(SUPPORT_SOURCE))))
SUPPORT_MODINC  = $(patsubst %,-I%,$(dir $(SUPPORT_OBJECTS)))

SUBJECT_SOURCE  = $(ROOT)/src/dynamo
SUBJECT_DIR     = $(DYNAMO_BUILD_ROOT)/dynamo
SUBJECT_SUBDIRS = $(shell find $(SUBJECT_DIR) -type d)
SUBJECT_MODINC  = $(patsubst %,-I%,$(SUBJECT_SUBDIRS))

export PFUNIT_INSTALL_DIR \
       PFUNIT_EXTRA_USAGE PFUNIT_EXTRA_INITIALIZE PFUNIT_EXTRA_FINALIZE

.PHONY: test
test: $(OBJ_DIR)/test | $(OBJ_DIR)/data
	@echo -e $(VT_BOLD)Running$(VT_RESET) $@
	$(Q)$(OBJ_DIR)/test -robust \
	                    -max-timeout-duration $(UNIT_TEST_RUN_TIME) \
	                    -max-launch-duration  $(UNIT_TEST_LAUNCH_TIME)

.PRECIOUS: $(OBJ_DIR)/%.F90
$(OBJ_DIR)/%.F90: %.pf | $(OBJ_DIR) $(OBJ_SUBDIRS) $(DRIVER_OBJ) \
                         $(SUPPORT_OBJECTS)
	@echo -e $(VT_BOLD)Generating$(VT_RESET) $@
	$(Q)$(PFUNIT_INSTALL_DIR)/bin/pFUnitParser.py $< $@

$(OBJ_DIR)/%.o: $(OBJ_DIR)/%.F90 $(SUPPORT_OBJECTS)
	@echo -e $(VT_BOLD)Compiling$(VT_RESET) $<
	$(Q)$(FC) $(FPPFLAGS) $(FFLAGS) \
	          $(F_MOD_DESTINATION_ARG)$(dir $@) \
	          $(SUBJECT_MODINC) $(SUPPORT_MODINC) \
	          -I$(PFUNIT_INSTALL_DIR)/mod -c -o $@ $<

.SECONDARY: $(OBJ_DIR)/test
$(OBJ_DIR)/test: $(SUPPORT_OBJECTS) $(TEST_OBJECTS) $(DRIVER_OBJ) \
                 $(SUBJECT_DIR)/modules.a
	@echo -e $(VT_BOLD)Linking$(VT_RESET) $@
	$(Q)$(LDMPI) $(LDFLAGS) $(LDFLAGS_FORTRAN) \
	             -L$(PFUNIT_INSTALL_DIR)/lib \
	             -o $@ $^ \
	             -lpfunit $(patsubst %,-l%,$(EXTERNAL_LIBRARIES)) \
	             $(patsubst %,-l%,$(LD_COMPILER_LIBRARIES))

$(OBJ_DIR)/support/feign_config_mod.f90: $(wildcard $(SUBJECT_SOURCE)/configuration/*.nld) | $(OBJ_DIR)/support
	$(TOOL_DIR)/GenerateFeigns -output $@ $^

# Replace OBJ_SUBDIR with mangling of % ?
$(OBJ_DIR)/%.o: %.F90 | $(OBJ_SUBDIRS)
	$(Q)$(FC) $(FPPFLAGS) $(FFLAGS) \
	          $(F_MOD_DESTINATION_ARG)$(dir $@) \
	          $(SUBJECT_MODINC) $(SUPPORT_MODINC) -c -o $@ $<

$(OBJ_DIR)/%.o: %.f90 | $(OBJ_DIR)/data $(OBJ_SUBDIRS)
	$(Q)$(FC) $(FPPFLAGS) $(FFLAGS) \
	          $(F_MOD_DESTINATION_ARG)$(dir $@) \
	          $(SUBJECT_MODINC) $(SUPPORT_MODINC) -c -o $@ $<

$(OBJ_DIR)/%.o: $(OBJ_DIR)/%.f90
	$(Q)$(FC) $(FPPFLAGS) $(FFLAGS) \
	          $(F_MOD_DESTINATION_ARG)$(dir $@) \
	          $(SUBJECT_MODINC) $(SUPPORT_MODINC) -c -o $@ $<

$(DRIVER_OBJ): | $(dir $(DRIVER_OBJ))
	$(MAKE) -f pfunit.mk DRIVER_OBJ=$@

$(OBJ_DIR)/data: data | $(OBJ_DIR)
	@echo Copying data directory
	# We use rsync so we can ignore hidden objects.
	$(Q)rsync -a --exclude='.*' $</ $@/

$(OBJ_DIR) $(OBJ_SUBDIRS) $(dir $(DRIVER_OBJ)) $(OBJ_DIR)/support:
	@echo -e $(VT_BOLD)Creating$(VT_RESET) $@
	$(Q)mkdir -p $@

.PHONY: ALWAYS
ALWAYS:

.PHONY: clean
clean:
	-rm -rf $(OBJ_DIR)
ifdef ALL
	$(MAKE) -f pfunit.mk clean
endif
