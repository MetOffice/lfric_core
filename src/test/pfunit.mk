##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

CMAKE ?= cmake

PFUNIT_SOURCE_DIR = $(abspath ../pfunit)
PFUNIT_BUILD_DIR = $(DYNAMO_BUILD_ROOT)/pfunit

include $(ROOT)/make/include.mk
include $(MAKE_DIR)/compiler.mk
include $(MAKE_DIR)/fortran/$(FORTRAN_COMPILER).mk

ifeq '$(FORTRAN_COMPILER)' 'ifort'
  PFUNIT_COMPILER_ID = Intel
  ifeq ($(shell test $(IFORT_VERSION) -lt 0130000; echo $$?), 0)
    $(error pFUnit will only compile with ifort v13 or later.)
  else ifeq ($(shell test $(IFORT_VERSION) -lt 0140000; echo $$?), 0)
    export CPPFLAGS += -DINTEL_13
    export FPPFLAGS += -DINTEL_13
  endif
else ifeq '$(FORTRAN_COMPILER)' 'gfortran'
  PFUNIT_COMPILER_ID = GNU
  ifeq ($(shell test $(GFORTRAN_VERSION) -lt 040900; echo $$?), 0)
    $(error pFUnit will only compile with gfortran v4.9.0 or later.)
  endif
else ifeq '$(FORTRAN_COMPILER)' 'nagfor'
  PFUNIT_COMPILER_ID = NAG
else ifeq '$(FORTRAN_COMPILER)' 'xlf'
  PFUNIT_COMPILER_ID = XL
else ifeq '$(FORTRAN_COMPILER)' 'crayftn'
  PFUNIT_COMPILER_ID = CRAY
else
  $(error Unrecognised compiler "$(FORTRAN_COMPILER)")
endif

ifdef CRAY_ENVIRONMENT
  CMAKE_FLAGS = -DCMAKE_TOOLCHAIN_FILE=$(ROOT)/src/test/cmake/cray.$(FORTRAN_COMPILER).cmake
endif

DRIVER_DIR = $(dir $(DRIVER_OBJ))

$(DRIVER_OBJ): $(PFUNIT_INSTALL_DIR)/include/driver.F90 \
               $(DRIVER_DIR)/testSuites.inc | $(DRIVER_DIR)
	@echo -e $(VT_BOLD)Compiling$(VT_RESET) $@
	$(Q)$(FC) $(FFLAGS) -c -I$(PFUNIT_INSTALL_DIR)/mod -I$(DRIVER_DIR) \
	          -DBUILD_ROBUST \
	          -DMPI=YES \
	          -DPFUNIT_EXTRA_USAGE=$(PFUNIT_EXTRA_USAGE) \
	          -DPFUNIT_EXTRA_INITIALIZE=$(PFUNIT_EXTRA_INITIALIZE) \
	          -DPFUNIT_EXTRA_FINALIZE=$(PFUNIT_EXTRA_FINALIZE) -o $@ $<

$(PFUNIT_INSTALL_DIR)/include/driver.F90: $(PFUNIT_BUILD_DIR)/Makefile
	@echo -e $(VT_BOLD)Building$(VT_RESET) pFUnit
	$(Q)$(MAKE) -C $(PFUNIT_BUILD_DIR)
	$(Q)$(MAKE) -C $(PFUNIT_BUILD_DIR) tests install

$(PFUNIT_BUILD_DIR)/Makefile: | $(PFUNIT_BUILD_DIR)
	@echo -e $(VT_BOLD)Configuring$(VT_RESET) pFUnit
	$(Q)cd $(PFUNIT_BUILD_DIR); $(CMAKE) $(CMAKE_FLAGS) \
                                          -DMPI=YES \
                                          $(MPI_RUN_SCRIPT) $(MPI_RUN_SCRIPT_COMPILER) \
                                          -DINSTALL_PATH=$(PFUNIT_INSTALL_DIR) \
                                          -DCMAKE_INSTALL_PREFIX=$(PFUNIT_INSTALL_DIR) \
	                                  $(PFUNIT_SOURCE_DIR)

$(PFUNIT_BUILD_DIR) $(dir $(DRIVER_OBJ)):
	@echo -e $(VT_BOLD)Creating$(VT_RESET) $@
	$(Q)mkdir $@

ALL_TESTS = $(shell find . -name "*.pf")

$(DRIVER_DIR)/testSuites.inc: $(DRIVER_DIR)/testSuites.inc.new $(ALL_TESTS)
	@echo -e $(VT_BOLD)Replacing$(VT_RESET) $@ with $<
	$(Q)mv -f $< $@

$(DRIVER_DIR)/testSuites.inc.new: ALWAYS | $(DRIVER_DIR)
	@echo -e $(VT_BOLD)Clearing$(VT_RESET) $@
	@echo > $@

$(ALL_TESTS): ALWAYS
	@echo -e $(VT_BOLD)Adding$(VT_RESET) $@
	@echo ADD_TEST_SUITE\($(notdir $(basename $@))_suite\) >> $(DRIVER_DIR)testSuites.inc.new

.PHONY: ALWAYS
ALWAYS:

.PHONY: clean
clean:
	-rm -rf $(PFUNIT_BUILD_DIR)
	-rm -rf $(PFUNIT_INSTALL_DIR)
