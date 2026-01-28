##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

MPIFORT_VN_STR := $(shell $(FC) --version)
MPIFORT_COMPILER := $(shell echo "$(MPIFORT_VN_STR)" | awk '{print $$1}')
$(info ** Chosen MPI Fortran compiler: $(MPIFORT_COMPILER))

ifeq '$(MPIFORT_COMPILER)' 'GNU'
  FORTRAN_COMPILER = gfortran
else ifeq '$(MPIFORT_COMPILER)' 'ifort'
  FORTRAN_COMPILER = ifort
else ifeq '$(MPIFORT_COMPILER)' 'Cray'
  FORTRAN_COMPILER = crayftn
else ifeq '$(MPIFORT_COMPILER)' 'nvfortran'
  FORTRAN_COMPILER = nvfortran
else
  $(error Unrecognised mpifort compiler option: "$(MPIFORT_COMPILER)")
endif

include $(LFRIC_BUILD)/fortran/$(FORTRAN_COMPILER).mk
