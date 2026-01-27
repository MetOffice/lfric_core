##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

MPIF90_VN_STR := $(shell $(FC) --version)
MPIF90_COMPILER := $(shell echo "$(MPIF90_VN_STR)" | awk '{print $$1}')
$(info ** Chosen MPI Fortran compiler: $(MPIF90_COMPILER))

ifeq '$(MPIF90_COMPILER)' 'GNU'
  FORTRAN_COMPILER = gfortran
else ifeq '$(MPIF90_COMPILER)' 'ifort'
  FORTRAN_COMPILER = ifort
else ifeq '$(MPIF90_COMPILER)' 'Cray'
  FORTRAN_COMPILER = crayftn
else ifeq '$(MPIF90_COMPILER)' 'nvfortran'
  FORTRAN_COMPILER = nvfortran
else
  $(error Unrecognised mpif90 compiler option: "$(MPIF90_COMPILER)")
endif

include $(LFRIC_BUILD)/fortran/$(FORTRAN_COMPILER).mk
