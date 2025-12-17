##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
# Various things specific to the Intel Fortran compiler.
##############################################################################
#
# This macro is evaluated now (:= syntax) so it may be used as many times as
# desired without wasting time rerunning it.
#
IFORT_VERSION := $(shell ifort -v 2>&1 \
                       | cut -d' ' -f3 \
                       | awk -F "." '/[0-9]\.[0-9]/ { yy = $$1 % 100; printf "%03i%02i%02i\n", yy,$$2,$$3}' )

$(info ** Chosen Intel Fortran compiler version $(IFORT_VERSION))

ifeq ($(shell test $(IFORT_VERSION) -lt 0150001; echo $$?), 0)
  $(error IFort is too old to build dynamo. Must be at least 15.0.1)
endif

F_MOD_DESTINATION_ARG = -module$(SPACE)
F_MOD_SOURCE_ARG      = -I
FORTRAN_RUNTIME       =

FFLAGS_OPENMP  = -qopenmp
LDFLAGS_OPENMP = -qopenmp

FFLAGS_NO_OPTIMISATION    = -O0
FFLAGS_SAFE_OPTIMISATION  = -O2 -fp-model strict
FFLAGS_RISKY_OPTIMISATION = -O3 -xhost
FFLAGS_DEBUG              = -g -traceback
#
# By default turning interface warnings on causes "genmod" files to be
# created. This adds unecessary files to the build so we disable that
# behaviour.
#
FFLAGS_WARNINGS           = -warn all -warn errors -gen-interfaces nosource
FFLAGS_UNIT_WARNINGS      = -warn all -gen-interfaces nosource
FFLAGS_INIT               = -ftrapuv

ifeq ($(shell test $(IFORT_VERSION) -ge 0190000 -a $(IFORT_VERSION) -lt 0190100; echo $$?), 0)
# It looks like there's a bad interaction between array shape checking and
# the "matmul" intrinsic in at least some iterations of v19.
#
FFLAGS_RUNTIME            = -check all,noshape -fpe0
else
FFLAGS_RUNTIME            = -check all -fpe0
endif

# Certain compile options cause XIOS failures on the Cray xc40 in
# those fast-debug jobs that write diagnostic. Therefore, we remove
# them for that platform. Note: the full-debug test can still use
# these options as it avoids such XIOS use.
ifdef CRAY_ENVIRONMENT
# On the Cray xc40 plaforms these options are switched off for fast-debug
FFLAGS_FASTD_INIT         =
FFLAGS_FASTD_RUNTIME      =
else
# Otherwise, use the same as the default full-debug settings
FFLAGS_FASTD_INIT         = $(FFLAGS_INIT)
FFLAGS_FASTD_RUNTIME      = $(FFLAGS_RUNTIME)
endif

# Option for checking code meets Fortran standard - currently 2008
FFLAGS_FORTRAN_STANDARD   = -stand f08

#########################################################################
# Application and file-specific options referenced in
# build/compile_options.mk files
#
# These variables need explanatory comments and need to be exported
#
# -qoverride-limits applied to PSy-layer code due to Intel compiler bug
# ref #1486
# When the Intel bug is fixed, this option will be removed by #1490
export FFLAGS_INTEL_FIX_ARG         = -qoverride-limits
#
# -warn noexternals applied to code that imports lfric_mpi_mod to avoid
# a warning-turned-error about missing interfaces for MPI calls in
# mpi.mod, such as MPI_Allreduce - switching to mpi_f08.mod resolves
# this via polymorphic interface declarations. Some SOCRATES functions
# do not currently declare interfaces either. Flag was introduced in
# Intel Fortran v19.1.0 according to Intel release notes.
ifeq ($(shell test "$(IFORT_VERSION)" -ge 0190100; echo $$?), 0)
  $(info ** Activating externals warning override for selected source files)
  export FFLAGS_INTEL_EXTERNALS = -warn noexternals
else
  export FFLAGS_INTEL_EXTERNALS =
endif
########################################################################

# The "-assume realloc-lhs" switch causes Intel Fortran prior to v17 to
# actually implement the Fortran2003 standard. At version 17 it becomes the
# default behaviour.
ifeq ($(shell test "$(IFORT_VERSION)" -lt 0170000; echo $$?), 0)
  $(info ** Activating Intel "Make it work" switch for version earlier than 17)
  FFLAGS_COMPILER += -assume realloc-lhs
endif

# There seems to be a recurring bug in IFort which causes it to trip over
# type checking functions which return types. This is how the pFUnit driver
# works so we need to sledge hammer this problem.
#
ifeq ($(shell test "$(IFORT_VERSION)" -le 0190000; echo $$?), 0)
  $(info ** Deactivating interface checking on unit tests)
  FFLAGS_UNIT_WARNINGS += -warn nointerfaces
endif

LDFLAGS_COMPILER =

FPPFLAGS = -P
