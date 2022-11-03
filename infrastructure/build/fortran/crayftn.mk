##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
# Various things specific to the Cray Fortran compiler.
##############################################################################
#
# This macro is evaluated now (:= syntax) so it may be used as many times as
# desired without wasting time rerunning it.
#
CRAYFTN_VERSION := $(shell ftn -V 2>&1 \
                     | awk -F "[. ]" '/[0-9]\.[0-9]\.[0-9]/ { printf "%03i%03i%03i", $$5,$$6,$$7}' )

$(info ** Chosen Cray Fortran compiler version $(CRAYFTN_VERSION))

ifeq ($(shell test $(CRAYFTN_VERSION) -lt 008007000; echo $$?), 0)
  $(error CrayFTN is too old. It must be at least 8.7.0)
endif

F_MOD_DESTINATION_ARG     = -J
F_MOD_SOURCE_ARG          = -p
OPENMP_ARG            = -h omp

FFLAGS_COMPILER           =
FFLAGS_NO_OPTIMISATION    = -O0
FFLAGS_SAFE_OPTIMISATION  = -O2
FFLAGS_RISKY_OPTIMISATION = -O3
FFLAGS_DEBUG              = -Gfast
FFLAGS_WARNINGS           = -m 0
FFLAGS_UNIT_WARNINGS      = -m 0
FFLAGS_RUNTIME            = -R bcdps -Ktrap=fp
# Option for checking code meets Fortran standards
FFLAGS_FORTRAN_STANDARD   = -en

LDFLAGS_COMPILER =

FPPFLAGS = -P

DEPRULE_FLAGS = -moduleobjects
