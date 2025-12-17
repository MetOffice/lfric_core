##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
# Various things specific to the Portland Fortran compiler.
##############################################################################
#
# This macro is evaluated now (:= syntax) so it may be used as many times as
# desired without wasting time rerunning it.
#
PGFORTRAN_VERSION := $(shell pgfortran --version \
9	                             | awk -F "[. -]" '/[0-9]+\.[0-9]+-[0-9]+/ { printf "%03i%03i%03i", $$(2), $$(3), $$(4) }' )

$(info ** Chosen Portland Fortran compiler version $(PGFORTRAN_VERSION))

ifeq ($(shell test $(PGFORTRAN_VERSION) -lt 015007000; echo $$?), 0)
  $(error PGFortran is too old. It must be at least 15.7-0)
endif

F_MOD_DESTINATION_ARG = -module$(SPACE)

FFLAGS_OPENMP  = -mp
LDFLAGS_OPENMP = -mp

FFLAGS_COMPILER           =
FFLAGS_NO_OPTIMISATION    = -O0
FFLAGS_SAFE_OPTIMISATION  = -O2
FFLAGS_RISKY_OPTIMISATION = -O4
FFLAGS_DEBUG              = -g -traceback
FFLAGS_RUNTIME            = -Mbounds -Mchkptr -Mchkstk
# Option for checking code meets Fortran standard (not available for PGI)
FFLAGS_FORTRAN_STANDARD   =

LDFLAGS_COMPILER = -g

FPPFLAGS = -P
