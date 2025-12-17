##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
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

FFLAGS_OPENMP  = -h omp
LDFLAGS_OPENMP = -h omp

FFLAGS_COMPILER           =
FFLAGS_NO_OPTIMISATION    = -O0
FFLAGS_SAFE_OPTIMISATION  = -O2 -hflex_mp=strict
FFLAGS_RISKY_OPTIMISATION = -O3 -hipa2

#Cray has debug levels tied to optimisation levels
ifeq ($(shell expr ${CRAYFTN_VERSION} \>= 015000000), 1)
  ifeq "$(PROFILE)" "full-debug"
    FFLAGS_DEBUG    = -G0
  else ifeq "$(PROFILE)" "fast-debug"
    FFLAGS_DEBUG    = -G2
  else ifeq "$(PROFILE)" "production"
    FFLAGS_DEBUG    =
  endif
endif

# Warnings
# ftn-664 : WARNING:  Actual argument "%s" has the PROTECTED attribute. The
#                     associateddummy argument "%s" does not have the INTENT(IN)
#                     attribute.
# ftn-7208 : WARNING: Privatized version of variable "var" is used before it is
#                     defined.
# ftn-7212 : WARNING: Variable "var" is used before it is defined.
ifeq "$(PROFILE)" "production"
    # Set warning output for production runs (not-debug)
    # -m 3 : Error, Warning (default)
    FFLAGS_WARNINGS = -m 3 -M E664,E7208,E7212
else
    # More verbose option
    # -m 0 : Error, Warning, Caution, Note, and Comment
    FFLAGS_WARNINGS = -m 0 -M E664,E7208,E7212
endif

FFLAGS_UNIT_WARNINGS      = -m 0
FFLAGS_RUNTIME            = -R bcdps
# fast-debug flags set separately as Intel compiler needs platform-specific control on them.
# Though, Cray will not set them to anything
FFLAGS_FASTD_RUNTIME      =
FFLAGS_FASTD_INIT         =

# Option for checking code meets Fortran standards
FFLAGS_FORTRAN_STANDARD   = -en

# Floating point checking with CCE causes XIOS failures. To allow
# flexibility for testing, do not apply above full-debug
ifeq "$(PROFILE)" "full-debug"
  # -Ktrap=fp : Trap on divz, inv, or ovf exceptions
  LDFLAGS_COMPILER = -Ktrap=fp
else
  LDFLAGS_COMPILER =
endif

DEPRULE_FLAGS = -moduleobjects

FPPFLAGS    = -P