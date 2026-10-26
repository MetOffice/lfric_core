##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Various things specific to the amdflang compiler.
##############################################################################

$(info ** Chosen AMD Flang)

F_MOD_DESTINATION_ARG     = -J
F_MOD_SOURCE_ARG          = -I

FFLAGS_OPENMP  = -fopenmp
LDFLAGS_OPENMP = -fopenmp

FFLAGS_COMPILER           =
FFLAGS_NO_OPTIMISATION    = -O0
FFLAGS_SAFE_OPTIMISATION  = -O2
FFLAGS_RISKY_OPTIMISATION = -Ofast
FFLAGS_DEBUG              = -g
FFLAGS_WARNINGS           =
FFLAGS_UNIT_WARNINGS      =
FFLAGS_INIT               =
FFLAGS_RUNTIME            =
# fast-debug flags set separately as Intel compiler needs platform-specific control on them
FFLAGS_FASTD_INIT         = $(FFLAGS_INIT)
FFLAGS_FASTD_RUNTIME      = $(FFLAGS_RUNTIME)

# Option for checking code meets Fortran standard (flang only supports 2018)
FFLAGS_FORTRAN_STANDARD   = -std=f2018

LDFLAGS_COMPILER =

FPPFLAGS = -P
