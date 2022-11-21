##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Various things specific to the Intel Fortran compiler.
##############################################################################

$(info Project specials for Intel compiler)

export FFLAGS_UM_PHYSICS = -r8

# Set fast-debug and full-debug options for lfric_atm Intel that
# differ from the setup in infrastructure used by other LFRic apps.
# First, the -check all option is removed from all configurations as
# it results in excessive numbers of warnings from UM code running
# within LFRic
FFLAGS_RUNTIME            = -fpe0
# Second, certain fast-debug options cause XIOS failures on the Cray
# in jobs that write diagnostic so they are removed. Note: the
# full-debug test can still use these options as it avoids such XIOS
# use.
ifdef CRAY_ENVIRONMENT
# On the Cray these options are switched off for fast-debug
FFLAGS_FASTD_INIT               = 
FFLAGS_FASTD_RUNTIME            =
else
# Otherwise, use the same as the default full-debug settings
FFLAGS_FASTD_INIT         = $(FFLAGS_INIT) 
FFLAGS_FASTD_RUNTIME      = $(FFLAGS_RUNTIME)
endif

$(info LFRic compile options required for files with OpenMP - see Ticket 1490)
%psy.o %psy.mod:   private FFLAGS_EXTRA = $(FFLAGS_INTEL_FIX_ARG)
# NOTE: The -qoverride-limits option contained in $(FFLAGS_INTEL_FIX_ARG) is
# not currently applied here. This is a temporary workaround for #3205
# which it was found to be inadvertently preventing compilation
# psy/%.o psy/%.mod: private FFLAGS_EXTRA = $(FFLAGS_INTEL_FIX_ARG)





