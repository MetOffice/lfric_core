##############################################################################
# (C) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Various things specific to the Nvidia Fortran compiler.
##############################################################################
#
# This macro is evaluated now (:= syntax) so it may be used as many times as
# desired without wasting time rerunning it.
#
NVFORT_VERSION := $(shell nvfortran -V | awk '/^nvfortran +[0-9]+\.[0-9]+/ { split($$2, a, "[.-]"); printf "%03i%02i%02i\n", a[1],a[2],a[3] }')
$(info ** Chosen Nvidia Fortran compiler version $(NVFORT_VERSION))
ifeq ($(shell test $(NVFORT_VERSION) -lt 0241100; echo $$?), 0)
  $(error nvFort is too old to build LFRic. Must be at least 24.11)
endif

F_MOD_DESTINATION_ARG = -module$(SPACE)

FFLAGS_COMPILER           =
FFLAGS_COMPILER          += -Mfree -Mpreprocess
FFLAGS_NO_OPTIMISATION    = -O0
FFLAGS_SAFE_OPTIMISATION  = -O2 -Mnofma -Mnovect
FFLAGS_RISKY_OPTIMISATION = -Ofast
FFLAGS_DEBUG              = -g -traceback
FFLAGS_RUNTIME            =
# Option for checking code meets Fortran standard (not available for PGI)
FFLAGS_FORTRAN_STANDARD   =

LDFLAGS_COMPILER = -g

# Flags for OpenMP threading / OpenMP offloading / OpenACC Offloading
# The LFRIC_OFFLOAD_DIRECTIVES env_variable is also queried in the PSyclone
# script to generate matching directives
ifeq ("$(LFRIC_OFFLOAD_DIRECTIVES)", "omp")
	FFLAGS_OPENMP  = -mp=gpu -gpu=mem:managed
	LDFLAGS_OPENMP = -mp=gpu -gpu=mem:managed -cuda
else ifeq ("$(LFRIC_OFFLOAD_DIRECTIVES)", "acc")
	FFLAGS_OPENMP  = -acc=gpu -gpu=mem:managed -mp=multicore
	LDFLAGS_OPENMP = -acc=gpu -gpu=mem:managed -mp=multicore -cuda
# Flags for the same targets but numerically matching the cpu results
else ifeq ("$(LFRIC_OFFLOAD_DIRECTIVES)", "omp_uniform")
	FFLAGS_OPENMP  = -mp=gpu -gpu=mem:managed,math_uniform
	LDFLAGS_OPENMP = -mp=gpu -gpu=mem:managed,math_uniform -cuda
else ifeq ("$(LFRIC_OFFLOAD_DIRECTIVES)", "acc_uniform")
	FFLAGS_OPENMP  = -acc=gpu -gpu=mem:managed,math_uniform -mp=multicore
	LDFLAGS_OPENMP = -acc=gpu -gpu=mem:managed,math_uniform -mp=multicore -cuda
else
	FFLAGS_OPENMP  = -mp
	LDFLAGS_OPENMP = -mp
endif

FPP = nvfortran -E
FPPFLAGS = -P -D__NVCOMPILER
