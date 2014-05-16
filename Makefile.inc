##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################
# Variables defined here will will be included in all Makefiles for this
# project.
# Expects ROOT to point to the root of the project tree.

ROOT ?= .
BUILD_DIR = $(ROOT)/build
TOOL_DIR = $(ROOT)/tools

ifndef COMPILER
  OS = $(shell uname -s)
  MACHINE = $(shell uname -m)
  ifeq '$(OS)' 'Linux'
    ifeq '$(MACHINE)' 'x86_64'
      COMPILER = ifort
    else
      $(error Unrecognised processor: $(MACHINE))
    endif
  else ifeq '$(OS)' 'AIX'
    COMPILER = xlf
  else
    $(error Unrecognised operating system: $(OS))
  endif
endif

ifeq '$(COMPILER)' 'ifort'
  # Intel compiler
  $(info ** Chosen Intel Fortran compiler)
  FC     = ifort
  FFLAGS += -g -O0 -xhost
  IFORT_VERSION = $(shell ifort -v 2>&1 | awk '{ print $$NF }' \
                   | awk -F . '{ printf "%03i%02i%02i", $$1, $$2, $$3 }')
  $(info ** Version $(IFORT_VERSION))

else ifeq '$(COMPILER)' 'xlf'
  # IBM compiler
  $(info ** Chosen IBM XL Fortran compiler)
  FC     = xlf
  FFLAGS +=
  F_MOD_DESTINATION_ARG = -qmoddir=$(OBJ_DIR)

else ifeq '$(COMPILER)' 'gfortran'
  # GNU compiler
  $(info ** Chosen GNU Fortran compiler)
  FC     = gfortran
  FFLAGS += -g -O0
  FWARNINGFLAGS = -Wall -Werror
  F_MOD_DESTINATION_ARG = -J $(OBJ_DIR)
  GFORTRAN_VERSION = $(shell $(FC95) --version 2>&1 \
                       | sed -n 's/.*[[:space:]]\{1,\}\([[:digit:]]\{1,\}\.[[:digit:]]\{1,\}\.[[:digit:]]\{1,\}\)[[:space:]]\{,\}.*/\1/p' \
                       | awk -F . '{ printf "%02i%02i%02i", $$1, $$2, $$3 }')
  $(info ** Version $(GFORTRAN_VERSION))

endif
