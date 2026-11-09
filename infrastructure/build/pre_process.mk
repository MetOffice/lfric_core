##############################################################################
# Copyright (c) 2026,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
#
# Locate and compile all Fortran .F90 and .X90 files and preprocess them.
#
# The following variables may be specified to modify the compile process...
#
# PRE_PROCESS_INCLUDE_DIRS: Space separated list of directories to search for
#                           inclusions.
# PRE_PROCESS_MACROS: Space separated list of macro definitions in the form
#                     NAME[=MACRO] to be passed to the compiler.


# ------
# ROOT: Project directory
# BIN_DIR: Path to directory for resulting executables.
#          Default: $(ROOT)/bin
# FFLAG_GROUPS: Space separated list of FFLAG_<group name> variables to use in
#               building up the FFLAGS variable. Only the group name is
#               specified.
# LINK_TYPE: 'static' or 'dynamic' linking.
#            Default: dynamic, except on Crays where it's static
# PROGRAMS: Names of programs to compile.
#           Default: Everything listed in programs.mk
# PROJECT_MAKE_DIR: Used to locate project specific targets modifiers and
#                   such.
# COMPILE_OPTIONS: Name of an optional file that can be included to list
#                  project specific compile options
#
##############################################################################

include $(LFRIC_BUILD)/fortran.mk

.SECONDEXPANSION:

# Build a set of "-I" arguments to seach the whole object tree:
INCLUDE_ARGS := $(addprefix -I, $(PRE_PROCESS_INCLUDE_DIRS))

# Build a set of "-D" argument for any pre-processor macros
#
MACRO_ARGS := $(addprefix -D,$(PRE_PROCESS_MACROS))

# Add a special compiler-specific define (in case that the compiler is not used
# as preprocesser, e.g. nvidia). This way, compiler-specific preprocessing
# can be done in all cases using consistent macros:
# IFORT_COMPILER, IFX_COMPILER, GNU_COMPILER, NVFORTRAN_COMPILER

UPPER_NAME := $(shell echo $(FORTRAN_COMPILER) | tr '[:lower:]' '[:upper:]')
MACRO_ARGS += -D$(UPPER_NAME)_COMPILER

PP_FILES := $(shell find $(WORKING_DIR) -name '*.[FX]90' -print)
# Convert to lower case
LOWER_PP_FILES := $(patsubst %.X90,%.x90,$(patsubst %.F90,%.f90,$(PP_FILES)))

pre_process: $(LOWER_PP_FILES)

# nvfortran as preprocessor needs different options
PREPROCESSOR = $(firstword $(subst -, ,$(notdir $(FPP))))
$(info Preprocessor $(FPP)  $(PREPROCESSOR) XXX $(FPPFLAGS))

ifeq ("$(PREPROCESSOR)", "nvfortran")
# nvfortran with -P will write the output file to the current directory, NOT
# to the directory specified with -o (and not to stdout either).
# So, instead of -P, we pipe the files through grep to remove linemarkers
%.x90: %.X90
	$(call MESSAGE,Preprocessing $(FPP) $(PREPROCESSOR), $<)
	$Q$(FPP) $(FPPFLAGS) $(INCLUDE_ARGS) $(MACRO_ARGS) $< | grep -v "^#" > $@
	$Qrm $<
%.f90: %.F90
	$(call MESSAGE,Preprocessing $(FPP) $(PREPROCESSOR), $<)
	$(call MESSAGE,Preprocessing, $<)
	$Q$(FPP) $(FPPFLAGS) $(INCLUDE_ARGS) $(MACRO_ARGS) $< | grep -v "^#" > $@
	$Qrm $<

else

%.x90: %.X90
	$(call MESSAGE,Preprocessing, $<)
	# Some compilers (e.g. gfortran) do not allow .X90 files to
	# be preprocessed. So, rename the .X90 into .F90 before preprocessing:
	mv $< $(basename $<).F90
	$Q$(FPP) $(FPPFLAGS) $(INCLUDE_ARGS) $(MACRO_ARGS) $(basename $<).F90 > $@
	$Qrm $(basename $<).F90
%.f90: %.F90
	$(call MESSAGE,Preprocessing, $<)
	$Q$(FPP) $(FPPFLAGS) $(INCLUDE_ARGS) $(MACRO_ARGS) $< > $@
	$Qrm $<

endif

# Get MESSAGE function from lfric.mk
include $(LFRIC_BUILD)/lfric.mk
