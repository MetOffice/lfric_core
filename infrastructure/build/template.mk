###############################################################################
# (c) Crown copyright 2023 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
###############################################################################
# Generate compilable Fortran source from template files.
#
# SOURCE_DIR    Directory in which to locate templates.
# SUBSTITUTIONS Name and list of values.
# TEMPLATES     Template file relative to source root.
# WORKING_DIR   Directory to take generated files.
#
###############################################################################
TOOL = $(LFRIC_BUILD)/tools/Templaterator

FORTRAN_FILES = $(subst .t90,.f90,$(subst .T90,.F90, $(TEMPLATES)))

SARGS = $(foreach key, $(SUBSTITUTIONS), $(foreach value, $(wordlist 2, 10, $(subst :, ,$(key))), -s $(firstword $(subst :, ,$(key)))=$(value)))

generate-from-template: $(addprefix $(WORKING_DIR)/, $(FORTRAN_FILES))

$(WORKING_DIR)/%.f90: $(SOURCE_DIR)/%.t90
	$(call MESSAGE,Templating, $<)
	$Q$(TOOL) $< -o $@ $(SARGS)

$(WORKING_DIR)/%.F90: $(SOURCE_DIR)/%.T90
	$(call MESSAGE,Templating, $<)
	$Q$(TOOL) $< -o $@ $(SARGS)

#include $(LFRIC_BUILD)/lfric.mk
