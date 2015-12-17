##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

DOXYGEN  ?= doxygen

DOCUMENT_DIRECTORIES = $(patsubst %/,%,$(dir $(wildcard */*.latex)))

.PHONY: most
most: documents api/infrastructure api/gungho

.PHONY: api/%
api/%: Doxyfile.% | api
	echo Building $* API documentation
	$(DOXYGEN) $<

api:
	echo Creating directory "$@"
	mkdir $@

.PHONY: documents
documents: $(DOCUMENT_DIRECTORIES)

.PHONY: $(DOCUMENT_DIRECTORIES)
$(DOCUMENT_DIRECTORIES):
	$(MAKE) -f ../make/latex.mk -C $@

.PHONY: gource
gource:
	$(MAKE) -C visualisation

clean: $(addprefix clean_,$(DOCUMENT_DIRECTORIES))
	-rm -rf api/*

clean_%: ALWAYS
	$(MAKE) -f ../make/latex.mk -C $* clean

.PHONY: ALWAYS
ALWAYS:

