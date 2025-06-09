################################################################################
# (c) Crown copyright 2023 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
################################################################################
# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build


BRANCH := $(shell git rev-parse --abbrev-ref HEAD)
DEPLOYDIR	  = $(HOME)/public_html

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile deploy

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

deploy:
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(DEPLOYDIR)/$(BRANCH)" $(SPHINXOPTS) $(O)

.PHONY: clean
clean:
	@echo "***** Removing, $(BUILDDIR) directory *****"
	@rm -rf $(BUILDDIR)

