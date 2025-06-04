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

.PHONY: help Makefile
.SECONDEXPANSION:

# Specific rule for HTML which generates dynamic content prior to
# documentation.
#
html: $(patsubst source/_plotly/%.py, \
                 generated/html/plots/%.html, \
                 $(wildcard source/_plotly/*))
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
#
# %: Makefile
# 	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
#
#deploy:
# 	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(DEPLOYDIR)/$(BRANCH)" $(SPHINXOPTS) $(O)

generated/html/plots/%.html: source/_plotly/%.py build/html/_static/plotly.min.js | $$(dir $$@)
	python3 $<

# This is a slightly alarming bodge to get hold of Plotly's JavaScript
# component from the Python component without making assumptions about which
# version it is.
#
build/html/_static/plotly.min.js: | $$(dir $$@) generated/html/js/
	python3 -c "from plotly.graph_objects import Figure; \
	            null_fig = Figure(); \
	            null_fig.write_html('generated/html/js/null.html', \
	                                include_plotlyjs='directory', \
	                                full_html=True)"
	mv generated/html/js/plotly.min.js $@
	rm generated/html/js/null.html

generated/html/plots/ generated/html/js/ build/html/_static/:
	mkdir -p $@

.PHONY: clean
clean:
	@echo "***** Removing, $(BUILDDIR) directory *****"
	@rm -rf $(BUILDDIR)

