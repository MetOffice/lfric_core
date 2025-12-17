##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################

HERE_DIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

export PYTHONPATH := ${PYTHONPATH}:$(HERE_DIR):$(HERE_DIR)/diagnostics_metadata_reconfigurator/src

ifdef VERBOSE
    VERBOSE_ARG = -v
endif
ifdef DEBUG
    DEBUG_ARG = -vv
endif

ifeq "x$(shell which flake8)" "x"
	TEST_PREREQ = pytest
else
	TEST_PREREQ = flake8
endif

.PHONY: test
test: $(TEST_PREREQ)

.PHONY: flake8
flake8: pytest
	flake8 $(VERBOSE_ARG) $(HERE_DIR)

.PHONY: pytest
pytest: mypy
	pytest $(VERBOSE_ARG) $(DEBUG_ARG) $(HERE_DIR)

.PHONY: mypy
mypy:
	mypy $(VERBOSE_ARG) $(HERE_DIR)
