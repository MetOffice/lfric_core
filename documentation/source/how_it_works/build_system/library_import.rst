.. ------------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------
.. _library_import:

Library Import
==============

The build system can automatically import libraries using their
``build/import.mk`` fragment.

In the ``build`` target of your project's top-level ``Makefile`` you need to
recursively call Make on this fragment::

    $(MAKE) $(QUIET_ARG) -f $LIBRARY/build/import.mk

This is usually done in a loop in order to easily import multiple libraries::

    $(Q)for SUBPROJECT in $(INTERNAL_DEPENDENCIES) ; do \
        $(MAKE) $(QUIET_ARG) -f $$SUBPROJECT/build/import.mk ; done

In which case any new library need only be added to the list.

Unit and integration testing require less manual handle-turning since they
provide this mechanism automatically::

    unit-tests/%: export IMPORT_PARTS   = $(CORE_ROOT_DIR)/infrastructure \
                                          $(CORE_ROOT_DIR)/components/science

What's more this mechanism can import them even if they don't have an
``import.mk`` fragment. This default import will do a source extract and
PSyclone any algorithms.
