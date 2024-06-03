.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section make:

Use of make
===========

The file ``Makefile`` that can be found in each application or code
library can be used to build or run aspects of the application or
library. Running ``make`` will build application (in an application
directory), build and run any unit tests, and build and run any
integration tests. The Makefile supports options to execute any of
these choices individually.

The Makefile can be used to run the Rose Stem system test
suite: run ``make test-suite``.

Each Makefile takes other optional arguments, for example to
obtain verbose output or to choose different sets of compile
options. Each application or component will describe the options
pertaining to its own build.
