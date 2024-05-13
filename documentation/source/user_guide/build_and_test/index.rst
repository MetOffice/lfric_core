.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section build and test:

Build and Test Systems
======================

The repository includes systems to build and test the infrastructure
code and applications from within the repository. The system is
supported by a number of tools that are each described.

To summarise, each directory of code including applications,
components and infrastructure, can include unit tests and integration
tests. Makefiles exist in each directory which build and run the tests
and, for applications, build the application.

Directories also include Rose Stem configurations to manage the build
and running of the unit and integration tests as well as building and
testing application with one or more configuration.

.. toctree::
   :caption: Contents:
   :maxdepth: 1

   use_of_make
