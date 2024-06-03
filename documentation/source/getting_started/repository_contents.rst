.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _repository_contents:

Contents of the Repository
--------------------------

The repository contains the following directories:

- The ``infrastructure`` directory contains code that supports the use
  of the core data model of LFRic, such as for creating, manipulating
  and grouping fields of data; supporting the generation of meshes and
  supporting the management of application configuration options.
- The ``components`` directory contains several libraries of code each
  of which can be used to support an LFRic model requirement. For
  example, the **lfric-xios** component provides an API to allow LFRic
  applications to use the XIOS IO library, and the **driver**
  component contains code that can support the construction of the
  top-level calling tree and data structures of a typical LFRic application.
- The ``mesh_tools`` directory contains application code to generate
  meshes used by LFRic applications.
- The ``apps`` directory contains application code.

Many of the directories contain directories of unit and integration
tests, directories of Rose metadata and Rose Stem tasks for supporting
running of tests on other host machines, and Makefiles for building
applications, building and running of local tests, or running Rose
suites.