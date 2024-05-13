.. ------------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------
.. _getting_started_index:

Getting Started
===============

Guide to the documentation
--------------------------

For initial orientation, a :ref:`quick overview <section repository
contents>` of the main contents of the LFRic core repository is
given. As noted, the Momentum atmosphere model, and related
applications, are developed in separate repositories. The LFRic core
repository does include some small LFRic applications for
training purposes or for developing and testing particular
technical capabilities.

Before giving an overview of the core infrastructure, an overview of
the :ref:`structure of a typical LFRic application <section
application structure>` is given. It briefly references several LFRic
core capabilities. These capabilities are then described in later
sections.

Any developer of the LFRic core or an LFRic application should have a
good understanding of the underlying principles behind LFRic and the
core data model, and the scientific model architecture known as
**PSyKAl**. The :ref:`LFRic data model and PSyclone <section psykal
and datamodel>` documentation describes key aspects of LFRic and of
PSyclone, the code autogeneration tool that LFRic applications depend
upon.

The :ref:`Application Documentation section<section applications>` provides
links to documentation for each application developed within the LFRic
core repository, describing the role of the application and including
pointers to the features of the LFRic core that it depends upon or
tests.

The :ref:`Meshes section<section meshes and tools>` describes the
LFRic mesh generator and LFRic meshes, including discussion of mesh
partitioning, mesh hierarchies and mesh maps, the LFRic mesh
generator.

The :ref:`Components section <section components>` describes
components which are code libraries delivering specific capabilities
required by some LFRic applications.

The :ref:`build and test system section <section build and test>` describes
the build and test system, and includes descriptions of the tools that
underpin it.

The :ref:`Technical Articles section <section technical articles>`
includes articles on several topics including the distributed memory
strategy and implementation, the clock and calendar.

Finally, detailed API documentation is given for the LFRic
infrastructure and each LFRic component.

.. _section repository contents:

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