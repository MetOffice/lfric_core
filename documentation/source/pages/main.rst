.. ------------------------------------------------------------------------------
     (c) Crown copyright 2023 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

##################
Introduction
##################
.. index::
   Introduction

The `LFRic Core <https://github.com/MetOffice/lfric_core>`_ project
aims to develop a software infrastructure primarily to support the
development of a replacement for the Unified Model but also to provide
a common library that underpins a range of modelling requirements and
related tools. The LFRic Core Project is being led by the `Core
Capability Development Team
<CoreCapabilityDevelopmentTeam@metoffice.gov.uk>`_ within the Science
IT group.

Development of the LFRic infrastructure and the new atmosphere model
are being done within the Momentum Partnership. Key initial aims for
the new model are as follows:

- The model will be scientifically as good as the UM atmosphere.
- The model will scale better on future exascale platforms.
- The infrastructure will be flexible enough to support future
  evolutions of the science.

Broadly speaking, the LFRic core aims to support similar earth system
modelling requirements. The detailed requirements documentation will
describe the current and intended capabilities in greater detail.

Active users and users at the Met Office are requested to work from
the head of main. New external users should contact us for advice on
how to get going.

Guide to the documentation
==========================

For initial orientation, a :ref:`quick overview <section repository
contents>` of the main contents of the LFRic core repository is
given. While LFRic core does include some small LFRic applications,
one should be aware that major applications, including the Momentum
atmosphere model application, are developed in other repositories and
will have their own application-specific documentation.

Before giving an overview of the core infrastructure, an overview of
the `structure of a typical LFRic application <section model
application structure>` is given. It briefly references several LFRic
infrastructure capabilities. These capabilities are then described in
later sections.

Any developer of the LFRic core or an LFRic application should have a
good understanding of the underlying principles behind LFRic and the
core data model, and the scientific model architecture known as
**PSyKAl**. The `LFRic data model and PSyclone <section psykal and
datamodel>` documentation describes key aspects of LFRic and of
PSyclone, the code autogeneration tool that LFRic applications depend
upon.

The `Application Documentation <section applications>` provides links
to documentation for each application developed within the LFRic core
repository, describing the role of the application and including
pointers to the features of the LFRic core that it depends upon or
tests.

Some technical features required by applications have been implemented
as LFRic "components". Each of these is described. TODO Add the section.

Several technical sections are given on major technical topics
including the distributed memory strategy and implementation, the
clock and calendar.

Finally, detailed API documentation is given for the LFRic
infrastructure and each LFRic component.

.. _section repository contents:

Contents of the Repository
==========================

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
