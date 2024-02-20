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
core data model and model architecture. The `LFRic data model and
PSyclone <section psykal and datamodel>` documentation describes key
aspects of LFRic and of PSyclone, the code autogeneration tool that
LFRic applications depend upon.

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
- The ``components`` directory...TODO
- The ``mesh_tools`` directory...TODO
- The ``apps`` directory...TODO

Many of the directories contain directories of unit and integration
tests, directories of Rose metadata and Rose Stem tasks for supporting
running of tests on other host machines, and Makefiles for building
applications, building and running of local tests, or running Rose
suites.

Overview of Key Aspects
=======================

Applications
------------

An LFRic application is an application that depends on the LFRic core
infrastructure, including, potentially, :ref:`LFRic components <section
components overview>`. The application will follow the ``PSyKAl``
architecture, described in the next section, to give the required
separation of concerns between the algorithmic descriptions of
high-level scientific processes and the lower-level computation that
implements these processes.

Applications can call model components, for example, to run an
atmosphere model. Applications can also be utilities; for example, to
regrid data passed from one model configuration to another.

The LFRic core repository includes just a few applications that are
used for training, testing or development of new features. The
Momentum atmosphere model exists as an application called `lfric_atm`
in a separate FCM repository alongside several other applications,
libraries of science code and science :ref:`components <section components
overview>`

.. _section psykal overview:

The PSyKAl Architecture
-----------------------

The LFRic infrastructure has been written to support a separation of
concerns between scientific and technical code, aiming to support
**portable performance**, where scientific models can be ported to new
machines and configurations with a small amount of work. A key part of
delivering such capability is the use of a tool called **PSyclone** to
generate technical code to manage parallelism and optimisation of the
scientific codes. PSyclone can enable the same science code to be
built and run with platform-specific technical code that can enable an
application ported to a new platform more easily, or may allow
different optimisation choices to be applied to the same model running
on the same platform but with different configuration options.

To delivering such capability PSyclone and the LFRic infrastructure
have been developed together. The LFRic Infrastructure is designed to
support scientific models written according to a software architecture
called **PSyKAl** that stands for the three-layers of a science model
design: High-level scientific Algorithms at the top, low-level
scientific Kernels at the bottom and the PSyclone-generated Parallel
Systems, or **PSy layer**, code in the middle.

The LFRic Infrastructure provides a parallel data model which at its
core supports storage and manipulation of model fields governed by
strict rules. The algorithms represent the basic mathematics of the
model operating on full fields only. The kernels implement individual
operations on a chunk of the model field data. Each kernel must be
written to the PSyKAl standard, which requires comprehensive metadata
descriptions of its input variables and of its basic structure. The
metadata enables PSyclone to generate the appropriate PSy layer code
that breaks up the field into appropriate chunks and calls the kernel
with each chunk.

As well as generating code to call kernels, PSyclone can use the
metadata to identify dependencies between kernels, can compute
appropriate loop bounds that take account of halos used in distributed
memory deployments, can insert halo swap calls at appropriate
points, and can apply shared memory optimisations, for example, to
parallelise the kernel calls.

A more :ref:`comprehensive description <section psykal and datamodel>` can
be found that details how the LFRic infrastructure and PSyclone work
together.

Application data structures
---------------------------

The PSyKAl architecture describes how individual fields are used
within a scientific model. To support more complex applications such
as the Momentum atmosphere model, other data structures exist:

- Fields can be stored in linked-list structures called field
  collections, meaning that complex algorithms that take a lot of
  fields (including optional fields) can have more compact argument
  lists of field collections rather than very long lists of fields.
- The ``modeldb`` data structure aims to fully-encapsulate all the
  data required to run a model such that it would be possible to run
  two such models in the same executable without conflict.
- Certain model-independent data can be held as singletons and so
  shared between instances of a model. Data structures exist that
  support holding of mesh-dependent data that can safely be shared
  between models that are using the same mesh.
- Support for ``runtime constants`` (in development) ...TODO

.. _section components overview:

Components
----------

The components directory contains packages of code delivering
self-contained infrastructure support or libraries of
functionality. Components may be dependent on LFRic infrastructure,
but LFRic infrastructure should never depend on a
component. Components may call out to other independent libraries.

Some examples of components held in LFRic core can illustrate their
role:

- The ``lfric-xios`` component provides an API that can read and write
  LFRic data using the XIOS package. It has been developed to minimise
  the dependency between LFRic applications and the IO infrastructure,
  assuming that in the future, support for other IO packages may be
  required.
- The ``inventory`` component is a self-contained component which
  supports applications that need to create key-value pair datasets
  based on a bespoke key.
- The ``driver`` component provides data structures and code that can
  be used to construct applications in a common way, for example to
  impose the `init`, `step`, `finalise` paradigm that is common for
  numerical modelling systems. The benefits of using a particular part
  of the driver component can be reduced develpment and maintenance
  cost. However, it is possible to pick and choose useful parts of the
  driver component; bespoke driver code can be used for application
  requirements that are not supported.

An application repository may choose to create components for holding
commonly-used science algorithms or interface code for calling out to
external science libraries.
