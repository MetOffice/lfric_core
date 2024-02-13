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

For initial orientation, a `quick overview <section repository
contents>` of the main contents of the LFRic core repository is
given. While LFRic core does include some small LFRic applications,
one should be aware that major applications, including the Momentum
atmosphere model, are developed in other repositories and will have
additional application-specific documentation.

Any developer of the LFRic core or an LFRic application should have a
good understanding of the underlying principles behind LFRic and the
core data model and model architecture. An `introduction <section
introduction to lfric infrastructure>` describes key aspects of LFRic
and of PSyclone, the code autogeneration tool that LFRic applications
depend upon.

While the `introduction <section introduction to lfric
infrastructure>` provides all that is needed for the simplest LFRic
application, more realistic model applications commonly need more
complex data structures to manage large numbers of fields, and control
flow and IO to manage the running and output of a numerical model. An
overview of aspects of an `LFRic model structure <section model
application structure>` introduces key aspects of LFRic core that
support these requirements.

Separate documentation exists for each application and component
developed within the LFRic core repository, describing the role of the
application and including pointers to the features of the LFRic core
that it depends upon or tests.

Some technical features required by applications have been implemented
as LFRic "components". Each of these is described.

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
infrastructure, including, potentially, `LFRic components <section
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
libraries of science code and science `components <section components
overview>`

The PSyKAl Architecture
-----------------------

The PSyKAl architecture describes the architecture of the core of a
scientific model...TODO: two to three summary paragraphs and the
picture.

Application data structures
---------------------------

The PSyKAl architecture describes how individual fields are used
within a scientific model. To support more complex applications such
as the Momentum atmosphere model, other data structures exist:

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
