.. ------------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------
.. _guide_to_documentation:

Guide to the documentation
==========================

The LFRic Core documentation is split into a number of main sections, this
:ref:`Getting Started<getting_started_index>` section contains information on
the dependencies required to run any LFRic based applications including those
found in LFRic Core. The :ref:`How to use it<how_to_use_it_index>` section
details how LFRic Core can be used and is aimed at users developing
scientific applications. The :ref:`How it works<how_it_works_index>` section
contains details of the implementation of LFRic Core, this is aimed at
experienced developers and those doing development of LFRic Core itself.
:ref:`How to contribute<how_to_contribute_index>` contains guidance for those
doing development within LFRic Core itself, particularly the procedures and
form development should take.

For initial orientation, a :ref:`quick overview<repository_contents>` of the
main contents of the LFRic core repository is given. As noted, the Momentum
atmosphere model, and related applications, are developed in separate
repositories. The LFRic core repository does include some small LFRic
applications for training purposes or for developing and testing particular
technical capabilities.

Before giving an overview of the core infrastructure, an overview of the
:ref:`structure of a typical LFRic application<application structure>` is
given. It briefly references several LFRic core capabilities. These
capabilities are then described in later sections.

Any developer of the LFRic core or an LFRic application should have a good
understanding of the underlying principles behind LFRic and the core data
model, and the scientific model architecture known as
**PSyKAl**. The :ref:`LFRic data model and PSyclone <psykal and datamodel>`
documentation describes key aspects of LFRic and of PSyclone, the code
autogeneration tool that LFRic applications depend upon.

The :ref:`Application Documentation section<core applications>` provides links
to documentation for each application developed within the LFRic core
repository, describing the role of the application and including pointers to
the features of the LFRic core that it depends upon or tests.

The :ref:`Meshes section<section mesh generation>` describes the LFRic mesh
generator and LFRic meshes, including discussion of mesh partitioning, mesh
hierarchies and mesh maps, the LFRic mesh generator.

The :ref:`Components section <components>` describes components which are code
libraries delivering specific capabilities required by some LFRic applications.

The :ref:`build and test system section <build and test>` describes the build
and test system, and includes descriptions of the tools that underpin it.

The :ref:`Technical Articles section <technical articles>`
includes articles on several topics including the distributed memory strategy
and implementation, the clock and calendar.

A:ref:`Glossary of Terms <glossary_of_terms>` defines commonly used terms and
is available from the right-hand side panel of each page.

Finally, detailed API documentation is given for the LFRic infrastructure and
each LFRic component.
