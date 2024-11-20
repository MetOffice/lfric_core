.. ------------------------------------------------------------------------------
     (c) Crown copyright 2023 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------
.. _software dependencies:

Software dependencies of LFRic
==============================

LFRic Applications
------------------

Given that most users of the LFRic core code are running applications such as
``lfric_atm`` that are stored in the separate "LFRic applications" repository,
`lfric_apps <https://code.metoffice.gov.uk/trac/lfric_apps/wiki>`_, it is worth
describing the relation between the two repositories.

Currently, the development of the two code bases is done
hand-in-hand. Therefore, certain revisions of the core code are tagged with the
version number of the relevant ``lfric_apps`` release. For example, revision
51381 is tagged ``core1.2`` and works with the 1.2 LFRic apps release.

Note, also, that any given revision of ``lfric_apps`` includes a
``dependencies.sh`` file in its top-level directory which references a specific
revision of the LFRic core code against which it will be built.

Compiler versions
-----------------

LFRic uses several features of modern Fortran, including object-oriented
features introduced at F2003 and later. Not all compilers have correctly
implemented all features, and therefore there may be problems when compiling and
running LFRic applications with some compilers.

The following compilers are routinely tested at the Met Office:

 * The Gnu compiler (versions 11.2 and 12.1)
 * The Intel ifort compiler (version 19.0)

Software Stack
--------------

To build and run typical LFRic applications, the following software will be
required. The numbers in parenthesis identify versions in use at the Met Office
for the revision of LFRic tagged ``core1.2``.

Common software which may already be installed on some HPC and research
platforms:

 * Python version 3 (3.7)
 * HDF5 (1.12.1)
 * NetCDF (4.8.1)
 * MPI (mpich 3.4.1)

More specialist software for developing, building and running LFRic
applications:

 * FCM (2021.05.0). LFRic code is held in a Subversion repository. `FCM
   <https://metomi.github.io/fcm/doc/user_guide>`_ is an application that wraps
   Subversion commands to help impose standard development workflows for LFRic
   development.
 * PSyclone (2.5.0), a code generation library used by LFRic for generating
   portable performance code. The `PSyclone documentation
   <https://psyclone.readthedocs.io/en/stable/>`_ list its own software
   dependencies, which include some Python packages and the following Fortran
   parser.
 * fparser (0.1.4), a Fortran parser used by PSyclone.
 * YAXT 0.10.0), an `MPI wrapper
   <https://swprojects.dkrz.de/redmine/projects/yaxt>`_ which supports MPI data
   exchange in LFRic application.
 * XIOS (r2252.2) an `IO server library
   <https://forge.ipsl.jussieu.fr/ioserver>`_ to support input and output of
   data to UGRID NetCDF files.
 * blitz (1.0.2), a `support library <https://github.com/blitzpp/blitz>`_
   required by XIOS.
 * rose-picker (2.0.0) available from the `GPL-utils project
   <https://code.metoffice.gov.uk/trac/lfric/browser/GPL-utilities>`_ in the
   core LFRic repository, supports parsing Rose metadata when generating
   namelist-reading code.

Additional specialist software essential for running the LFRic infrastructure
and application tests:

 * `Rose (2.3.1) <https://metomi.github.io/rose/doc/html/index.html>`_ and `Cylc
   (8+) <https://cylc.github.io/>`_ for running the full test-suite that
   includes application configurations, integration tests, unit tests, style
   checker and metadata validation checks.
 * PFUnit (3.3.3) A `Fortran unit testing
   <https://github.com/Goddard-Fortran-Ecosystem/pFUnit>`_ framework.
 * stylist (0.4.1a): A `code style-checker
   <https://github.com/MetOffice/stylist>`_.
 * plantuml (1.2021.7): the LFRic repository holds formal descriptions of key
   LFRic classes using a format that can be rendered into UML diagrams by
   `plantuml <https://plantuml.com>`_.


Future releases
---------------

About three releases of ``lfric_apps`` are planned to take place each year.  The
LFRic core code will be appropriately tagged against each such release.

Late 2024, a version 2.0 release is planned:
  * The primary aim of the release is to upgrade to PSyclone version 3.
  * This release is expected to include an upgrade the Rose testing of LFRic
    core to Cylc 8 to match ``lfric_apps``.
  * The version of pfUnit will hopefully be upgraded to version 4.
