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
``lfric_apps``, it is worth describing the relation between the two
repositories.

Currently, the development of the two code bases is done
hand-in-hand. Therefore, certain revisions of the core code are tagged with the
version number of the relevant ``lfric_apps`` release. For example, revision
50658 is tagged ``core1.1`` and works with the 1.1 LFRic apps release.

Note, also, that any given revision of ``lfric_apps`` includes a
``dependencies.sh`` file which references a specific revision of the LFRic core
code against which it was tested.

Compiler versions
-----------------

LFRic uses several features of modern Fortran including object oriented
features. Not all compilers have correctly implemented all features. The
following compilers are routinely tested at the Met Office:

 * The Gnu compiler (versions 11.2 and 12.1)
 * The Intel ifort compiler (version 19.0)

Software Stack
--------------

To build and run typical LFRic applications, the following software will be
required. The numbers in parenthesis identify versions in use at the Met Office
for the revision tagged ``core1.2``.

Common software which may already be installed on some HPC and research
platforms:

 * Python version 3 (3.7)
 * HDF5 (1.12.1)
 * NetCDF (4.8.1)
 * MPI (mpich 3.4.1)

More specialist software for building and running LFRic applications:

 * FCM. LFRic code is held in a Subversion repository. FCM is an application
   that wraps Subversion commands to help impose standard development workflows
   for LFRic development. https://metomi.github.io/fcm/doc/user_guide/
 * PSyclone (2.5.0), a code generation library used by LFRic for generating
   portable performance code. The PSyclone documentation
   https://psyclone.readthedocs.io/en/stable/ list its own software
   dependencies, which include some Python packages and the following Fortran
   parser.
 * fparser (0.1.4), a Fortran parser used by PSyclone.
 * YAXT 0.10.0), a library which supports MPI data exchange in LFRic
   applications: https://swprojects.dkrz.de/redmine/projects/yaxt
 * XIOS (r2252.2) to support input and output of data to UGRID NetCDF files.
 * blitz (1.0.2), required by XIOS.
 * rose-picker (2.0.0) available from the LFRic repository.

Additional specialist software essential for running the LFRic infrastructure
and application tests:

 * Rose and Cylc used for running the full test-suite that includes
   application configurations, integration tests, unit tests, style checker and
   metadata validation checks.
 * PFUnit (3.3.3) A Fortran unit testing framework.
 * stylist (0.4.1a): A code style-checker.
 * plantuml: used to generate UML diagrams from the UML descriptions that are
   maintained.

Future releases
---------------

It is expected that 3 releases of ``lfric_apps`` will take place each year, and
the LFRic core code will be appropriately tagged against each such release.

Late 2024, a version 2.0 release is planned:
  * The primary aim of the release is to upgrade to PSyclone version 3.
  * This release is expected to include an upgrade the Rose testing of LFRic
    core to Cylc 8 to match ``lfric_apps``.
  * The version of pfUnit will hopefully be upgraded to version 4.
