.. ------------------------------------------------------------------------------
     (c) Crown copyright 2023 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _build_and_run:

Quick Start Guide for LFRic
===========================

The LFRic repository includes some small applications for testing and training
purposes. In particular, the ``simple_diffusion`` app is used in LFRic training
courses and the ``skeleton`` app is designed to be an extremely simple and
minimal LFRic app.

A separate ``lfric_apps`` repository holds key science applications, such as
``lfric_atm``. The ``lfric_apps`` repository has a different set of instructions
for building and testing apps, but largely they do the same thing.

Check-out a working copy
========================

To checkout a working copy of the code to a new directory, named 'trunk' in this
example, run this command:

.. code-block::

  fcm co https://code.metoffice.gov.uk/svn/lfric/LFRic/trunk trunk

.. topic:: FCM Keywords

   Users of FCM may take advantage of its support for keywords to shorten this.

   Met Office developers should find they can use a site-wide keyword ``lfric.x``

   .. code-block::

      fcm co fcm:lfric.x-tr lfric_core_trunk

   Those without these site-wide keywords can set up their own locally. In
   ``~/.metomi/fcm/keyword.cfg`` just add the following lines:

   .. code-block::

     # LFRic repository
     location{primary}[lfric] = https://code.metoffice.gov.uk/svn/lfric/LFRic

   You may now use the shortened URL:

   .. code-block::

     fcm co fcm:lfric-tr

To create a branch and check it out, run the following:

.. code-block::

   fcm bc MyBranchName fcm:lfric.x-tr

The command will create a branch with your chosen name prefixed by the revision
number of the head of trunk. Assuming the current trunk revision is ``1234``, to
check-out the branch run the following, which checks the branch out to a new
directory that is named after the branch if the optional working copy name is
not supplied.

.. code-block::

   fcm co r1234_MyBranchName [working_copy_name]

Building an application
=======================

The applications can be found in directories within the ``application``
directory. An application can be built by running the ``make`` command from
within its directory.

.. code-block::

  cd r1234_MyBranchName/applications/simple_diffusion
  make

.. warning::

   Instructions for building an application from the ``lfric_apps`` repository,
   such as the ``lfric_atm`` atmosphere model are different. This is because
   these applications need to include steps to import external code, including
   the LFRic core code. Refer to the ``lfric_apps`` documentation for
   instructions.

The ``make`` command will build the ``simple_diffusion`` executable and place it
in the ``simple_diffusion/bin`` directory. It will also build and run any
integration and unit tests that the application has. The ``make`` command has a
number of optional arguments:

  * ``make build`` builds just the application executable.
  * ``make unit-tests`` builds and runs any the unit tests.
  * ``make integration-tests`` builds and runs any integration tests.

The ``Makefile`` has a number of optional variable settings that can be
overridden. To see these, look in the ``Makefile``. More than one option can be
supplied in a given invokation of ``make``.

For example,

 * ``make PROFILE=full-debug`` applies the set of ``full-debug`` compile flags
   rather than the default ``fast-debug`` settings.
 * ``make PSYCLONE_TRANSFORMATION=my_transforms`` applies the ``my_transforms``
   transformations to PSyclone-generated code instead of the default settings
   (which are, typically, either ``minimum`` or ``none``). The ``my_transforms``
   transformation scripts must be defined in the
   ``applications/[app-name]/optimisations`` directory. See PSyclone
   documentation for an explanation.
 * ``make VERBOSE=1``: The verbose option causes output of information
   from the dependency analyser, the compile command for each compilation
   process, and timings of these processes.

In addition to the ``bin`` directory that is created to hold the application
executable, the build process creates a ``working`` directory to hold the
products of the build system process.

Running ``make clean`` will remove the ``working`` and ``bin`` directories.

For most applications, it is common practice for applications in the ``lfric``
and ``lfric_apps`` repository to include a test configuration of the application
in an ``example`` directory within the application directory. Having built the
application it should be possible to ``cd`` to the ``example`` directory and run
the application, as follows:

.. code-block::

   ../bin/simple_diffusion configuration.nml

The ``configuration.nml`` contains a set of namelist configurations. Depending
on the application, the ``example`` directory may contain other files required
to run it such as a mesh definition file, a start dump or other input files.

Running the Cylc test suite
---------------------------

The test suite requires Cylc and Rose to run. Once they are installed and run,
tests are run using the ``rose stem`` command to select the group of tests to
run and ``cylc`` to initiate the suite. For example, to run the developer tests
for all the applications and components within an LFRic working copy run the
following from the top-level of the working copy:

.. code-block::

   rose stem --group=developer
   cylc play <working_copy_name>
