.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _build_and_run:

Quick Start Guide for LFRic
===========================

The LFRic repository includes some small applications for testing and training
purposes. In particular, the :ref:`simple_diffusion <simple diffusion
application>` app is used in LFRic training courses and the :ref:`skeleton
<skeleton application>` app is designed to be an extremely simple and minimal
LFRic app.

A separate ``lfric_apps`` repository holds key science applications, such as
``lfric_atm``. The lfric_apps repository has a different set of instructions for
building and testing apps, but largely they do the same thing.

Obtain the code
---------------

Obtain the code from the Github repository and install it on your machine.

The examples below assume you have cloned the repository or your fork of the
repository and placed it into a directory called ``lfric_core``.

Building an application
-----------------------

The applications can be found in directories within the ``application``
directory. An application can be built by running the ``make`` command from
within its directory.

.. code-block::

  cd lfric_core/applications/simple_diffusion
  make

The ``make`` command will build the simple_diffusion executable and place it in
the ``simple_diffusion/bin`` directory. It will also build and run any
integration and unit tests that the application has.

The make command uses the Makefile in the same directory as the application. The
Makefile has a number of optional arguments:

  * ``make build`` builds just the application executable.
  * ``make unit-tests`` builds and runs any the unit tests.
  * ``make integration-tests`` builds and runs any integration tests.

.. warning::

   Instructions for command-line building of an application from the
   lfric_apps repository, such as the lfric_atm atmosphere model, are
   different:

   https://metoffice.github.io/lfric_apps/developer_guide/local_builds.html

   The method is different because it includes steps to import external code,
   including the LFRic core code.

The Makefile has a number of optional variable settings that can be
overridden. To see these, look in the file. More than one option can be supplied
in a given invokation of ``make``.

A few of the more commonly used variables are:

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
products of an executable build and a ``test`` directory to hold products of a
build of the unit and integration tests.

If you are developing a change and testing builds within a branch, from time to
time you will want to commit changes to the upstream repository. A
``.gitignore`` file, found at the top-level of the directory tree, should
prevent you inadvertently including build artefacts in your commit if you ever
run a ``git add .`` command. But do use ``git status`` to be sure of what your
changeset includes. Note that it is always safer to use ``git add <filename>`` on the specific files you want to add/change.

Running ``make clean`` will remove the working, test and bin directories.

After building an application from the command line, it can be useful to do a
quick test to ensure it can run. Most applications in the lfric and lfric_apps
repository hold a simple example configuration in their ``example`` directory
(which is also included in the ``.gitignore`` file, so if you want to add or change files in the ``example`` directory, you will need to use ``git add -f ...``).

After building, go into to the example directory and run the application, as
follows:

.. code-block::

   ../bin/simple_diffusion configuration.nml

The ``configuration.nml`` file contains a set of namelist
configurations. Depending on the application, the ``example`` directory may
contain other files required to run it such as a mesh definition file, a start
dump or other input files.

Running the Cylc test suite
---------------------------

To run the test suites, Cylc and Rose need to be installed.

The ``rose stem`` command selects the group of tests to run. Then ``cylc play``
starts the suite running. For example, to run the developer tests for all the
applications and components within an LFRic working copy, run the following from
the top-level of the working copy:

.. code-block::

   rose stem --group=developer
   cylc play <working_copy_name>

Building the documentation
--------------------------

Documentation for the code exists as RST files in the documentation
directory. The directory includes a Makefile that runs Sphinx (currently using
v8.1.0 with the PyData Sphinx Theme 0.16.1) to generate html web pages. To build
the documentation and then view it with the Firefox browser:

.. code-block::

   cd lfric_core/documentation
   make html
   firefox build/html/index.html
