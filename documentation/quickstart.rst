Quick Start Guide for Dynamo
============================

So you want to play with Dynamo? These instructions should get you up and
running quickly and easily.

.. contents:: Table of Contents

.. note::
   The canonical version of this document is held as reStructured text in
   the repository at `source:Dynamo/trunk/documentation/quickstart.rst`:trac:.
   Any changes in a branch which render this document inaccurate should also
   include updates to the version of this document on that branch. The version
   displayed on the wiki is generated from the head of trunk.

Build Environment
-----------------

Within the Met Office (and on certain external machines) we make use of the 
`Environment Modules system <http://modules.sourceforge.net/>`_ to manage our
development environment. If you do not have access to this skip to the
`Outwith the Met Office`_ section for more details.

Met Office Dynamo Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a terminal window you need to source the module setup script and load the
Dynamo environment to set up libraries and compiler:

 1. Source the setup script
 
      +--------------------+----------------------------------------------+
      | Met Office Desktop | ``. /data/users/lfric/modules/setup``        |
      +--------------------+----------------------------------------------+
      | Met Office XC-40   | ``. /data/d03/lfric/modules/setup``          |
      +--------------------+----------------------------------------------+
      | MONSooN            | ``. /projects/umadmin/mhambl/modules/setup`` |
      +--------------------+----------------------------------------------+
      | ARCHER             | ``. /fs2/n02/n02/mhambley/modules/setup``    |
      +--------------------+----------------------------------------------+

 #. Load the modules

      +--------------------+------------------------------------------------+
      | Met Office Desktop | ``module load environment/dynamo/intel``       |
      |                    +------------------------------------------------+
      |                    | ``module load environment/dynamo/gnu``         |
      +--------------------+------------------------------------------------+
      | Met Office XC-40   | ``module load meto-environment/dynamo/cce``    |
      |                    +------------------------------------------------+
      |                    | ``module load meto-environment/dynamo/intel``  |
      +--------------------+------------------------------------------------+
      | MONSooN            | ``module load meto-environment/dynamo/intel``  |
      +--------------------+------------------------------------------------+
      | ARCHER             | ``module load lfric-environment/dynamo/intel`` |
      +--------------------+------------------------------------------------+

Render your own video using the "gource" documentation target is only available
on the Met Office desktop. You will need to do the following::

  module load base-environment/gource

Outwith the Met Office
~~~~~~~~~~~~~~~~~~~~~~

The "extra" directory contains a script which will construct a self-contained
development environment on a Debian Linux system. It is run as a script to
build the environment and then sourced to set up the required environment
variables. It is not guaranteed to work although it is tested. If nothing else
it should codify the requirements outlined below.

The build system makes use of environment variables to understand the system on
which it finds itself.

+---------------+----------------------+------------------------------------------------------------------------------+
| Variable Name | Tool                 | Notes                                                                        |
+===============+======================+==============================================================================+
| FC            | Fortran Compiler     | Supported compilers: Intel, GNU                                              |
+---------------+----------------------+------------------------------------------------------------------------------+
| FPP           | Fortran preprocessor | This will be ``cpp -traditional-cpp`` for GNU Fortran and ``fpp`` for Intel. |
+---------------+----------------------+------------------------------------------------------------------------------+
| LDMPI         | MPI linker           | The automatic choice used if this is not defined is usually correct.         |
+---------------+----------------------+------------------------------------------------------------------------------+

The build system expects paths to libraries to be presented in compiler
environment variables.

+---------------+-----------------------------------+--------------------------+--------------------------------------------------------+
| Variable Name | Prepend with                      | Purpose                  | Example                                                |
+===============+===================================+==========================+========================================================+
| ``FFLAGS``    | ``-I<path to module directory>``  | Point to ``.mod`` files. | ``-I/opt/esmf/mod/mod0/Linux.ifort.64.mpich2.default`` |
+---------------+-----------------------------------+--------------------------+--------------------------------------------------------+
| ``LDFLAGS``   | ``-L<path to library directory>`` | Point to ``.so`` files.  | ``-L/opt/esmf/lib/lib0/Linux.ifort.64.mpich2.default`` |
+---------------+-----------------------------------+--------------------------+--------------------------------------------------------+

Since Dynamo is dynamically linked the directory in which the ``.so`` files are
found must be included in the ``LD_LIBRARY_PATH`` variable.

Normally the build system suppresses much output in order to reduce clutter.
While getting an initial build to work some of the suppressed output can be
invaluable in making sure the correct paths are being searched. To get this
output use ``make VERBOSE=1``

Currently the following compilers are supported:

+---------------+-----------------+
| Compiler      | Minimum Version |
+===============+=================+
| Intel Fortran | 15.0.0          |
+---------------+-----------------+
| GNU Fortran   | 4.9.2           |
+---------------+-----------------+

The requirements for building Dynamo are:

+---------+----------+-----------------+---------------------------------------+------------------------------------------+
| Package | Version  | Dependency Type | Purpose                               | Pointed to by...                         |
+=========+==========+=================+=======================================+==========================================+
| CMake   | 2.8.12.2 | Build           | pFUnit builds are prepared using this | PATH                                     |
+---------+----------+-----------------+---------------------------------------+------------------------------------------+
| ESMF    | 7.0.0    | Build, Run      | Framework & Infrastructure library    | FFLAGS, LDFLAGS, LD_LIBRARY_PATH         |
+---------+----------+-----------------+---------------------------------------+------------------------------------------+
| GMake   | 3.81     | Build           | GNU's version of "make"               | PATH                                     |
+---------+----------+-----------------+---------------------------------------+------------------------------------------+
| HDF5    | 1.8.12   | Build, Run      | Parallel filesystem used by NetCDF    | FFLAGS, LDFLAGS, LD_LIBRARY_PATH         |
+---------+----------+-----------------+---------------------------------------+------------------------------------------+
| Jinja2  | 2.6      | Build           | Required by configuration handler     | PYTHONPATH                               |
+---------+----------+-----------------+---------------------------------------+------------------------------------------+
| MPICH   | 3.1      | Build, Run      | MPI library used by NetCDF            | FFLAGS, LDFLAGS, PATH, LD_LIBRARY_PATH   |
+---------+----------+-----------------+---------------------------------------+------------------------------------------+
| NetCDF  | 4.3.1.1  | Build, Run      | File I/O                              | FFLAGS, LDFLAGS, LD_LIBRARY_PATH         |
+---------+----------+-----------------+---------------------------------------+------------------------------------------+
| Python  | 2.7      | Build           | Required by the dependency analyser   | PATH                                     |
+---------+----------+-----------------+---------------------------------------+------------------------------------------+
| libxml2 | 2.9.1    | Build, Run      | Required by pFUnit                    | CPPFLAGS, LDFLAGS, PATH, LD_LIBRARY_PATH |
+---------+----------+-----------------+---------------------------------------+------------------------------------------+
| Zlib    | 1.2.8    | Build, Run      | Compression library used by NetCDF    | CPPFLAGS, LDFLAGS, LD_LIBRARY_PATH       |
+---------+----------+-----------------+---------------------------------------+------------------------------------------+

To build the documentation you will need:

+---------+----------+-----------------------------------------------------------+
| Package | Version  | Purpose                                                   |
+=========+==========+===========================================================+
| Doxygen | 1.8.7    | API documentation                                         |
+---------+----------+-----------------------------------------------------------+
| LyX     | 2.0.7    | Some of the developer documentation is prepared using LyX |
+---------+----------+-----------------------------------------------------------+

Gource is just a bit of fun so it gets a separate section:

+-----------+---------+---------------------------------+-----------------------------------------------------+
| Package   | Version | Note                            | Purpose                                             |
+===========+=========+=================================+=====================================================+
| Boost     | 1.55.0  |                                 | C++ template library used by Gource                 |
+-----------+---------+---------------------------------+-----------------------------------------------------+
| GLEW      | 1.10.0  |                                 | OpenGL Extension Wrangler used by SDL               |
+-----------+---------+---------------------------------+-----------------------------------------------------+
| SDL       | 2.0.3   |                                 | Simple Directmedia Layer used by Gource             |
+-----------+---------+---------------------------------+-----------------------------------------------------+
| SDL Image | 2.0.0   |                                 | Image loading component of SDl, used by Gource      |
+-----------+---------+---------------------------------+-----------------------------------------------------+
| FFMPEG    | 2.1.4   |                                 | Video transcoding tool, used to create video output |
+-----------+---------+---------------------------------+-----------------------------------------------------+
| Gource    | 0.40    | http://code.google.com/p/gource | Repository visualisation tool                       |
+-----------+---------+---------------------------------+-----------------------------------------------------+

Checkout Dynamo
---------------

To checkout a working copy of the code to a new directory, named 'trunk' in
this example, run one of the following commands::

  svn co https://code.metoffice.gov.uk/svn/lfric/Dynamo/trunk trunk
  svn co $REPO/trunk trunk # Alternative once environment/dynamo/base module is loaded

or::

  fcm co https://code.metoffice.gov.uk/svn/lfric/Dynamo/trunk trunk

FCM users, particularly developers, may want to set up a suitable keyword to
make this easier. In ``~/.metomi/fcm/keyword.cfg`` add the following lines::

  # Gung-Ho / LFRic : Dynamo project
  location{primary}[Dynamo] = https://code.metoffice.gov.uk/svn/lfric/Dynamo

This would allow the URL for the trunk in an fcm command to be specified as
``fcm:Dynamo_tr`` and the root of the path to all branches to be
``fcm:Dynamo_br``

To work in the working copy, change directory to the working copy you've just
checked out. ('trunk' in our example)::

  cd trunk


Running Make
------------

The current build system used for Dynamo is "Make". Makefiles have been set up
such that by running ``make`` in top level of a working copy will build the
executable, build the unit tests and execute the unit tests::

  make

It must be GNU make and a sufficiently modern version must be used. The build
system will check and complain if the version isn't up to scratch but will just
fail if a non-GNU version is used.

Three targets are offered:

+------------+---------------------------------------+---------+
| Target     | Result                                | Default |
+============+=======================================+=========+
| full-debug | No optimisation and run-time checking |         |
+------------+---------------------------------------+---------+
| fast-debug | Safe optimisation only                | Yes     |
+------------+---------------------------------------+---------+
| production | Risky optimisation                    |         |
+------------+---------------------------------------+---------+

All targets include debug symbols into the executable code.

In order for PSyclone to select the correct optimisation script it must know
the platform you are building on. This is achieved by setting the
`DYNAMO_BUILD_TARGET` environment variable to a single target platform in the
same form as used for `DYNAMO_TEST_SUITE_TARGETS`, described below.

This is done for Met Office users by the LFRic module system.

Run ``make clean`` to remove all compiled application and unit test output
should make fail to perform a rebuild and try again. If you want to delete the
compiled pFUnit framework as well use ``make clean-all``

The binary for Dynamo can be found in the ``bin`` directory in the top level of
your working copy.

Makefiles exist in several subdirectories so individual parts can be built by
running make in those subdirectories:

 * ``documentation/Makefile`` - to build the documentation
 * ``src/dynamo/Makefile`` - to just build the Dynamo binary
 * ``src/test/Makefile`` - to build and run the unit tests

.. NOTE::
   At present the unit tests are only know to compile with GNU Fortran, Intel
   Fortran and CMake. If not using these running ``make`` in the top level of
   the working copy will produce errors which can be ignored if you're not
   interested in the unit tests.

Running The Dynamo Binary
-------------------------

The binary for Dynamo can be found in the ``bin`` directory in the top level of
your working copy. It expects to find a grid file in the current working
directory. The quickest way to execute it is as follows::

  cd data
  ../bin/dynamo

Explicitly Running The Unit Tests
---------------------------------

The unit tests can be built and run by from the  top level directory ('trunk')
with the following::

  make test

Running the Test Suite
----------------------

The test suite requires `Cylc <https://github.com/cylc/cylc>`_ and `Rose
<https://github.com/metomi/rose>`_ to run. It makes use of the "Rose Stem"
test launch tool.

Once they are installed and working, set the environment variable
`DYNAMO_TEST_SUITE_TARGETS` to a space separated list of target platforms taken
from `rose-stem/opt/rose-suite-<target>.conf`.

The `test-suite` target may then be used thus::

  export DYNAMO_TEST_SUITE_TARGETS=place-machine
  make test-suite

The Met Office environment module sets up this variable for local platforms.

Using the command ``rose stem`` will launch the suite against the Met Office
SPICE server farm. This is useful during development.

Building The Documentation
--------------------------

The Documentation can be built and run by from the top level directory with the
following::

  make docs

To view Doxygen documentation for the !GungHo science code point a browser at::
``documentation/api/gungho/html/index.html``

For the software infrastructure, point a browser at::
``documentation/api/infrastructure/html/index.html``

A pdf of the scientific formulation is found at::
``documentation/formulation/dynamo_formulation.pdf``

A pdf that provides an introduction to the data model is found at::
``documentation/datamodel/dynamo_datamodel.pdf``


Possible Issues
---------------

Slow builds
~~~~~~~~~~~

You may find that builds stall around dependency analysis. If this is the case
refer to `Dynamo/BuildSystem#RelocateBuildArtifacts`:trac:.
