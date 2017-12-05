Quick Start Guide for LFRic
===========================

So you want to play with LFRic models? These instructions should get you up and
running quickly and easily.

.. contents:: Table of Contents

.. note::
   The canonical version of this document is held as reStructured text in
   the repository at
   `source:LFRic/trunk/gungho/documentation/quickstart.rst`:trac:.
   Any changes in a branch which render this document inaccurate should also
   include updates to the version of this document on that branch. The version
   displayed on the wiki is generated from the head of trunk.

Build Environment
-----------------

LFRic uses very new versions of compilers and applications to get the
best support for Fortran 2003. Typically, these do not come as
standard and a bespoke environment needs to be set up.

For Met Office users and users of Monsoon and Archer refer to `DevelopmentEnvironment`:trac:
wiki:DevelopmentEnvironment

For external users or users of non-Met Office Linux installations,
refer to `LFRicTechnical/LFRicBuildEnvironment`:trac:

Normally the build system suppresses much output in order to reduce clutter.
While getting an initial build to work some of the suppressed output can be
invaluable in making sure the correct paths are being searched. To get this
output use::

  make VERBOSE=1

Checkout a Working Copy
-----------------------

To checkout a working copy of the code to a new directory, named 'trunk' in
this example, run this command::

  svn co https://code.metoffice.gov.uk/svn/lfric/LFRic/trunk trunk

Users of FCM may take advantage of its support for keywords to shorten this.

Met Office developers should find they can use a site-wide keyword "lfric.x"::

  fcm co fcm:lfric.x-tr

Those without these site-wide keywords can set up their own locally. In
``~/.metomi/fcm/keyword.cfg`` just add the following lines::

  # LFRic repository
  location{primary}[lfric] = https://code.metoffice.gov.uk/svn/lfric/LFRic

You may now use the shortened URL::

  fcm co fcm:lfric-tr

The ``-br`` suffix may be used to access the branches directory::

  fmc co fcm:lfric.x-br/dev/joebloggs/r1234_MyBranch

To work in a working copy just change into the directory::

  cd r1234_MyBranch


Running Make
------------

The current LFRic build system uses "Make". Everyting has been set up such
that running ``make`` in top level of a working copy will build the
executables and build and execute the unit tests::

  cd workingcopy
  make

It must be GNU make and a sufficiently modern version is needed. The build
system will check and complain if the version isn't up to scratch but will just
fail if a non-GNU version is used.

Three build profiles are offered:

+------------+---------------------------------------+---------+
| Profile    | Result                                | Default |
+============+=======================================+=========+
| full-debug | No optimisation and run-time checking |         |
+------------+---------------------------------------+---------+
| fast-debug | Safe optimisation only                | Yes     |
+------------+---------------------------------------+---------+
| production | Risky optimisation                    |         |
+------------+---------------------------------------+---------+

All profiles include debug symbols into the executable code.

Pass the ``PROFILE`` variable to make in order to select one of these profies::

  make PROFILE=production

In order for PSyclone to select the correct optimisation script it must know
the platform you are building on. This is achieved by setting the
``LFRIC_TARGET_PLATFORM`` environment variable to a single platform identifier,
as defined above.

Use ``make clean`` to remove all compiled application and unit test output.

LFRic model binaries may be found in the ``bin`` directory in the top level of
your working copy.

.. NOTE::
   At present the unit tests are only know to compile with Intel and GNU
   Fortrans. If you are not using these, running ``make`` in the top level of
   the working copy will produce errors which can be ignored if you're not
   interested in the unit tests.

Running The Gung Ho Binary
--------------------------

The binary for Gung Ho can be found in the ``bin`` directory in the top level
of your working copy. It expects to find a grid file in the current working
directory. The quickest way to execute it is as follows::

  cd gungho/data
  ../../bin/gungho

Explicitly Running The Unit Tests
---------------------------------

The unit tests can be built and run from the top level directory with the
following::

  make test

Running the Test Suite
----------------------

The test suite requires `Cylc <https://github.com/cylc/cylc>`_ and `Rose
<https://github.com/metomi/rose>`_ to run. It makes use of the "Rose Stem"
test launch tool.

Once they are installed and working, set the environment variable
`TEST_SUITE_TARGETS` to a space separated list of target platforms taken
from `rose-stem/opt/rose-suite-<target>.conf`.

The `test-suite` make target may then be used thus::

  export TEST_SUITE_TARGETS=place-machine
  make test-suite

The Met Office environment module sets up this variable for local platforms.

Using the command ``rose stem`` will launch the suite against the Met Office
SPICE server farm. This is useful during development.

Building The Documentation
--------------------------

The Documentation can be built and run from the top level directory with the
following::

  make docs

To view Doxygen documentation for the !GungHo science code point a browser at::
``documentation/gungho/api/index.html``

For the software infrastructure, point a browser at:
``documentation/infrastructure/api/index.html``

A pdf of the scientific formulation is found at:
``documentation/gungho/formulation.pdf``

A pdf that provides an introduction to the data model is found at:
``documentation/gungho/datamodel.pdf``

To build the documentation you will need:

+----------+----------+-------------------------------------------------------------+
| Package  | Version  | Purpose                                                     |
+==========+==========+=============================================================+
| Doxygen  | 1.8.12   | API documentation                                           |
+----------+----------+-------------------------------------------------------------+
| Inkscape | 0.47     | Used to convert vector graphic formats                      |
+----------+----------+-------------------------------------------------------------+
| LaTeX    | 2e       | Technical and science documentation is prepared using LaTeX |
+----------+----------+-------------------------------------------------------------+

Possible Issues
---------------

Slow builds
~~~~~~~~~~~

You may find that builds stall around dependency analysis. If this is the case
refer to the section titled "Relocate Build Artifacts" in `LFRicTechnical/BuildSystem`:trac:.
