The LFRic Build System
=======================

This is a quick introduction to the LFRic Subversion project's build system.

.. note::
   The canonical version of this document exists as a reStructured text file
   in the repository at
   `source:LFRic/trunk/infrastructure/documentation/buildsystem.rst`:trac:.
   Branches which render it inaccurate should update the version of this
   document on that branch. The version displayed in the wiki is generated
   from the head of trunk.

Platform Identifiers
--------------------

In a number of places it is necessary to uniquely identify a target system.
This is achieved using a platform identifier.

Platform identifiers take the form <site>-<machine>. For example the Met
Office's XC40 is referred to as meto-xc40 while the current MONSooN machine
is known as monsoon-xc40.

Targets
-------

The project is split into a number of sub-projects. Each of these has a
Makefile offering targets "build" and "test-suite". Targets "unit-tests",
"integration-tests" and "documentation" are offered where appropriate. Other
targets specific to a particular sub-project may also be provided.

There is also a top level Makefile. This provides targets which
operate on a selection of sub-projects. By default it operates on
``infrastructure``, ``mesh_tools`` and ``gungho`` but it may be
changed by setting the ``OPERATE_ON`` variable to a space-separated
list of required sub-projects.

Building
~~~~~~~~

To perform a full build of ``gungho`` simply change to the sub-project
directory in the working copy and issue a ``make``::

  cd r1234_MyWorkingCopy/gungho
  make

Three build profiles are offered: ``fast-debug``, ``full-debug`` and
``production``. They are specified using the ``PROFILE`` variable passed to
make. For instance, to produce a build with all debugging enabled simply use
``make PROFILE=full-debug``. The default is ``fast-debug`` if none is
specified.

All three profiles insert debug symbols to the executable which provides
source file and line numbers in backtraces. They have no run-time impact other
than to increase the size of the executable.

``full-debug`` turns optimisation off and run-time checking, of e.g. array
out-of-bounds, on.

For those who need faster execution but also want bit reproducibility
``fast-debug`` limits optimisation to a safe level and turns off run-time
checking.

``production`` increases the optimisation level but does not guarantee
bit reproducibility. This is due to potential run-time optimisations,
particularly in the order of shared and distributed memory processing. We would
expect inter-run reproducibility but it can not be guaranteed.

Relocate Build Artifacts
^^^^^^^^^^^^^^^^^^^^^^^^

By default all build related files end up in the directory ``working``
within the sub-project directory. This can be changed by exporting the
``WORKING_DIR`` environment variable. It should contain an absolute
path to a suitable temporary space in which to perform the compile.

On some systems (such as some HPC systems) parts of the file system
may be optimised for large files and may perform poorly with the many
serial accesses to small files involved in compiling. Such systems may
provide other file systems (including temporary file systems) which
can be pointed to so as to speed up compilation.

Verbose Building
^^^^^^^^^^^^^^^^

By default the build system will suppress as much output as possible to reduce
clutter. When resolving problems it may be useful to see that clutter. If you
want to see the actual commands being run by the build set the ``VERBOSE``
variable::

  make VERBOSE=1

Linking
^^^^^^^

By default the build system will link dynamically. If you want a statically
linked binary then you should pass the ``LINK_TYPE`` variable::

  make LINK_TYPE=static

Static linking is the default on Cray systems as they do not seem to play well
with dynamic linking.

Cleaning
~~~~~~~~

As with many build systems it is possible to ``make clean`` to delete all build
artefacts. These include working files and complete executable binaries.

Testing
~~~~~~~

Unit tests and integration tests will be built and run as part of a
normal build. 

Each sub-project includes a Rose Stem test suite which can run on a
range of compute platforms. The Rose Stem test suite is manually
invoked from the top-level directory using::

  make test-suite 

This command will invoke tests for each sub-project on each of the
platforms listed in the ``TEST_SUITE_TARGETS`` as a space-separated
list. Each test on each platform is run as a separate Rose Stem suite.

Configurations are held in ``rose-stem/opt``. Each filename has the form
``rose-suite-<platform id>.conf`` using the platform identifiers described
above.

For those using a Met Office module collection the core module will set this up
for you. e.g. On the desktop do::

  module load common-environment/lfric

By default rose stem tests are run for ``infrastructure``,
``mesh_tools`` and ``gungho`` but this may be changed by setting the
``OPERATE_ON`` variable to a space-separated list of required
sub-projects.

For further information an testing see `LFRicTechnical/Testing`:trac:.

Choosing a Compiler
-------------------

A number of compilers are supported by the build system. However
currently only Intel and Gnu can be used to both build and run
applications. The following lists the environment variable settings
required to select them:

Intel::
  FC=ifort      # Fortran compiler
  FPP=fpp       # Fortran preprocessor
  LDMPI=mpif90  # MPI Linker

Gnu::
  FC=gfortran                  # Fortran compiler
  FPP=cpp -traditional-cpp     # Fortran preprocessor
  LDMPI=mpif90                 # MPI Linker



Each compiler needs a different set of arguments. These are defined by the
LFRic build system in ``infrastructure/build/fortran/<compiler>.mk``. Each
file is named after the command used to launch the compiler.

A number of universal compiler options and build system variables are defined
first.

=========================  =================================================================================
Variable                   Purpose
=========================  =================================================================================
``<compiler>_VERSION``     Holds the version of the running compiler. This is held as a zero padded integer. e.g. 4.9.2 becomes 040902
``F_MOD_DESTINATION_ARG``  The argument to be used to tell the compiler where to put ``.mod`` files.
``OPENMP_ARG``             The argument passed to enable interpretation of OpenMP directives.
``DEPRULE_FLAGS``          Additional arguments to pass to the dependency generator script.
=========================  =================================================================================

Notice that the ``FFLAGS`` variable is not modified directly. Instead a number
of argument groups are defined.

=============================  ===========================================================================================================
Argument Group                 Purpose
=============================  ===========================================================================================================
``FFLAGS_COMPILER``            Control the way the compiler works. For instance mdoifying unexpected behavior.
``FFLAGS_NO_OPTIMISATION``     Turn off all optimisation.
``FFLAGS_SAFE_OPTIMISATION``   Enable only bit reproducibility safe optimisation.
``FFLAGS_RISKY_OPTIMISATION``  Enable optimisations which may affect bit reproducibility.
``FFLAGS_DEBUG``               Enable the generation of debugging symbols, tracebacks and such.
``FFLAGS_WARNINGS``            Control warning handling. For development we like to receive all warnings and have them act as errors.
``FFLAGS_INIT``                Turn on any declaration initialisation the compiler may support to fill declared memory with a known value.
``FFLAGS_RUNTIME``             Turn on any run-time checking the compiler may support. Array bound checking for instance.
=============================  ===========================================================================================================

The same is done with ``LDFLAGS``.

=========================  =======================================================================================================
Variable                   Purpose
=========================  =======================================================================================================
``LDFLAGS_COMPILER``       Some linkers need to be told to retain debug symbols. That could be done here.
``LD_COMPILER_LIBRARIES``  Specify any special libraries needed at link stage as a space separated list.
                           These are not filenames, drop the "lib" from the front.
=========================  =======================================================================================================


Automatic Dependency Analysis
-----------------------------

A knowledge of dependencies between source files allows only those things
affected by a change to be recompiled. Maintaining such a dependency list by
hand is tedious and error prone. This is why the build system will generate
and maintain dependency information automatically.

The dependency analyser scans Fortran source looking for "use" statements
to determine the prerequisites of a source file. It relies on source files
having the same name as the module they contain, which implies that they
contain only one module.

e.g. if the module is called "my_special_stuff_mod" then the file name
should be "my_special_stuff_mod.f90". If this is not the case these modules will
be rebuilt every time.

The current dependency analyser caches its findings in a dependency
database which allows reruns of make to run quicker. However, a re-run
will fail for some operations such as renaming of files. If you get a
failure, run ``make clean`` prior to rebuilding. This will force the
dependency database to be rebuilt.


External Libraries
------------------

If you want to make use of an external library then some simple edits to
``Makefile`` will be required.

There are two critical variables: ``IGNORE_DEPENDENCIES`` and
``EXTERNAL_*_LIBRARIES``.

``IGNORE_DEPENDENCIES`` is a space-separated list of modules which should be
ignored by the dependency analyser when discovered in ``use`` statements. This
prevents the dependency analysis trying to rebuild your library. For example,
the infrastructure makes use of the ESMF library: the name to add would be
"esmf" since we ``use esmf``.

``EXTERNAL_*_LIBRARIES`` is a space-separated list of library names to be passed
to the linker using "little l" arguments. Libraries which are available as a
".so" file should be listed in ``EXTERNAL_DYNAMIC_LIBRARIES`` while those only
available as a ".a" file are listed in ``EXTERNAL_STATIC_LIBRARIES``. Returning
to the ESMF example, the library file is ``libesmf.so`` so the string ``esmf``
would be added to ``EXTERNAL_DYNAMIC_LIBRARIES``.

In addition to adding knowledge of the library to the build system you have to
make it findable. If you are using a module from the Met Office module
collection then this is handled for you automatically. If you are not then you
will need to do a little additional work.

The path to the directory containing the ``mod`` files should be added to
``FFLAGS`` in "big I" notation. e.g. ``export FFLAGS="$FFLAGS -I/path/to/mods"``

The path to the directory containing the library files should be added to
``LDFLAGS`` in "big L" notation.
e.g. ``export LDFLAGS="$LDFLAGS -L/path/to/libs"``

This will get you compiling. To run the resulting executable you have to make
sure the run-time linker can find the library files. To achieve this modify
``LD_LIBRARY_PATH``.
e.g. ``export LD_LIBRARY_PATH=/path/to/libs:$LD_LIBRARY_PATH``


Automatically Generated Code
----------------------------

The LFRic infrastructure makes use of the PSyclone tool, created as part of the
GungHo project, to automatically generate some of the code. PSyclone is
designed to generate the Parallel Systems layer (PSy-layer) code, but it also
modifies the algorithm code as written by science developers.

Algorithm routines in ``gungho/source/algorithm`` are written with a ``.x90``
extension. From these both compilable Fortran source for the algorithm and
corresponding PSy-layer are generated. Normally this process is handled for
you by the build system.

Where PSyclone does not yet support a feature that you wish to use, you can
override the automatically generated code with a manually written PSy-layer
module.

The simplest way to achieve this is to start by building GungHo with a stock
checkout. This will generate all the PSy-layer modules as normal. Then move the
one you wish to override from "working/gungho/algorithm" to
``gungho/source/psy`` and modify it according to your needs. The file will have
the form ``*_psy.f90``. Algorithm files are automatically generated in all cases
so should not be moved in this way.

.. attention::
   You should ensure that the LFRic and PSyclone developers are aware of your
   need to modify the PSy-layer code, to ensure that your changes fit in with
   the ongoing development of LFRic and PSyclone.

Optimisation
~~~~~~~~~~~~

PSyclone is able to make use of scripts to apply optimisations to the code as
it is generating it.

It is necessary for PSyclone to know which platform you are building on in
order to select the correct optimisation scripts. This is achieved through the
``LFRIC_TARGET_PLATFORM`` environment variable. It should contain a single
platform identifier following the convention outlined above.

The build system will look for ``optimisation/<platform id>/global.py`` which
will be applied to all algorithm files. Platform identifiers are specified
above.

If a file ``optimisation/<platform id>/<algorithm>.py`` exists it will be used
in preference to the global script. The algorithm name is taken from
``gungho/source/algorithms/<algorithm>.x90``.

UM physics codes
~~~~~~~~~~~~~~~~

If you wish to compile the GungHo core with UM physics you should use the
"um_physics" sub-project. This is not built by default by the top level
make file due to the time it takes to compile. It can be added to the
``OPERATE_ON`` environment variable should you want it to build from the top
level.

In order to build the UM code with GungHo, fcm make is used to extract and
preprocess it; this is then rsync'd with the working build directory such that
the LFRic build system can proceed with analysing and building the whole code.

The ``LFRIC_TARGET_PLATFORM`` environment variable is used to determine how to
extract the UM source. A lookup table in ``um_physics/fcm-make/target-map.txt``
is used to map LFRic target platform to UM target platform.

The UM code is a considerable size so filters are specified in
``um_physics/fcm-make/extract.cfg`` to restrict the number of files extracted.

The make procedure will first carry out the FCM extract invocation, then
subsequently rsync the extracted and preprocessed code to the ``working``
directory tree. It is the intention that the UM code for the build is kept
separate from the main LFRic source and any modifications on the UM side should
be made through the branches incorporated at the FCM extract stage (these could
be a separate working copy). To change the UM branches incorporated into the
build, modify the ``um_sources`` environment variable in
``um_physics/fcm-make/parameters.sh``.

Since the UM code will continue to evolve and we will want to source difference
versions/branches from the UM repository, the environment variables needed by
FCM extract (including the paths to the repository/branches/working) are set
in ``um_physics/fcm_make/parameters.sh``.  Typically, this script will be
maintained in the LFRic repository but can be overridden if a different UM
source/configuration is required.

The um_physics sub-project offers a special make target ``partial-clean`` or
``pclean`` which deletes only the working Gung Ho copy, not the UM copy.

Creating a Sub-project
----------------------

If you need a new sub-project it should be relatively easy to create. Start
with a project directory in the working copy root. This should have the name
of your new effort.

Below this you are free to do whatever you like but in order to be useable by
the top-level Makefile your project must include a Makefile with a default
target which will generally be used to build the project. It will also need a
"clean" target and one for "test-suite".

You may well want to take advantage of the existing LFRic build infrastructure.
This is held in ``infrastructure/build`` and consists of a number of make files
and tools.

There is a skeleton project, called ``skeleton_project`` which you can use as
a guide and potentially a template.

If your project has a ``build`` directory at its top level then the LFRic
build system will look in it for project specifics. In particular a
``build/fortran`` directory containing makefiles named after different
compilers in which you can put compiler specifics. For general project related
changes use the ``build/project.mk`` file.
