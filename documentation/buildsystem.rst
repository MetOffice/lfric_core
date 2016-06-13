The Dynamo Build System
=======================

This is a quick introduction to Dynamo's build system. Such a system is
intended to make developers' lives easier but to do so it requires a few
conventions to be followed.

.. note::
   The canonical version of this document exists as a reStructured text file
   in the repository. Branches which render it inaccurate should
   update the version of this document on that branch. The version
   displayed in the wiki is generated from the head of trunk.

Targets
-------

A number of options, referred to as targets, are offered by the top
level ``make`` file, as described in the following sections.

Building
~~~~~~~~

It is necessary for PSyclone to know which platform you are building on in
order that it be able to select the correct optimisation scripts. This is
achieved through the `DYNAMO_BUILD_TARGET` environment variable. It should
contain a single platform identifier following the same convention as outlined
for test suite targets, detailed in the section on testing.

This variable will be set up for Met Office developers using the LFRic module
system.

Three build targets are offered: ``full-debug``, ``fast-debug`` and
``production``. The default is ``fast-debug`` if no target is specified.

All three insert debug symbols to the executable which provides source file and
line numbers in backtraces. They have no run-time impact other than to
increase the size of the executable.

``Full-debug`` turns optimisation off and run-time checking, of e.g. array
out-of-bounds, on.

For those who need faster execution but also want bit reproducibility
"fast-debug" limits optimisation to a safe level and turns off run-time
checking.

``Production`` increases the optimisation level but does not guarantee
bit reproducibility. This is due to potential run-time optimisations,
particularly in the order of shared and distributed memory processing. We would
expect inter-run reproducibility but it can not be guaranteed.

If you want fine control over how the build is performed a number of "knobs"
are provided in the form of variables which may be passed to the make command.

  ``make [OPTIMISATION=(NONE|SAFE|RISKY)] [SYMBOLS=(YES|NO)] [CHECKS=(YES|NO)]``

These may be mixed with the targets so if you wanted a safe build without
symbols but with run-time checking you could achieve it with:

  ``make fast-debug SYMBOLS=NO CHECKS=YES``

Relocate Build Artifacts
^^^^^^^^^^^^^^^^^^^^^^^^

By default all build related files end up in the directory "build" within the
working copy. This is convenient for interactive development as it allows the
generated source to be examined. Unfortunately there are a number of situations
where this is less than ideal.

The database engine used to store dependency information interacts badly with
NFS exported mounts. If your working copy is stored on such a mount you may
find that builds are very slow. They will tend to stall during the construction
of ``programs.mk``.

Another filesystem which can have problems is Lustre, often used on
supercomputers. It is tuned for large parallel file access to a few files, not
the many serial accesses to small files involved in compiling. Thus a working
copy held on such a file system will be slow to compile and place unreasonable
demands on the file server.

Both these problems can be aleviated by relocating the build directory.

This is achieved by simply exporting the ``DYNAMO_BUILD_ROOT`` environment
variable. It should contain an absolute path to a suitible temporary space in
which to perform the compile

Verbose Building
^^^^^^^^^^^^^^^^

By default the build system will suppress as much output as possible to reduce
clutter. When resolving problems it may be useful to see the output. If you want
to see the actual commands being run by the build use ``make VERBOSE=1``.

Building Without a Dependency Check
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the moment building the dependency information is much slower than it should
be. This can hinder interactive development.

Until the problem is fixed properly you can set the ``SKIPDEP`` knob to bypass
the dependency check. Only do this if you know for a fact you have not changed
the dependencies.

  ``make SKIPDEP=1``

MPI
^^^

If you are attempting to build with MPI support do *not* set ``FC`` to
``mpif90``. For a start some compilers (Cray, IBM) deal with MPI themselves.

Instead define the variable ``LDMPI`` with the command to link against the MPI
libraries.

Linking
^^^^^^^

By default the build system will link dynamically. If you want a statically
linked binary then you should pass the ``LINK`` variable: ``make LINK=STATIC``.

Static linking is the default on Cray systems as they do not seem to play well
with dynamic linking.

Cleaning
^^^^^^^^

As with many build systems it is possible to ``make clean`` to delete Dynamo
and all its build artefacts. You should, however, be aware that this does not
delete the pFUnit installation made as part of the build process. This is by
design as people do not generally want to rebuild pFUnit every time they clean
out their project build.

When you want to rebuild your code with a different compiler, use
``make clean-all``. This will clean pFUnit as well.

Testing
^^^^^^^

Unit tests will be built and run as part of a normal build, so there is no
need to worry about that. The test suite, on the other hand, must be
manually invoked.

When ``make test-suite`` is used a number of instances of Rose Stem will be
launched. The environment variable ``DYNAMO_TEST_SUITE_TARGETS`` holds a space
separated list of the different test suite configurations to use.

Note that each suite is named after the test suite configuration and
not the branch. Therefore, you cannot run ``make test-suite`` on two
branches at the same time.

Configurations are held in ``rose-stem/opt``. Each filename has the form
``rose-suite-<name>.conf`` and it is these names which should be used in the
list.

For those using a Met Office module collection the core dynamo module will set
this up for you. e.g. On the desktop do
``module load base-environment/dynamo/base``.

For further information an testing see `DynamoTesting`:trac:.


Choosing a Compiler
-------------------

A number of compilers are supported by the build system. To choose between them
simply set the ``FC`` environment variable. e.g. ``export FC=ifort`` This is
often done by an automated system, as it is for users of the Met Office Dynamo
module system.

Each compiler needs a different set of arguments. These are defined in an
hierarchical fashion.

make/fortran/<compiler>.mk
~~~~~~~~~~~~~~~~~~~~~~~~~~

Universal compiler options and build system variables are put in this file.

=========================  =================================================================================
Variable                   Purpose
=========================  =================================================================================
``<compiler>_VERSION``     Holds the version of the running compiler. This is held as a zero padded integer. e.g. 4.9.2 becomes 040902
``F_MOD_DESTINATION_ARG``  The argument to be used to tell the compiler where to put ``.mod`` files.
=========================  =================================================================================

src/dynamo/make/fortran/<compiler>.mk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes it is necessary to pass additional arguments to the dependency
generator script. This is achieved using the ``DEPRULE_FLAGS`` variable.

The various arguments have been split into groups to aid clarity but are all
brought together in ``FFLAGS`` in the end.

=============================  ===========================================================================================================
Argument Group                 Purpose
=============================  ===========================================================================================================
``FFLAGS_COMPILER``            Control the way the compiler works. Examples include turning on OpenMP.
``FFLAGS_NO_OPTIMISATION``     Turn off all optimisation.
``FFLAGS_SAFE_OPTIMISATION``   Enable only bit reproducibility safe optimisation.
``FFLAGS_RISKY_OPTIMISATION``  Enable optimisations which may affect bit reproducibility.
``FFLAGS_DEBUG``               Enable the generation of debugging symbols, tracebacks and such.
``FFLAGS_WARNINGS``            Control warning handling. For development we like to receive all warnings and have them act as errors.
``FFLAGS_INIT``                Turn on any declaration initialisation the compiler may support to fill declared memory with a known value.
``FFLAGS_RUNTIME``             Turn on any run-time checking the compiler may support. Array bound checking for instance.
=============================  ===========================================================================================================

The same is done with ``LDFLAGS``.

====================  =======================================================================================================
Argument Group        Purpose
====================  =======================================================================================================
``LDFLAGS_COMPILER``  Some linkers need to be told to link an OpenMP executable. This is where to specify things of that ilk.
====================  =======================================================================================================

Sometimes the link stage needs special libraries. To specify these use
``LD_COMPILER_LIBRARIES``. It should hold a space separated list of additional
libraries. These should be library names, not filenames. i.e. drop the "lib"
from the front.

Automatic Dependency Analysis
-----------------------------

A knowledge of dependencies between source files allows only those things
affected by a change to be recompiled. Maintaining such a dependency list by
hand is tedious and error prone. This is why the build system will generate
and maintain dependency information automatically.

The dependency analyser scans Fortran source looking for "use" statements
to determine the prerequisites of a source file. It relies on source files
having the same name as the module they contain, which implies that they contain
only one module.

e.g. if the module is called "my_special_stuff_mod" then the file name
should be "my_special_stuff_mod.f90". If this is not the case these modules will
be rebuilt every time.

The current dependency analyser works well for most common uses, but it will
fail when you rename files. If you rename a file, then run ``make clean`` prior
to rebuilding. This will force the dependency database to be rebuilt.


External Libraries
------------------

If you want to make use of an external library then some simple edits to
``src/dynamo/Makefile`` will be required.

There are two critical variables: ``IGNORE_DEPENDENCIES`` and
``EXTERNAL_*_LIBRARIES``.

``IGNORE_DEPENDENCIES`` is a space-separated list of modules which should be
ignored by the dependency analyser when discovered in ``use`` statements. This
prevents the dependency analysis trying to rebuild your library. For example,
Dynamo uses the ESMF library: the name to add would be "esmf" since we
``use esmf``.

``EXTERNAL_*_LIBRARIES`` is a space-separated list of library names to be passed
to the linker using "little ell" arguments. Libraries which are available as a
".so" file should be listed in ``EXTERNAL_DYNAMIC_LIBRARIES`` while those only
available as a ".a" file are listed in ``EXTERNAL_STATIC_LIBRARIES``. Returning
to the ESMF example, the library file is ``libesmf.so`` so the string ``esmf``
would be added to ``EXTERNAL_DYNAMIC_LIBRARIES``.

In addition to adding knowledge of the library to the build system you have to
make it findable. If you are using a package from the Met Office's Dynamo
module collection then this is handled for you automatically. If you are not
then you will need to do a little additional work.

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

Dynamo uses the PSyclone tool, created as part of the GungHo project,
to automatically generate some of the code. PSyclone is designed to generate the
Parallel Systems layer (PSy-layer) code, but it also modifies the
algorithm code as written by science developers.

Algorithm routines in "src/dynamo/algorithm" are written with a ``.x90``
extension. From these both compilable Fortran source for the algorithm and
corresponding PSy-layer source are generated. Normally this process is handled
for you by the build system.

Where PSyclone does not yet support a feature that you wish to use, you can
override the automatically generated code with a manually written PSy-layer
module.

The simplest way to achieve this is to start by building Dynamo with a stock
checkout. This will generate all the PSy-layer modules as normal. Then move the
one you wish to override from "build/dynamo/psy" to "src/dynamo/psy" and modify
it according to your need. The file will have the form "psy_*.f90". Algorithm
files are automatically generated in all cases so should not be moved in this
way.

.. attention::
   You should ensure that the Dynamo and PSyclone developers are aware of your
   need to modify the PSy-layer code, to ensure that your changes fit in with
   the ongoing development of Dynamo and PSyclone.

Optimisation
~~~~~~~~~~~~

PSyclone is able to make use of scripts to apply optimisations to the code as
it is generating it.

The build system will look for ``optimisation/<platform>/global.py`` which will
be applied to all algorithm files. The "platform" specifier follows the same
form as that used for test suite targetting.

If a file ``optimisation/<platform>/<algorithm>.py`` exists it will be used in
preference to the global script. The algorithm name is taken from
``src/dynamo/algorithms/<algorithm>.x90``.
