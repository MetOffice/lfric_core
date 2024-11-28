Build System Component: Dependerator
====================================

If you need to work on the LFRic build system you will need some grounding in
the "Dependerator" tool which is an integral part of it.

Dependency Analysis
~~~~~~~~~~~~~~~~~~~

In order to effectively and efficiently build a Fortran source tree it is
necessary to understand how the various source files relate to each other. This
dependency analysis allows them to be built in the correct order and is
performed by the "Dependerator" tool.

There are two stages:

#. Analyse ever source file
#. Generate dependency information

The reason every source file is analysed is that ahead-of-time it is not
possible to know which files are needed and which are not. That is part of the
purpose of dependency analysis.

Performance is gained by the rest of the build system which ensures only those
files which have changed are analysed.

The information harvested in the first stage is collected and transmitted to the
second stage through the use of a database file.

Examine the Source
~~~~~~~~~~~~~~~~~~

To harvest pertinent information the ``DependencyAnalyser`` tool is used::

    infrastructure/build/tools/DependencyAnalyser <database file> <source file>

The database is used to maintain the harvested information in preparation for the
second stage. This is what is referred to by the first argument and is created
by the tool.

There are additional optional arguments which may be useful:

``-ignore <module name>``
    May be specified multiple times.
    
    Usage of this module will not be included in the database. This is
    essential for 3rd party libraries which will not be found as part of the
    build.

``-include <directory>``
    May be specified multiple times.
    
    The preprocessor will use this directory for inclusions. It is passed
    through to an ``-I`` argument.

``-macro <name>[=<value>]``
    May be specified multiple times.
    
    The preprocessor will use this macro definition. If the value is not
    specified then the macro is defined but uninitialised.
    
    It is passed through to a ``-D`` argument.

You may also use the standard ``-help``, ``-version``, ``-verbose`` and
``-debug`` for their normal purposes.

Analyse Dependencies
~~~~~~~~~~~~~~~~~~~~

Once all the source is examined there are a couple of tools which take that
information, perform analysis and create fragments of Make file. These can
then be included by the build system to obtain the dependency information
needed.

Build Order
-----------

Creating a tree of dependencies which determine build order is performed by the
``DependencyRules`` tool::

    infrastructure/build/tools/DependencyRules <output file>

If not specified the tool will look for a database file called
``dependencies.db`` in the same directory the output file will be created in.
If it is somewhere else use the ``-database`` argument to specify it.

Two arguments, ``-moduledir`` and ``-objectdir`` are used while generating the
Make file fragment. If they are not specified they default to the directory
where the fragment will be created.

It is the responsibility of the build system to make sure these files actually
end up in these locations.

Also available is ``-moduleobjects`` to support compilers such as those from
Cray which stores module information in the object files rather than their own
separate file.

Also supported are the common ``-help``, ``-version``, ``-verbose`` and
``-debug`` arguments.

Resulting Make Fragment
.......................

What ends up in the resulting Make fragment (the build system calls this
``dependencies.mk``) are a number of rules for each file.

Matters are complicated by the fact that Make operates on a file-by-file basis
while Fortran operates partly file-by-file and partly module-by-module.

Each source file produces a single object file but potentially several module
files.

Thus an object file depends on a single source file (not captured by this
tool but inferred in the build system) but multiple module files, one for each
modules used by any program unit contained in the source file::

    path/to/third_mod.o: $(MOD_DIR)/first_mod.mod $(MOD_DIR)/second_mod.mod

You will notice the ``MOD_DIR`` variable. This is a place holder used in the
generated output but not defined there. The wider build system defines this to
point to the directory which holds the modules.

For compilers which generate module files (i.e. the ``-moduleobjects`` argument
has not been used) lines are also generated for each module the source file
holds::

    $(MOD_DIR)/third_mod.mod: path/to/third_mod.o
    $(MOD_DIR)/fourth_mod.mod: path/to/third_mod.o

Build Composition
-----------------

The second analysis tool, ``ProgramObjects`` generates a list of all the object
files needed to link the programs embodied by the source:

    infrastructure/build/tools/ProgramObjects <output file>

As per ``DependencyRules`` this tool expects to find a ``dependencies.db`` in
the same directory as the output file unless you use the ``-database``
argument.

It also takes an ``-objectdir`` argument if objects are not found in the
output file directory.

Again the common ``-help``, ``-version``, ``-verbose`` and
``-debug`` arguments may be used.

Resulting Make Fragment
.......................

The build system puts the output of this tool into ``programs.mk``.

The file defines a variable for each program listing the module objects required
to link that program::

    SOME_PROGRAM_OBJS = special/first_mod.o normal/second_mod.o path/to/third_mod.o

It also creates a variable ``PROG_OBJS`` while lists all the program objects::

    PROG_OBJS = some_program.o

The names match so the module objects needed to link with ``some_program.o`` are
held in ``SOME_PROGRAM_OBJS``.
