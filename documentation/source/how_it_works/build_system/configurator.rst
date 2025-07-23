.. -----------------------------------------------------------------------------
    (c) Crown copyright 2024 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

Configurator
============

The Configurator is a tool that generates Fortran source based on an
:ref:`extended form of Rose metadata <extended-rose-metadata>`. The
Fortran code reads namelist configuration files aligned with the
metadata and stores the configuration choices in generated data
structures and functions with meaningful names. Applications can use
these structures and functions to access the configuration choices. To
support parallel applications, the generated code manages the
distribution of choices to all MPI ranks.

Usage
-----

The Configurator calls three commands which may be found in
``infrastructure/build/tools`` and a separate tool
:ref:`rose_picker<Rose Picker>` which
converts the extended Rose metadata file into a JSON file.

The first command takes the JSON file created by ``rose_picker`` and
creates a module for each namelist. Each module has procedures to read
a namelist configuration file for the namelist, to MPI broadcast
configuration choices and to access configuration choices::

    GenerateNamelist [-help] [-version] [-directory PATH] FILE

The ``-help`` and ``-version`` arguments cause the tool to tell you about
itself, then exit.

The ``FILE`` argument points to the metadata JSON file to
use. Generated source is put into the current working directory, or
into ``PATH`` if specified.

The second command generates the code that calls procedures from the
previously generated namelist loading modules to actually read a
namelist configuration file::

    GenerateLoader [-help] [-version] [-verbose] FILE NAMELISTS...

As before, ``-help`` and ``-version`` options reveal details about
the tool before exiting.

The ``FILE`` is that of the resulting generated source file. Finally,
the ``NAMELISTS`` are a space-separated list of one or more namelist
names that the code will read.

The final command generates a module which provides procedures to
directly configuring the contents of a namelist. This module ought not
be used within a normal application. Instead, it is to allow test
systems to feign the reading of a namelist so they can control the
test environment::

    GenerateFeigns [-help] [-version] [-output FILE1] PATH2

Once again, ``-help`` and ``-version`` cause the command to exit after
giving its details.

The resulting source file is written to ``FILE1``, or to
``feign_config_mod.f90`` in the current working directory if ``FILE1``
is not specified. The ``FILE2`` argument points to the JSON metadata
file created by ``rose-picker``.
