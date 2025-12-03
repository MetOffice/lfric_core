.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _configurator:

Configurator
============

The Configurator is a tool that generates Fortran source based on an
:ref:`extended form of Rose metadata <extended rose metadata>`. The
Fortran code reads namelist configuration files aligned with the
metadata and stores the configuration choices in generated data
structures and functions with meaningful names. Applications can use
these structures and functions to access the configuration choices. To
support parallel applications, the generated code manages the
distribution of choices to all MPI ranks.

The Configurator provides several python scripts found in
``infrastructure/build/tools``. Each of these scripts generate
Fortran source code that is specific to an applications metadata.

.. dropdown:: **GenerateNamelistLoader**

  Takes JSON file created by ``rose_picker`` from an applications ``rose-meta.conf`` file.
  For each namelist described in the JSON file, a Fortran module is generated.
  Each module has procedures to read a the specifc namelist from the configuration file,
  MPI broadcast the configuration choices *and to access configuration choices*:

  .. admonition:: Usage

    GenerateNamelistLoader
      *[-help] [-version] [-directory PATH]* FILE

  ``-help`` | ``-version``:
    Returns script information, then exits.
  ``-directory PATH``:
    Location of generated source. Defaults to the current working directory.
  ``FILE``:
    JSON file containing the application metadata.

.. dropdown:: **GenerateConfigLoader**

  The second command generates the code that calls procedures from the
  previously generated namelist loading modules to actually read a
  namelist configuration file:

  .. admonition:: Usage

    GenerateConfigLoader
      *[-help] [-version] [-verbose] [-duplicate LISTNAME]* FILE NAMELISTS...

  ``-help`` | ``-version``:
    Returns the script help information, then exits.
  ``-duplicate LISTNAME``:
    Optional argument to add LISTNAME to the set of namelists allowed
    to have duplicate instances.
  ``FILE``:
    The generated Fortran source code.
  ``NAMELISTS``:
    Space-separated list of one or more namelist names that
    the code will read.

.. dropdown:: **GenerateAppConfigType**

  .. admonition:: Usage

    GenerateConfigType
      *[-help] [-version] [-directory PATH]* FILE

.. dropdown:: **GenerateExtenedNmlType**

  .. admonition:: Usage

     GenerateExtendedNamelistType
       *[-help] [-version] [-directory PATH]* FILE

.. dropdown:: **GenerateFeigns**

  The final command generates a module which provides procedures to
  directly configuring the contents of a namelist. This module ought not
  be used within a normal application. Instead, it is to allow test
  systems to :ref:`feign <feigning configuration>` the reading of a
  namelist so they can control the test environment:

  .. admonition:: Usage

    GenerateFeigns
      *[-help] [-version] [-output FILE1]* FILE2

  ``-help`` | ``-version``:
    Caused the tool to tell you about itself, then exit.
  ``-output FILE1``:
    Generated source file is written FILE1, defaults to ``feign_config_mod.f90``
    in the current working directory,
  ``FILE2``:
    JSON metadata file created by ``rose-picker``.

Ultimately, these scripts require an applications
extended Rose metadata in the form of a JSON file.

For convienence, a separate tool, (:ref:`rose_picker<Rose Picker>`)
is used to convert the extended Rose metadata file into a JSON file.





