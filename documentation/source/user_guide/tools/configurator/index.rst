.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

Configurator
============

Used by the build system this tool takes an extended Rose metadata file and
produces Fortran source to manage the namelists produced from that metadata.

.. toctree::
   :caption: Further Details
   :maxdepth: 1

   picker
   metafile

Usage
~~~~~

The Configurator consists of three commands which may be found in
``infrastructure/build/tools`` and a separate tool :ref:`Rose Picker` which
converts the extended Rose metadata file into a JSON file.

The first takes the metadata and creates namelist loading modules::

    GenerateNamelist [-help] [-version] [-directory PATH] PATH

The ``-help`` and ``-version`` arguments cause the tool to tell you about
itself, then exit.

The final path is the metadata JSON file to use and the optional ``-directory``
tells the tool where to put the generated source. If it is not specified the
current working directory is used.

The next command generates an orchestration which makes use of the previously
generated namelist loading modules to load a namelist file::

    GenerateLoader [-help] [-version] [-verbose] PATH NAMELIST ...

As before ``-help`` and ``-version`` reveal details about the tool and exit.

The ``PATH`` is that of the resulting generated source file. Finally the
namelists expected to appear in the file are presented as a space separated
list of ``NAMELIST`` names.

The final command generates a module which which fakes the loading of a
namelist. This is useful for controlling a test environment::

    GenerateFeigns [-help] [-version] [-output PATH] PATH

.. warning::

    This tool is deprecated and will be removed once it is no longer needed.

Once again ``-help`` and ``-version`` exit after giving their information.

The resuling source file will be named by the ``-output`` argument and defaults
to ``feign_config_mod.f90`` in the current directory. The final ``PATH`` is,
once again, the JSON metadata file.
