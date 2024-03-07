.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _rose picker:

Rose Picker
===========

Since Rose already implements code to read its own configuration files we would
ideally make use of it to save effort and prevent diversion of two
implementations.

Unfortunately Rose is licenced under the GPL which is incompaible with the BSD
licence we use for LFRic. As such we chose to produce a GPL stub program which
makes use of the Rose code to read the metadata file and produce a JSON
representation from it. This file can be used by our tools without fear of
licence infringement.

In order to maintain the separation of the two licences ``rose_picker`` is
held in a separate repository within the LFRic collection.

Usage
~~~~~

The tool is invoked like this::

    rose_picker [-help] [-directory PATH] [-include_dirs PATH] PATH

As ever, usage information is available using ``-help``.

The location where generated files are placed is provided by ``-directory`` and
defaults to the current working directory.

Metadata files which use the "include" directive to import additional files
require those additional files to be located. They will be searched for in the
directories specified by ``-include_dirs``. This argument may be specified
multiple times, once for each directories.

The final ``PATH`` is that of the metadata file.

The reason the output path is specified as a directory is that as well as the
JSON file a text file containing a list of namelists defined is also generated.
