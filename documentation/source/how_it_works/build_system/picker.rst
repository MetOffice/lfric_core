.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _rose picker:

Rose Picker
===========

The ``rose_picker`` tool is derived from Rose and maintained
separately from LFRic code for software licencing reasons.

Usage
-----

The tool is invoked like this::

    rose_picker [-help] [-directory PATH1] [-include_dirs PATH2] FILE

Usage information is provided if the ``-help`` option is specified.

The tool generates a JSON file and a text file containing a list of
namelists from the input metadata file specified by ``FILE``. If
specified, output files are written to the ``PATH1``
directory. Otherwise, the files are written to the current working
directory.

If a metadata file uses the ``include`` directive to import additional
files an ``-include_dirs`` option must be used to specify the
``rose-meta`` directories containing them. Where two (or more)
``rose-meta`` directories need to be searched the ``-include_dirs``
argument is specified once for each. For example::

    rose_picker -include_dirs rose-meta                           \
                -include_dirs ../lfric_core/rose-meta             \
                rose-meta/lfric-gungho_model/HEAD/rose-meta.conf
