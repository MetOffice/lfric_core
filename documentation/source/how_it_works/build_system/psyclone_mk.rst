.. -----------------------------------------------------------------------------
     (c) Crown copyright 2025 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _psyclone_mk:

PSyclone Transformation Scripts
-------------------------------

The LFRic build system searches for PSyclone transformation scripts prior to
applying any PSyclone transformations. Initially, the build system will search
for a Python (``.py``) module with a matching filename to the Fortran module
currently being operated on. If a match is made, the build system proceeds to
call the ``psyclone`` command on the current module, specifying the path to the
transformation script in the command.
If a match is not made, the build system proceeds to search for a more
generic transformation script, ``global.py``, and applies it to the current
module via PSyclone. The transformations called in ``global.py`` will be applied
to any module that does not have a transformation script of its own. This
principle is true for both LFRic and non-LFRic source. Finally, if no
optimisation script is present, global or otherwise, the module is passed to
PSyclone with no transformation script specification. This will also remove any
existing directives in the source file.

.. note::
    LFRic is written in the PSyKAl format and therefore uses PSyclone's
    ``psykal`` method of operation to preprocess modules. There is no support
    for other PSyclone methods of operation in the LFRic Core library.


To direct the build system to find a transformation script, certain GNU
Make variables must be set. Namely, the ``OPTIMISATION_PATH`` and ``DSL``
variables.

- The ``OPTIMISATION_PATH`` variable should be set to the common path between
  all transformation scripts within an application.
- The ``DSL`` variable directs the build to search in specified sub-directories
  of the ``OPTIMISATION_PATH``.

The purpose of ``DSL`` is to separate transformations that would otherwise
conflict. For example, in the ``skeleton`` miniapp, transformation scripts are
grouped by the PSyclone method of operation they require (psykal,  transmute)::

    optimisation/
    └── platform/ (ex1a, minimum, xc40, etc.)
        ├── psykal
        │   ├── global.py
        │   └── sub_directory/
        └── transmute/
            ├── global.py
            └── sub_directory/
