.. ------------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _sphinx style guide:

Style Guide for LFRic Core Sphinx Documentation
===============================================

Documentation of LFRic Core is written using restructured text and is
typically rendered by Sphinx.

This section recommends practices to be followed to ensure the
documentation follows a consistent style.

Copyright
---------

Include a copy of the standard Copyright text in all files (with an
appropriate year):

Text Formatting
---------------

Lines of text must wrap at a maximum of 80 characters. Other lines
(code snippets and so forth) should not exceed 80 characters unless
there is an exceptional reason.

Headings
--------

Headings used in restructured text documents have no defined order
other than the order within a file. LFRic core will use the following
hierarchy.

-  '#' with overline, for parts
-  '*' with overline, for chapters
-  '=' for sections
-  '-' for subsections
-  '^' for subsubsections
-  '"' for paragraphs

Not all levels of the hierarchy need to be used in a file: many files
will start at the "sections" level of this hierarchy.

Links
-----

Include links to key headings and to figures so other documentation
writers can easily reference them.

Sphinx can automatically infer links to headings from the title of the
heading. Do not rely on this feature as heading titles can change
regularly.  When a document needs to refer to a heading elsewhere in
the documentation, provide an explict link preceding the section
header:

::

    .. _sphinx style guide:

    Style Guide for LFRic Core Sphinx Documentation
    ===============================================

Code snippets
-------------

When writing documentation, it often useful to include code snippets. When this is done:

-    use the correct syntax highlighting
-    define key variables used in the snippet.

For example, the following will use syntax highlighting appropriate to
Fortran code:

::

   .. code-block:: fortran

     type( field_collection_type ), pointer :: collection
     type( field_type ), pointer :: field

     collection => modeldb%fields%get_field_collection("my_collection")
     call my_collection%get_field("my_field", field)


Note that the above example may include some default Sphinx
highlighting. The expected rendering of the above text using the
Fortran style is as follows:

.. code-block:: fortran

     type( field_collection_type ), pointer :: collection
     type( field_type ), pointer :: field

     collection => modeldb%fields%get_field_collection("my_collection")
     call my_collection%get_field("my_field", field)
