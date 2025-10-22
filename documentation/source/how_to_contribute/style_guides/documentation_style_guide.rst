.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _sphinx style guide:

Style Guide for LFRic Core Sphinx Documentation
===============================================

Documentation of LFRic Core is written using restructured text and is
rendered by Sphinx, typically into html.

This section recommends practices to be followed to ensure the
documentation follows a consistent style.

Copyright
---------

Include a copy of the standard Copyright text in all files.

.. code-block:: rst

  .. -----------------------------------------------------------------------------
      (c) Crown copyright Met Office. All rights reserved.
      The file LICENCE, distributed with this code, contains details of the terms
      under which the code may be used.
     -----------------------------------------------------------------------------

Purposes of each section
------------------------

There are four main sections to the documentation: "Getting Started",
"How to use it", "How it works" and "LFRic Core API".

At the time of writing, the How to use it is in the process of being
written, but the other sections are at a very early stage of
development, so a summary of the purpose of each section is given
here.

Getting Started
~~~~~~~~~~~~~~~

The :ref:`Getting Started <getting_started_index>` section should
contain the most important information for someone who has never used
LFRic Core before, details on how to get started building the included
applications, the general layout of the repository and its contents as
well as any more detailed tutorials.

How to use it
~~~~~~~~~~~~~

The :ref:`How to use it<how_to_use_it_index>` section is for developers of
LFRic applications. It provides an introduction to the structure of a
typical model application and how to use them.

How it works
~~~~~~~~~~~~

The :ref:`How it works<how_it_works_index>` section is for developers of the
LFRic core software and experienced developers who may need to know the
implementation details of the LFRic Core components they are using.

How to contribute
~~~~~~~~~~~~~~~~~

The :ref:`How to contribute<how_to_contribute_index>` section is for
developers of the core LFRic code. It describes coding and documentation
standards.

Text Formatting
---------------

Lines of text must wrap at a maximum of 80 characters. Other lines
(code snippets and so forth) should not exceed 80 characters unless
there is an exceptional reason.

Headings
--------

LFRic core documentation will use the following hierarchy for marking
headings (based on Sphinx recommendations) for most moderately-sized
files:

-  '=' for sections
-  '-' for subsections
-  '^' for subsubsections
-  '"' for paragraphs

If a file is large, and it is necessary to go to a deeper level of
headings, one or both of the following can be used before the above.

-  '#' with overline, for parts
-  '*' with overline, for chapters

.. caution::

    Headings used in restructured text documents have no defined order
    other than the order of occurrence within a file. Be sure to check
    the order adheres to the working practices.

Links
-----

Include a reference label to allow other text to include links to
headings and to figures. Add the reference label before a heading as
in the following example:

.. code-block:: rst

    .. _sphinx style guide:

    Style Guide for LFRic Core Sphinx Documentation
    ===============================================

While Sphinx can automatically infer links to headings from the title
of the heading, do not rely on this feature as heading titles can
change regularly.

It is acceptable to add a reference to a heading or figure even when
no other text currently refers to it as it makes it easier should a
future writer want to add a link to the item.

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

Note that Spinx may have applied some default highlighting to the
above text. But, when using the Fortran style, the expected rendering
of the above text is as follows:

.. code-block:: fortran

     type( field_collection_type ), pointer :: collection
     type( field_type ), pointer :: field

     collection => modeldb%fields%get_field_collection("my_collection")
     call my_collection%get_field("my_field", field)
