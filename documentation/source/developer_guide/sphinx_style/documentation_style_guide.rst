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
appropriate year).

Purposes of each section
------------------------

There are four main sections to the documentation: "Getting Started",
"User Guide", "Developer Guide" and "LFRic Core API".

Getting Started
~~~~~~~~~~~~~~~

The :ref:`Getting Started <getting_started_index>` section provides an
overview of the code repository and the documentation.

At the time of writing, the User Guide is in the process of being
written, but the other sections are at a very early stage of
development.

User Guide
~~~~~~~~~~

The :ref:`User Guide <user_guide_index>` section is for developers of
LFRic applications. It provides an introduction to the structure of a
typical model application, an overview of the LFRic data model, and
usage examples for key aspects of the LFRic infrastructure.

Developer Guide
~~~~~~~~~~~~~~~

The :ref:`Developer Guide <developer_guide_index>` section is for
developers of the core LFRic code. It describes coding and
documentation standards, gives detailed technical overviews of
aspects of the infrastructure and components, describes the build
system including the use of templated code and code generation, and
describes the testing strategy.

LFRic Core API
~~~~~~~~~~~~~~

The :ref:`LFRic Core API <API_index>` section describes the technical
API for LFRic core software.

Text Formatting
---------------

Lines of text must wrap at a maximum of 80 characters. Other lines
(code snippets and so forth) should not exceed 80 characters unless
there is an exceptional reason.

Headings
--------

Headings used in restructured text documents have no defined order
other than the order within a file. LFRic core will use the following
hierarchy (based on Sphinx recommendations) for most moderately-sized
files:

-  '=' for sections
-  '-' for subsections
-  '^' for subsubsections
-  '"' for paragraphs

If a file is large, and it is necessary to go to a deeper level of
headings, one or both of the following can be used before the above.

-  '#' with overline, for parts
-  '*' with overline, for chapters

Links
-----

Include links to key headings and to figures so other documentation
writers can easily reference them.

While Sphinx can automatically infer links to headings from the title
of the heading, do not rely on this feature as heading titles can
change regularly.

When a document needs to refer to a heading elsewhere in the
documentation, provide an explicit link preceding the section header:

::

    .. _sphinx style guide:

    Style Guide for LFRic Core Sphinx Documentation
    ===============================================

Anticipate the need for other issues to refer to sections by including
links to documentation that is likely to be referenced, even if no
text has a link yet.

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


Note that rendering of the above example may include some default
Sphinx highlighting. When using the Fortran style, the expected
rendering of the above text is as follows:

.. code-block:: fortran

     type( field_collection_type ), pointer :: collection
     type( field_type ), pointer :: field

     collection => modeldb%fields%get_field_collection("my_collection")
     call my_collection%get_field("my_field", field)
