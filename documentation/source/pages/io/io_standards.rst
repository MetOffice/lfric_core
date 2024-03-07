.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

LFRic IO Standards
------------------

Standards for XIOS Iodef Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**File Definitions**


.. admonition:: Rules

   #. Global file settings (for instance, the file type) are to be specified precisely once: in the “Local file definitions” section of every top-level iodef file.
   #. Non-global file settings are to be specified at the level of individual files.

**Rationale**

.. warning::

    When XIOS aggregates file definitions - possibly spread across several XML files -, the last set of attributes processed "wins", in the sense that all the others are silently ignored.

    This can create considerable confusion and frustration for users.

With the current layout of the LFRic top-level iodef files, the last set of attributes to be processed are the ones in the local file definitions section.
