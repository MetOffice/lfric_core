.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

LFRic IO Standards
==================

Standards for XIOS Iodef Files
------------------------------

**File Definitions**


.. admonition:: Rules

   #. Global file settings (for instance, the file type) are to be
      specified precisely once: in the “Local file definitions”
      section of every top-level iodef file.
   #. Non-global file settings are to be specified at the level of
      individual files.
   #. ``id`` values shall be unique for a given ``XIOS`` ``XML`` element
      type.

**Rationale**

.. warning::

    When XIOS aggregates definitions - possibly spread across
    several XML files -, then attributes are overriden by defined elements
    that use a ``_ref``.

    The overriding operates in a chain, so users may need to act cautiously
    to understand attribute inheritance and overrides.

    This is why reuse of ``id`` values for different elements is problematic,
    not allowed, and tested by a code inspection rule for ``field`` elements.

With the current layout of the LFRic top-level iodef files, the last
set of attributes to be processed are the ones in the local file
definitions section.
