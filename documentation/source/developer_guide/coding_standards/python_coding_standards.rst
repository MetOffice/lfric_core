.. ------------------------------------------------------------------------------
     (c) Crown copyright 2025 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

Python Coding Standards
=======================

Rules for coding standards are driven by the need for readability and
consistency (which aids readability). While some people are happy to read
inconsistent code, other people find inconsistency to be distracting and their
needs should be respected.

Python version
~~~~~~~~~~~~~~

The LFRic project uses only Python version 3. We track version deprecation so
as support ends for Python versions we stop supporting them too. At the time
of writing the oldest supported version was 3.9.

Style guide
~~~~~~~~~~~

There exists already a widely used style guide known as `PEP 8`_. All new or
modified Python code in LFRic trunk must be PEP8 compliant.

To check whether your code is PEP8 compliant, run the `flake8`_ tool style
checker on your source file, i.e::

    flake8 /path/to/file.py

Resolve *all* errors and warnings (but see :ref:`E402` exception below).

Flake8 is available in the LFRic software stack on SPICE Desktop and servers.

.. _PEP 8: https://www.python.org/dev/peps/pep-0008
.. _flake8: https://flake8.pycqa.org/en/latest/


Copyright
~~~~~~~~~

Every Python source file should have the following mast-head:

.. code-block:: python

    ##############################################################################
    # (c) Crown copyright 2025 Met Office. All rights reserved.
    # The file LICENCE, distributed with this code, contains details of the terms
    # under which the code may be used.
    ##############################################################################


Whitespace
~~~~~~~~~~

Unlike Fortran, Python supports the tab character but we should continue to use
only spaces for indentation. Spaces don't change size between editors and
platforms.

Labels
~~~~~~

Labels appear as names for variables, classes and functions. In order to match
our Fortran style we mandate pothole_case for these labels. However we do not
mandate the addition of suffixes like "_type". These exist in the Fortran
standard for technical reasons which do not apply to Python.

.. _E402:

E402
~~~~

When running plotting scripts in batch mode, it is common to set the backend
before importing pyplot:

.. code-block:: python

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

When ``flake8`` is run on this code, the following error is raised::

    test.py:3:1: E402 module level import not at top of file

This deviation from PEP8 is considered to be acceptable for plotting scripts
which require batch processing only. In these cases please add the
``# noqa: E402`` directive to the ``matplotlib.pyplot`` import line.


Parenthesis
~~~~~~~~~~~

.. deprecated::
    2020

    We mandated a space between parenthesis and argument lists (as we do for
    Fortran), however this is not PEP8 compliant (the existing code shall be
    updated when it is modified).
