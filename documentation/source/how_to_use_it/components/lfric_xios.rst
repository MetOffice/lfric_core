.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _lfric xios component:

LFRic-XIOS Component
====================

This component contains a classes and routines which can be used to interface with XIOS from within an LFRic model.
These building blocks can be used to create a bespoke I/O solution for any application.

**Running the tests**

Building the code and running the tests requires a Fortran compiler as well as software environment with builds of XIOS, NetCDF and HDF5. The integration tests also require Python >3.9.
To run the unit tests run the command ``make unit-tests`` from this directory.

The integration tests can be run with ``make integration-tests``

    To run the integration tests with plotting enabled, run the command

    ``make integration-tests PLOT_TEST_OUTPUT=true``

    This requires a Python environment with the matplotlib and xarray packages installed.
