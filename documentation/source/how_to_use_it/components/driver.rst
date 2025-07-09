.. -----------------------------------------------------------------------------
    (c) Crown copyright 2024 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _driver component:

Driver Component
================

The driver component contains some convenience modules that can be used to help
construct the driver layer of an application. Applications can use some or all
of the driver component, mixing and matching as necessary.

Bear in mind that the driver component only supports typical usage. If your specific
needs are not covered, then your application will need its own specific driver code.
The driver component may act as a starting point for application specific driver code,
but it should not be modified with any application specific logic.


.. _driver logging:

Logging
-------

The logging component ``driver_log_mod`` provides procedures for
initialising and finalising the :ref:`logging system<logging>` using
the namelist configuration options of an application.

The ``init_logger`` procedure takes two arguments that are passed
directly to the logging system: an MPI communicator and a name that
will be used when constructing the name of the output log files.

The procedure reads the options that define whether the application is
configured to write to rank zero only, and the logging level for the
application. These are passed to the relevant calls of the logging
system.

.. warning::

   Several model applications configure the logging system based on
   namelist input. That means that namelists must be read prior to the
   logging system being initialised. The namelist code does include
   some calls to output log messages if there are errors in the
   namelist. Since the logger will not have been initialised, these
   error messages will appear in the standard output rather than in
   the PET output files.
