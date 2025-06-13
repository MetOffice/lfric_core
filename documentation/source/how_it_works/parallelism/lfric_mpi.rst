.. -----------------------------------------------------------------------------
     (c) Crown copyright 2025 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _lfric_mpi:

LFRic Message Passing Interface
===============================

All global MPI operations are wrapped in the lfric_mpi object. Processor to
processor (i.e. halo exchanges) are handled by ``halo_comms``. Most of the
functions in the `lfric_mpi` object are simple wrappers to MPI library
functions. They provide a way to encapsulate the MPI functionality and insulate
the rest of the code from changes to the external library.

For details of the API to the ``lfric_mpi`` wrappers, see 
:ref:`How to use: lfric_mpi <lfric_message_passing_interface>`

Different versions
------------------

The ``lfric_mpi`` wrapper provides three versions of the object that provide
different levels of functionality. These are selected through preprocessor
directives. The build system sets these directives based on how environment
variables are set.

No MPI
^^^^^^

If for some reason, it is not possible or not desirable to link to the system
MPI library, it is possible to preprocess the source so all calls to the system
MPI library are eliminated. Single processor serial running of the model is
still supported in this mode. This option is chosen by setting the preprocessor
directive ``NO_MPI``, which is set by the build system in response to the
environment variable ``NO_MPI``. The logic is negative, because the default (in
the absence of any environment variables or preprocessor directive) should be to
link to the MPI library.

MPI_F08/Legacy
^^^^^^^^^^^^^^

As an alternative to not using MPI at all, the lfric MPI wrapper can
use one of two interfaces to the Fortran MPI library. The original "legacy" interface stores many of its
internals (such as the handle used to identify communicators) as integers. This
means that it impossible to implement effective type checking - one integer
looks just like any other to the compiler. A second interface was implemented
where each these integers is wrapped in its own object (e.g. the communicator is
held in an ``mpi_comm`` object) - allowing type checking. This interface is
called "mpi_f08". The ``lfric_mpi`` wrapper object supports both interfaces.
When directly building the code, select the legacy interface by setting the
proprocessor directive ``LEGACY_MPI``, Not setting any preprocessor directives
will result in code that uses the mpi_f08 interface. However, the logic is
reversed in the current build system. Here, setting the environment variable
``USE_MPI_F08`` makes the code use the "mpi_f08" interface. The default, in the
absence of any environment variable being set is to use the legacy mpi
interface.

Moving between the different interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To insulate LFRic users from how the system library holds the communicator,
LFRic users should only pass around the LFRic construct that contains the
communicator: the ``lfric_comm_type`` object.

But even with this abstraction, the code has to interact with external libraries
that may still be using the legacy interface (such as XIOS, Yaxt etc.) For these
cases, the ``lfric_comm_type`` object can provide the legacy, integer version of
the communicator by calling the function:
``lfric_comm_type%get_comm_mpi_val()``.

When the Oasis coupler is used, the Oasis3-MCT library insists on initialising
MPI and returning the integer version of the communicator for LFRic to use. For
this reason, it is also possible to initialise the ``lfric_comm_type`` from the
integer version of the communicator by calling the subroutine:
``lfric_comm_type%set_comm_mpi_val(comm)``.
