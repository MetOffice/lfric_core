.. -----------------------------------------------------------------------------
     (c) Crown copyright 2025 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _lfric_message_passing_interface:

LFRic Message Passing Interface
===============================

LFRic supports running in parallel over a distributed memory parallel system by
passing messages between the parallel components, For this, it uses a system
library that conform to the Message Passing Interface (MPI) standard.

Most communication is between neighbouring components and is made via halo
exchanges. This is handled by the ``halo_comms`` (and associated) objects, which
subsequently call into the Yaxt library (which then, in turn, calls into the MPI
library) to perform the halo exchanges.

That leaves a small selection of global communication functionality that needs
to be supported, such as global reductions (min/max/sum) and broadcast of
information from one process to all the others. For this, the functionality
provided directly by the MPI library is used.

Calls into the system MPI library are encapsulated into a wrapper object. This
insulates the rest of the model from external changes to the system MPI
library. The object that wraps around the MPI library in the LFRic
infrastructure is called ``lfric_mpi``.

When passing messages between parallel components, the syatem MPI library
provides a way to describe which components should be involved in the
communication, through use of a concept called "communicators". This
communicator is a construct of the system MPI library, so even though it can
sometimes look like a simple integer variable (in plder versions of MPI), it
shouldn't be used in user code. For this reason the lfric_mpi wrapper includes a
wrapper object that sits around the communicator called ``lfric_comm``

The wrapper
module provides access to the communicator being used.

Application Programming Interface
---------------------------------

Helper funcitions
^^^^^^^^^^^^^^^^^

These functions are not part of the ``lfric_mpi_type`` but are usuful helper
functions that are used around use of the object.

* ``subroutine create_comm(comm)`` : Returns a "world" communticator (of type
  ``lfric_comm_type``) by initialsing the system MPI library.
* ``subroutine destroy_comm()`` : Finalises the system MPI library and releases
  the "world" communicator.
* ``function get_lfric_datatype(fortran_type, fortran_kind)
  result(mpi_datatype)`` : Converts a Fortran type/kind into a datatype
  enumerator from the system MPI library.

The ``lfric_comm_type``
^^^^^^^^^^^^^^^^^^^^^^^

This is a wrapper that sits around the system MPI communicator and can be passed
around user code.

* ``function get_comm_mpi_val()``
* ``subroutine set_comm_mpi_val(comm)``


The ``lfric_mpi_type``
^^^^^^^^^^^^^^^^^^^^^^

This is a wrapper that sits around the system MPI library and povides message
passing functinonality.

* ``subroutine initialise(in_comm)`` : Initialises the ``lfric_mpi``
  object based on the given communicator.
* ``subroutine finalise()`` : Finalises the ``lfric_mpi`` object
* ``function get_comm() result(communicator)`` : Reurns the communicator.
* ``function is_comm_set() result(comm_state)`` : Returns whether a communicator
  has been set - i.e. whether the ``lfric_mpi`` object has been initialised
* ``subroutine global_sum(l_sum, g_sum)`` : All parallel tasks provide a local
  value in ``l_sum`` and will be returned the global sum of these values in
  ``g_sum``. This subroutine can be used with 32-bit integers, 32-bit reals and
  64-bit reals.
* ``subroutine global_min(l_min, g_min)`` : All parallel tasks provide a local
  value in ``l_min`` and will be returned the global minimum of these values in
  ``g_min``. This subroutine can be used with 32-bit integers, 32-bit reals and
  64-bit reals.
* ``subroutine global_max(l_max, g_nax)`` : All parallel tasks provide a local
  value in ``l_max`` and will be returned the global maximum of these values in
  ``g_max``. This subroutine can be used with 32-bit integers, 32-bit reals and
  64-bit reals.
* ``subroutine all_gather(send_buffer, recv_buffer, count)`` : Gather integer
  data from all MPI tasks into a single array in all MPI tasks. The data in
  ``send_buffer`` from the jth process is received by every process and placed
  in the jth block of the ``recv_buffer``.
* ``subroutine broadcast(buffer, count, root)`` : Broadcasts the information
  held in buffer on processor root to all other parallel tasks. The variable
  count gives the size of buffer (this should be ommitted if buffer is a scalar).
  Broadcast can be used on logicals, 32-bit integers, 32-bit reals and 64-bit
  reals. Scalars, 1d, 2d and 3d arrays are supported. Broadcast can also be used
  on simple string variables.
* ``function get_comm_size()`` : Returns the number of parallel tasks in the
  current communicator.
* ``function get_comm_rank()`` : Returns the number of the current rank.
