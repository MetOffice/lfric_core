.. -----------------------------------------------------------------------------
     (c) Crown copyright 2025 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _lfric_message_passing_interface:

LFRic Message Passing Interface
===============================

LFRic supports running in parallel over a distributed memory parallel system by
passing messages between the parallel components, For this, it uses the system
library that conforms to the Message Passing Interface (MPI) standard
(see the `MPI-Forum <https://www.mpi-forum.org/>`_).

Most communication is between neighbouring components and is made via halo
exchanges. This is handled by the ``halo_comms`` (and associated) objects, which
subsequently call into the Yaxt library
(see the `Yaxt docs <https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/>`_).
which then, in turn, calls into the
MPI library to perform the halo exchanges.

That leaves a small selection of global communication functionality that needs
to be supported, such as global reductions (min/max/sum) and broadcast of
information from one process to all the others. 

This functionality is provided by calls to the system MPI library. These are
encapsulated in a wrapper object. This insulates the rest of the model code from
external changes to the system MPI library. The object that wraps around the MPI
library in the LFRic infrastructure is called ``lfric_mpi``.

For more details about the inner workings of the wrapper, see
:ref:`How it works: lfric_mpi <lfric_mpi>`

Application Programming Interface
---------------------------------

Some of the functionality in the programming interface will only be used by
advanced users, but the full API is listed here for completeness.

Helper functions
^^^^^^^^^^^^^^^^^

These functions are not part of the ``lfric_mpi_type`` but are useful helper
functions that are used around use of the object.

* ``subroutine create_comm(comm)`` : Returns a "world" communicator (of type
  ``lfric_comm_type``) by initialising the system MPI library.
* ``subroutine destroy_comm()`` : Finalises the system MPI library and releases
  the "world" communicator.
* ``function get_lfric_datatype(fortran_type, fortran_kind)
  result(mpi_datatype)`` : Converts a Fortran type/kind into a datatype
  enumerator from the system MPI library.

The ``lfric_comm_type``
^^^^^^^^^^^^^^^^^^^^^^^

This is a wrapper around the system MPI communicator and can be passed around
user code and contains two procedures:

* ``function get_comm_mpi_val()`` : Returns an integer implementation of the
  communicator that can be required when interfacing with external software.
* ``subroutine set_comm_mpi_val(comm)`` : This sets the ``lfric_comm`` object to
  be the same communicator that the given integer points at.


The ``lfric_mpi_type``
^^^^^^^^^^^^^^^^^^^^^^

This is a wrapper that sits around the system MPI library and provides message
passing functionality.

* ``subroutine initialise(in_comm)`` : Initialises the ``lfric_mpi``
  object based on the given communicator.
* ``subroutine finalise()`` : Finalises the ``lfric_mpi`` object
* ``function get_comm() result(communicator)`` : Returns the lfric communicator
  object.
* ``function is_comm_set() result(comm_state)`` : Returns whether a communicator
  has been set - i.e. whether the ``lfric_mpi`` object has been initialised
* ``subroutine global_sum(l_sum, g_sum)`` : All parallel tasks provide a local
  value in ``l_sum`` and the global sum of these values will be returned in
  ``g_sum``. This subroutine can be used with 32-bit integers, 32-bit reals or
  64-bit reals.
* ``subroutine global_min(l_min, g_min)`` : All parallel tasks provide a local
  value in ``l_min`` and the global minimum of these values will be returned in
  ``g_min``. This subroutine can be used with 32-bit integers, 32-bit reals or
  64-bit reals.
* ``subroutine global_max(l_max, g_max)`` : All parallel tasks provide a local
  value in ``l_max`` and the global maximum of these values will be returned in
  ``g_max``. This subroutine can be used with 32-bit integers, 32-bit reals or
  64-bit reals.
* ``subroutine all_gather(send_buffer, recv_buffer, count)`` : Gather integer
  data from all MPI tasks into a single array on all MPI tasks. The data in
  ``send_buffer`` from the jth process is received by every process and placed
  in the jth element of the ``recv_buffer``.
* ``subroutine broadcast(buffer, count, root)`` : Broadcasts the information
  held in buffer on processor root to all other parallel tasks. The variable
  count gives the size of buffer (this should be omitted if buffer is a scalar).
  Broadcast can be used on logicals, 32-bit integers, 32-bit reals or 64-bit
  reals. Scalars, 1d, 2d and 3d arrays are supported. Broadcast can also be used
  on simple string variables.
* ``function get_comm_size()`` : Returns the number of parallel tasks in the
  current communicator.
* ``function get_comm_rank()`` : Returns the number of the current rank.
