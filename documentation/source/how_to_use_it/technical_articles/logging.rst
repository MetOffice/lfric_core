.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _logging:

Logging
=======

A simple logging framework is provided by the ``log_mod`` module found in the
``utilities`` source directory. It should be used exclusively in favour of
``print`` or ``write`` to ``output_unit`` or ``error_unit``.

All logged messages will be collated into files in the current directory with
the pattern ``PETxxx.Log``` where `xxx` is the MPI rank which generated the
message.

Event Levels
~~~~~~~~~~~~

Like many logging systems the concept of a "level" is used. Each message is
assigned a level at the time it is raised. The logger compares this level with
its own level and if the message is more important it will be emmitted.
Otherwise it will be silently ignored.

Supported levels are:

* Error
* Warning
* Info
* Debug
* Trace

The first, "Error" is special as it will not only emmit a message but will also
terminate execution. This is used when an unrecoverable issue has occurred.

The logger's level is set using ``log_set_level()``. If this is not called it
will default to "Info" level.

Timestep
~~~~~~~~

As it is intended for use with scientific models which work iteratively on time
the logging system understand the concept of "time steps." At the beginning of
each timestep the logger is updated and will display the step as part of any
message logged.

This is handled automatically by the ``model_clock_type`` object.


Log An Event
~~~~~~~~~~~~

It may come as no surprise that logging an event is achieved using the
``log_event`` procedure. This takes a message and the level of the message.


Output Streams
~~~~~~~~~~~~~~

By default the logger will send low level messages to the standard output and
high level messages (warnings and errors) to the standard error stream.

If you have alternative needs an I/O unit number can be passed to
``log_set_info_stream()`` or ``log_set_alert_stream()``.
