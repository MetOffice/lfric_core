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

When running in parallel, all logged messages will be collated into
files in the current directory with the pattern ``PETxxx.name.Log``
where ``xxx`` is the number of the MPI rank which generated the message
and ``name`` is a name defined by the application. When running a
serial job, log messages are sent to the standard output units.

Configuration
-------------

Important aspects of the logging system are configured when it is
initialised by an application. The :ref:`driver component<driver
component>` includes :ref:`support<driver logging>` for initialising
the logger from the application namelist
configuration.

The initialisation routine sets up the MPI communicator, sets the log
level, and opens the PET files for each rank. To reduce output, there
is an option to configure the logger to ignore messages except for
those written from rank 0.

If an application outputs a log messages before the logging system has
been initialised, the message is output to file but is prefixed by a
warning message that states that the system has not yet been
initialised.

As the uninitialised logging system has no knowledge of the
application's parallelism, log messages written by multiple ranks of a
parallel job will be written to the same output file and may overwrite
each other in a confusing manner. Therefore, it is recommended that an
application initialises the logging system as early as possible, and
that writing of log messages prior to its initialisation is kept to a
minimum, or only done on one rank.

Event Levels
------------

Like many logging systems the concept of a "level" is used. The
application defines a logging level, and each message is assigned a
level at the time it is raised. The logger compares the message level
with the level defined by the application. If the message level is
higher than the application level it will be output.  Otherwise it
will be silently ignored.

Supported levels from highest to lowest are:

* Always
* Error
* Warning
* Info
* Debug
* Trace

The "Error" level is special as it will not only emit a message but
will also terminate execution. This is used when an unrecoverable
issue has occurred.

The application logging level is set using ``log_set_level()``. If
this is not called it will default to "Info" level.

If the application logging level is "Debug" or "Trace", the output
unit will be flushed every time a log message is output. Flushing the
buffer can ensure that no or fewer messages go missing should the
model subsequently terminate with an error, but likely it also slows
down the progress of the application.

Typically, an application will have very few or no messages at the
"Always" level, which means that setting a log level of "Error" or
"Warning" results in a very small amount of output showing
job-terminating error messages or warning messages that indicate there
is something wrong happening.

An application should have few "Info" messages such that messages give
a flavour of the progress of a run without providing too much
output. The "Debug" and "Trace" can be used for steadily more in-depth
information about key processes that may be useful for monitoring the
progress of an application configuration that is giving incorrect
answers.

Timestep
--------

As it is intended for use with scientific models which work iteratively on time
the logging system understand the concept of "time steps." At the beginning of
each timestep the logger timestep can be updated so that it will
display the step as part of any message logged.

The ``model_clock_type`` objects include calls to update the logging
system timestep.


Log An Event
------------

Messages are output using the ``log_event`` procedure which
takes a message and the level of the message. The logger includes a
workspace character string that can be used to construct a formatted
message prior to calling ``log_event``.

.. code-block:: fortran

   call log_event("Outputting results", log_level_debug)
   write(log_scratch_space, '(A,I5,A,F12.5)') &
         "Cell ",cell_number,"   value: ", cell_value
   call log_event(log_scratch_space, log_level_debug)

Get the log level
-----------------

The ``log_level`` function returns the current log level. Where a log
message is outputting a value that requires some optional computation,
checking the log level of the message against the application log
level can be used to avoid the computation.

.. code-block:: fortran

   app_log_level = log_level()
   ! Avoid calculating average if the value is not going to be printed
   if (app_log_level <= log_level_debug) then
     average_value = calculate_average_value(field)
     write(log_scratch_space,'(A,F12.5)')'Average is ',average_value
     call log_event(log_scratch_space, log_level_debug)
   end if

Output Streams
--------------

By default the logger will send low level messages to the standard output and
high level messages (warnings and errors) to the standard error stream.

If you have alternative needs an I/O unit number can be passed to
``log_set_info_stream()`` or ``log_set_alert_stream()``.
