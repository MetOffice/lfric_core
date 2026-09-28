.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _profiling:

Stopwatch
=========

The stopwatch type is designed for timing each independent call to a
subroutine or function. Upon destruction or a manual call to stop timing, the
stopwatch will output the recorded time and the name of the stopwatch to the
:ref:`logger<logging>`.

The stopwatch type is designed to use one instance of the type for each
watch section. Overlapping watch sections cannot use the same stopwatch
instance since the watch is reset when calling ``start``.


Start and stop
~~~~~~~~~~~~~~

At the most basic level, a stopwatch is designed to start and stop, telling
the time between those two actions. This is how the stopwatch type is used.

``start`` will begin timing. ``stop`` will end timing and output the elapsed
time to the :ref:`logger<logging>`.

``start`` has a required argument ``name``, a character array that will be
printed with the elapsed time by the logger.

.. code-block:: fortran

    type(stopwatch_type) :: stopwatch_1

    call stopwatch_1%start("My Example Watch")

    ! Code to be timed

    call stopwatch_1%stop()


Pause and resume
~~~~~~~~~~~~~~~~

A pair of procedures are provided for suspending the stopwatch temporarily.

``pause`` instructs the stopwatch instance to start recording a new ``paused``
region, which will be subtracted from the total time recorded by the stopwatch.

``resume`` causes the stopwatch to stop timing the ``paused`` region. Since
the paused time is cumulative, a stopwatch instance can be paused and resumed
multiple times.

.. code-block:: fortran

    type(stopwatch_type) :: stopwatch_1

    call stopwatch_1%start("My Example Watch")

    ! Code to be timed

    call stopwatch_1%pause()

    ! Code not timed

    call stopwatch_1%resume()

    call stopwatch_1%stop()

The stopwatch can be stopped without calling ``resume``, as this will be
handled by the stopwatch anyway.

.. code-block:: fortran

    ...

    call stopwatch_1%pause()

    ! Code not timed

    ! call stopwatch_1%resume() - not called

    call stopwatch_1%stop()

Output
~~~~~~

The stopwatch output always uses the following format:

.. code-block:: text

    (STOPWATCH) Total time for <stopwatch_name> : <elapsed_time> (s)


Reusing Watches
~~~~~~~~~~~~~~~~

If you intend to use a stopwatch more than once (not including pausing), ensure
that the two regions do not overlap. When calling ``start``, the stopwatch
is reset to the default values, which also calls ``stop`` on the watch if
it is running.

If you do use the same stopwatch for multiple regions, good usage looks like
this:

.. code-block:: fortran

    type(stopwatch_type) :: stopwatch_1

    call stopwatch_1%start("Region One")

    ! Region One stops
    call stopwatch_1%stop()

    ! Gap

    ! Region Two starts
    call stopwatch_1%start("Region Two")

    call stopwatch_1%stop()
