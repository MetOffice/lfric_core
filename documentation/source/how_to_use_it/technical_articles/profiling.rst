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


Start and stop
~~~~~~~~~~~~~~

At the most basic level, a stopwatch is designed to start and stop, telling
the time between those two actions. This is how the stopwatch type is used.

``start`` will begin timing. ``stop`` will end timing and output the elapsed
time to the :ref:`logger<logging>`.

``start`` has a required argument ``name``, a character array that will be
printed with the elapsed time by the logger.


Pause and resume
~~~~~~~~~~~~~~~~

A pair of procedures are provided for suspending the stopwatch temporarily.

``pause`` instructs the stopwatch instance to start recording a new ``paused``
region, which will be subtracted from the total time recorded by the stopwatch.

``resume`` causes the stopwatch to stop timing the ``paused`` region. Since
the paused time is cumulative, a stopwatch instance can be paused and resumed
multiple times.


Output
~~~~~~~~~~~~~~~~~~

The stopwatch output always uses the following format:

.. code-block
    (STOPWATCH) Time taken for <stopwatch_name> : <elapsed_time> (s)


