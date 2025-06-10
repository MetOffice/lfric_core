.. -----------------------------------------------------------------------------
    (c) Crown copyright 2024 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _simple io:

Simple I/O
==========

Normally there will never be a need for applications to manage I/O at the
lowest levels. Instead they will make use of the Configurator to access values
from namelists and :ref:`XIOS integration <lfric xios component>` to obtain
field data.

If you feel you might like to access files directly take into account that the
existing methods encapsulate a lot of development and thought. In particular
they are optimised for parallel running, e.g. read once and broadcast to all
ranks. Doing it yourself can lead to problems.

If you do need to do it yourself a module ``io_utility_mod`` is provided which
manages I/O unit numbers. This is critical to avoiding collisions with unit
numbers already opened by other parts of the code.


Catch and Release
~~~~~~~~~~~~~~~~~

A pair of procedures are provided for accessing I/O units.

``claim_io_unit`` is a function which returns a unit number and removes it from
the pool of potential unit numbers.

``release_io_unit`` gives the unit number back to the pool once you are done
with it. Do not use the number after you have released it.


Open and Shut Case
~~~~~~~~~~~~~~~~~~

A couple of cenvenience procedures are provided which not only claim an I/O
unit but use it with a file.

``open_file`` will claim a unit number or use the one provided and open a file
on that unit for reading. All associated error checking is performed for you.

``close_file`` closes the file on a given unit and releases the unit back to
the pool along with all associated error checking.


Convenient Reading
~~~~~~~~~~~~~~~~~~

There is even a function ``read_line`` which reads one text line from an I/O
unit into a buffer. It performs error checking and returns a logical value
indicating whether the end of the file has been reached which allows the call
to be easily included in a loop over all lines in the file.

