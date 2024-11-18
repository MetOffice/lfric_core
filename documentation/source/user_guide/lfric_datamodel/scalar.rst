.. ------------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _section scalar:

LFRic scalars
=============

An LFRic scalar is a data structure that holds a scalar value but is
aware of the distributed memory domain decomposition.

Unlike fields, the value of a scalar is public.

To set the value of a scalar field:

.. code-block:: fortran

   type(scalar_type)       :: my_value

   my_scalar%value = my_value

Assuming a kernel is called to calculate the value of a scalar, the
resulting value will apply to the local domain. The LFRic scalar has
methods to compute sums, and to obtain the minimum or maximum value of
the scalar across all the domains of the model.

.. code-block:: fortran

   sum = my_scalar%get_sum()
   min_value = my_scalar%get_min()
   max_value = my_scalar%get_max()

Scalars are supported for 32-bit and 64-bit real values and for 32-bit
integer values. As with fields, an application can support
configurable precision by defining a ``scalar_type`` whose setting
depends on compile defs:

.. code-block:: fortran

   #if (RDEF_PRECISION == 32)
      use scalar_real32_mod, only: scalar_type       => scalar_real32_type
   #else
      use scalar_real64_mod, only: scalar_type       => scalar_real64_type
   #endif
