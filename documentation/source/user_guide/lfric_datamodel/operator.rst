.. ------------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _section operator:

LFRic operators
===============

An LFRic operator is a data structure that maps fields from one
function space to another. Like fields, operators for 64-bit and
32-bit reals are supported.

To initialise an operator that maps from function space ``fs_source``
to ``fs_target``:

.. code-block:: fortran

   type(operator_real64_type) :: my_operator

   call my_operator%initialise(fs_target, fs_source)

Operators can be copied by using the ``deep_copy`` method:

.. code-block:: fortran

   operator_copy = my_operator%deep_copy()

All application of operators to fields is done by kernel.

Operator proxies
----------------

Like fields, to impose the separation of concerns between algorithms
and kernels, the data in an operator is hidden from the operator
object and must be accessed using a proxy object.

The principle of operator proxies is largely the same as for field
proxies. But while there are limited circumstances where a developer
`may` need to access field data outside of the PSy layer, the
likelihood of doing this for operators is minuscule. For this reason,
this section on user documentation will not go into further details.
