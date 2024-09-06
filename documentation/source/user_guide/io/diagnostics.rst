.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section diagnostics:

Diagnostics
===========

Diagnostics are outputs from a model application used to analyse the
scientific progress of its run. Commonly, a model can compute very
many diagnostics, but a user can choose just a subset of them to
output.

The LFRic infrastructure provides a framework for linking the model to
a diagnostic processing system. Currently, LFRic applications usually
send the diagnostic data to the XIOS library. In principle, a
different output library could be used without making changes to a
model.

LFRic diagnostic support
------------------------

The LFRic infrastructure aims to support the following principles for
outputting diagnostics:

#. Any field can be output as a diagnostic.
#. A field can be output from anywhere within the model code.
#. The model should not have to compute a diagnostic if it has not
   been requested.
#. It should be possible to choose a different diagnostic processing
   system without significant changes to a model.

The LFRic infrastructure supports these principles in the following way.

#. Each LFRic field has a single generic output method.
#. As the output method is a method of the field, it can be called
   from anywhere where the field is in scope.
#. The output method of fields calls a procedure pointer, permitting
   an application to choose which diagnostic processing application to
   use for each diagnostic. Details on creating output methods are
   found in the :ref:`write_interface procedures <section write
   interface procedures>` section.

Writing out a diagnostic
++++++++++++++++++++++++

An LFRic field has a ``write_field`` method which sends the field to
the diagnostic system via a ``write_interface`` procedure. The
``write_interface`` can be chosen by the application. Its role is to
package the data from the field and send it to the chosen diagnostic
output system or library.

The following links the write method of a particular LFRic field
``my_diagnostic_field`` with a ``write_interface`` procedure called
``send_diagnostic`` which would link to a particular diagnostic
system.

.. code-block:: fortran

   use field_parent_mod,                only: write_interface
   use send_diagnostic_mod,             only: send_diagnostic
   ! <snip>

   procedure(write_interface), pointer :: write_behaviour

   write_behaviour => send_diagnostic
   call my_diagnostic_field%set_write_behaviour(write_behaviour)

Once the diagnostic field has been computed, the ``write_field``
method is called, triggering the sending of the field to the chosen
diagnostic system:

.. code-block:: fortran

   call my_diagnostic_field%write_field('my_diagnostic_field_name')

The ``write_field`` method calls ``send_diagnostic`` with the name of
the field and the :ref:`field proxy <section field proxy>`. The field
proxy holds pointers to the data and to other information about the
field. The field proxy information should be sufficient for
``send_diagnostic`` to work out how to pass the diagnostic data to the
diagnostic processing library it links to.

The field name is an optional argument of ``write_field``: if a string
is not supplied, the name of the field will be used.

.. _section optional diagnostics:

Optional diagnostics
++++++++++++++++++++

If a diagnostic is not requested and is not otherwise used by the
model, then to save memory and time it is beneficial to avoid
initialising the field or computing the data.

The code assumes that ``write_behaviour`` has already been assigned to
a ``write_interface``.

The first example avoids any code relating to the field if the
diagnostic is not required:

.. code-block:: fortran

   type(field_type) :: my_diagnostic_field
   type(function_space_type), pointer :: vector_space

   if (my_diagnostic_flag) then
     ! Diagnostic has been requested
     vector_space => function_space_collection%get_fs( mesh, element_order, W3 )
     call my_diagnostic_field%initialise(vector_space, name='my_diagnostic_name')
     call my_diagnostic_field%set_write_behaviour(write_behaviour)

     call invoke(compute_my_diagnostic_type(my_diagnostic_field))

     call my_diagnostic_field%write_field()
  end if

In the second example, it is not possible to avoid calling code that
references the diagnostic when it is not requested. This scenario can
occur when a complex science kernel computes diagnostics alongside
computing prognostic fields.

When algorithms call kernels, the PSy layer code requires that all
fields are initialised. To save memory, the LFRic infrastructure
allows fields to be initialised without any field data, meaning that
the field takes up a minimal amount of memory. The following example
illustrates the approach:

.. code-block:: fortran

   use empty_data, only, : empty_real_data

   ! Function space for the diagnostic field
   vector_space => function_space_collection%get_fs( mesh, element_order, W3 )
   if (my_diagnostic_flag) then
     ! Diagnostic has been requested
     call my_diagnostic_field%initialise(vector_space, name='my_diagnostic_name')
     call my_diagnostic_field%set_write_behaviour(write_behaviour)
   else
     ! Diagnostic is not required
     call my_diagnostic_field%initialise(vector_space, name='my_diagnostic_name' &
                                         override_data = empty_real_data)
   end if

   call invoke(big_science_kernel_type([lots of fields], my_diagnostic_field))

   if (my_diagnostic_flag) then
     call my_diagnostic_field%write_field()
   end if

When a field is initialised with override data, then the override data
array takes the place of the field data array, and no memory is
allocated to hold any field data, thus saving memory.

The above example uses an array from a module. This allows the kernel
to use the same module, which enables it to check whether fields being
passed in are pointing to the override data array. If a field is
pointing to the override array then it does not need to be, and should
not be, computed:

.. code-block:: fortran

   subroutine big_science_kernel_type([lots of fields], my_diagnostic_data)
     use empty_data, only, : empty_real_data

     ! <snip>

     if ( .not. associated(my_diagnostic_data, empty_real_data) ) then
       ! Field is properly allocated so will be computed
       call compute_my_diagnostic(my_diagnostic_data)
     end if

This approach saves having to pass an extra logical into the kernel.

Diagnostics from existing fields
++++++++++++++++++++++++++++++++

For reasons described above, the same field `name` should not be
written out as a diagnostic twice in one time-step, but the same
`field` can be written out as a diagnostic as long as a different name
is used in each case.

This may occur when interim values of a field need to be written out
from different parts of the model.

In the code examples above, it is implicitly assumed that the
underlying function uses the field name to identify the field. But the
``write_field`` method takes a field name as an optional argument,
which can override the field name.

The following code block illustrates a situation where
one might want to output a diagnostic from the same field before and
after a kernel has processed it:

.. code-block:: fortran

   subroutine science_algorithm(my_diagnostic_field)

     type(field_type), intent(inout) :: my_diagnostic_field

     call my_diagnostic_field%write_field('my_diagnostic_at_start')

     call invoke(compute_my_diagnostic_type(my_diagnostic_field))

     call my_diagnostic_field%write_field('my_diagnostic_at_end')

   end subroutine science_algorithm


Enhanced approach
+++++++++++++++++

The above code demonstrate the LFRic diagnostic system using
simple examples where fields are initialised and named with hard-wired
choices. The LFRic infrastructure includes an interface to the XIOS
system, and features of this system can enable diagnostic-writing code
to be simplified from the perspective of a science model developer.

For example, the diagnostic configuration supplied to XIOS by a model
run provides information about which diagnostics are requested on
which time-steps, and what the output format of the diagnostics will
be.

Knowledge about which time-step a diagnostic is output can be used to
set the ``my_diagnostic_flag`` used in the :ref:`Optional
Diagnostics<section optional diagnostics>` section above.

Knowledge about the function space that the field lives on can be
inferred from the output format of the diagnostic.

These two aspects can be combined into a single generic function,
illustrated by a rewrite of the second code example in the
:ref:`Optional Diagnostics<section optional diagnostics>` section, as follows:

.. code-block:: fortran

   my_diagnostic_flag = init_diag(my_diagnostic_field, 'my_diagnostic_name')

   call invoke(big_science_kernel_type([lots of fields], my_diagnostic_field))

   if (my_diagnostic_flag) then
     call my_diagnostic_field%write_field()
   end if

Such a function has been written to support the LFRic atmosphere. See
the LFRic atmosphere documentation for more details.

.. _section write interface procedures:

write_interface procedures
++++++++++++++++++++++++++

The ``write_interface`` procedure acts as the interface between LFRic
and a diagnostic system. Its abstract interface is defined as follows.


.. code-block:: fortran

   subroutine write_interface(field_name, field_proxy)
     import field_parent_proxy_type
     character(len=*), optional,      intent(in) :: field_name
     class(field_parent_proxy_type ), intent(in) :: field_proxy
   end subroutine write_interface

Sending a field to a diagnostic processing system involves sending
various pieces of information about the field. The diagnostic name and
the data are clearly required, but also the size and dimensions of the
field may be required. Furthermore, the data order of the data in the
field may not be the same as the data order expected by the diagnostic
system, so the data may need to be reordered in some way.

Providing the ``write_interface`` with the field proxy provides it
with access to all information about the field, which should be
sufficient to pass it on to a diagnostic system.


Example: XIOS integration
-------------------------

Several LFRic applications use the XIOS library as a diagnostic
processing system by integrating to the :ref:`lfric_xios component
<section lfric xios>`. See the `diagnostics documentation
<https://metoffice.github.io/simulation-systems/WorkingPractices/Diagnostics/lfric_diagnostics.html>`_
in the lfric_apps repository for more usage examples.

Initialising fields for lfric_xios
++++++++++++++++++++++++++++++++++

For diagnostic fields supported by the ``lfric_xios`` component it is
possible to infer the field type (its function space) from the XIOS
metadata, and to query XIOS to determine whether a diagnostic is
required for a given time-step.

For the Momentum(R) atmosphere, an ``init_diag`` procedure has been
written to support initialisation of diagnostics. To add and compute a
new diagnostic, one can write the following:

.. code-block:: fortran

    type(field_type), allocatable :: optional_diagnostic
    logical(l_def)                :: optional_diagnostic_flag

    optional_diagnostic_flag = init_diag(optional_diagnostic, &
                                         'optional_diagnostic_id')

    if (optional_diagnostic_flag) then
      ! Compute diagnostic
      !<snip>
    end if

Writing fields for lfric_xios
+++++++++++++++++++++++++++++

As the field type of a field can be inferred from its metadata, so can
the method required to send the field to the XIOS library. Therefore,
all diagnostic fields that use the ``lfric_xios`` component are
assigned the same ``write_interface``: the ``init_diag`` contains:

.. code-block:: fortran

   write_behaviour => write_field_generic
   call field%set_write_behaviour(write_behaviour)

After the diagnostic is computed, the``write_field`` method is called,
but only if the diagnostic flag was set to ``.true.``.

.. code-block:: fortran

    if (optional_diagnostic_flag) then
      call optional_diagnostic%write_field()
    end if

As with the previous example, the ``write_field`` procedure calls
``write_field_generic`` with the field name and the field proxy. The
``write_field_generic`` procedure reformats the data based on the
field properties. For example, XIOS has been configured to expect
data in level-first data order whereas LFRic fields are column-first
data order, so the procedure reorders the data appropriately.
