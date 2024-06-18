.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section diagnostics:

Diagnostics
===========

Diagnostics are outputs from a model used to analyse the scientific
progress of its run. Commonly, a model has the capability to compute
very many diagnostics, but a user can choose just a subset of them to
output.

The LFRic infrastructure provides a framework for adding diagnostics
to a model. Currently, the framework is usually interfaced to the XIOS
library. In principle, a different output library could be used.

LFRic diagnostic support
------------------------

The LFRic infrastructure aims to support the follow principles for
outputting diagnostics:

#. Any field can be output as a diagnostic.
#. Each diagnostic output in a time-step should be uniquely
   identifiable.
#. A field can be output from anywhere within the time-step.
#. The model should not have to compute a diagnostic if it has not
   been requested.

The LFRic infrastructure supports these principles in the following way.

#. Each LFRic field has a generic output method that can be overloaded
   with a function to send the field name and data to the chosen
   output library.
#. Each field can be given a string name which can be the diagnostic
   name.
#. As the output method is a method of the field, it can be called
   from anywhere where the field is in scope.

The LFRic infrastructure cannot enforce the principle that each output
in a time-step is unique; the model design must ensure uniqueness. If
two fields are output that have the same name and the same time-step,
then the underlying output system cannot disambiguate them.

Writing out a diagnostic
------------------------

To output a field to a diagnostic, the field must be associated with
an output function. The following associateds the LFRic field
``my_diagnostic_field`` with the procedure ``my_diagnostic_method``.

.. code-block:: fortran

   use field_parent_mod,                only: write_interface
   use my_diagnostic_system_mod,        only: my_diagnostic_method
   ! <snip>

   procedure(write_interface), pointer :: write_behaviour

   write_behaviour => my_diagnostic_method
   call my_diagnostic_field%set_write_behaviour(write_behaviour)

Once the diagnostic field has been computed, the write procedure is
called to send it to the diagnostic system.

.. code-block:: fortran

   call my_diagnostic_field%write_field('my_diagnostic_field_name')

If the string is not supplied, the name of the field will be used.

The ``write_field`` method ``my_diagnostic_method`` will be called
with the name of the field and the field proxy. The field proxy holds
pointers to the data and to other information about the field. The
information is used to work out how to send the data to the output
system supported by the method.

.. _section optional diagnostics:

Optional diagnostics
--------------------

If a diagnostic is not requested and is not otherwise used by the
model, then to save memory and time it is beneficial to avoid
initialising the field or computing the data.

The code assumes that the ``write_behaviour`` has been defined
previously.

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

Sometimes, a complex science kernel may compute diagnostics alongside
computing prognostic fields, such that the kernel is always called
even when the diagnostics are not requested. When algorithms call
kernels the PSy layer code requires that all fields are
initialised. To save memory, the LFRic infrastructure allows fields to
be initialised without any field data. The following example
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

Any array can be used to override the field data. The above example
uses an array from a module with one element from a module. Because
the array is in a module, the kernel can use it to check whether a
field has been properly initialised, and can avoid computing fields
that are not.

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
--------------------------------

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

To illustrate, the following code block illustrates a situation where
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
-----------------

The above code examples demonstrate the LFRic diagnostic system using
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



[DELETE BELOW???]




A diagnostic field is a field used for sending data for output to a
file. Depending on the construction of an application, a user can
often select which diagnostics need to be output, and can choose times
to output the diagnostics. To save compute resources, the LFRic
infrastructure supports the capability to initialise fields and
compute diagnostics only if they are required on a given
time-step.

Diagnostic values are often computed local to a science routine. Therefore, as
long as the diagnostic output system can support output of diagnostics
from anywhere in the code, the diagnostic field to hold the value can
be allocated locally. The benefit of allocating a field locally is
that the application does not need to hold so many diagnostic fields
in memory at the same time.

In an ideal world, if a diagnostic is not requested (and if it is not
also required as a by-product of computing another diagnostic that is
requested) then it should neither be initialised nor be computed.
Code referencing the diagnostic should not be visited when running the
configuration. However, some scientific routines have APIs containing
a mix of diagnostics and prognostics, or a number of diagnostics where
it may be that only some are requested. Therefore, the model must be
able to correctly run code that references, but does not compute,
diagnostics that are not required.

Diagnostic fields computed by Kernels
-------------------------------------

Assuming an algorithm has correct logic it is possible to avoid
initialising a field for a diagnostic that has not been requested,
then avoid computing or outputting it. It is even possible, if not
recommended, to pass an uninitialised field to another algorithm as
long as the algorithm being called also avoids computing the
field. However, it is not possible to pass an uninitialised field to a
kernel: all fields must be initialised before passing to a kernel as
the PSy layer code inserted between the algorithm and kernel will
attempt to dereference information from the field.

Therefore, if a kernel that computes an optional diagnostic must be
called (either to evolve the model or to compute other optional
diagnostics) then all the diagnostic fields passed to the kernel must
be initialised.

To support memory efficiency, LFRic supports initialisation of a field
without allocating an associated data field.

::

    use empty_data_mod, only :: empty_array

    type(function_space_type), pointer :: vector_space
    type(field_type), pointer :: mandatory_diagnostic
    type(field_type), pointer :: optional_diagnostic

    ! Get pointer to function space required for the diagnostic field
    vector_space => function_space_collection%get_fs( twod_mesh, 0, W3, 1)
    ! Initialise the diagnostic fields
    call mandatory_diagnostic%initialise(vector_space)

    if (optional_diagnostic_requested) then
      call optional_diagnostic%initialise(vector_space)
    else
      ! Intialise the diagnostic without a data array
      call optional_diagnostic%initialise(vector_space, &
                                          override_data = empty_data )
    end if

    call invoke( compute_diagnostics_type( mandatory_diagnostic, &
                                           optional_diagnostic )

Instead of passing a pointer to a field array to the kernel, the PSy
layer code will pass a pointer to the empty data array (an array with
one element should be used as arrays with zero elements behave
unpredictably). Therefore, the kernel must not attempt to write to it!
If the ``empty_data`` module used by the algorithm is also used by the
kernel, then it is easy to add a check to the kernel code to see
whether the data array is to be computed:

::

    module compute_diagnostic_kernel_mod
    use empty_data_mod, only : empty_data

    <snip>

    if ( .not. associated(optional_diagnostic_data, empty_data) ) then
      ! Compute optional diagnostic
      <snip>
    end if

Alternatively, identical logical tests can be used in the algorithm
and the kernel to prevent diagnostics being computed when they are
initialised with dummy data.

Example: XIOS integration
-------------------------

Several LFRic applications use the XIOS library by integrating to the
:ref:`lfric-xios component <section lfric xios>`. For
diagnostic fields supported by the ``lfric_xios`` component it is
possible to infer the field type from the XIOS metadata, and to
query XIOS to determine whether a diagnostic is required for a given
time-step.

For the Momentum(R) atmosphere, an ``init_diag`` procedure has been
written to support initialisation of diagnostics. To add, compute and
output a new diagnostic, one can write the following:

::

    type(field_type), allocatable :: optional_diagnostic
    logical(l_def)                :: optional_diagnostic_flag

    optional_diagnostic_flag = init_diag(optional_diagnostic, &
                                         'optional_diagnostic_id')

    if (optional_diagnostic_flag) then
      ! Compute diagnostic
      <snip>
    end if

    if (optional_diagnostic_flag) then
      call optional_diagnostic%write_field('optional_diagnostic_id')
    end if

Where a field is not required for a given time-step, the ``init_diag``
procedure will initialise it with an empty data array as described
above.

The ``init_diag`` procedure should support many applications that use
the same meshes and physics spaces as used by the Momentum atmosphere,
with some modifications to support different types of multidata fields.

Clearly, for a large algorithm that computes a lot of diagnostics
alongside other processing, the initialise and writing processes can
be modularised into separate procedures.
