.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section diagnostics:

Diagnostic Fields
=================

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
