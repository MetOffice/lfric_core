.. ------------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _section field:

LFRic Fields
------------

This section provides an overview of fields and how they should be
initialised and used in an application.

An LFRic field holds data over the horizontal domain of a mesh. Its
design supports the LFRic separation concern by preventing direct
access to the data. Model manipulation of data should only be done by
passing the field to a kernel or PSyclone built-in.

The data held by a field represents one of a number of :ref:`function
spaces <section function space>`. The function space describes the
layout of data points on each 3-dimensional cell. A field can
represent 32-bit real, 64-bit real or integer data, and can hold more
than one set of data of same function space and data type.

Initialising new fields
=======================

To create a field, first construct a :ref:`function space <section
function space>` and a 3-dimensional mesh (note that a 3D mesh with a
single level may be referred to in the code as a 2-dimensional
mesh). The code for creating a new field based on an existing mesh
``mesh_id`` and function space ``W2`` is:

.. code-block:: fortran

   type(function_space_type), pointer :: vector_space
   type(field_type)                   :: wind_field

   ! Get a reference to a lowest order W2 function space
   vector_space => function_space_collection%get_fs(mesh_id, 0, W2)

   ! Create a field to hold wind data
   call wind_field%initialise(vector_space, name = "wind")

The name is an optional argument. A name would be required if passing
the field to other parts of the infrastructure. For example, fields
added to field collections would require a name so that they can be
referenced later.

Once created, a field can be passed to a call to an ``invoke`` for
processing by a kernel or a PSyclone built-in.

A field can be constructed from another field as follows.

.. code-block:: fortran

   call wind_field%copy_field_properties(new_field, name = "wind_copy")

This call initialises the new field with the same mesh and function
space, but it does not initialise the data. If no name argument is
supplied, the new field will be unnamed rather than adopting the name
of the original field.

.. warning::

   There is an LFRic function called ``copy_field_serial`` which
   copies the field properties `and` the field data. However, as the
   name suggests, any such copy would be done serially and would not
   take advantage of any shared memory parallelism. Therefore, use of
   ``copy_field_serial`` is not advised. If the data needs to be
   copied, then the ``setval_x`` built-in can be invoked after the
   field is initialised. Initialising new fields with ``setval_x``
   allows PSyclone to optimise the copy.

   .. code-block:: fortran

      call wind_field%copy_field_properties(new_field, name = "wind_copy")
      call invoke( setval_x(new_field, wind_field) )

The function space and mesh used to initialise a field have a
particular halo depth. By default, a field is initialised with the
same halo depth. Optionally, a smaller halo depth can be requested:

.. code-block:: fortran

   ! Create a field to hold wind data
   call wind_field%initialise(vector_space, name = "wind", halo_depth = 1)

The function will fail if the requested halo depth is larger than the
function space halo.

The field_proxy object
----------------------

The data held in a field is private, meaning it cannot be accessed
using field methods. Clearly, the data does need to be accessed
somewhere in the code and the field proxy provides the methods for
doing so. The field proxy object must be used with care to maintain
the integrity of the application's data.

Keeping the data private within the field is a way of enforcing the
PSyKAl design that underpins key LFRic applications. The application
needs to monitor the status of halos: whether or not they are "dirty"
or out of date with the corresponding owned data points on the
neighbouring ranks. PSyclone generates code that does this monitoring
correctly. If additional code is using and modifying data without
PSyclone's knowledge, the data can become inconsistent.

The field proxy object may be used by application writers in the
following limited circumstances:

 #. For writing PSyKAl-lite code. PSyKAl-lite code represents
    hand-written PSy layer code where PSyclone does not support your
    requirement. The PSy layer accesses field information using the
    field proxy so it can be passed to kernels. Ideally, PSykal-lite
    code should be written in a style that, plausibly, PSyclone
    `could` generate if it were extended to support the new
    requirement.
 #. For writing an :ref:`external field <section external field>`
    interface to copy data between the LFRic application and another
    application.
 #. For debugging purposes, or within unit or integration tests.

If the field proxy is to be used, a good understanding of the
:ref:`distributed memory design <section distributed memory>` is
required so that code maintains the integrity of the data and its
halos. For example, if the data is updated in such a way that the
halos may be inconsistent with data on neighbouring partitions, then
either a halo swap needs to be performed or the halos needs to be
marked as dirty.

Data can be accessed using the proxy as follows:

.. code-block:: fortran

   real(r_def), pointer :: wind_field_data(:)
   type(field_proxy_type) :: wind_field_proxy

   wind_field_proxy = wind_field%get_proxy()
   wind_field_data => wind_field_proxy%data

.. _section mixed precision field:

Mixed precision fields
======================

The ``field_type`` object referenced in a lot of code examples found
in the documentation is either a 32-bit or a 64-bit field. The choice
of precision is made at build-time: the default is 64-bit, but 32-bit
can be chosen by setting compile def ``RDEF_PRECISION`` to 32. See the
``field_mod`` module for how the ``field_type``, and the matching
``field_proxy_type``, precision are defined. Key parts are shown here:

.. code-block fortran

   module field_mod

   #if (RDEF_PRECISION == 32)
   use field_real32_mod, only: field_type         => field_real32_type, &
                               field_proxy_type   => field_real32_proxy_type
   #else
   use field_real64_mod, only: field_type         => field_real64_type, &
                               field_proxy_type   => field_real64_proxy_type
   #endif

   implicit none
   private

   public :: field_type, &
             field_proxy_type

   end module field_mod

The choice of compile def will point ``field_type`` fields to one of
two concrete implementations of the field object:
``field_real32_type`` or ``field_real64_type``. Similarly, there are
32-bit and 64-bit versions of the ``field_proxy_type``.

The choice made at build-time applies to all ``field_type`` variables.
Where an application requires a combination of 32-bit and 64-bit
fields an application can define additional field types that are
controlled by separate compile defs. The science code has to be
written such that code in one part of the application uses the
different field type definitions.

Different field type objects are made available by taking a copy of
the ``field_mod`` module, changing the name of the public types and
ensuring they are configured by a different compile def. With the
appropriate configuration of three individual compile defs, and
definitions of two additional field types, each of the fields declared
in the following code can be either 32-bit or 64-bit:

.. code-block fortran

   type(field_type)            :: wind_field
   type(r_tran_field_type)     :: dry_mass
   type(r_solver_field_type)   :: theta_advection_term

.. attention::

   The ``field_mod`` module also declares a ``field_pointer_type``
   which points to a field pointer of the chosen default precision.
   The ``field_pointer_type`` is used in ``select type`` operations
   within the infrastructure code for field collections as field
   collections can hold both fields and pointers to fields.

   Field collections support the ability to access an individual named
   field and also the ability to iterate over all the fields in the
   field collection. Support of iterators, however, is only available
   for fields with the default precision as defined in the
   ``field_mod`` module.

Integer fields
==============

The infrastructure supports 32-bit integer fields:
``integer_field_type``. There creation and usage is essentially the
same as for real fields. One key difference is that real fields and
integer fields have their own set of PSyclone built-ins.

Currently, there is no known requirement for 64-bit integer fields, so
a 64-bit integer field has not been made available.

Multidata fields
================

Multidata fields hold more than one quantity on the same mesh and
function space. The number and list of quantities is defined by the
application. Illustrative examples used in the Momentum\ :sup:`®`
atmosphere model are multidata fields that contain fields for several
different vegetation or surface types. But the Momentum\ :sup:`®`
model also uses multidata fields to store data at different soil
levels (rather than have a 3D mesh representing soil layers) and a
12-item multidata field to store monthly climatology data.

Multidata fields are created by first obtaining a multidata function
space.

.. code-block:: fortran

   integer, intent(in)                :: surface_tiles
   type(function_space_type), pointer :: fspace_surface_tiles
   type(field_type)                   :: canopy_water

   ! Get a reference to a lowest order W3 multidata function space
   fspace_surface_tiles =>
        function_space_collection%get_fs(mesh2D, 0, W3, surface_tiles)

   ! Create a field to hold wind data
   call canopy_water%initialise(fspace_surface_tiles, name = "canopy_water")

Multidata fields can be used when an array of fields needs to be
passed into a kernel. Within the kernel, data on each point in the
mesh is contiguous in memory.

To illustrate, the following kernel code loops over the field types
within the inner-most loop, where the kernel is called to operate on a
single column of the mesh.

The ``map(dof))`` reference points to the data point for one of the
dofs of the cell at the first level of the 3D mesh. For each level
up, the code steps ``nfield_types`` points.

.. code-block:: fortran

   do levels = 0, nlayers-1
     do dof = 1, ndofs_per_cell
       do field_type = 0, nfield_types-1
          data(map(dof) + levels * nfield_types + field_type) = ...
       end do
     end do
   end do
