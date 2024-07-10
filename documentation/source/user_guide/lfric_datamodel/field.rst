.. ------------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _section field:

LFRic Fields
------------

An LFRic field holds data over the horizontal domain of a mesh. Its
design supports the LFRic separation concern by preventing direct
access to the data. Model manipulation of data should only be done by
passing the field to a kernel or PSyclone built-in.

[TODO Grab stuff from the lfric data model to flesh this out].

Initialising new fields
=======================

To create a field, first construct a :ref:`function space <section
function space>` and a 3-dimensional mesh (note that a 3D mesh with a
single level may be referred to in the code as a 2-dimensional
mesh). The method for creating a new field based on an existing mesh
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
supplied, the new field will be unnamed.

There is an LFRic function called ``copy_field_serial`` which copies
the field properties `and` the field data. However, as the name
suggests, any such copy would be done serially and would not take
advantage of any shared memory parallelism. Therefore, use of
``copy_field_serial`` is not advised. If the data needs to be copied,
then the ``setval_x`` built-in can be invoked after the field is
initialised. Initialising new fields with ``setval_x`` allows PSyclone
to optimise the copy.

.. code-block:: fortran

   call wind_field%copy_field_properties(new_field, name = "wind_copy")
   call invoke( setval_x(new_field, wind_field) )

The function space and mesh used to initialise a field have a
particular halo depth. By default, a field is initialised with the
same halo depth, but optionally, a smaller halo depth can be
requested:

.. code-block:: fortran

   ! Create a field to hold wind data
   call wind_field%initialise(vector_space, name = "wind", halo_depth = 1)

The function will fail if the requested halo depth is larger than the
function space halo.

The field_proxy object
----------------------

The data held in a field is private meaning it cannot be accessed
using field methods. Clearly the data does need to be accessed
somewhere in the code: the field proxy provides the methods for doing
so. But the field proxy object must be used with care to maintain the
integrity of the application's data.



Keeping the data private within the field is a way of enforcing the
PSyKAl design that underpins key LFRic applications. Critically, the
application needs to understand when and how data is modified so that
it can implement strategies correctly for applying redundant
computation into halos and calling halo swaps. If data is meddled with
outside of the PSyKAl design then errors may be introduced into the
application.

The field proxy object may be used in the following limited
circumstances:

 #. For writing PSyKAl-lite code. PSyKAl-lite code represents
    hand-written PSy layer code where PSyclone does not support your
    requirement. The PSy layer does access field information using the
    field proxy.
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
