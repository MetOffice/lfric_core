.. ------------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _section field:

LFRic Fields
------------

An LFRic field holds data over the horizontal domain of a mesh. [Grab
stuff from the lfric data model to flesh this out].

Key commands and functions
==========================

To create a field, one first needs to construct a :ref:`function space
<section function space>` and a 3-dimensional mesh (note that a 3D
mesh with a single level may be referred to in the code as a
2-dimensional mesh). Optionally, the field can be named. The method
for creating a new field based on an existing mesh ``mesh_id`` and
function space ``W2`` is:

.. code-block:: fortran

   type(function_space_type), pointer :: vector_space
   type(field_type)                   :: wind_field

   ! Get a reference to a lowest order W2 function space
   vector_space => function_space_collection%get_fs(mesh_id, 0, W2)

   ! Create a field to hold wind data
   call wind_field%initialise(vector_space, name = "wind")

Once created, a field can be passed to a call to an ``invoke`` for
processing by a kernel or a PSyclone built-in.

A field can be constructed from another field as follows.

.. code-block:: fortran

   call old_field%copy_field_properties(new_field, name)

The name argument is optional. If no name is supplied, the new field
will be unnamed.

Note that ``copy_field_properties`` copies everything about the field
except the data and the name. If the data needs to be copied, then the
``setval_x`` built-in can be invoked. Initialising new fields like
this allows PSyclone to optimise the copy.

.. code-block:: fortran

   call invoke( setval_x(new_field, old_field) )
