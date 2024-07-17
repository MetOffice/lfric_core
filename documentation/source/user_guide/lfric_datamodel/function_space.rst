.. ------------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _section field:

LFRic Function spaces
---------------------

LFRic function spaces are used to map data points and finite element
basis functions onto a domain represented by a 3D mesh. A high level
introduction into :ref:`LFRic function spaces<section function space
intro>` describes this mapping in more detail. This section is focused
on how to create and use function spaces based on existing function
space types. If new fundamental function spaces types are required,
refer to the developer documentation.

Initialising a function space is a prerequisite to initialising a
field: field initialisation allocates the correct amount of data to
represent data on all the cells of the mesh.

A new LFRic function space is initialised as follows:

.. code-block:: fortran

   type(function_space_type), pointer :: vector_space

   ! Get a reference to a lowest order W2 function space
   vector_space => function_space_collection%get_fs(mesh_id, element_order, fs_type)

The ``mesh_id`` identifies the 3D mesh that the function space must
cover. The ``element_order`` is typically ``0`` for lowest order or
``1`` for next-to-lowest order. The ``fs_type`` refers to one of a
number of available function space types. A `high-level description
<https://psyclone.readthedocs.io/en/stable/dynamo0p3.html#supported-function-spaces>`_
of supported function spaces can be found in the PSyclone
documentation. Refer to GungHo documentation for a comprehensive
description of each function space.

The function space constructor takes two optional arguments:

#. Specifying an integer ``ndata`` allows the creation of
   :ref:`multidata fields<section multidata field> with ``ndata`` dofs
   per dof location. The default value is `1`.
#. Spacifying ``ndata_first = .true`` allows creation of fields that
   with data ordered layer-by-layer instead of column-by-column.
