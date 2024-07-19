.. ------------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _section function space:

LFRic Function spaces
---------------------

LFRic function spaces map data points and finite element basis
functions onto a domain represented by a 3D mesh. A high level
introduction into :ref:`LFRic function spaces<section function space
intro>` describes the mapping in more detail. This section is focused
on how to create and use function spaces based on function space types
supported by the infrastructure. If new fundamental function space
types are required, refer to the developer documentation.

Initialising an LFRic function space is a prerequisite to initialising
a field: field initialisation allocates the correct amount of data to
represent data for the function space type on all the cells of the
mesh.

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

The function space constructor can take two optional arguments:

#. Specifying an integer ``ndata`` allows the creation of
   :ref:`multidata fields<section multidata field>` with ``ndata`` dofs
   per dof location. The default value is `1`.
#. Spacifying ``ndata_first = .true`` allows creation of fields that
   with data ordered layer-by-layer instead of column-by-column.

Quadrature rules
----------------

Quadrature rules define how function space basis functions will be
integrated when kernels are called. This section does not describe the
mathematical process of using quadrature rules to apply basis
functions. Suffice to say, the infrastructure includes a number of
different quadrature rule objects. An appropriate rule is selected
from the ones supported by the infrastructure. A quadrature type is
called with the rule. Using the rule, the quadrature type calculates
weights for a defined set of points in a reference cell.

PSyclone will generate the code that calls the chosen quadrature type
to compute the basis functions at the required set of points defined
by the type, applying the weights defined by the quadrature rule.

A kernel can request quadrature for several function spaces and for
either or both basis functions or differential basis functions.

For example, the following creates a type based on a Gaussian
rule. The ``quadrature_xyoz_type`` function defines a set of points in
the cell based on the rule and the number ``nqp`` passed to the
function.

.. code-block:: fortran

   type( quadrature_rule_gaussian_type ) :: quadrature_rule

   qr = quadrature_xyoz_type(nqp, quadrature_rule)

The rule ``qr`` can be passed to a kernel whose metadata defines that it
accepts this form of quadrature rule:

.. code-block:: fortran

   type, public, extends(kernel_type) :: compute_total_pv_kernel_type
     private
     type(arg_type) :: meta_args(8) = (/              &
        arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),  &
     <snip>
     type(func_type) :: meta_funcs(4) = (/                &
        func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS),  &
        func_type(W0,          GH_DIFF_BASIS),            &
        func_type(W1,          GH_BASIS),                 &
        func_type(W3,          GH_BASIS)                  &
        /)
     integer :: operates_on = CELL_COLUMN
     integer :: gh_shape = GH_QUADRATURE_XYoZ
   contains
      procedure, nopass :: compute_total_pv_code
   end type

Prior to calling the kernel, PSy layer code will compute the basis and
differential basis functions for each of the function spaces defined
in the ``meta_funcs`` type, and pass the values to the kernels.

The following rules are supported:

+---------------------+--------------------------------------+
| Rule name           | Object                               |
+=====================+======================================+
| Gaussian            | quadrature_rule_gaussian_type        |
+---------------------+--------------------------------------+
| Gaussian-Lobatto    | quadrature_rule_gauss_lobatto_type   |
+---------------------+--------------------------------------+
| Newton-Cotes        | quadrature_rule_newton_cotes_type    |
+---------------------+--------------------------------------+

Functions defining locations within a cell include

+-----------------------+--------------------------------------+
| Object name           | Location description                 |
+=======================+======================================+
| quadrature_edge_type  | Cell edges                           |
+-----------------------+--------------------------------------+
| quadrature_face_type  | Cell faces                           |
+-----------------------+--------------------------------------+
| quadrature_xyz_type   | 3D array in the cell volume          |
+-----------------------+--------------------------------------+
| quadrature_xyoz_type  | 2D horizontal and 1D vertical        |
+-----------------------+--------------------------------------+
| quadrature_xoyoz_type | 1D in each direction                 |
+-----------------------+--------------------------------------+
