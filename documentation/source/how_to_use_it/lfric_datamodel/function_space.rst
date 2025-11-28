.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _function space:

LFRic function spaces
=====================

LFRic function spaces map data points and finite element basis
functions onto a domain represented by a 3D mesh. A high level
introduction into :ref:`LFRic function spaces<function space
intro>` describes the mapping in more detail. This section is focused
on how to create and use function spaces based on function space types
supported by the infrastructure. If new fundamental function space
types are required, expert advice will be needed.

Initialising an LFRic function space is a prerequisite to initialising
a field: field initialisation allocates the correct amount of data to
represent degrees of freedom (dofs) for the function space type on all
the cells of the mesh that was used to create the LFRic function
space.

A new LFRic function space is initialised as follows:

.. code-block:: fortran

   type(function_space_type), pointer :: vector_space

   ! Get a reference to a lowest order W2 function space
   vector_space => function_space_collection%get_fs(mesh_id,         &
                                                    element_order_h, &
                                                    element_order_v, &
                                                    fs_type)

The ``mesh_id`` identifies the 3D mesh that the LFRic function space
must cover. The ``element_order`` options dictate the polynomial
order of the function space in the horizontal and vertical directions
respectively. They are typically ``0`` for lowest order or ``1`` for
next-to-lowest order. The GungHo dynamical core is tested
with both of these settings whereas the Momentum\ :sup:`®` atmosphere
model runs with lowest order function spaces.

The ``element_order_h`` and ``element_order_v`` options are ordered
horizontal order first, and vertical order second. This is the
convention adopted throughout the code, so functions will always
take the horizontal version of an entity, and then the vertical
version, and finally the next entity / argument.

The ``fs_type`` refers to one of a number of available function space
types. A `high-level description
<https://psyclone.readthedocs.io/en/stable/dynamo0p3.html#supported-function-spaces>`_
of supported function spaces can be found in the PSyclone
documentation. Refer to GungHo documentation for a comprehensive
description of each function space.

The function space constructor can take two other optional arguments:

#. Specifying an integer ``ndata`` allows the creation of
   :ref:`multidata fields<multidata field>` with ``ndata`` dofs
   per dof location. The default value is `1`.
#. Spacifying ``ndata_first = .true.`` allows creation of fields with
   data ordered layer-by-layer instead of column-by-column.

Quadrature rules
----------------

This section summarises how quadrature rules are applied to model data
without attempting to describe the details of the mathematical
processes.

As described in the overview of :ref:`LFRic function spaces<function
space intro>`, a function space includes a set of `basis
functions` that are used, along with dofs, to describe a field
throughout a cell. To manage the cost of running kernels, the
function spaces are computed at only a selection of locations in the
cell. Quadrature rules are created to ensure the computation can be
done sufficiently accurately and efficiently.

Simply put, the infrastructure supports several quadrature `rules`,
which are types that define `how` function space basis functions are
integrated, and several quadrature types which define a list of
locations in a cell `where` the rules will be computed. The choice of
rule and location governs the accuracy of the method but also
influences the cost (for example, more points require more compute
time).

A finite element kernel includes metadata that describes which method
for assigning quadrature locations is appropriate. Prior to calling
the kernel, the algorithm initialises the quadrature type with
arguments that include a quadrature rule:

.. code-block:: fortran

   type( quadrature_rule_gaussian_type ) :: quadrature_rule
   type( quadrature_xyoz_type )          :: qr

   qr = quadrature_xyoz_type(nqp, quadrature_rule)

Here, the `nqp` argument is an integer that is used by the
``quadrature_xyoz_type`` initialisation procedure to define the number
and location of the points in a reference cell at which to compute the
quadrature.

PSyclone will generate the code that calls the computation method of
the chosen quadrature type for each function space that is specified
by the kernel. The computation method uses the quadrature rule to
compute the basis functions at the required set of points defined by
the type, returning an array of numbers representing the evaluated
basis functions.

As well as defining the relevant quadrature type, a kernel will
request quadrature for one or more function spaces and for either or
both basis functions or differential basis functions.

For example, the following kernel requires both basis and differential
basis functions to be computed for the argument associated with the
``ANY_SPACE_9`` field, differential basis functions for the ``W0``
field and basis functions for the ``W1`` and ``W3`` fields. The
quadrature type is ``quadrature_XYoZ``. This type contains points and
weights stored in 2D (xy represents the horizontal direction) and 1D
(z represents the vertical direction).

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

The subroutine argument list includes arguments alligned to each of
these components. Namely, the five arguments suffixed ``_basis`` or
``_diff_basis`` and the last four arguments that describe,
respectively, the number of points in the horizontal and vertical
direction and the weights in the horizontal and vertical direction.

.. code-block:: fortran

   subroutine compute_total_pv_code(                                         &
                     nlayers,                                                &
                     pv,                                                     &
                     xi,                                                     &
                     theta,                                                  &
                     rho,                                                    &
                     chi1, chi2, chi3, panel_id,                             &
                     omega, f_lat,                                           &
                     ndf_w3, undf_w3, map_w3,    w3_basis,                   &
                     ndf_w1, undf_w1, map_w1,    w1_basis,                   &
                     ndf_w0, undf_w0, map_w0,    w0_diff_basis,              &
                     ndf_chi, undf_chi, map_chi, chi_basis, chi_diff_basis,  &
                     ndf_pid, undf_pid, map_pid,                             &
                     nqp_h, nqp_v, wqp_h, wqp_v )


The following quadrature rules are supported by the infrastructure:

+---------------------+--------------------------------------+
| Rule name           | Object                               |
+=====================+======================================+
| Gaussian            | quadrature_rule_gaussian_type        |
+---------------------+--------------------------------------+
| Gaussian-Lobatto    | quadrature_rule_gauss_lobatto_type   |
+---------------------+--------------------------------------+
| Newton-Cotes        | quadrature_rule_newton_cotes_type    |
+---------------------+--------------------------------------+

Functions defining locations within a cell include the following.

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
