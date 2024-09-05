.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. role:: raw-html(raw)
    :format: html

.. _unit_test_canned_data_routines:

Unit test canned-data support routines
======================================

There are a number of modules available that provide different sets of canned data. click on a module name, below to see the data it provides and how to call the routines that provide that data.

.. dropdown:: ``get_unit_test_sizes_mod.f90``

   This module holds functions that return a set of sizes that often need to be passed into kernels. The functions are named as follows:

   * ``get_{function_space_name}_m3x3_q3x3x3_size``

   So a call to ``get_w0_m3x3_q3x3x3_size(ndf, undf, ncells, dim_space, dim_space_diff, nqp_h, nqp_v, nlayers=nlayers)`` will return the following sizes for a field on a w0 function space using a 3x3 mesh and a 3x3x3 quadrature: :raw-html:`<br />` ``ndf`` - the number of dofs in any cell :raw-html:`<br />` ``undf`` - the total number of unique dofs in the 3x3 domain :raw-html:`<br />` ``ncells`` - the number of cells in the 3x3 domain (9) :raw-html:`<br />` ``dim_space`` - the number of dimensions in this function space :raw-html:`<br />` ``dim_space_diff`` - the number of dimensions in this function space when differentiated :raw-html:`<br />` ``nqp_h`` - the number of  quadrature points in the horizontal for a 3x3x3 quadrature (9) :raw-html:`<br />` ``nqp_v`` - the number of quadrature points in the vertical for a 3x3x3 quadrature (3) :raw-html:`<br />` The following can be optionally provided: :raw-html:`<br />` ``nlayers``- the number of cells in the mesh in the vertical direction (if no value is provided the routines will assume a value of 3 for nlayers)

.. dropdown:: ``get_unit_test_dofmap_mod.f90``

   This module holds functions that return dofmaps and stencil dofmaps. The functions are named as follows:

   * ``get_{function_space_name}_m3x3_dofmap`` :raw-html:`<br />` So a call to ``get_w2_m3x3_dofmap(map_w2)`` will return the dofmap (``map_w2(:,:)``) for a field on a w2 function space on a 3x3 mesh. :raw-html:`<br />`
   * ``get_m3x3_stencil_dofmap_{stencil_shape}`` :raw-html:`<br />` where ``{stencil_shape}`` is one of ``point``, ``cross``, ``1dx`` or ``1dy``. :raw-html:`<br />` The routine takes a dofmap as input and returns the appropriate stencil dofmap, :raw-html:`<br />` so ``get_m3x3_stencil_dofmap_cross(stencil_map, w2_dofmap)`` returns ``stencil_map(:,:,:)`` for a cross-shaped stencil dofmap which will be based on the input w2 dofmap.

.. dropdown:: ``get_unit_test_basis_mod.f90``

   This module holds functions that return basis functions and differential basis functions for lowest order function spaces with a 3x3x3 quadrature. These routines are named as follows:

   * ``get_{function_space_name}_q3x3x3_basis`` :raw-html:`<br />` A call to ``get_w1_q3x3x3_basis(basis_w1)`` returns the basis functions for a w1 field in ``basis_w1(:,:,:,:)``
   * ``get_{function_space_name}_q3x3x3_diff_basis`` :raw-html:`<br />` A call to ``get_w1_q3x3x3_diff_basis(diff_basis_w1)`` returns the differential basis functions (``diff_basis_w1(:,:,:,:)``) for a w1 field.

.. dropdown:: ``get_unit_test_quadrature_mod.f90``

   This module holds routines that return quadrature information for a simple 3x3x3 Gaussian quadrature.

   * ``get_gaussian_q3x3x3_quadrature_points_xy`` and
   * ``get_gaussian_q3x3x3_quadrature_points_z`` :raw-html:`<br />` A call to ``get_gaussian_q3x3x3_quadrature_points_xy(points_xy)`` will return the coordinates of the quadrature points for a 3x3x3 Gaussian quadrature in the horizontal (``points_xy(:,:)``) and ``get_gaussian_q3x3x3_quadrature_points_z(points_z)`` will return the vertical points (``points_z(:)``).

   * ``get_gaussian_q3x3x3_quadrature_weights_xy`` and
   * ``get_gaussian_q3x3x3_quadrature_weights_z`` :raw-html:`<br />` A call to ``get_gaussian_q3x3x3_quadrature_weights_xy(weights_xy)`` will return the weights to apply to the quadrature points for a 3x3x3 Gaussian quadrature in the horizontal (``weights_xy(:)``) and ``get_gaussian_q3x3x3_quadrature_weights_z(weights_z)`` will return the vertical weights (``weights_z(:)``).

.. dropdown:: ``get_unit_test_entity_map_mod.f90``

   This module holds routines that return the list of cell entities on which each dof of a function space lives. It does so for element order 0 or 1.

   * ``get_{function_space_name}_{order}_entity_map(entity_map)`` :raw-html:`<br />` Where ``order`` id either ``o0`` or ``o1``. ``get_w2_o1_entity_map(entity_map)`` gets the list for the `W2` function space for element order 1. See the routine comments for the entity numbering which is derived from the reference cell.

.. dropdown:: ``get_unit_test_qfaces_mod.f90``

   This module holds routines that return information about computing quadrature on the faces of a cube cell in a 3x3 mesh. Some of the functions can return data for horizontal faces only, vertical faces only, or all faces - replace ``{hv}`` in the names with a `_h_`, `_v_` or `_hv_`. The functions for returning dofmaps are named as follows:

   * ``get_number_quadrature_points_per_face(nqp)`` :raw-html:`<br />` returns the number of quadrature points per face.
   * ``get_quadrature_faces_{hv}_weights(wqp)`` :raw-html:`<br />` returns the weights of each quadrature point for each face.
   * ``get_{function_space_name}_qfaces_{hv}_basis(basis_w2_face)`` :raw-html:`<br />` returns the basis functions for each face. The module is not fully complete. Refer to the list of public functions to find out what data has been canned. If other data is needed, it can be computed and added.

.. dropdown:: ``get_unit_test_nodal_coords_mod.f90``

   This module holds routines which return co-ordinate information for the function space DoFs. These routines are named as follows:

   * ``get_{function_space_name}_nodal_coords`` :raw-html:`<br />` A call to ``get_w0_nodal_coords`` will return the !DoF co-ordinates for the W0 functionspace in ``coords_w0(:,:)``.

.. dropdown:: ``get_unit_test_planar_mesh_mod.f90``

   This module holds routines to set input arrays with some basic data about the regular 3x3x3 biperiodic mesh used in most of the unit tests. As the mesh is regular, the array inputs are fixed and need to be allocated by the test routine.

   * ``get_m3x3_adjacent_face`` :raw-html:`<br />` allocates and initialises a 4-by-9 integer array with information about the 4 faces of the 9 (nominally 3x3) horizontal neighbour cells. Namely, if a cell `cell_number` numbers a face `face_number`, then the number of the face according to the cell that shares this face is `adjacent_face(face_number,cell_number)`
   * ``get_out_face_normal`` :raw-html:`<br />` allocates and initialises a 3-by-6 real array with vectors normal to the face, facing outward from the cell.
   * ``get_normals_to_faces`` :raw-html:`<br />` allocates and initialises a 6-by-3 real array with vectors normal to the face in the mesh coordinate system.

.. dropdown:: ``get_unit_3x3x3_chi_mod.f90``

   This module holds routines to create coordinate fields (chi1, chi2, chi3) for the W0 and WChi coordinate spaces. These routines utilise dofmaps generated by the ``get_unit_test_dofmap_mod.f90`` module.

   * ``get_w0_3x3x3_field`` :raw-html:`<br />` returns a set of three coordinate fields (chi1, chi2, chi3) for the lowest order, W0 space on a planar domain.
   * ``get_wchi_3x3x3_field`` :raw-html:`<br />` returns a set of three coordinate fields (chi1, chi2, chi3) for the W_Chi (W3, k=1) space on a planar domain.
   * ``get_wchi_3x3x3_latlon_field`` :raw-html:`<br />` returns a set of three coordinate fields (chi1, chi2, chi3) for the W_Chi (W3, k=1) space on a portion of a spherical domain with lat-lon variation.

.. dropdown:: ``get_unit_m3x3_cma_data_mod.f90``

   This module holds routines to compute the scalars and maps for Column Matrix Assembly (CMA) operators. These routines utilise the dofmaps for the spaces they are on.

   * ``get_cma_size`` :raw-html:`<br />` computes the CMA scalars from inputs of `nlayers` `ndf` and the number of dofs that are face or interior valued for the ''to'' and ''from'' spaces.
   * ``get_cma_prod_size`` :raw-html:`<br />` Computes the CMA scalars from the input CMA scalars of the product operands.
   * ``get_wtheta_m3x3_cma_data`` :raw-html:`<br />` Computes the CMA scalars, allocates and populates the maps for a Wtheta-Wtheta CMA operator.
   * ``get_wtheta_w3_m3x3_cma_data`` :raw-html:`<br />` Computes the CMA scalars, allocates and populates the maps for a Wtheta-W3 CMA operator.
   * ``get_w3_w3_m3x3_cma_data`` :raw-html:`<br />` Computes the CMA scalars, allocates and populates the maps for a W3-W3 CMA operator.
   * ``get_w2v_w2v_m3x3_cma_data`` :raw-html:`<br />` Computes the CMA scalars, allocates and populates the maps for a W2v-W2v CMA operator.
