.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section stretch transform namelist:

.. _stretch_transform_nml:
=======================
``&stretch_transform``
=======================
Optional control namelist for stretched grid transformation. The stretched
grid transform modifies the base planar mesh resulting in two regions of
differing cell size with a stretch region separating them. The intention is to
provide (for a given axis) an inner cell size (inner region) which gradually
transitions (stretch region) to the outer cell size (outer region). Only valid
for planar meshes, this transform is enabled by the
:ref:`apply_stretch_transform<apply_stretch_transform>` logical.

* ``cell_size_inner``: **<real>,<real>**
    Cell size for domain inner region of stretch grid along ``x-axis,y-axis``,
    units in `m` or `degrees` depending on planar mesh type.
* ``cell_size_outer``: **<real>,<real>**
    Cell size for domain outer region of stretch grid along ``x-axis,y-axis``,
    units in `m` or `degrees` depending on planar mesh type.
* ``n_cells_outer``: **<integer>,<integer>**
    Depth (in cells) of domain outer region along ``x-axis,y-axis``.
* ``n_cells_stretch``:  **<integer>,<integer>**
    Depth (in cells) of domain stretch region along ``x-axis,y-axis``.
* ``stretching_on``: **'<string>'**
    Features to use as anchor points for stretch transform.

    ``cell_centres``
      Use cell centres as anchor points.

    ``cell_nodes``
      Use cell nodes as anchor points.

    ``p_points``
      Use p-points as anchor points.

* ``transform_mesh``: **'<string>'**
    Principle mesh to apply stretch transform to. Any meshes connected to this
    mesh via InterMesh maps (:ref:`mesh_maps<mesh_maps>`) will have their node
    locations updated accordingly.
