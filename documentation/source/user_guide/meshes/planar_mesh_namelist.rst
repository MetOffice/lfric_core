.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section planar_mesh_nml:

.. _planar_mesh_nml:

=============================
``&planar_mesh``
=============================

Control namelist for planar mesh generation, required for use with the
`planar_mesh_generator`. This creates a mesh which uses the planar mesh
base strategy and requires ``edge_cells`` along `both` axes aswell as the
domain boundary periodicity in order to define the mesh connectivity. The
remaining options are mesh transformations or logical triggers. Mesh
tranformations are applied to the node coordinates after the base mesh is
generated; the connectively of the mesh elements are not altered by
transformations.

.. _apply_stretch_transform:

* ``apply_stretch_transform``: **<logical>**
    Apply the stretch transform to base planar mesh node coordinates.
* ``create_lbc_mesh``: **<logical>**
    Generate a rim mesh which is derived one of the priciple meshes
    (:ref:`mesh_names<mesh_names>`).
* ``domain_centre``: **<real>,<real>**
    Location of domain centre for all principle meshes. Coordinates aligned
    with the :ref:`coord_sys<coord_sys>` with units of `[m]` or `[degrees]` as
    appropriate, `i.e.` `x,y` or `lon,lat` coordinates.
* ``domain_size``: **<real>,<real>**
    Domain size for all principle meshes.  Domain sizes aligned with the
    :ref:`coord_sys<coord_sys>` with units of `[m]` or `[degrees]` as
    appropriate, `i.e.` `x,y` or `lon,lat` domain extents.
* ``edge_cells_x``: **<integer>, ...**
    Number of cells along x-axis of domain. Order of integers map to
    :ref:`mesh_names<mesh_names>`.
* ``edge_cells_y``: **<integer>, ...**
    As above, but for y-axis.
* ``lbc_parent_mesh``: **'<string>'**
    Name of the principle mesh to create a rim mesh from, this name should
    exist in the :ref:`mesh_names<mesh_names>` variable.
* ``lbc_rim_depth``: **<integer>**
    Depth (in cells) of rim mesh. This is the number of cells radially across
    the rim mesh from the domain centre to a domain boundary.

.. _periodic_x:

* ``periodic_x``: **<logical>**
    Periodicity across pair of domain boundaries in x-axis.

.. _periodic_y:

* ``periodic_y``: **<logical>**
    Periodicity across pair of domain boundaries in y-axis.
