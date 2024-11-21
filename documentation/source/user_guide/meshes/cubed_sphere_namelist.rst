.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section cubedsphere_mesh nml:

.. _cubedsphere_mesh_nml:

======================
``&cubedsphere_mesh``
======================
Control namelist for cubed-sphere mesh generation, required for use with the
`cubedsphere_mesh_generator`. This creates a mesh which uses the
cubded-sphere base strategy and requires ``edge_cells`` to define the mesh
connectivity. The remaining options are mesh transformations that are
subsequently applied to the node coordinates; the connectively of the mesh
elements are not altered by transformations.

* ``edge_cells``: **<integer>, …**
    Number of cells along edge of each mesh panel. The sequence of integers
    will map to entries given by :ref:`mesh_names<mesh_names>`.
* ``equatorial_latitude``: **<real>**
    Real world latitude (DegreesN) of cubed-sphere mesh equator after applying
    Schmit transform. The `top` (or `bottom`) panels of the cubed-sphere are
    reduced in size while maintaining the same connectivity. This has the
    effect of a localised increase in resolution over a panel of the
    cubed-sphere without increasing the overall number of cells in the mesh.
* ``smooth_passes``: **<integer>**
    Number of interations of smoothing function applied to mesh node
    locations.

