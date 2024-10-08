.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section mesh generation:

******************************
Mesh generators
******************************

The mesh generators create meshes for use by LFRic model applications.
They generate 2D-mesh quadrilateral cell topologies [#f1]_ comprised
of faces, edges and nodes .

Two separate generators are provided for the creation of either
''CubedSphere'' or ''Planar'' meshes. Development of these tools
are driven by application requirements on the LFRic Core infrastructure.

Additional mesh attributes not covered by the UGRID convention are added
via NetCDF attributes in order to support specific use cases arising from
application requirements.

============================
Cubed-Sphere mesh generator
============================

Generates a base Cubed-Sphere 2D mesh the follow strategy.
 * Six panels (1 per face of the cube), each of side ``n X n`` cells.
 * Panels 1:4 band the equator, with Panel-1 centred on the null island (0DegreesN,0DegreesE).
 * Panels 5 & 6 are centred on the North and South poles respectively.

============================
Planar mesh generator
============================
Generates a base gridded 2D mesh the follow strategy.
 * Single panel of side ``n X m`` cells.
 * Location using domain centre/extents.
 * Axes aligned with `<longitude>,<latitude>` or `<x>,<y>`.
 * Allows various configurations of periodicity at domain boundaries.

.. NOTE:: All (lon,lat) coordinates are with respect to a real world frame of reference. While mesh reference features (i.e. poles, equator, null-island) are aligned to their real world equivalents in the base meshes, Subsequent transforms (e.g. rotation, stretcthing) may move the real world location of the mesh feature.

.. rubric:: Footnotes

.. [#f1] NetCDF(.nc) file compliant with UGRID v1.0 convention.
