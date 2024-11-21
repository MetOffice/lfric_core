.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section mesh generators:

******************************
Mesh generators
******************************

The mesh generators output mesh input files intended for use by LFRic
applications. They generate 2D-mesh topologies comprised of faces
(quadrilateral cells), edges and nodes elements.

Two separate generators are provided for the creation of either
''Cubed-Sphere'' or ''Planar'' meshes. Development of these tools is driven by
application requirements on LFRic core infrastructure.

Output files are in NetCDF format following UGRID conventions [#f1]_. While
UGRID may describe the mesh topologies, this may not be sufficient for some
use cases in LFRic applications. Where required, additional mesh information
may be added using variable attributes on the mesh.


============================
Cubed-Sphere mesh generator
============================

Generates one (or more) 2D Cubed-Sphere meshes that each follow the base
strategy.

* Six panels (1 per face of the cube), each of side ``n X n`` cells.
* Panels 1:4 band the equator, with Panel-1 centred on the null island
  [#f2]_.
* Panels 5 & 6 are centred on the North (90N,0E) and South (90S,0E) poles
  respectively.

--`usage:`

    | >> **cubedsphere_mesh_generator <configuration.nml>**

============================
Planar mesh generator [#f3]_
============================
Generates one (or more) 2D gridded meshes that each follow the base strategy.
 * Single panel of side ``n X m`` cells.
 * Located using combination of the specified domain centre and extents.
 * Axes aligned with `<longitude>,<latitude>` or `<x>,<y>`.
 * Non-periodic at domain boundaries.

--`usage:`

    | >> **planar_mesh_generator <configuration.nml>**


.. rubric:: Footnotes

.. [#f1] NetCDF(.nc) file compliant with UGRID v1.0 convention.
.. [#f2] `Null Island`: where the prime medridian and the equator intersect,
	 ``i.e.`` 0N, 0E.
.. [#f3] For historical reasons, the planar mesh generator is not restricted
	 to a planar surface domain. References to planar meshes with respect
	 the following sections should be taken as a mesh produced by the
	 `planar_mesh_generator` program.
