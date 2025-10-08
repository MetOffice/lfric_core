.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _section ugrid compliance:

UGRID compliance
========================

Mesh topologies are produced by the mesh generators following the UGRID
convention [#f1]_. The convention specifies mandatory and optional data.
The following sections document variables/data which are output for a given
mesh topology in a :ref:`mesh output file<section mesh file content>`:

* UGRID mandatory data

* UGRID optional data

* Non-UGRID data for application support

.. _ugrid_madatory:

Mandatory UGRID data
---------------------------

* | ``node_coordinates``
  | Node coordinates, [lon,lat] for spherical domains, [x,y] for
  | flat cartesian domains.

* | ``face_node_connectivity``
  | Integer identifiers of nodes that construct a given face.

.. _ugrid_optional:

Optional UGRID data
---------------------------
* | ``face_coordinates``
  | These are derived from the node coordinates as follows:

  * For *Cubed-Sphere meshes*, the faces are not necessarily of similar size
    along the [lon,lat] axes.

    1. For the nodes connected to a given face, the node coordinates are
       converted from (lon, lat, radius=1) to cartesian [x,y,z] coordinates.

    2. The face centre is calculated as the mean of vector sum of the node
       coordinates, allowing for a radius ratio.

    3. The face centre coordinates are converted back to (lon, lat).

  * For *Planar meshes*, the domains are constructed to be aligned with the
    coordinate axes, whether [lon,lat] or [x,y]. The face centre coordinates
    are obtained from an offset to the cell's NW node location., *i.e.* (x,y) :sub:`NW node` + 0.5(:math:`\Delta` x, :math:`-\Delta` y) :sub:`cell`.

* | ``edge_node_connectivity``
  | Integer identifiers of nodes that construct a given edge.

* | ``face_face_connectivity``
  | Integer identifiers of faces adjacent to a given face.

* | ``face_edge_connectivity``
  | Integer identifiers of edges that construct a given face.

Non-UGRID data
------------------------
* | **Geometry**
  | Indicates the geometrical shape the mesh domain forms, *i.e.* spherical.
  | *Related attribute:* ``geometry``

* | **Coordinate system**
  | Indicates the coordianate system used to interpret coordinate values.
  | *Related attributes:* ``coord_sys``

* | **Periodicity**
  | Indicates the periodic state of the domain boundaries on a given axis.
  | *Related attribute:* ``topology``, ``periodic_x``, ``periodic_y``

* | **InterMesh maps**
  | Indicates the `targets` of InterMesh maps that are present in the file,
    where the given mesh is the `source`.
  | *Related attributes:* ``n_mesh_maps``, ``maps_to``

* | **Domain extents**
  | Domain extents along the coordinate system axes.
  | *Related attributes:* ``domain_extents``

* | **Domain feature locations**
  | Coordinate locations of significant features on the mesh, *e.g.* North pole.
  | *Related attributes:* ``north_pole``, ``null_island``,
    ``equatorial_latitude``

.. rubric:: Footnotes

.. [#f1] NetCDF(.nc) file compliant with UGRID v1.0 convention.
