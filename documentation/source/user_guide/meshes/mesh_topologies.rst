.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section meshes:

******************************
Mesh topologies
******************************

=================
Cubed-Sphere mesh
=================
A Cubed-Sphere mesh an unstructured 2D topology similar to the
surface of a gridded cube (Fig.Xa,Cubed-mesh), which has the
spatial geometry of a sphere (Fig.Xb Cubed-Sphere-Mesh).

Use case: Global Circulation Model (GCM)

===========
Planar mesh
===========
The planar mesh topology is connected as a grid with nodes located on a:

Curved surface
    A gridded regional domain on a spherical surface used for Limited Area Models (LAMs). This mesh topology implicitly uses a spherical coordinate system (longitude-latitude) for the node locations. Full periodicity at the domain boundaries is not supported.

Flat surface
    A gridded regional domain on a flat surface, generally used for idealised modelling cases. This mesh topology implicity uses a cartesian coordinate system (xyz) for node locations. Full-periodicity at the domain boundaries is permitted.

=============================
Eave mesh (under development)
=============================
Eave meshes are a set of planar mesh topologies derived from a
Cubed-Sphere parent. An eave mesh overlays a given panel (side)
of the cubed-sphere, with "eaves" which extends past the associated
panel boundaries. On the cubed-sphere panels, nodes of cubed-sphere
and eave meshes are co-located. Nodes in the "eave" sections follow
the distribution pattern of nodes on the cubed-sphere panel the
given eave mesh is mapped to. An eave mesh allows activities normally
performed on a structured grid to be applied to the unstructured
Cubed-Sphere parent by addressing each panel separately and extracting
the information at the co-located nodes.

Use case: Ancillary field generation.

========
Rim mesh
========
The LBC mesh topology is provided for the use case of storing data for Lateral Boundary Conditions (LBCs) for a regional model. This mesh topology is made up cells around the boundaries (rim) of a regional domain. This rim has a depth (in cells) extending inwards from the domain boundaries. Due to nature of the use case, a given LBC mesh is only valid for use with the parent planar mesh it was generated from. The all LBC mesh cells are geographically co-located with a mapped cell on the parent planar mesh.
