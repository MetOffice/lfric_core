.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section meshes:

Meshes
======

LFRic applications base their data on three-dimensional meshes. These
3D meshes are constructed by reading in two-dimensional meshes and
extruding them into a vertical coordinate.

******************
Mesh Tools
******************
Programs for the generation of 2D mesh topologies. Output aims to be formatted as UGRID (v1.0) compliant NetCDF(.nc) files. Development of these tools are driven by application requirements on the LFRic Core infrastructure.

Generated meshes are two-dimensional, using quadrilateral cells arranged for a given use case. Additional mesh attributes not covered by the UGRID convention are added via NetCDF attributes in order to support specific use cases arising from application requirements.

******************
Cubed-Sphere
******************
The Cubed-Sphere mesh topology's primary use case is for Global General Circulation Models (GCMs) as an alternative to the regular lon-lat grid. Cells in this mesh topology are connected as a gridded-cube as the name suggests, though the geographical location of the the nodes are such that the mesh geometry is spherical.


******************
Planar mesh
******************
The planar mesh topology's use case is for Regional GCMs. The planar mesh topology is connected as a regular grid with nodes located on a:

Curved surface
    A gridded regional domain on a spherical surface used for Limited Area Models (LAMs). This mesh topology implicitly uses a spherical coordinate system (longitude-latitude) for the node locations. Full periodicity at the domain boundaries is not supported.

Flat surface
    A gridded regional domain on a flat surface, generally used for idealised modelling cases. This mesh topology implicity uses a cartesian coordinate system (xyz) for node locations. Full-periodicity at the domain boundaries is permitted.

******************
Supporting mesh topologies
******************
LBC mesh (Planar mesh support)
******************
The LBC mesh topology is provided for the use case of storing data for Lateral Boundary Conditions (LBCs) for a regional model. This mesh topology is made up cells around the boundaries (rim) of a regional domain. This rim has a depth (in cells) extending inwards from the domain boundaries. Due to nature of the use case, a given LBC mesh is only valid for use with the parent planar mesh it was generated from. The all LBC mesh cells are geographically co-located with a mapped cell on the parent planar mesh.

******************
Eave mesh(es) (Cubed-Sphere mesh support, under development)
******************
The eave meshes are topologically the same as a planar mesh on a curved surface, although, like LBC meshes, they are generated from a parent as supporting meshes. They are under development to support processes which are complicated by the geometric poles (corners) of the cubed-sphere mesh topology. Nodes forming cells on an eave mesh will be geographically co-located with the nodes forming a panel (side) of the cubed-sphere parent. In addition, the "eaves" of an eave mesh are a rim of cells which extend past the extents of the corresponding cubed-sphere panel. Similar to LBC meshes the rim may be a specified number of cells deep.  
