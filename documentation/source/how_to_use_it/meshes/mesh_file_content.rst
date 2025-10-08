.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _section mesh file content:

Mesh file contents
========================

The mesh generators produce one or more mesh files containing the
following entities:

* One or more principal meshes, specified by namelist
  (:ref:`mesh names<mesh_names>`).

* (Optional) Derived meshes from a specified principal mesh, *e.g.* Rim mesh,
  Eave mesh.

* (Optional) InterMesh maps between two specified principal meshes.

* (Local mesh) Information on partition (halos etc.) and non-partitioned global
  mesh.

.. _principal_meshes:

Principal meshes
------------------------
These meshes are explicitly defined from the
:ref:`mesh configuration namelists<section configuration namelists>`

.. _derived_meshes:

Derived meshes
------------------------
These mesh types may require additional input, though crucially, are
derived from one of the specified principal meshes. They cannot be
separately generated from the parent principal mesh. Derived meshes
are implicitly linked to their parent mesh via a one-way
:ref:`InterMesh map<intermesh_maps>` (Derived -> Principal
InterMesh map).

.. _intermesh_maps:

InterMesh maps
------------------------
InterMesh maps provide a means of identifying co-located (or partially
co-located) cells between two mesh topologies. InterMesh maps are
attached to the ``source`` mesh, and record which cells on the
``target`` mesh are co-located to a given cell on the ``source`` mesh.
Explicit InterMesh maps are requested via
:ref:`mesh configuration namelists<section configuration namelists>`
with the following restrictions:

* Requested InterMesh maps are only available between two specified
  principal meshes during a mesh generation event.

* The resolutions of the principal mesh pairing are such that one mesh
  is a sub-division of the other.

For explicitly requested InterMesh maps, the generator will output an
InterMesh map for each direction between the meshes, *i.e.*
``source -> target``, ``target -> source``. The InterMesh map will
always be attached to the ``source`` *`relative`* to the mapping.

.. _local_meshes:

Local meshes
-----------------------
When partitioning global meshes for a specific domain decomposition,
additional information is provided in order to populate the
applications ``local_mesh_type`` objects. Partitioning results in
multiple output files, one for each partition. Each file contains
topologies for the principal meshes, although a local mesh only
contains a subset of the cells from the global mesh. The subset
of cells encompass the partition itself plus a halo region.

Additional information relating to the partition halos and the
non-partitioned global mesh is output to file when the generators
are requested to partition the principal meshes.
