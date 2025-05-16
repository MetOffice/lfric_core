.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _distributed memory:

Distributed Memory Documentation
################################

Introduction
============

All cores within a node of the supercomputer share main memory. Shared
memory parallelism can be used over these cores. Between the nodes there
is no such shared memory. Information has to be specifically passed
between nodes using messages. This is known as distributed memory
parallelism.

Distributed memory parallelism is achieved by partitioning the full
(global) domain into smaller sections (partitioned or local domains).
Information is passed between the tasks through the use of halo regions.
Halos are an area that surrounds a partitioned domain. Data in these
areas are provided from other partitions that own that section of the
domain. The filling of halos with the current data involves the exchange
of information from adjacent partitions, so they are often known as
"halo exchanges" or "halo swaps".

None of the above issues are unique to LFRic. Running large model codes
over a number of parallel tasks has been done for many years. The LFRic
infrastructure has a data model with significant differences to other
codes like the UM which will affect the way distributed memory
parallelism is implemented:

#. LFRic is designed to support horizontally unstructured meshes, which
   means it must support non-uniform shaped partitions with complex
   shaped halo regions.

#. Because the domain is horizontally unstructured, moving from one data
   point to a horizontally adjacent data point requires using a lookup
   table, which is slower than the direct addressing used in the UM. In
   order to recover some code performance, fields are laid out in memory
   so that vertically adjacent data points are next to each other, which
   means that, for optimisation reasons, the looping structure is such
   that the innermost loop (the fastest-changing loop) is over these
   vertical columns of data. This is often called k-innermost or k-first
   looping, as k is used to denote the vertical direction.

#. LFRic supports a finite-element formulation which means that data is
   held as "degrees of freedom" (dofs) on all the different entities
   within the mesh: cell volumes, faces, edges and vertices.
   Partitioning the domain along cell boundaries will always lead to
   partitions which need to share data values. Sharing data makes
   assigning ownership of that data, and therefore, halo swapping more
   complex.

There are a couple of key requirements that the distributed memory
implementation needs to support. It is important to consider the
ordering of cells in the cell list separately to the order of the data
in memory:

#. Looping over specific ranges of cells. In the design used for LFRic,
   the looping over the horizontal range of cells is performed in the
   PSy layer. There are good reasons to want to loop over specific
   ranges of horizontal cells (e.g. for overlapping communications with
   compute), so it is important to order the cells to support this.

#. Working on specific parts of the data held in memory. It is important
   to allow operations on specific sections of the field data. The
   layout of the data should be in an order that supports this. For
   example, it is advantageous to be able to swap halos to a specific
   depth or to be able to perform redundant calculations on a specific
   depth of halo data to reduce the cost of communications. To achieve
   these requirements the data within a field must be ordered such that
   the data owned by the partition comes first, followed by increasingly
   deep halo data.

For full details of the data model supported by the LFRic infrastructure
see the LFRic design documentation [LFRic_Design]_.

Parallel implementation
=======================

Partitioning
------------

Field data in a finite-element representation are held as "degrees of
freedom" (dofs). These dofs can be associated with any mesh entity (cell
volume, face, edge or vertex). If the dof is located on an entity that
is shared by more than one cell (i.e. a face, edge or vertex dof), then
with any cell-based partitioning, it could be shared by more than one
partition. The solution is to represent those dofs on both partitions
(with the necessary surrounding data in "halos" of cells around each
partition). The calculations for the duplicated dofs is then performed
on each partition. :numref:`partitioning` shows how this
technique works for the case where the dofs located on
vertices are calculated from dofs located in the cell volume.

.. _partitioning:
.. figure:: images/partitioning.svg

  Partitioning method shown here for a 1-d line of cells. a) The domain is
  partitioned along a cell boundary. If the calculation of dofs on vertices
  (dots) in some way depends on the dofs in the volume of the cells surrounding
  them (crosses), then the calculation of the vertex dofs at the edge of a
  partition (i.e. Vert34 and Vert35) will be difficult, as they depend on the
  volume dof from the other partition. b) The solution is to introduce halos at
  the boundaries and halo exchange the dofs in the volume of the cells, then
  the calculation of the vertex dofs at the edge of the partition will now be
  correct. There will be redundant calculation - both partitions calculate the
  values for Vert34 and Vert35. The dofs at the outer edges of the of the halo
  will be incorrect, but if required, these can be corrected with a halo
  exchange of the vertex dofs.

If the value at a particular dof depends on other dofs from a stencil
bigger than one cell, then deeper halos will be required, but the
redundant calculation still only needs to run over the first depth of
halo to get all the dofs on owned cells correct.

To ensure the computation of data in owned cells is correct, the data in
halos need to be correct. Correct values are obtained from neighbouring
partitions by executing halo swaps using the MPI library. For
performance and ease of use reasons, the MPI library is not called
directly, instead a wrapper is used. The YAXT (Yet Another eXchange
Tool) library written by DKRZ is used for this wrapper.

YAXT halo exchanges require that a routing table is pre-generated at
initialisation time, so each rank needs to supply the following to the
YAXT call:

-  A list of all the global dof ids of all the locally owned dofs

-  A list of all the global dof ids of all the dofs in the increasing
   levels of halo depth

YAXT then generates a routing table for each depth of halo exchange.
These are stored until needed to perform a halo exchange.

Mesh Objects
------------

The global mesh object holds all the information that describes the
topology and geometry of the two-dimensional (surface) mesh that the
model will use. The data to fill the object is read from a file (the
file uses a convention for holding unstructured data in netCDF format
called UGRID) [UGRID_Convention]_. This is the mesh information for the
whole un-partitioned domain. It, therefore, provides a global indexing
of cells.

Mesh entities such as vertices or edges are located at the boundary
between two cells. When the mesh is partitioned, the information held on
these entities should have a definitive owner for future halo exchanges,
so it is useful to assign which cell "owns" each mesh entity. The owning
cell will hold the data for any dofs on the entities and the other cell
that shares the entity can hold its own data value or get the reference
value from the owning cell. The global mesh object holds this entity
"ownership" information. The ownership is decided by assigning each
entity to the neighbouring cell with the highest global index.

The mesh object (separate from the global mesh object) holds a full,
three-dimensional representation of the mesh, for only the cells owned
the local partition (and the halo cells that surround it). It is
generated by applying the partition object and the method for generating
the vertical mesh above each cell (known as extrusion) to the global
mesh object.

Partition Object
----------------

The purpose of the partition object is to provide each partition with
the list of cells on that partition. This list of cells is passed with
the global mesh object (and vertical mesh information) to generate the
local mesh. The partition object works by by applying a partitioning
algorithm to the global mesh object. The partitioning algorithm is
supplied to the partitioner through the use of a function passed in via
a function pointer.

For efficiency reasons, the algorithm for partitioning a cubed-sphere
mesh is a specific algorithm that will only work for a cubed-sphere mesh
(each panel of the cubed-sphere is split into rectangular partitions).
This specific partitioner could be replaced with a more general
partitioner for use on fully unstructured meshes.

Computation within the science parts of the code is performed by
iterating over cells. The order of the computation is determined by
order the cells appear in the partition object, so it should be defined
to support a chosen iteration strategy.

By default, the cells in the partition are ordered using the mesh-ids
from the global mesh. This makes cells that are geographically
nearby close to each other in the list. This can aid with
performance, through better cache use.

Alternatively, the cells can be ordered so that they support the
optimisation known as "overlapping communications with compute"

.. _overlap comms compute:

Overlapping Communications with Compute
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A potential iteration strategy involves iterating over all cells that
aren't affected by halo swaps at the same time as performing the
halo swap, then completing the calculations on the (hopefully few)
remaining cells. This is known as overlapping communications with
compute. In order to implement this, it must be possible to iterate over
the cells that aren't affected by halo exchanges for all different
sizes of stencil operations, so all those sets of cells must be held
contiguously in their array. The cells, therefore, are ordered as
follows:

-  Inner n+1 cells

-  Inner n cells ...

-  Inner 1 cells

-  Edge cells

-  Halo 1 cells ...

-  Halo n-1 cells

-  Halo n cells

-  Halo n+1 cells

In the above, a halo cell refers to a cell that is outside the local
domain, but is still required to perform calculations on the local
domain. Larger stencil operations require a greater depth of halo. An
inner halo cell is a cell on which calculations can be performed without
needing information from halo cells. Obviously, larger stencil
operations reach further, so more inner halos will be required for
greater stencil sizes. The layout of halo and inner halo cells can
be seen in the :numref:`reversehalos`.

.. _reversehalos:
.. figure:: images/ReverseHalos.svg

   Description of the way halos and inner halos are organised.

The ordering of cells is achieved by building a linked list. The cells
around the perimeter of the domain are defined and placed in the linked
list. These are the "Edge" cells. A depth=1 stencil is applied to all
these edge cells and any new cells found that are on the outside of the
domain are the "Halo(1)" cells. These are added to the list after the
"Edge" cells. This is repeated, passing the stencil over all the
"Halo(1)" cells to find "Halo(2)" and so on. The stencil is then
reapplied to the "Edge" cells and all the new cells that are inside the
domain are the "Inner(1)" cells. These are added to the linked list
before the "Edge" cells. Again, this is repeated to get the multiple
depths of inner halos. Finally, an array is allocated to hold all the
cells in the list, and they are copied into the array, before the linked
list is deleted.

Function Space Object
---------------------

The function space object defines how dofs are placed on the mesh. It
has the following features to support running in distributed memory
parallel:

-  The dofs held in the field data array have to be in a specific order.
   Any set of dofs to be swapped in a halo exchange needs to be held
   contiguously in memory to allow for efficient message passing. The
   order that dofs are held within a field is defined by the "dofmap".
   The dofmap is a lookup table that describes (for each cell) the
   location within the large one-dimensional field data array of the
   bottom dof of every vertical column of dofs. The data within a column
   of dofs is held contiguously, so other dofs in the column can be
   found by iterating upwards from the bottom dof. The LFRic
   infrastructure also supports the concept of multiple data values at
   each dof location. These are called multi-data fields and the
   different values for a multi-data field sit above the entries in the
   dofmap in the same way the different vertical levels of data do. The
   dofmap needs to be constructed to reflect the order required to make
   halo data contiguous.

-  A globally consistent index for each dof needs to be provided, so
   that halo dofs can be correctly associated with their owned dof.

.. _dofmap generation:

Generating the dofmap to support halo exchanges
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to be able to exchange halo data, the dofs have to be ordered
so that the columns of dofs that need to be exchanged are held
contiguously in the field arrays. To achieve this, dofs are split into
three categories: "Owned dofs" (these are dofs that are on owned cells
and are actually owned by the local rank), "Annexed dofs" (these are
dofs that are on owned cells but are actually owned by another rank that
shares them) and finally "Halo dofs" (these are dofs on halo cells, that
are owned by other ranks). :numref:`dofownership` describes
the three types of dof.

Ownership of dofs that are shared between cells follows the ownership of
the entity they lie on (see the global mesh object, above) and ownership
is therefore granted to the cell (and therefore that cell's owning
partition) with the highest global cell id.

.. _dofownership:
.. figure:: images/DofOwnership.svg

   Examples of the three different types of dof. Four
   partitions, P0-P3, of a distributed, lowest order W0 field are
   shown. Only the halo cells for partition P1 are shown - other halos
   are omitted from the diagram for clarity. The cells within each
   partition are shown with their global cell id. Each dof in P1 is
   shown as a dot.

Because columns of dofs are held contiguously, the ordering of field
data can be done by considering the bottom level of dofs only: the dofs
above will simply follow the lowest level dof.

The ordering of the field array is done in two stages. Each cell entity
(for every cell in the bottom level) is processed and all the columns of
dofs above that entity are categorised as either owned, annexed or halo.
The columns of dofs are stored by entity type, ordered using the global
id of that entity which means that shared columns of dofs only appear
once in the arrays (shared columns are on the same entity with the same
global id). In the second stage, the cell by cell dofmap is generated by
combining the entity arrays so that shared columns of dofs point to the
same space in the field array.

The code that constructs the dofmap as described above is held in the
``dofmap_setup`` function which can be found in the source file:
``function_space_constructor_helper_functions_mod.F90``.

The entity arrays are generated in the following way:

#. Initialise two counters - ``id_owned`` for the index into the field
   array for the base of columns of owned dofs (counts upwards from 1)
   and ``id_halo`` for all the others (counts downwards from -1)

#. Loop over all the cells in the bottom level of the local domain
   including the first level of halo.

#. Loop over all the columns of dofs in the volume of the cells and add
   the index (``id_owned`` or ``id_halo`` depending on whether the cell
   is owned or not) to a store of volume dofs. Increase (or decrease)
   the appropriate index to leave room in the field data for a whole
   column of dof data.

#. Loop over all the columns of dofs on the faces of the cells (start
   with the vertically oriented faces, then do the top and bottom) and
   add the index (depending on whether the face is owned by the local
   rank or not) to a store of face dofs.

#. Loop over all the columns of dofs on the edges of the cells (start
   with the top and bottom edges, then do the vertical edges) and add
   the index (depending on whether the edge is owned by the local rank
   or not) to a store of edge dofs.

#. Loop over all the columns of dofs on the vertices of the cells and
   add the index (depending on whether the vertex is owned by the local
   rank or not) to a store of vertex dofs.

#. Repeat process over the next depth of halo cells until the maximum
   halo depth is reached.

The sequence above results in four arrays, one for each type of cell
entity (volume, face, edge and vertex) with dof indices. There are no
duplicated columns of dofs in these arrays and the columns have indices
that are positive for owned dofs and increasingly negative for
annexed/halo dofs as the depth of halo increases.


The final dofmap needs to be ordered by cell. :numref:`dofmap` shows a graphical
representation of the field ordering process. All the cells are looped over
again, and the indices for the columns of dofs for each cell are extracted. The
ordering by entity is retained within each cell, which can be seen in
:numref:`dofmap`-a for cell number 1. Once all cells in a horizontal level have
been looped over, the data is ordered by entity first, then by cell number
(:numref:`dofmap`-b). If two different cells share a column of dofs, both will
point to the same column of the field data.

The annexed/halo negative indices are flipped and added to the end of
the owned dofs (:numref:`dofmap`-c), so that annexed/halo dofs
appear contiguously at the end of the field data :numref:`dofmap`-d.

.. _dofmap:
.. figure:: images/dofmap.svg

  The ordering of data in a field data array. a) The order of data shown in
  detail for Cell number 1, only (the dofs for other cells follow similar
  ordering). b) The order of data once all the cells have been looped over.
  Owned and annexed dofs can only appear on owned cells and halo dofs can only
  appear on Halo cells. Dofs that are shared between cells appear with the
  first cell processed that shares them. c) Once all cells have been processed
  the order of annexed and halo dof indices are flipped and the dof indices
  are placed after the owned dofs. d) The final ordering of dofs in the
  field array.

Generating a unique global dof index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each dof must be given a global index so it can be uniquely identified.
If four partitions share a dof, each partition knows how to find the
same global index and its "owner"

Four qualities uniquely identify every dof:

#. The global index of the "owning" cell

#. The number of the dof within the cell

#. The number of the vertical level of the cell that contains the dof

#. The number of the entry within a multi-data field (this will be 1 for
   a standard single-valued field)

So, the equation for setting a zero-based global dof index is

::

        global_dof_id = (global_cell_id-1) * ndofs*nlayers*ndata +
                        (idof-1) * nlayers*ndata +
                        (k-1) * ndata +
                        (m-1)

| Where:
| ``global_cell_id`` is the global cell index of the owning cell,
| ``ndofs`` is the maximum number of dofs in a cell,
| ``nlayers`` is the number of layers,
| ``ndata`` is the total number of entries in a multi-data field,
| ``idof`` is the number of the specific dof within the cell,
| ``k`` is the layer number,
| ``m`` is the multi-data entry number.

The difficulty, here, is making sure the global cell index of the owning
cell is correct. As mentioned above, the owning cell is the cell with
the highest global cell index that shares a dof. The global index can be
calculated locally on each partition easily for interior cells, but any
dof on the boundary (i.e. at the edge of the outermost halo) can't
know the cell index of all the cells that surround it (as they are not
all known to the local partition). An extra halo known as the "ghost"
halo is added only for the purpose of assigning ownership of these edge
dofs. The ghost halos are not used in any other part of LFRic - they are
not iterated over and they are not involved in halo exchanges.

Application Programming Interface
=================================

The distributed memory functionality is only applicable to field and
scalar objects and is only available within the Psy layer. To prevent
its use elsewhere, it has been made available from the field proxy - an
object that is only available within the Psy layer. Scalar objects are
not used outside the Psy layer, so no proxy is required.

The following functions provide distributed memory support on field
data:

-  ``field_proxy%halo_exchange( depth )``

   | ``depth`` (in) The depth to which the halos should be exchanged
   | Performs a blocking halo exchange operation on the field.

-  ``field_proxy%halo_exchange_start( depth )``

   | ``depth`` (in) The depth to which the halos should be exchanged.
   | Starts a halo exchange operation on the field. The halo exchange is
     non-blocking, so this call only starts the process. On Return,
     outbound data will have been transferred, but no guarantees are
     made for in-bound data elements at this stage.

-  ``field_proxy%halo_exchange_finish( depth )``

   | ``depth`` (in) The depth to which the halos have been exchanged.
   | Wait (i.e. block) until the transfer of data in a halo exchange
     (which has been initiated by a call to ``halo_exchange_start``) has
     completed.

-  | ``field_proxy%get_sum()``
   | Perform a global sum operation on the field.

-  | ``field_proxy%get_min()``
   | Calculate the global minimum of the field.

-  | ``field_proxy%get_max()``
   | Calculate the global maximum of the field.

-  ``field_proxy%is_dirty( depth )``

   | ``depth`` (in) The depth at which to check the halos
   | Returns whether the halos at the given depth are dirty or clean.

-  | ``field_proxy%set_dirty()``
   | Flags all halos as being dirty.

-  ``field_proxy%set_clean( depth )``

   | ``depth`` (in) The depth up to which to set the halo to clean
   | Flags all the halos up the given depth as clean.

The following functions provide distributed memory support for scalar
data:

-  | ``scalar%get_sum()``
   | Perform a global sum operation on the scalar.

-  | ``scalar%get_min()``
   | Calculate the global minimum of the scalar.

-  | ``scalar%get_max()``
   | Calculate the global maximum of the scalar.


References
==========

.. [LFRic_Design] The LFRic Team, *The design, data model and implementation of
   the GungHo dynamical core for the Met Office.*, LFRic Documentation, 2016.

.. [UGRID_Convention] Jagers, B., *et al.*, *UGRID Conventions Version 1.0*,
   http://ugrid-conventions.github.io/ugrid-conventions, 2016.
