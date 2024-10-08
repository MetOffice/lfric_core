.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section meshes and tools:

Meshes
=====================

LFRic fields require an association with a 3D-mesh. Such a mesh requires formal definition of its entities (nodes,edges,faces,volumes) and how they are interconnected. Field 3D-meshes are constructed internally via the extrusion of pre-generated 2D-meshes dimensioned to the models surface extents.

This section aims to describe supported 2D-meshes and how to generate them. 

.. toctree::
   :caption: Contents
   :maxdepth: 2

   quick_start
   mesh_generator_namelists
   mesh_generation
   mesh_topologies
   transformations
   intermesh_mappings
   prepartitioning
   mesh_metadata
