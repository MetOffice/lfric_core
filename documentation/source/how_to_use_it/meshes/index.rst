.. -----------------------------------------------------------------------------
    (c) Crown copyright 2024 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _section mesh generation:

Mesh generation
=====================

LFRic fields require an association with a 3D-mesh object. Such a mesh object
requires formal definition of its entities (nodes,edges,faces,volumes) and how
they are interconnected. LFRic 3D-mesh objects are constructed internally by
the LFRic application via extrusion of 2D-meshes which are read in from file.

This section aims to describe supported 2D-meshes and how to generate them for
use in an LFRic application.

.. toctree::
   :caption: Contents
   :maxdepth: 4

   mesh_generators
   mesh_configuration_namelists
   mesh_file_content
   mesh_ugrid_compliance
