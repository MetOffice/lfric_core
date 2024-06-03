.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section mesh generator:

The Mesh Generator
==================

The mesh generator creates the meshes required by LFRic model
applications. It generates single-layer meshes comprising faces, edges
and nodes.

It supports generation of the following types of meshes:

 - Cubed-sphere meshes of arbitrary number of cells. A `Cn` mesh
   describes a cubed-sphere mesh where each of the six panels has `n`
   by `n` cells. The coordinates of the mesh map onto a sphere. By
   default the mesh is TODO with six identically-sized panels and with
   top and bottom panels centred on the North and South poles. Various
   options can be applied to smooth, stretch or rotate the mesh
   coordinates
 - Planar meshes which are rectangular and may be periodic on one or
   both pairs of edges. The coordinates can be on a flat surface with
   cell sizes defined, or can be mapped onto a sphere to support
   limited area models (LAMs). Planar meshes can be stretched to make
   certain rows or columns larger. LAM meshes are constructed such that,
   if rotated around the sphere such that it is centred on the
   equator, the rows and columns will align with lines of latitude and
   longitude, respectively.
 - LBC meshes comprising cells that overlay the rim of one or more
   cells deep of a planar LAM mesh.

The mesh generator can construct hierarchical meshes, where each mesh
in the hierarchy is created by subdividing all the cells in a coarser
mesh into the same number of subcells.

Where meshes are related, the mesh generator can output mesh maps that
map the cells of one mesh onto the cells of another. Mesh maps are
relevant both for mapping hierarchical meshes (planar or cubed-sphere)
or for mapping LAM and LBC meshes.

The mesh generator can partition a set of meshes. Partitioned meshes
are important for high resolution models running on large numbers of
processes as the run-time cost of partitioning meshes in such
configurations can be high.

Mesh output is aligned with the UGRID standard.


TODO: Fully-document the mesh generator.
