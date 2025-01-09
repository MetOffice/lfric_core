.. ------------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _lfric requirements:

Requirements summary
====================

The following lists some of the major requirements that need to be met by the
LFRic Infrastructure. It should be noted that some of these requirements
underpin other complex requirements that are not described in detail, such as
the need for a comprehensive diagnostic system.

Where implementation is substantially incomplete, or where there are
notable omissions, the current status of the LFRic implementation in
meeting these requirements is also described. Otherwise, it can be
assumed that the requirement has been met to a certain degree, noting
that further enhancements in some capabilities are anticipated, but not
documented here.

-  Support for fields on global meshes, rectangular limited area and
   bi-periodic meshes, and lateral boundary condition (LBC) meshes. In
   principle, LFRic can support any mesh. In practice, the global mesh
   support is currently limited to the cubed sphere. Support for LBC
   meshes is in development.

-  Support for fields on meshes made up of columns of cells which are
   quadrilateral in the horizontal.

-  Support for finite element function spaces with higher order
   functions. This means the library must allow configuration of and
   access to the basis functions for the field, and the ability to store
   field data comprising more than one data point on more than one mesh
   entity (vertices, edges, faces or cell volumes).

-  Support for finite difference scalar and vector fields. Currently
   this support is provided by using a subset of the finite element
   implementation.

-  Support for organising and grouping fields as required by
   applications. For example, the LFRic atmosphere requires support for
   collections of fields of different types and support for tiled
   fields.

-  Distributed memory support is required based on horizontal domain
   decomposition with support for halos and halo exchange operations.
   Currently, LFRic includes support for partitioning a global
   cubed-sphere mesh, a limited area mesh and a bi-periodic mesh into
   rectangular meshes.

-  Support is required for different looping strategies to enable an
   operation over all, or a subset, of data points or cells in a field.
   LFRic provides the ability for kernels to operate on columns of cells
   (and the data points they contain), or individual levels (for certain
   field types). PSyclone support for the latter is in development.

-  Fields require support for halos of any depth. Support is required
   for looping over data points or on data on cells up to different halo
   depths. Currently this is supported, but all fields in an application
   are partitioned with the same depth of halo even if a smaller halo is
   needed for some fields.

-  To allow for concurrent communication and computation, support for
   looping over data on cells within or outside so-called "inner halos"
   is required. This requirement and implementation is described in
   detail in the :ref:`distributed memory
   <overlap comms compute>` section.

-  For kernels that operate on continuous fields (those in which
   data points are shared between neighbouring columns), support is
   required for looping strategies (such as "colouring" or "tiling") that
   ensure concurrent shared-memory operations on the same shared data point
   are prevented.

-  Support for the GungHo dynamics is required. GungHo is a mixed finite
   element scheme, therefore support for placing of fields with
   different function space types on the same mesh is required.

-  Support for fields on hierarchical meshes. Each mesh in the hierarchy
   covers the same geographical domain, but the resolution of the lower resolution
   meshes is an integer multiple of the higher resolution mesh. For example, if
   the resolution of the higher resolution mesh is 10km then the lower
   resolution mesh may have a resolution of 20km or 30km. If the two meshes are
   aligned, then each cell in the lower resolution mesh would contain
   :math:`2 * 2` or :math:`3 * 3` cells of the higher
   resolution mesh. The hierarchy is intended to support multigrid solvers, and will
   enable easier transformation of data between different science
   schemes running at different resolutions, where interpolation can be done
   local to a distributed memory domain.

-  Support for mapping stencils of various shapes and depths to allow
   operations on cell data that are dependent on data in nearby cells.

-  Support for output and input of field data to and from files, to
   underpin a requirement for reading and writing dumps, reading
   ancillary files and writing diagnostic data.

-  Support for basic infrastructure such as clocks and calendars for
   managing time-stepping and long runs, and for log messages to provide
   progress messages to output files.

-  Support for interfacing to external model couplers such as OASIS.
