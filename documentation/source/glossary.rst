.. -----------------------------------------------------------------------------
    (c) Crown copyright 2024 Met Office. All rights reserved.
    The file LICENSE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

:html_theme.sidebar_secondary.remove: true
:orphan:

.. _glossary_of_terms:

Glossary
========

.. glossary::
  Annexed DoF
    DoFs on owned cells that are actually owned by another rank that
    shares them.

  Data Model
    LFRic core infrastructure which provides the framework for working
    with data objects supporting LFRic ScienceModel(s).

  DoF
    Degree of Freedom. DoFs describe the state of a field. They typically
    correspond to the value of the field at specific locations in the mesh
    (although they don't have to), so are equivalent to data values in
    finite difference models

  Edge Cell
    These are the cells that lie on the boundary of a partition. They contain
    DoFs that are shared with halo cells. Thery are not to be confused with
    "edge dofs" which are the DoFs that lie oon the edges of a cell.

  Ghost Cell
    In an LFRic context, ghost cells are not the same as halo cells. They
    are a special subset of cells used during initialisation to determine
    ownership of shared degrees of freedom (DoFs) on the outer faces,
    edges, or vertices of the outermost halo cells. They are not part
    of the active computation and are discarded after the initialisation
    phase.

  Global DoF Index
    A unique identifier for every DoF. Every partition uses the same
    global DoF index for the same DoF. If four partitions share a DoF,
    each partition knows how to identify that same DoF and its owner.


  Global Mesh
    An object describing a 2D-mesh which encompasses the entirety of
    an LFRic application Model Domain.

  Halo
    A halo is made from halo cells and refers to the extra layers of cellsa
    surrounding a processor's local domain. These halo cells are used to
    store data from neighboring domains, enabling each processor to
    perform calculations that depend on adjacent values without needing
    constant inter-processor communication.

  Halo DoF
    DoFs on halo cells, that are owned by other ranks.

  Inner Halo
    A subset of owned cells that are not affected by halo exchanges, and
    so, can be used for computation while communication is still ongoing.

  InterMesh Map
    Mapping between 2D-mesh cell IDs from source-to-target meshes. The
    mapping lists the cell IDs (local to the target mesh) that
    geographically overlap a given cell id in the source mesh (ID
    local to the source mesh).

  LFRic application
    A program that uses the LFRic infrastructure.

  Local Mesh
    An object describing a 2D-mesh which encompasses a sub-section
    (partition) of a Global Mesh. The Local Mesh object contains
    related mesh information similar to that of a Global Mesh. In
    addition, it also stores partition information unique to its
    intended process rank.

  Mesh
    An object describing a 3D-mesh which encompasses a sub-section
    (partition) of the Global Mesh. A mesh object is derived from the
    corresponding Local Mesh object for a given process rank and an
    extrusion object configured by an LFRic application.

  Mesh Entity
    A fundamental geometric component of the computational mesh. The mesh
    entities comprise of the 3d volume, 2d faces, 1d edges and 0d vertices.

  Model Domain
    The geographical domain extents to which an LFRic application is
    configured, *e.g.* GCM, LAM.

  Modeldb
    An object that encapsulates all the data required to describe both
    the scientific and technical state of a model.

  Owned DoF
    DoFs on owned cells and are owned by the local rank.

  Partitioning
    A way of dividing a system's computational domain (including data and
    workload) into smaller sections. These sections (called partitions)
    can be processed in parallel to improve performance, scalability, and
    maintainability.

  Science Model
     A library of code that simulates a particular science process,
     *e.g.* a radiation or land surface model. **Note:** An LFRic
     application may be written to access any number of Science
     Models.
