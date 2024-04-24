.. ------------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENSE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

##################
Glossary
##################
.. glossary::
  DataModel
    LFRic core infrastructure which provides the framework for working with data objects supporting LFRic ScienceModel(s).

  Global Mesh
    An object describing a 2D-mesh which encompasses the entirety of an LFRic application Model Domain.

  InterMesh Map
    Mapping between 2D-mesh cell IDs from source-to-target meshes. The mapping lists the cell IDs (local to the target mesh) that geographically overlap a given cell id in the source mesh (ID local to the source mesh).

  LFRic application
    Top level program executable. An application is not part of LFRic core, though will employ the LFRic core DataModel to work with the ScienceModel(s) accessed by the application.

  Local Mesh
    An object describing a 2D-mesh which encompasses a sub-section (partition) of a Global Mesh. The local mesh object contains related mesh information similar to that of a Global Mesh. In addition, it also stores partition information unique to it's intended process rank.

  Mesh
    An object describing a 3D-mesh which encompasses a sub-section (partition) of the Global Mesh. A mesh object is derived from the corresponding Local Mesh object for a given process rank and an extrusion object configured by an LFRic application.

  Model Domain
    The geographical domain extents to which an LFRic application is configured, *e.g.* GCM, LAM.

  ScienceModel
     Model representing a real world scientific process(es).










