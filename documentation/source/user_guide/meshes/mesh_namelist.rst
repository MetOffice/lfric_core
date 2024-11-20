.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section mesh_nml:

.. _mesh_nml:
============================
``&mesh``
============================
This is the main controlling namelist for configuring the principle meshes
[#f1]_ in the output file. Variables in this namelist are common to both mesh
generators and are applied to each mesh topology independently. This namelist
is required in the configuration file, `e.g.` configuration.nml.

.. °C °E

.. _coord_sys:

* ``coord_sys``: **'<string>'**
    The coordinate system used to locate mesh nodes/features in the output
    file. Valid options:

    ``ll``
      Spherical coordinates (lon,lat) in (Degrees East, Degrees North).

    ``xyz``
      Cartesian coordinates, This is only supported for a flat planar mesh at
      z=0, as a result features are located with only 2-coordinates, ``(x,y)``
      in metres.

* ``geometry``: **'<string>'**
    Geometrical shape of the surface domain. Valid options:

    ``planar``
      Planar surface geometry.

    ``spherical``
      Curved surface geometry (spherical).

.. _mesh_names:

* ``mesh_names``: **'<string>', …**
    Names applied to principle mesh topologies.  The number of names
    should match the value given by :ref:`n_meshes<n_meshes>`.
    The order of appearance has no effect on the generation of each
    requested mesh.

.. _mesh_maps:

* ``mesh_maps``: **'<string>', …**
    Each listing requests generation of intermesh mappings between
    the specified meshes in the output file. Intermesh maps provide
    a link from cells on one mesh that spatially overlap cells on
    another mesh. These maps are restricted to pairs of meshes where
    the cells of one mesh are a sub-division of cells in the other.

    A mapping is given as a string in the form
    ``'<mesh_name_A>:<mesh_name_B>'``, where these names should
    appear in the variable :ref:`mesh_names<mesh_names>`.

.. _n_meshes:

* ``n_meshes``: **<integer>**
    Number of principle mesh topologies output by the generator.

.. _partition_mesh:

* ``partition_mesh``: **<logical>**
    Partition principle meshes according to configuration given
    by :ref:`&partitions<partitions_nml>` namelist.

.. _rotate_mesh:

* ``rotate_mesh``: **<logical>**
    Transform principle meshes according to configuration given
    by :ref:`&rotation<rotation_nml>` namelist.

* ``topology``:  **'<string>'**
    Describes periodicity type for opposing domain bounds. Valid options:

    ``non_periodic``
      (`Planar meshes only`) All domain boundaries are non-periodic, as a
      result there is no connectivity information at the domain
      boundaries. Crossing a non-periodic boundary enters a void [#f2]_
      region.

    ``channel``
      (`Planar meshes only`) Of the four domain boundaries, two are linked as
      a periodic pair, the remaining boundaries are non-periodic. Crossing one
      of the periodic domain boundaries re-enters the domain at a point on
      other periodic domain boundary.

    ``periodic``
      Opposing domain boundaries are linked as periodic pairs.



.. rubric:: Footnotes

.. [#f1] Principle meshes are those meshes explicity requested
   and named in the configuration file, `i.e.` |nbsp|
   :ref:`mesh_names<mesh_names>`. These meshes are generated from
   the base strategy as opposed to being derived from another mesh
   topology `e.g.` rim mesh.

.. |degree| unicode:: U+00B0 .. degree symbol
.. |nbsp|   unicode:: U+00A0 .. no-break space symbol
   :trim:
