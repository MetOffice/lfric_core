.. -----------------------------------------------------------------------------
    (c) Crown copyright 2024 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _section configuration namelists:

========================
Configuration namelists
========================

The mesh generators are controlled via a configuration file containing Fortran
namelists. The following namelists are mandatory (depending on which generator
is used):

* :ref:`&mesh<mesh_nml>`: Required
* :ref:`&cubedsphere_mesh<cubedsphere_mesh_nml>`: Cubed-Sphere mesh generator
  only.
* :ref:`&planar_mesh<planar_mesh_nml>`: Planar mesh generator only.

All other namelists are optional, although may still only be applicable
depending on generator.

.. _mesh_nml:

.. dropdown:: ``&mesh``

  This is the main controlling namelist for configuring the principle meshes
  [#f1]_ in the output file. Variables in this namelist are common to both mesh
  generators and are applied to each mesh topology independently. This namelist
  is required in the configuration file, `e.g.` configuration.nml.

  .. _coord_sys:

  * ``coord_sys``: **'<string>'**
      The coordinate system used to locate mesh nodes/features in the output
      file. Valid options:

      ``ll``
        Spherical coordinates (lon,lat) in (°E,°N).

      ``xyz``
        Cartesian coordinates, This is only supported for a flat planar mesh
	at z=0, as a result, features are located with only 2-coordinates,
	``(x,y)`` in metres.

  * ``file_prefix``: **'<string>'**
      File prefix used to construct output filename(s). Depending on the logical,
      :ref:`partition_mesh<partition_mesh>`, the output filename(s) are constructed as:

      **<file_prefix>.nc**
        Single file output for non-partitioned mesh topologies. Requested meshes
	encompass an entire model domain.

      **<file_prefix>_<partition_id>_<n_partitions>.nc**
        Multiple files output for pre-partitioned mesh topologies. Each file
	contains the cells of the requested meshes on a given partition of
	the model domain. Pre-partitioned meshes allow a simple load operation
	to populate an application's local mesh objects on a mpi rank.

	**Note:** Pre-partitioned meshes are only suitable for applications
	running with mpi ranks equal to the number of partitions.

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
        (`Planar meshes only`) Of the four domain boundaries, two are linked
	as a periodic pair, the remaining boundaries are
	non-periodic. Crossing one of the periodic domain boundaries re-enters
	the domain at a point on other periodic domain boundary.

      ``periodic``
        Opposing domain boundaries are linked as periodic pairs.

.. _cubedsphere_mesh_nml:

.. dropdown:: ``&cubedsphere_mesh``

  Control namelist for cubed-sphere mesh generation, required for use with the
  `cubedsphere_mesh_generator`. This creates a mesh which uses the
  cubed-sphere base strategy and requires ``edge_cells`` to define the mesh
  connectivity. The remaining options are mesh transformations that are
  subsequently applied to the node coordinates; the connectively of the mesh
  elements are not altered by transformations.

  * ``edge_cells``: **<integer>, …**
      Number of cells along edge of each mesh panel. The sequence of integers
      will map to entries given by :ref:`mesh_names<mesh_names>`.
  * ``equatorial_latitude``: **<real>**
      Real world latitude (°N) of cubed-sphere mesh equator after
      applying Schmit transform. The `top` (or `bottom`) panels of the
      cubed-sphere are reduced in size while maintaining the same
      connectivity. This has the effect of a localised increase in resolution
      over a panel of the cubed-sphere without increasing the overall number
      of cells in the mesh.
  * ``smooth_passes``: **<integer>**
      Number of interations of smoothing function applied to mesh node
      locations.

.. _planar_mesh_nml:

.. dropdown:: ``&planar_mesh``

  Control namelist for planar mesh generation, required for use with the
  `planar_mesh_generator`. This creates a mesh which uses the planar mesh
  base strategy and requires ``edge_cells`` along `both` axes aswell as the
  domain boundary periodicity in order to define the mesh connectivity. The
  remaining options are mesh transformations or logical triggers. Mesh
  tranformations are applied to the node coordinates after the base mesh is
  generated; the connectively of the mesh elements are not altered by
  transformations.

  .. _apply_stretch_transform:

  * ``apply_stretch_transform``: **<logical>**
      Apply the stretch transform to base planar mesh node coordinates.
  * ``create_lbc_mesh``: **<logical>**
      Generate a rim mesh which is derived one of the priciple meshes
      (:ref:`mesh_names<mesh_names>`).
  * ``domain_centre``: **<real>,<real>**
      Location of domain centre for all principle meshes. Coordinates aligned
      with the :ref:`coord_sys<coord_sys>` with units of `[m|°]`
      as appropriate, `i.e.` `x,y` or `lon,lat` coordinates.
  * ``domain_size``: **<real>,<real>**
      Domain size for all principle meshes.  Domain sizes aligned with the
      :ref:`coord_sys<coord_sys>` with units of `[m|°]` as
      appropriate, `i.e.` `x,y` or `lon,lat` domain extents.
  * ``edge_cells_x``: **<integer>, ...**
      Number of cells along x-axis of domain. Order of integers map to
      :ref:`mesh_names<mesh_names>`.
  * ``edge_cells_y``: **<integer>, ...**
      As above, but for y-axis.
  * ``lbc_parent_mesh``: **'<string>'**
      Name of the principle mesh to create a rim mesh from, this name should
      exist in the :ref:`mesh_names<mesh_names>` variable.
  * ``lbc_rim_depth``: **<integer>**
      Depth (in cells) of rim mesh. This is the number of cells radially
      across the rim mesh from the domain centre to a domain boundary.

  .. _periodic_x:

  * ``periodic_x``: **<logical>**
      Periodicity across pair of domain boundaries in x-axis.

  .. _periodic_y:

  * ``periodic_y``: **<logical>**
      Periodicity across pair of domain boundaries in y-axis.

.. _partitions_nml:

.. dropdown:: ``&partitions``

  Optional control namelist for partitioning of mesh domains. Use of this
  namelist allows meshes to be partitioned by the mesh generators rather than
  by an application at runtime. Principle meshes are partitioned and written
  to file as 1 partition per file.

  This functionality results in multiple mesh files with the corresponding
  portion of each principle mesh. In addition to the mesh tolopologies,
  partition information is also written to each file. In essence, each file
  provides the information required to load and instantiate 2D-mesh objects
  which are `local` to a given process rank. This namelist is enabled if
  triggered by the :ref:`partition_mesh<partition_mesh>` logical.

  .. NOTE:: Care should be taken when specifying partition
	    configurations. Partitions flag mesh cell ids as being members of
	    that parition, there is no restriction that limits the shape or
	    continuity of a partition.

  * ``max_stencil_depth``: **<integer>**
      For communication across partitions, the `local` mesh domains need to be
      extended with a `halo` region applied around each partition. The
      required depth of the halo region depends on a specific LFRic
      application's runtime configuration. The `max_stencil_depth` (in cells)
      specifies the anticipated stencil size that the partitioned mesh should
      support.

  .. _n_partitions:

  * ``n_partitions``: **<integer>**
      The total number of requested partitions. For `cubed-sphere` meshes,
      this shoud be restricted to 1 or a multiple of 6.

  .. _panel_decomposition:

  * ``panel_decomposition``: **'<string>'**
      Specifies panel partition strategy applied to principle meshes. The
      generators use the partitioning module support from LFRic core
      infrastruture. Valid options:

      ``auto``
        The infrastructure code will attempt to group cells into partitions
        which are uniform and square as possible.

      ``row``
        Forces the mesh to be divided as a single row (in x-axis) of
        :ref:`n_partitions<n_partitions>`.

      ``column``
        Forces the mesh to be divided as a single column (in y-axis) of
        :ref:`n_partitions<n_partitions>`.

      ``custom``
        Forces the partitioner to attempt to configure the panel into
	partitions given by :ref:`panel_xproc<panel_xproc>` and
	:ref:`panel_yproc<panel_yproc>`.

  .. _panel_xproc:

  * ``panel_xproc``: **<integer>**
      Number of partitions in local x-direction of mesh panel, this variable
      is only valid when requesting a ``custom`` decomposition
      (:ref:`panel_xproc<panel_xproc>`).

  .. _panel_yproc:

  * ``panel_yproc``:  **<integer>**
      Number of partitions in local y-direction of mesh panel, this variable
      is only valid when requesting a ``custom`` decomposition
      (:ref:`panel_yproc<panel_yproc>`).

  * ``partition_range``: **<integer>,<integer>**
      Specifies the ``start``, ``end`` partition ids to output, valid ids range from
      [0: :ref:`n_partitions<n_partitions>`-1]. The generators will produce 1
      file per requested partition, with the partition id tagged to the output
      filename.

.. _rotation_nml:

.. dropdown:: ``&rotation``

  Optional control namelist for rotation of node coordinates. Rotation in
  [longitude, latitude] is specified by using a reference location, which is
  rotated such that it arrives at the specified target location. All nodes
  then undergo the same coordinate transformation.

  Enabled by :ref:`rotate_mesh<rotate_mesh>` (cartesian planar meshes not
  supported).

  .. NOTE:: When considering meshes with rotation applied, all (lon,lat)
	    coordinates which are referenced in the configuration/output files
	    are with respect to a real world frame of reference.



  * ``rotation_target``: **'<string>'**
      Feature to use as reference for rotation.

      ``north_pole``
        Use North Pole (90°N,0°E) as reference location.

      ``null_island``
        Use Null Island (0°N,0°E) as reference location.

  * ``target_north_pole``: **<real>,<real>**
      Target location of North Pole after rotation.
  * ``target_null_island``: **<real>,<real>**
      Target location of Null Island after rotation.

.. _stretch_transform_nml:

.. dropdown:: ``&stretch_transform``

  Optional control namelist for stretched grid transformation. The stretched
  grid transform modifies the base planar mesh resulting in two regions of
  differing cell size with a stretch region separating them. The intention is
  to provide (for a given axis) an inner cell size (inner region) which
  gradually transitions (stretch region) to the outer cell size (outer
  region). Only valid for planar meshes, this transform is enabled by the
  :ref:`apply_stretch_transform<apply_stretch_transform>` logical.

  * ``cell_size_inner``: **<real>,<real>**
      Cell size for domain inner region of stretch grid along
      ``x-axis,y-axis``, units in `[m|°]` depending on planar mesh
      type.
  * ``cell_size_outer``: **<real>,<real>**
      Cell size for domain outer region of stretch grid along
      ``x-axis,y-axis``, units in `[m|°]` depending on planar mesh
      type.
  * ``n_cells_outer``: **<integer>,<integer>**
      Depth (in cells) of domain outer region along ``x-axis,y-axis``.
  * ``n_cells_stretch``:  **<integer>,<integer>**
      Depth (in cells) of domain stretch region along ``x-axis,y-axis``.
  * ``stretching_on``: **'<string>'**
      Features to use as anchor points for stretch transform.

      ``cell_centres``
        Use cell centres as anchor points.

      ``cell_nodes``
        Use cell nodes as anchor points.

      ``p_points``
        Use p-points as anchor points.

  * ``transform_mesh``: **'<string>'**
      Principle mesh to apply stretch transform to. Any meshes connected to
      this mesh via InterMesh maps (:ref:`mesh_maps<mesh_maps>`) will have
      their node locations updated accordingly.


.. rubric:: Footnotes

.. [#f1] Principle meshes are those meshes explicity requested
   and named in the configuration file, `i.e.` |nbsp|
   :ref:`mesh_names<mesh_names>`. These meshes are generated from
   the base strategy as opposed to being derived from another mesh
   topology `e.g.` rim mesh.

.. [#f2] The void region is where there is no connected mesh after crossing the bounds of the domain.

.. |degree| unicode:: U+00B0 .. degree symbol
.. |nbsp|   unicode:: U+00A0 .. no-break space symbol
   :trim:
