.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section Configuration file:


===================
Configuration file
===================

The mesh generators are controlled via a configuration file containing FOTRAN namelists. Some namelists are only relevant
depending on which generator program is used or are triggered by settings in the configuration file.

.. _mesh_nml:

&mesh (Cubed-Sphere|Planar)
  Control namelist for common variables of principle meshes requested (Cubed-Sphere|Planar). Principle meshes are those which are not derived from another mesh. 

.. _coord_sys:

* ``coord_sys``:``['xyz'|'ll']``
    Coordinate system used to locate mesh nodes/features. As the mesh topologies are 2D, there are only 2-coordinates per node. Cartesian coordinates (x,y) in metres or Spherical coordinates (lon,lat) in degrees. In practice, cartesian coordinates are generally reserved for idealised cases. For cubedsphere meshes, set to 'll'
* ``geometry``:``['planar'|'spherical']``
    Geometrical shape of domain.
* ``topology``:``['non_periodic'|'channel'|'periodic']``
    Describes periodicity types at domain boundaries. For cubedsphere meshes, set to 'periodic'

.. _n_meshes:

* ``n_meshes``: ``<integer>``
    Number of principle mesh topologies to be requested by the generator.

.. _mesh_names:

* ``mesh_names``: ``<comma-separated strings>``
    Names to apply to principle mesh topologies, only the first :ref:`n_meshes<n_meshes>` are applied.

.. _mesh_maps:

* ``mesh_maps``: ``<comma-separated strings>``
    Requested pairs of meshes which require creation of intermesh maps. Format of each string is ``'<mesh_name_A>:<mesh_name_B>'``, where these names will appear in :ref:`mesh_names<mesh_names>`.

.. _rotate_mesh:

* ``rotate_mesh``:``[.true.|.false.]``
    Transform principle meshes according to :ref:`&rotation<rotation_nml>` namelist.

.. _partition_mesh:

* ``partition_mesh``:``[.true.|.false.]``
    Partition principle meshes according :ref:`&partitions<partitions_nml>` namelist.

.. _rotation_nml:

&rotation (Cubed-Sphere|Planar, :ref:`rotate_mesh<rotate_mesh>` =.true.)
  Control namelist for rotation of node coordinates (Cubed-Sphere|Planar).

* ``rotation_target``:``['north_pole'|'null_island']``.
    Feature to use as fixed reference for rotation.
* ``target_north_pole``:``<longitude>,<latitude>``.
    Real world location of mesh north pole after rotation.
* ``target_null_island``:``<longitude>,<latitude>``.
    Real world location of mesh null island after rotation.

.. _stretch_transform_nml:

&stretch_transform (Planar, :ref:`apply_stretch_transform<apply_stretch_transform>` =.true.)
  Control namelist for stretched grid transformation (Planar).

* ``cell_size_inner``:``<x_cell_size>,<y_cell_size>``
    Domain inner region cell sizes along x/y axes.
* ``cell_size_outer``:``<x_cell_size>,<y_cell_size>``
    Domain outer region cell sizes along x/y axes.
* ``n_cells_outer``:``<x_n_cells>,<y_n_cells>``
    Depth (in cells) of domain outer region along x/y axes.
* ``n_cells_stretch``:``<x_n_cells>,<y_n_cells>``
    Depth (in cells) of domain stretch region along x/y axes.
* ``stretching_on``:``['cell_centres'|'cell_nodes'|'p_points']``
    Features to use as anchor points for stretch transform.
* ``transform_mesh``:``<mesh_name>``
    Principle mesh to apply stretch transform to. Any meshes connected to this mesh via InterMesh maps (:ref:`mesh_maps<mesh_maps>`) will have their node locations updated accordingly.

.. _partitions_nml:

&partitions (Cubed-Sphere|Planar, :ref:`partition_mesh<partition_mesh>` =.true.)
  Control namelist for partitioning of mesh domains. This results in multiple mesh files containing meshes local to a given processor rank (Cubed-Sphere|Planar).

* ``max_stencil_depth``:``<integer>``
    The maximum stencil depth supported by the resulting local meshes. In effect, this controls the halo depth around each mesh partition. 

.. _n_partitions:

* ``n_partitions``:``<integer>``
    The total number of requested partitions. For cubed-spheres this shoud be 1 or a multiple of 6.
* ``panel_decomposition``: ``'auto'|'row'|'column'|'custom'``
    Specifies partition strategy to apply to cubedsphere panel/planar mesh.
* ``panel_xproc`` (``panel_decomposition='custom'``)
    Number of partitions in local x-direction of panel
* ``panel_yproc`` (``panel_decomposition='custom'``)
    Number of partitions in local y-direction of panel
* ``partition_range`` (``<start partition id>,<end partition id>``)
    Specifies start/end partition ids to output, valid ids range from 0 to :ref:`n_partitions<n_partitions>`-1. Will produce 1 file per requested partition.

&cubedsphere_mesh (Cubed-Sphere)
  Control namelist for Cubed-Sphere mesh generation.

* ``edge_cells``:``<comma separated integers>``
    Number of cells along edge of mesh panel. Only the first n-entries given be :ref:`n_meshes<n_meshes>` are used. 
* ``smooth_passes``:``<integer>``
    Number of interations of smoothing function applied to mesh node locations. This has most impact at the corners of the cube topology.
* ``equatorial_latitude``:``<real>``
    Real world latitude (degrees) of Cubed-Sphere mesh equator after applying Schmit transform.
    The `top`(or `bottom`) panels of the cubed-sphere are reduced in size while maintaining the
    same connectivity. This has the effect of a localised increase in resolution over a panel of
    the cubed-sphere without increasing the overall number of cells in the mesh.

&planar_mesh(Planar)
  Control namelist for Planar mesh generation.

.. _apply_stretch_transform:

* ``apply_stretch_transform``:``<logical>``
    Transform grid to a stretched grid with inner, outer and transition regions.
* ``create_lbc_mesh``:``<logical>``
    Generate a rim mesh from a priciple mesh.
* ``domain_centre``:``[x,y]``
    Location of domain centre for all requested meshes, [m/degrees determined by :ref:`coord_sys<coord_sys>`].
* ``domain_size``:``[x,y]``
    Domain size for all requested meshes [m/degrees determined by :ref:`coord_sys<coord_sys>`].
* ``edge_cells_x``: ``<integer>, ...``
    Number of cells along x-axis of domain. List of integers corresponding to :ref:`mesh_names<mesh_names>`. Only the first :ref:`n_meshes<n_meshes>` entries used. (>1 for :ref:`periodic_x<periodic_x>` =.true.)
* ``edge_cells_y``: ``<integer>, ...``
    As above, but for y-axis. (>1 for :ref:`periodic_x<periodic_y>` =.true.``)
* ``lbc_parent_mesh``:``<string>``
    Name of the principle mesh to create a rim mesh from, this name should exist in the :ref:`mesh_names<mesh_names>` variable.
* ``lbc_rim_depth``:``<integer>``
    Depth (in cells) of rim mesh.

.. _periodic_x:

* ``periodic_x``:``<logical>``
    Periodicity across domain boundaries in x-direction.

.. _periodic_y:

* ``periodic_y``:``<logical>``
    Periodicity across domain boundaries in y-direction.

.. note:: Any configuration variable of polar coordinates are to be specified in ``<Degrees Longitude>, <Degrees Latitude>`` on real world frame of reference.
