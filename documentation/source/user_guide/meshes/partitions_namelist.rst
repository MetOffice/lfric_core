.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section partitions namelist:

.. _partitions_nml:

=================
``&partitions``
=================

Optional control namelist for partitioning of mesh domains. Use of this
namelist allows meshes to be partitioned by the mesh generators rather than by
an application at runtime. Principle meshes are partitioned and written to
file as 1 partition per file.

This functionality results in multiple mesh files with the corresponding
portion of each principle mesh. In addition to the mesh tolopologies,
partition information is also written to each file. In essence, each file
provides the information required to load and insantiate 2D-mesh objects which
are `local` to a given process rank. This namelist is enabled if triggered by
the :ref:`partition_mesh<partition_mesh>` logical.

.. NOTE:: Care should be taken when specifying partition
	  configurations. Partitions flag mesh cell ids as being members of
	  that parition, there is no restriction that limits the shape or
	  continuity of a partition.

* ``max_stencil_depth``: **<integer>**
    The maximum stencil depth (in cells) supported by the partitioned
    principle meshes.

.. _n_partitions:

* ``n_partitions``: **<integer>**
    The total number of requested partitions. For `cubed-sphere` meshes, this
    shoud be restricted to 1 or a multiple of 6.

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
      Forces the partitioner to attempt to configure the panel into partitions
      given by :ref:`panel_xproc<panel_xproc>` and
      :ref:`panel_yproc<panel_yproc>`.

.. _panel_xproc:

* ``panel_xproc``: **<integer>**
    Number of partitions in local x-direction of mesh panel, this variable is
    only valid when requesting a ``custom`` decomposition
    (:ref:`panel_xproc<panel_xproc>`).

.. _panel_yproc:

* ``panel_yproc``:  **<integer>**
    Number of partitions in local y-direction of mesh panel, this variable is
    only valid when requesting a ``custom`` decomposition
    (:ref:`panel_yproc<panel_yproc>`).

* ``partition_range``: **<integer>,<integer>**
    Specifies `start`,`end` partition ids to output, valid ids range from
    [0: :ref:`n_partitions<n_partitions>`-1]. The generators will produce 1
    file per requested partition, with the partition id tagged to the output
    filename.
