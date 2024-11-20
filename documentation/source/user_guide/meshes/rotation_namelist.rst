.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section rotate_mesh_nml:

.. _rotation_nml:
============================
``&rotation``
============================
Optional control namelist for rotation of node coordinates. Rotation in
[longitude, latitude] is specified by using a reference location, which is
rotated such that it arrives at the specified target location. All nodes then
undergo the same coordinate transformation. Enabled by
:ref:`mesh:rotate_mesh<rotate_mesh>` (cartesian planar meshes not supported).

* ``rotation_target``: **'<string>'**
    Feature to use as reference for rotation.

    ``north_pole``
      Use North Pole (90N,0E) as reference location.

    ``null_island``
      Use Null Island (0N,0E) as reference location.

* ``target_north_pole``: **<real>,<real>**
    Target location of North Pole after rotation.
* ``target_null_island``: **<real>,<real>**
    Target location of Null Island after rotation.
