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

.. NOTE:: When considering meshes with rotation applied, all (lon,lat)
	  coordinates which are referenced in the configuration/output files
	  are with respect to a real world frame of reference.
