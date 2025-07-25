.. -----------------------------------------------------------------------------
    (c) Crown copyright 2024 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _using configuration:

Configuration code generation
=============================

Applications that use the :ref:`extended Rose metadata <extended rose
metadata>` can run the LFRic :ref:`Configurator <configurator>` as
part of the application build process. The Configurator generates all
the code required to read namelist configuration files. The code also
makes the configuration information available to application coders in
a user-friendly format.

This section describes how to load an application configuration into
the application, and how code can use the various types of application
configuration.

Loading the configuration
=========================

The Configurator generates a procedure ``read_configuration`` to read
a namelist configuration file aligned with the Rose metadata and to
store the configuration choices in a ``namelist_collection_type``

.. code-block:: fortran

  use configuration_mod, only: read_configuration

  ...
  type(namelist_collection_type) :: configuration

  call read_configuration( filename, configuration )

The LFRic infrastructure provides a :ref:`component<driver
configuration>` that orchestrates both reading of the namelist
configuration file and cross-checking the contents to ensure all
required namelists are present, and can be used instead of directly
calling the above procedure.

.. _configuration object:

Using the Configuration Object
==============================

The term "configuration object" refers to an object of type
``namelist_collection_type``. It is designed to hold comprehensive
information about an application configuration.

A configuration object holds a number of ``namelist_type`` objects. To
access a namelist object, call the ``get_namelist`` function:

.. code-block:: fortran

  use namelist_mod, only: namelist_type
  use namelist_collection_mod, only : namelist_collection_type

  type(namelist_collection_type) :: configuration

  type(namelist_type), pointer :: base_mesh_nml

  base_mesh_nml => configuration%get_namelist('base_mesh')

To get a configuration value:

.. code-block:: fortran

  character(str_def) :: mesh_name

  call base_mesh_nml%get_value('mesh_name', mesh_name)

To get, and to use, an enumeration, one has to get the value
representing the choice, but also one or more of the enumeration list
to check against. Enumerations are stored as ``i_def`` integers. The
enumeration options are parameters that can be obtained directly from
the Configurator-generated ``_config_mod`` file:

.. code-block:: fortran

   use base_mesh_config_mod, only: geometry_spherical, geometry_planar

   integer(i_def) :: geometry_choice
   real(r_def)    :: domain_bottom

   base_mesh_nml => configuration%get_namelist('base_mesh')
   call base_mesh_nml%get_value('geometry', geometry_choice)

   select case (geometry_choice)
   case (geometry_planar)
     domain_bottom = 0.0_r_def
   case (geometry_spherical)
     domain_bottom = earth_radius
   case default
     call log_event("Invalid geometry", LOG_LEVEL_ERROR)
   end select

.. admonition:: Hidden values

  Use of enumerations can be better than using numerical options or
  string variables.

  A parameter name is more meaningful and memorable than a numerical
  option, making code more readable. There is also a clearer link
  between the name and the metadata, as the metadata can be easily
  searched to find information about the option.

  Code that compares integer options and parameters is safer than
  code that compares string options and parameters. If there are
  spelling errors in the names, the former will fail at compile time
  whereas problems with the latter only arise at run-time.

.. _config_mod files:

Where namelists are duplicated, the possible values of the instance
variable can be used to distinguish between them. For example, for a
namelist ``partitioning`` with an instance key of ``mesh_choice``,
the relevant parts of the Rose metadata may look as follows::

  [namelist:partitioning]
  duplicate=true
  instance_key_member=mesh_type

  [namelist:partitioning=mesh_type]
  !enumeration=true
  values='source', 'destination'

As the instance key is ``mesh_type``, there may be two copies of the
namelist each with a different ``mesh_choice``::

  &partitioning
  mesh_choice           = 'destination',
  partitioner           = 'cubedsphere',
  panel_decomposition   = 'auto',
  /

  &partitioning
  mesh_choice           = 'source',
  partitioner           = 'planar',
  panel_decomposition   = 'auto',
  /

The different namelist options can be extracted with the following
code (noting that the possible ``mesh_choice`` values must be known in
the code):

.. code-block:: fortran

  ! Get namelist objects for the source and destination partitioning
  source_partitioning_nml =>                             &
              configuration%get_namelist('partitioning', &
                                         'source')
  destination_partitioning_nml =>                        &
              configuration%get_namelist('partitioning', &
                                         'destination')

  ! Extract information from the two different namelist objects
  call source_partitioning_nml%get_value('partitioner',       &
                                          source_partitioner)
  call destination_partitioning_nml%get_value('partitioner',  &
                                          destination_partitioner)


Using config_mod files
======================

In the examples above, the ``config_mod`` files were used only to
obtain the parameters that represent the options of an enumeration
variable.

It is normally possible to obtain any variables direct from the
``config_mod`` files rather than going through the configuration
object functions. However, this method cannot work where namelists are
duplicated. Namelists can be duplicated by metadata definition as
described above, but applications can also be hard-wired to read two
separate namelist configurations where the same namelist appears in
both. In these cases, the data in the configuration object will be
specific for each part of the application, but the values for the first
namelist in the ``config_mod`` file will be overwritten by the second
namelist.

In the following example, the same requirement as the example above is
met by directly using the value of the ``geometry`` option from the
``config_mod`` file:

.. code-block:: fortran

   use base_mesh_config_mod, only: geometry_spherical, &
                                   geometry_planar,    &
                                   geometry

   real(r_def)    :: domain_bottom

   select case (geometry)
   case (geometry_planar)
     domain_bottom = 0.0_r_def
   case (geometry_spherical)
     domain_bottom = earth_radius
   case default
     call log_event("Invalid geometry", LOG_LEVEL_ERROR)
   end select
