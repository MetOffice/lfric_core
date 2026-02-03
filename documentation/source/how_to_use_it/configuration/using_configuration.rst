.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _generating_configuration:

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

.. _loading_configuration:

Loading the configuration
=========================

The Configurator generates a procedure (``read_configuration``) to read
a configuration file based on an application's metadata file (``.json``).
The configuration entails one or more Fortran namelists, which
are each read and stored in a namelist specific type
(`<namelist_name>_nml_type`). These namelist objects are in turn stored
in a configuration object (`config_type`).

This allows for multiple configurations to be loaded into a given
application.

.. code-block:: fortran

   use config_loader_mod, only: read_configuration

   type(config_type) :: config_A, config_B
   ...

   call config_A%initialise( 'ConfigurationName_A' )
   call config_B%initialise( 'ConfigurationName_B' )

   call read_configuration( namelist_file_A, config=config_A )
   call read_configuration( namelist_file_B, config=config_B )

The LFRic infrastructure provides a :ref:`driver configuration
component<driver configuration>` that orchestrates both reading of the
namelist configuration file and cross-checking the contents to ensure
all required namelists are present. The driver configuration component
can be used instead of directly calling the above procedure.

.. _using_config_object:

Using the Config Object
==============================

The term "config object" refers to an object of type
``config_type``. It holds a number of extended ``namelist_type``
objects each of which holds the configuration choices for one of the
namelists. To access a configuration value, simply reference
the configuration hirarchy in the ``config_type``, `e.g.`

.. code-block:: fortran

   use config_mod, only: config_type

   type(config_type)  :: config
   character(str_def) :: name

   name = config%base_mesh%mesh_name()

.. _config_enumerations:

Enumerations
------------

An enumeration is a variable that can take one of a small number of
fixed values. In the namelist, the permitted values are strings, but
within the code, the option and each of the permitted values are
converted into integers.

To get, and to use, an enumeration, one has to get the value
representing the choice, but also one or more of the enumeration list
to check against. Enumerations are stored as ``i_def`` integers. The
enumeration options are parameters that can be obtained directly from
Configurator-generated ``_config_mod`` modules.

To illustrate, Rose metadata specifies that the value of the
namelist variable ``geometry`` can be either the string "spherical" or
the string "planar". In the following code, the namelist entry is
checked against two allowed choices of geometry: ``spherical`` and
``planar``, referenced by the two integer parameters in the
``base_mesh_config_mod`` module. The names of the parameters are
prefixed with the name of the variable to ensure there is no
duplication of parameter names with other enumeration variables:

.. code-block:: fortran

   use base_mesh_config_mod, only: geometry_spherical, geometry_planar

   integer(i_def) :: geometry_choice
   real(r_def)    :: domain_bottom

   geometry_choice = config%base_mesh%geometry()

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

   Code that compares integer options and parameters is safer than code
   that compares string options and parameters. If there are spelling
   errors in the names in the code, the former will fail at compile
   time whereas problems with the latter only arise at run-time.

.. _config_duplicate_namelists:

Duplicate namelists
---------------------

When a defined namelist is allowed to have multiple instances in a namelist
input file, the namelist is said to allow `duplicates` (with a given
configuration). This is indicated by the Rose metadata as ``duplicate=true``.
These namelist instances have the same variable names, though those varibles
may contain different values.

Instances of a ``duplicate`` namelist may be differentiated using one of the
namelists members as a key. The :ref:`extended Rose metadata <extended rose
metadata>`, ``!instance_key_member``
indicates to the configurator tool which namelist variable to use as the
instance key. For example, consider a ``partitioning`` namelist, with the
variable, ``mesh_choice`` used as the instance key member, the relevant parts
of the Rose metadata may look as follows::

  [namelist:partitioning]
  duplicate=true
  !instance_key_member=mesh_choice

  [namelist:partitioning=mesh_choice]
  !enumeration=true
  values='source', 'destination'

As the instance key is ``mesh_choice``, there may be two copies of the
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

To extract a *specific instance*, the possible ``mesh_choice`` string must be known:

.. code-block:: fortran

   type(partitioning_nml_iterator_type) :: iter
   type(partitioning_nml_type), pointer :: partitioning_nml

   call iter%initialise(config%partitioning)
   do while (iter%has_next())

     partitioning_nml => iter%next()
     mesh_choice = partitioning_nml%get_profile_name()

     select case (trim(mesh_choice))

     case('source')
         srce_partitioner = partitioning_nml%partitioner()

     case('destination')
         dest_partitioner = partitioning_nml%partitioner()

     end select

   end do

.. _config_mod_files:

Using configuration module files
================================

In the examples above, the ``config_mod`` modules were used **only** to
obtain the parameters that represent the options of an enumeration
variable. The preferred practice is to only use global scope
``config_mod`` modules to access fixed runtime parameters.

.. admonition:: Configuration access via ``use`` statements.

   While existing code will allow access to any variables direct from
   the ``config_mod`` files, this legacy practice is strongly
   **discouraged**. Access via ``use`` statements cannot work where
   namelists may have multiple instances or when an application wishes
   to load multiple configuration files. It is recommended that access
   to configuration namelists be limited to usage of the ``modeldb%config``
   item described above.

   Direct access from ``*_config_mod`` modules will be restricted to
   enumeration parameter values only by April 2026.
