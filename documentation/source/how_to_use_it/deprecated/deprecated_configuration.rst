.. _DeprecatedConfiguration:

Modeldb Configuration Item
==============================================

The ``configuration`` item (``namelist_collection_type``) within the ``modeldb``
object stores the input namelists used to configure an instance of modeldb.
Once the configuration has been populated, the configuration values are
immutable, unlike other components of ``modeldb``. This item is **deprecated**
and use of the ``config`` item in ``modedb`` is the preferred configuration
access method. The component ``modeldb%configuration`` is marked for removal
after April 2026.

.. Should provide a link to the namelist collection type (when it's written) #PR206

Initialising the configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The ``modeldb%configuration`` item is populated using a module generated
by the ``configurator`` tool. A namelist input file is simply read in and
any valid namelists are added the ``modeldb%configuration`` item.

The ``configuration`` item is first initialised before reading a namelist input
file.

.. code-block:: fortran

   use configuration_mod, only: read_configuration

   call modeldb%configuration%initialise()
   call read_configuration( filename, configuration=modeldb%configuration )

.. _access_configuration_data:

Accessing configuration data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To access configuration data from the ``configuration`` item, a pointer to
the required namelist is first requested before the namelist variable
member can be retrieved. This allows for different namelists to have
member variables with the same name.

.. code-block:: fortran

   type(namelist_type), pointer :: config_nml

   config_nml => modeldb%configuration%get_namelist('<namelist_name>')
   call config_nml%get_value( '<namelist_variable_name>', nml_var )

All namelists in the collection must have a unique `<namelist_name>`, unless
the namelist metadata specfies ``duplicate=.true.``. For namelists which may
appear multiple times in a namelist file, the ``profile_name`` must also be
specified.

.. code-block:: fortran

   type(namelist_type), pointer :: config_nml

   config_nml => modeldb%configuration%get_namelist( '<namelist_name>', &
                                                     profile_name='<profile_name>' )
   call config_nml%get_value( '<namelist_variable_name>', nml_var )
