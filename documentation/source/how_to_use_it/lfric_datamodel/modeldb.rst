.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _modeldb:

Modeldb
=======

The modeldb class is a data structure that is designed to encapsulate
all the data required to describe both the scientific and technical
state of a model.

The modeldb class can hold:

* Model field data such as all the standard field collections
  (e.g. the depository, prognostics and diagnostics) along with all
  the other model specific field collections and field bundles.
* Other model data such as single values, arrays or objects
* Configuration that is read in from namelists
* Clock/calendar objects
* MPI communicator
* I/O information

There is a requirement that the LFRic code be able to run more than
one instance of a "model" from the same executable. For example, you
might want to run a linear model and its adjoint in the same
executable. Or you might want to run a number of members of an
ensemble simultaneously from the same executable. In order to make
this possible, multiple versions of the model state can be stored in
multiple copies of modeldb.

How to use
----------

The modeldb class contains a collection of field collections. Fields
can be stored in field collections, then these field collections
stored in modeldb.

Other values, such as single values, arrays or even objects can be
stored in a collection of key-value pairs that is stored with modeldb

Items common to all models like clock and calendar and MPI information
can be accessed directly from modeldb

Fields
------

The modeldb object contains an item called ``fields``. This is just a
collection (actually a linked list) of field collections. These field
collections hold fields in the same way as any other field collection.

To put a new collection into "fields"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To add a collection to ``modeldb%fields`` use:

.. code-block:: fortran

   call modeldb%fields%add_empty_field_collection("my_collection", &
                                                  table_len = 100)

where ``"my_collection"`` is the name of the field collection you want
adding and the ``table_len`` is the length of the hash table that is
used to hold the fields in the collection. Small collections will be
more efficient with a small ``table_len``, but collections with many
fields will be more efficient with a larger ``table_len``.

Unlike other modeldb items, the ``fields`` item is initialised within the
call to ``add_empty_field_collection``.

To put a field into one of the collections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

   type( field_collection_type ), pointer :: my_collection
   type( field_type ),                    :: my_field

   my_collection => modeldb%fields%get_field_collection("my_collection")
   call my_collection%add_field(my_field)

This will put a copy of ``my_field`` into the collection. If you want
to use the version held in the collection, you will need to retrieve a
pointer to it from the collection.

To get a field out
^^^^^^^^^^^^^^^^^^

Assuming the field, ``my_field``, has the name "my_field", use:

.. code-block:: fortran

   type( field_collection_type ), pointer :: my_collection
   type( field_type ),            pointer :: my_field

   my_collection => modeldb%fields%get_field_collection("my_collection")
   call my_collection%get_field("my_field", my_field)

This returns a pointer to the actual field held in the collection. Any
changes to the field you have extracted will instantly change the
version in modeldb. They both refer to the same location in memory.

Values
------

The modeldb object contains an item called ``values``. This is a
collection of key-value pairs. Many "values" can be stored in this
collection and are accessed via their key (or name). The "values" that
can be stored in a key-value pair must be one of:

* Scalar (or arrays of) 32-bit real value(s) (``real(real32) ::``)
* Scalar (or arrays of) 64-bit real value(s) (``real(real64) ::``)
* Scalar (or arrays of) 32-bit integer value(s) (``integer(int32) ::``)
* Scalar (or arrays of) 64-bit integer value(s) (``integer(int64) ::``)
* Scalar (or arrays of) logical value(s) (``logical ::``)
* A text string or an array of text strings (``character(*) ::``)
* Any Fortran object that inherits from ``abstract_value_type`` from
  the module ``key_value_mod`` (``type, extends(abstract_value_type)
  ::``)

(the definitions used here for the variable "kinds" are from the
intrinsic module ``iso_fortran_env``)

The final item on the list ("any Fortran object") makes the collection
of values very powerful. Anything from the object that encapsulates a
time-stepping algorith to the object that stores all the
atmosphere-ocean coupling information can be stored in modeldb.

Prior to using the ``values`` collection, it must be have been initialised
as shown in the code below.

To put a value in
^^^^^^^^^^^^^^^^^

.. code-block:: fortran

   real(real64) :: my_value

   my_value = 7.0_real64
   call modeldb%values%initialise()
   call modeldb%values%add_key_value('my_value', my_value)

Again, this will put a copy of the value into the collection

To get a value out
^^^^^^^^^^^^^^^^^^

.. code-block:: fortran

   real(real64), pointer :: my_value

   call modeldb%values%get_value("my_value", my_value)

This returns a pointer to the value held in the collection. Any
subsequent maths performed on what is returned (the pointer) will
change the value held in the collection. They both refer to the same
location in memory.

Configuration
-------------

The ``config`` item (`config_type`) within the ``modeldb`` object stores
the input namelists used to configure an instance of modeldb. Once the
``config`` item has been populated the configuration values are immutable,
unlike other components of modeldb.

.. Should provide a link to the namelist collection type (when it's written)

Initialising the configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The ``config`` item is populated using a module generated by the
:ref:`Configurator<configurator>` tool. A namelist input file is simply read
in and any valid namelists are added to the ``config`` item.

As with the ``values`` item, the ``config`` item must be initialised
prior to its first use.

.. code-block:: fortran

   use config_loader_mod, only: read_configuration

   call modeldb%config%initialise()
   call read_configuration( filename, config=modeldb%config )

.. _access_config_data:

Accessing configuration data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To access configuration data from the ``config`` item, simply
reference the namelist member via its location in the configuration
hierachy.

.. code-block:: fortran

   MemberValue = modeldb%config%<NamelistName>%<MemberName>()

The access pattern for namelists which allow multiple instances
(with metadata ``duplicate=true``) is via an iterator which cycles
through instances of the namelist type.

.. code-block:: fortran
   :force:

   type( <NamelistName>_nml_iterator_type ) :: iter
   type( <NamelistName>_nml_type ), pointer :: config_nml

   call iter%initialise( modeldb%config%<NamelistName> )
   do while ( iter%has_next() )
     config_nml => iter%next()
     MemberValue = config_nml%<MemberName>()
   end do

I/O contexts
------------

An I/O context is used to describe how, when and where data are read
from or written to disk. Different groups of data can be read/written
in different circumstances, so there is a requirement to hold a number
of I/O contexts. The modeldb object contains an item called
``io_contexts`` for this purpose. It is simply a collection of the
different io_contexts that are required.

As with the ``configuration`` item, the ``io_context`` item must be
initialised prior to its first use.

To put an I/O context into the collection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: fortran

   type( lfric_xios_context_type ) :: my_io_context

   call modeldb%io_context%initialise()
   call modeldb%io_contexts%add_context(my_io_context)

This will put a copy of ``io_context`` into the collection. If you
want to use the version held in the collection, you will need to
retrieve a pointer to it from the collection.

To get an I/O context out of the collection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Assuming the context, ``my_io_context``, has the name "my_io_context", use:

.. code-block:: fortran

   type( lfric_xios_context_type ) :: my_io_context

   call modeldb%io_contexts%get_io_context("my_io_context", my_io_context)

This returns a pointer to the actual I/O context held in the collection.
