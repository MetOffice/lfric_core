.. -----------------------------------------------------------------------------
     (c) Crown copyright 2023 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

*******
Modeldb
*******

Introduction
============

It is required that the LFRic code be able to run more than one instance of a "model" from the same executable. For example you might want to run a linear model and its adjoint in the same executable. Or you might want to run a number of members of an ensemble simultaneously from the same executable.

In order to make this possible, multiple versions of the scientific and technical state of the model need to be held. This data will include:

* Model field data such as all the standard field collections (e.g. the depository, prognostics and diagnostics) along with all the other model specific field collections and field bundles.
* Other model data such as single values or arrays
* Configuration that is read in from namelists
* Clock/calendar objects
* MPI communicator
* I/O information

This data is held in a structure called "modeldb"

How to use
==========

Items like clock and calendar and the object holding MPI information can be accessed directly from modeldb: ``modeldb%clock``, ``modeldb%calendar`` and ``modeldb%mpi``, but the modeldb object can also hold fields and other values. 

Fields
~~~~~~

The modeldb object contains an item called ``fields``. This is just a collection (actually a linked list) of field collections. These field collections hold fields in the same way as any other field collection.

To put a new collection into "fields"
-------------------------------------

To add a collection to ``modeldb%fields`` use:
::

  call modeldb%fields%add_empty_field_collection("my_collection", table_len = 100)

where ``"my_collection"`` is the name of the field collection you want adding and the ``table_len`` is the length of the hash table that is used to hold the fields in the collection. Small collections will be more efficient with a small ``table_len``, but collections with many fields will be more efficient with a larger ``table_len``.

To put a field into one of the collections
------------------------------------------
::

  type( field_collection_type ), pointer :: my_collection
  type( field_type ),                    :: my_field
  
  my_collection => modeldb%fields%get_field_collection("my_collection")
  call my_collection%add_field(my_field)

This will put a copy of ``my_field`` into the collection. If you want to use the version held in the collection, you will need to retrieve a pointer to it from the collection.

To get a field out
------------------

Assuming the field, ``my_field``, has the name "my_field", use:
::

  type( field_collection_type ), pointer :: my_collection
  type( field_type ),            pointer :: my_field
  
  my_collection => modeldb%fields%get_field_collection("my_collection")
  call my_collection%get_field("my_field", my_field)

This returns a pointer to the actual field held in the collection. Any changes to the field you have extracted will instantly change the version in modeldb. They both refer to the same location in memory.

Values
~~~~~~

The modeldb object contains an item called ``values``. This is a collection of key-value pairs. Many "values" can be stored in this collection and are accessed via their key (or name). The "values" that can be stored in a key-value pair must be one of:

* Scalar (or arrays of) 32-bit real value(s) (``real(real32) ::``)
* Scalar (or arrays of) 64-bit real value(s) (``real(real64) ::``)
* Scalar (or arrays of) 32-bit integer value(s) (``integer(int32) ::``)
* Scalar (or arrays of) 64-bit integer value(s) (``integer(int64) ::``)
* Scalar (or arrays of) logical value(s) (``logical ::``)
* A text string or an array of text strings (``character(*) ::``)
* Any Fortran object that inherits from ``abstract_value_type`` from the module ``key_value_mod`` (``type, extends(abstract_value_type) ::``)

(the definitions used here for the variable "kinds" are from the intrinsic module ``iso_fortran_env``)

To put a value in
-----------------
::

  real(real64) :: my_value

  my_value = 7.0_real64
  call modeldb%values%add_key_value('my_value', my_value)

Again, this will put a copy of the value into the collection

To get a value out
------------------
::

  real(real64), pointer :: my_value

  call modeldb%values%get_value("my_value", my_value)

This returns a pointer the the value held in the collection. Any subsequent maths performed on what is returned (the pointer) will change the value held in the collection. They both refer to the same location in memory.

Configuration
~~~~~~~~~~~~~

See configuration documentation???

I/O contexts
~~~~~~~~~~~~

An I/O context is used to describe how and when data are read from or written to disk. Different groups of data can be read/written in different ways, so there is a requirement to hold a number of I/O contexts. The modeldb object contains an item called ``io_contexts`` for this purpose. It is simply a collection of io_contexts.

To put an I/O context into the collection
-----------------------------------------
::

  type( lfric_xios_context_type )        :: my_io_context
  
  call modeldb%io_contexts%add_context(my_io_context)

This will put a copy of ``io_context`` into the collection. If you want to use the version held in the collection, you will need to retrieve a pointer to it from the collection.

To get an I/O context out of the collection
-------------------------------------------

Assuming the context, ``my_io_context``, has the name "my_io_context", use:
::

  type( lfric_xios_context_type )        :: my_io_context
  
  call modeldb%io_contexts%get_io_context("my_io_context", my_io_context)

This returns a pointer to the actual I/O context held in the collection.
