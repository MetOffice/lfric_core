.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
   Some of the content of this file has been produced with the assistance of
   Met Office Github Copilot Enterprise.

.. _field_cache:

Field Caching
=============

The field caching mechanism provides a way to save and restore the state of
fields in a field collection. This is particularly useful in scenarios where
you need to temporarily modify fields during computation and later restore them
to their original state, such as during iterative algorithms, testing, or
recovery procedures.

Overview
--------

The ``field_cache_alg_mod`` module provides the ``field_cache_type`` derived
type, which manages a cache of field copies from a
``field_collection_type``. It supports caching and restoration of three field
types:

* ``field_type`` - Real-valued scalar fields
* ``integer_field_type`` - Integer-valued scalar fields
* ``field_array_type`` - Arrays of fields

Basic Workflow
--------------

The typical workflow for using the field cache involves three steps:

1. **initialise** the cache with a source field collection and the names of
   fields to cache
2. **Cache the fields** when you want to save their current state
3. **Reset the fields** to restore them to their cached state

Creating and initialising a Cache
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, declare and initialise a ``field_cache_type`` object:

.. code-block:: fortran

    use field_cache_alg_mod, only : field_cache_type
    use field_collection_mod, only : field_collection_type

    type(field_cache_type) :: my_cache
    type(field_collection_type), target :: fields
    character(len=:), allocatable :: field_names_to_cache(:)

    ! Define which fields to cache
    allocate(character(len=10) :: field_names_to_cache(3))
    field_names_to_cache(1) = "temperature"
    field_names_to_cache(2) = "pressure"
    field_names_to_cache(3) = "humidity"

    ! initialise the cache
    call my_cache%initialise(fields, field_names_to_cache)

The ``initialise`` method takes:

* ``source_field_collection`` - The target field collection that contains the
  fields to cache
* ``field_names`` - A list of field names to cache

All specified fields must exist in the source collection, otherwise an error is
logged.

Saving Field State
~~~~~~~~~~~~~~~~~~

Once the cache is initialised, you can save the current state of the fields:

.. code-block:: fortran

    ! Save the current state of all specified fields
    call my_cache%cache_fields()

This operation copies the current values of all specified fields into the
internal cached field collection. You can call this multiple times to update
the cached state.

Restoring Field State
~~~~~~~~~~~~~~~~~~~~~

When you need to restore the fields to a previously cached state:

.. code-block:: fortran

    ! Restore fields to the last cached state
    call my_cache%reset_fields()

This operation copies the cached values back to the fields in the source
collection.

Practical Example
-----------------

Here's a complete example demonstrating the field cache usage:

.. code-block:: fortran

    program field_cache_example
      use field_cache_alg_mod, only : field_cache_type
      use field_collection_mod, only : field_collection_type
      implicit none

      type(field_cache_type) :: cache
      type(field_collection_type), target :: field_collection
      character(len=:), allocatable :: fields_to_cache(:)

      ! Assume field_collection has been populated with fields
      ! initialise the field collection here...

      ! Specify fields to cache
      allocate(character(len=10) :: fields_to_cache(2))
      fields_to_cache(1) = "u_wind"
      fields_to_cache(2) = "v_wind"

      ! initialise cache
      call cache%initialise(field_collection, fields_to_cache)

      ! Cache the initial state
      call cache%cache_fields()

      ! Perform some computations that modify the fields
      ! (e.g., run a dynamical core step)

      ! Later, restore the fields to their cached state if needed
      call cache%reset_fields()

    end program field_cache_example

Important Considerations
------------------------

**Field Existence Validation**

During initialisation, the cache validates that all specified fields exist in
the source collection. If a field does not exist, an error message is logged.

**Memory Usage**

The cache maintains a complete copy of each cached field. Be mindful of memory
constraints when caching large fields or many fields.

**Thread Safety**

The field cache is not thread-safe. Do not use a single ``field_cache_type``
instance from multiple threads concurrently.

**Field Collection Lifetime**

The cache maintains a pointer to the source field collection. Ensure that the
source collection remains valid for the entire lifetime of the cache.

**Supported Field Types**

Only the following field types are supported:

* ``field_type`` (real-valued fields)
* ``integer_field_type`` (integer-valued fields)
* ``field_array_type`` (arrays of fields)

Attempting to cache other field types will result in an error message.
