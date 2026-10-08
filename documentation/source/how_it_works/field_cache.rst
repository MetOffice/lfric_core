.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
   Some of the content of this file has been produced with the assistance of
   Met Office Github Copilot Enterprise.

.. _how_it_works_field_cache:

Field Caching Mechanism
=======================

This section describes the internal design and operation of the field
caching mechanism provided by the ``field_cache_alg_mod`` module.

Architecture Overview
---------------------

The ``field_cache_type`` is a derived type that manages the caching and
restoration of fields. It consists of three key components:

1. **Source Field Collection Pointer** - A non-owning reference to the
   original field collection
2. **Cached Field Collection** - An internal field collection that stores
   copies of cached fields
3. **Field Names List** - An array of field names specifying which fields
   to cache

The cache operates independently of the source collection; modifications to
cached fields do not affect the source collection until an explicit reset is
performed.

initialisation Phase
--------------------

When ``initialise`` is called, the following operations occur:

1. The source field collection pointer is stored (the source collection is
   not copied)
2. The field names are stored for later use
3. An empty ``cached_fields`` collection is initialised with a hash table
   size of 100
4. Each specified field name is validated against the source collection

**Validation**: If a requested field does not exist in the source collection,
an error is logged. All fields to be cached must be present in the source
collection at the time of initialisation.

Caching Phase
-------------

The ``cache_fields`` method performs deep copies of specified fields into the
cache:

**Main Loop**

For each field name in the cache specification:

1. Retrieve the field from the source collection using
   ``get_abstract_field``
2. Perform type dispatch using Fortran's ``select type`` construct
3. Handle the field based on its specific type

**Type-Specific Handling**

Real-Valued Fields (``field_type``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each real-valued field:

1. Check if the field already exists in the cache collection
2. If not, create a cached field by copying field properties from the source
   field
3. Add the new cached field to the cache collection
4. Retrieve the cached field reference
5. Copy current data values from the source field to the cached field using
   PSyclone with ``invoke(setval_x(...))``

The ``copy_field_properties`` method copies metadata about the field (such as
mesh, function space, and other properties) but not the data values.

Integer-Valued Fields (``integer_field_type``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Integer fields follow the same process as real-valued fields but use:

* ``int_setval_x`` instead of ``setval_x`` for copying data values

Array of Fields (``field_array_type``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For field arrays:

1. Check if the field array already exists in the cache
2. If not, create a new field array with the same bundle size as the source
3. Add the field array to the cache collection
4. For each field in the array bundle:

   a. Get references to both the source and cached field
   b. Copy field properties if the cached field is not yet initialised
   c. Copy data values from source to cached field

**Unsupported Types**

If a field type is not recognized, an error is logged.

Reset Phase
-----------

The ``reset_fields`` method reverses the caching operation:

**Pre-Reset Validation**

First, the method checks if there are any cached fields. If the cache is
empty, an error is logged.

**Main Loop**

For each field name in the cache specification:

1. Retrieve the current field from the source collection
2. Perform type dispatch
3. Handle each type by copying cached values back to the source

**Finalization**

The ``finalise_field_cache`` method (called automatically when the cache object
goes out of scope) nullifies the pointer to the source field collection to avoid
dangling references.

Data Copy Mechanism
-------------------

The actual data copying is performed using PSyclone builtins (via
``invoke``):

* ``setval_x`` - Copies real field values
* ``int_setval_x`` - Copies integer field values

These builtins handle:

* Distributed memory operations (data synchronization across MPI ranks)
* OpenMP parallelisation for shared memory performance

This approach ensures that caching operations respect the parallel structure of
the fields and properly synchronise data across multiple processes.

Design Considerations
---------------------

**Non-Owning Source Reference**

The cache maintains a non-owning pointer to the source collection. This design
allows the cache to remain valid as long as the source collection exists, but
avoids circular references and simplifies memory management. The caller is
responsible for ensuring the source collection remains valid during the cache's
lifetime.

**Lazy Field initialisation in Cache**

Fields are only created in the cache collection when first cached. This
allows:

* Multiple calls to ``cache_fields`` to update cached values
* Efficient memory usage (only cached fields consume memory)

**Incremental Caching**

Each call to ``cache_fields`` updates the cached values. Calling it multiple
times overwrites previous cached data.

Performance Characteristics
---------------------------

**Time Complexity**

* initialisation: O(n) where n is the number of fields to cache
* Caching: O(n × m) where m is the average field size (linear in field data)
* Reset: O(n × m)

**Space Complexity**

The cache requires O(n × m) additional memory to store field copies.

**Parallel Efficiency**

The underlying ``invoke`` builtins handle all MPI communication and OpenMP
parallelisation, ensuring that caching and reset operations scale efficiently
in parallel environments.

Limitations and Future Extensions
----------------------------------

**Current Limitations**

* Only three field types are supported (``field_type``,
  ``integer_field_type``, ``field_array_type``)
* No selective field subsetting (must cache entire fields)
* No time series capability (only one snapshot at a time)

**Potential Future Extensions**

* Support for ``real32_field_type`` and ``real64_field_type`` (when PSyclone
  provides kernels)
* Multiple named snapshots for comparison
* Partial field caching (specific subsets of mesh regions)
