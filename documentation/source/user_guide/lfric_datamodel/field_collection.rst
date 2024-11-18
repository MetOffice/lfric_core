.. ------------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _section field collection:

Field collections
=================

A field collection is an object that can store several fields or
pointers to fields.

A field collection is initialised as follows:

.. code-block:: fortran

   type(field_collection_type)           :: my_field_collection
   character(len=*),                     :: my_collection_name

   call my_field_collection%initialise(name = my_collection_name, &
                                       table_len = table_len)

Both arguments are optional. The field collection name is required
only if the field collection is to be added to another object such as
a list of field collections.

The table length defines the length of a hash table that stores linked
lists of fields. A number roughly the order of the number of fields in
the table should be acceptable. The table length defaults to 100 which
should be acceptable even for small field collections.

The following illustrates two instantiated fields being added to a
collection. One of the fields is a pointer:

.. code-block:: fortran

   type(field_type), intent(in)               :: my_field
   type(field_type), intent(in),pointer       :: my_ptr_field

   call my_field_collection%add_field(my_field)
   call my_field_collection%add_reference_to_field(my_ptr_field)

The fields must have a name. Names must not be a duplicate of a name
of a field already in the collection.

Fields can be removed by name:

.. code-block:: fortran

   call my_field_collection%remove_field("my_field_name")

To check whether a field of a given name exists in a collection call
the following logical function:

.. code-block:: fortran

   if ( my_field_collection%field_exists("my_field_name") ) then
     ...

The following gets a pointer to a field in a collection. Ensure the
kind of the field pointer passed to ``get_field`` matches the kind of
the field being requested.

.. code-block:: fortran

   call my_field_collection%get_field("my_field_name", my_field)

The field collection object stores fields in a linked list. A separate
:ref:`field collection iterator<section field collection iterator>`
object uses the ``get_next_item`` function of the field collection to
step through field collections. Most likely, the ``get_next_item``
would not be used in isolation, but is described here for
completeness.

The following function gets the first and second linked list item in
the list (assuming there are at least two items):

.. code-block:: fortran

   type(linked_list_item_type), pointer :: start, first

   ! A null input moves to the start of the list
   nullify(start)

   first = my_field_collection%get_next_item(start)
   second = my_field_collection%get_next_item(first)

The function returns a linked list item with a payload that can be any
kind of field. Therefore, a ``select type`` construct is required to
disambiguate all the different possible field types.

Three self-explanatory functions return information about the field
collection: ``get_length`` returns the number of fields in the
collection; ``get_table_len`` returns the length of the internal hash
table; ``get_name`` returns the field collection name.

Finally, the ``clear`` function removes all items from the field
collection.

.. _section field collection iterator:

Field collection iterator
-------------------------

The field collection iterator supports the ability to loop through all
the fields in a field collection. For example:

.. code-block:: fortran

   type(field_collection_type), intent(inout) :: state

   type( field_collection_iterator_type) :: iter
   class( field_parent_type ), pointer :: field

   call iter%initialise(state)

   do

   if ( .not.iter%has_next() ) exit
     field => iter%next()

     ! Select type to disambiguate different field types
     select type(field)
       type is (field_real32_type)
         ! Do something
         ...
