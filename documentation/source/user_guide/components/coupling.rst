.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section coupling component:

Coupling Component
==================

The coupling component provides an interface to the external library that
provides coupling functionality. At the moment the coupling functionality is
provided by the OASIS3-MCT library.

.. attention::

   All coupling should be performed through the coupling component, we
   should no longer access Oasis functionality directly from elsewhere in the 
   model. Outside of the coupling component, there should be no lines
   containing ``use mod_oasis``

The coupling object
-------------------

This object holds all the information related to the coupling.
It is initialised with 

.. code-block:: fortran

   call coupling%initialise(cpl_name, comm_out, comm_in, comm_is_split)

where:

* ``cpl_name`` is the name that will be given to the coupling component (this
  is the name that the component will be referred to as in the coupling
  configuration file.
* ``comm_in`` is the MPI communicator that Oasis will split
* ``comm_out`` is the MPI communicator that Oasis will return for you to run
  your model in.
* ``comm_is_split`` is returned as true if Oasis has split the communicator.

Once it is initialised, the coupling object need to be populated with
information on:

1. The way the data is held in memory and how it is distributed over the
   parallel tasks. This is the "partition" definition used by the coupling.
2. The different fields of data that have to be coupled. These are called
   the "variable" definitions.

The routines used to set up this information are called on the coupling object
and have the following interface:

.. code-block:: fortran

  call coupling%define_partitions( mpi, mesh )

where:

* ``mpi`` is the model's MPI object - which includes the MPI communicatior,
* ``mesh`` is the mesh object for the fields being exchanged. Currently only 2d
  meshes are supported.

.. code-block:: fortran

  call coupling%define_variables( cpl_snd_2d, cpl_rcv_2d, cpl_snd_0d )

where:

* ``cpl_snd_2d`` is a field collection that holds all the 2d fields that will be
  sent to the coupler
* ``cpl_rcv_2d`` is a field collection that holds all the 2d fields that will be
  received from the coupler
* ``cpl_snd_0d`` is a collection of fields that will have some form of reduction
  operation performed on them, so the value passed through the coupler will be
  a (0d) scalar.

.. note::

   Clearly, this is not a complete interface. It is the interface that was
   required to support the coupling required by the ``lfric_atm`` application.
   Further types of field can easily be supported - they just need adding to the
   code-base 

When the coupling of a field is instigated, an array that specifies the order
of the data to be coupled is required. This can be extracted from the coupling
object:

.. code-block:: fortran

  local_index = coupling%get_local_index()

where:

* ``local_index`` is returned containing the index to sort data for sending or
  receiving.

Coupling exchange objects
-------------------------

In order to send data through the coupler, it has to be extracted from within a
field object so it can be passed to the coupler. This will break the
encapsulation enforced by the field object. The agreed method for doing this is
to use the "external field" mechanism. We create an object called
``coupler_exchange_2d_mod`` which is inherited from an external field. This is
initialised with:

.. code-block:: fortran

  call coupler_exchange_2d_mod%initialise(lfric_field_ptr, sorting_index)

where:

* ``lfric_field_ptr`` is a pointer to the lfric field that contains the data to
  be coupled.
* ``sorting_index`` is an array of indices that determines the order the data
  will be sent to the coupler. This would usually be extracted from the
  partition information in the coupling object.

Data can then be sent to or received from the coupler with calls to:

.. code-block:: fortran

  call coupler_exchange_2d_mod%copy_from_lfric(return_code)

where:

* ``return_code`` is an optional argument that can use used to access the return
  code from the coupling.

.. code-block:: fortran

  call coupler_exchange_2d_mod%copy_to_lfric(return_code)

where:

* ``return_code`` is an optional argument that can use used to access the return
  code from the coupling.

Obviously, sending a scalar to the coupler is much simpler. There is no
encapsulation to break, so it is simply a case of calling:

.. code-block:: fortran

  call coupler_send_0d( scalar, name, var_id, model_clock )

where:

* ``scalar`` should contain the scalar value to be passed to the coupler
* ``name`` is the name to associated with the field used to create the scalar
* ``var_id`` is the coupling identifier that is stored with (and can be
  retrieved from the field (e.g. ``var_id = field%get_cpl_id(1)``)
* ``model_clock`` contains the time within the model.

The scalar value should be generated by performing a reduction operation on a
field and it is this field that is held in the collection of fields used for
(0d) scalar coupling. The receive counterpart for coupling a scalar has not been
implemented as it doesn't make a great deal of sense: it is not clear how to
generate a full field from the reduced scalar passed from the coupler. 
