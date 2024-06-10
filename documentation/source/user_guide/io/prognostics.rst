.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section prognostics:

Prognostic Fields
=================

A typical earth system model will run several different science
components in each time-step, passing fields between different science
components as it runs. Some fields will be passed from the end of one
time-step to the start of the next. The LFRic infrastructure provides
mechanisms to allocate all these fields in one place, to ensure each
can be passed to the right science components and, if required, be
written to and read from checkpoint dumps. As these fields exist in
the model state throughout the run they are commonly referred to as
prognostic fields, though the definition here may differ from other
definitions of prognostic fields that refer more strictly to the key
physical quantities describing the system being modelled.

Where supported, the initialisation of a prognostic field in an LFRic
application can be integrated with the set-up of the IO subsystem such
that there is a single source of truth that describes both the
representation of the field in the model and its representation in
model output files. For example, the Momentum\ :sup:`®` atmosphere model uses
the field definition in the XIOS ``iodef.xml`` file to determine the
type of model prognostic field to create.

.. topic:: Illustrative examples

   #. Fields such as wind, potential temperature, density and moisture
      are physical quantities at the core of the numerical simulation
      of an atmosphere model. Such fields need to be initialised at
      the beginning of a model run, exist until the run completes, and
      be written to any checkpoint files.
   #. Some fields may be initialised from file data. It is often
      convenient to read files at the start of a run to keep the IO
      code separate from the scientific code. Hence, like prognostics,
      the fields need to be initialised at the beginning of the run. But
      they do not need to be written to a checkpoint file as the
      original file can be read again when the model restarts.
   #. Where a field is passed from one section of science to another
      (such as from radiation to boundary layer) it can be convenient
      to treat it as a prognostic field: to allocate it alongside the
      above prognostic fields at the start of the run rather than
      within the time-step code. If the field is to be passed across a
      time-step boundary then it will also need to be written to a
      checkpoint file.


Overview
--------

It is recommended that LFRic models are implemented with an initialise
stage, a run stage (which executes a single time-step) and a finalise
stage, using the model database structure referred to as ``modeldb``
to hold the model state and other data.

For models designed this way, prognostic fields can be initialised
during the initialise stage and stored in ``modeldb``. As well as
storing individual fields, ``modeldb`` can also store field
collections. A field collection is a Fortran derived type that can
store several fields or pointers to fields. Use of field collections
permits large numbers of fields to be passed down a calling tree
within a single argument. Fields can be added to a field collection
dynamically, meaning the contents of the field collection can depend
on the run-time configuration of the model.

During the initialisation stage, a field collection can be created
that includes references to all fields that need to be written to a
checkpoint dump (noting that a pointer to a field can be included in
more than one field collection). Maintaining such a list simplifies
the process for writing all such fields out to a checkpoint dump.

In summary, in the initialisation phase, for each field, the following
needs to be determined:

 #. Is the field required as a prognostic? If so, it needs to be
    initialised and passed to the science that uses or updates it; for
    example, by adding the field to the appropriate field collection.
 #. Does the field need to be checkpointed? If so, the field can be
    added to a field collection used by the checkpoint/restart system.
 #. What form of field should be initialised? How many levels? Which
    function space? Is it a multidata field? All these choices are
    dependent on a combination of the model configuration and the type
    of field.

Each of the above processes can be written as hard-wired
code. Alternatively, some aspects can be based on sources of metadata
so as to automate the process and to reduce the risk of error. For
example, the Momentum\ :sup:`®` atmosphere makes use of metadata provided as
part of its integration to the XIOS library to determin the form of
field to be initialised.

During the time-step running stage of the model, field collections and
fields set up at initialisation time can be extracted from the
``modeldb``. The following code extracts a field from a field
collection held by ``modeldb`` using its text label.

::

    ! Get the main data structure from the model database
    model_data => modeldb%model_data
    ! Get a pointer to the radiation field collection
    radiation_fields => model_data%radiation_fields
    ! Get a pointer to the long-wave downward radiation at the surface
    call radiation_fields%get_field('lw_down_surf', lw_down_surf)

Once all the model steps are completed, a checkpoint dump may need to
be written to enable the simulation to be extended in a separate run
of the application. The field_collection object supports an iterator
that allows a function to loop through all the fields in a field
collection. Assuming the initialisation phase has created a field
collection containing all the fields that need to be checkpointed, a
function can use the iterator to loop through the field collection's
contents, writing each field out to the checkpoint dump.

Example: The Momentum\ :sup:`®` Model
-------------------------------------

In the Momentum\ :sup:`®` atmosphere model, the set-up of the
prognostics field is closely aligned with the model's integration with
the XIOS library. All fields that the model may read or write,
including all prognostic fields, are registered in an ``iodef.xml``
file which is the file recognised by the XIOS IO library. The set-up
of prognostics fields uses some of the information in the
``iodef.xml`` file to determine how to initialise each field.

The record in the ``iodef.xml`` file includes a string ID for the
field, the formal name of the field (such as the CF name), the units
and information about the XIOS domain used to describe the format of
the data in the input or output file.

In the Momentum\ :sup:`®` atmosphere model, a routine called
`create_physics_prognostics` includes the code for creating most of
the prognostic fields, including dealing with the steps described in
the preceding section.

To initialise a prognostic field, a function is called that creates
the prognostic field based on a combination of input arguments and on
information obtained from XIOS.

The input arguments include the XIOS string ID, an optional flag to
indicate whether a field is to be checkpointed (determined by model
logic) and a pointer to the field collection that will transport the
field through the model call tree.

For example, the following will add two fields to the radiation field
collection. If certain logical conditions are met, the field will also
be marked for checkpointing: a pointer to the field will also be added
to the dedicated field collection used by the checkpoint routine.

::

    if (surface == surface_jules .and. albedo_obs) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call processor%apply(make_spec('albedo_obs_vis', main%radiation, &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('albedo_obs_nir', main%radiation, &
        ckp=checkpoint_flag))

.. topic:: The processor apply method

   Due to the way XIOS operates, Momentum's
   ``create_physics_prognostics`` routine is called twice. The first
   call and the second call use different ``processor`` types to
   ``apply`` different methods on each call. XIOS defines what are
   known as contexts to support a particular set of input or output
   requests. The first call is required to set up an XIOS "context"
   which involves XIOS reading the ``iodef.xml`` file. Once the
   context is set up, the model can use XIOS to query the field
   definitions that were read into the context from the ``iodef.xml``
   file.

   The field reference in the ``iodef.xml`` is sufficient to enable
   XIOS to register a field for potential diagnostic output, but
   additional work is required by the model to register fields for
   writing or reading via the checkpoint-restart system. The first
   call to ``create_physics_prognostics``, done early in the model
   initialisation phase, sets up of XIOS context, allowing the
   ``apply`` method to registers fields required for checkpointing
   with the context.

   Only once the XIOS context has been set up can the model read field
   definitions read from the ``iodef.xml`` file. These definitions are
   used to correctly initialise each field. In the second call to
   ``create_physics_prognostics``, a different ``processor%apply``
   method queries XIOS to find out the XIOS domain a field is on. The
   XIOS domain uniquely determines the particular internal model field
   type. For example, a field defined on the ``half_level_face_grid``
   in the ``iodef.xml`` file lives on the :math:`\mathbb{W}_{3}`
   function space within the model whereas a field on the
   ``half_level_edge_grid`` is on the :math:`\mathbb{W}_{2h}` function
   space. The ``apply`` routine initialises the field and adds it to
   the required field collections.

   In addition to the requested field collection, fields to be written
   to or read from checkpoint files may be included in a second
   checkpoint field collection. Fields can be included in more than
   one field collection by including the actual field in one
   collection and pointers to the field in other collections.  For
   consistency and simplicity, Momentum's ``apply`` method always adds
   fields to the requested collection and the checkpoint collection as
   pointers. The actual field is held in another field collection
   called the ``depository``. The ``depository`` just acts as a
   central hold-all for fields that need to remain in scope throughout
   the model run: while it exists in ``modeldb``, it is not intended
   to be used within the model.

During the finalise stage of a model, a procedure in the
:ref:`lfric-xios component <section lfric xios>` can be called and
passed the field collection containing the fields that need to be
checkpointed. After the fields are checkpointed, the collections can
be cleared and the fields can go out of scope.
