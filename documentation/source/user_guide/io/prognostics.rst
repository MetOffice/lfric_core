.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section prognostics:

Prognostic Fields
=================

In an earth system model, the fields that exist throughout a model run
are often referred to as the "prognostic fields". For some such
fields, it is essential that they exist during the whole run as they
contain the true model state. For other fields, implementing them as
prognostic fields can be a practical convenience where fields are
shared between more than one area of the model or where they are set
by reading in data from file.

The LFRic infrastructure provides support and recommendations for
creating and managing the model prognostic state for applications
aligned with the :ref:`standard LFRic model application structure
<section application structure>`. In such applications, fields are
created and initialised during the initialisation phase and stored in
the ``modeldb`` data structure. For a large model with many prognostic
fields, :ref:`field collections <section field collection>` can be
used to hold many fields.

In the time-step phase of the model, fields and field collections are
extracted from ``modeldb`` and passed to science code through the
argument list. Packaging fields in field collections is a useful way
of keeping argument lists manageable when a procedure takes a large
number of fields as arguments.

The following figure expands on the figure found in the
:ref:`application structure section <section application structure>`
to illustrate the role of each stage of a model in creating and using
prognostic fields.

.. figure:: images/prognostics.svg

   Illustrates the role of each stage of a model when dealing with
   prognostic fields.

.. topic:: Examples of prognostic fields

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
stage. The infrastructure provides a model database structure called
``modeldb`` which can hold the model state and other data through all
three stages.

For such models, prognostic fields are initialised during the
initialise stage and stored in ``modeldb``. As well as storing
individual fields, ``modeldb`` can also store :ref:`field collections
<section field collection>`. Use of field collections permits large
numbers of fields to be passed down a calling tree within a single
argument. Fields can be added to a field collection dynamically,
meaning the contents of the field collection can depend on the
run-time configuration of the model.

Some prognostic fields need to be written to, and then read from,
checkpoint dumps if a long run is to be broken into shorter segments.
As fields can be referenced by more than one field collection, storing
all fields required for checkpointing in a dedicated field collection
can simplify the process of writing a checkpoint dump.

Note, that the actual field can only be included in one field
collection. Other field collections can include a pointer to the same
field. To simplify matters, it is recommended that all of the actual
prognostic fields are added to a single large ``depository`` field
collection, and that other field collections hold pointers to the
prognostic fields they need.

In summary, in the initialisation phase, for each field, the following
steps are recommended.

 #. If the field is a prognostic it needs to be initialised and added
    to the ``depository`` field collection.
 #. To pass each field to the science sections that need it, a pointer
    to the field is added to relevant field collections passed down
    the calling tree to the science code.
 #. If the field needs to be checkpointed, a pointer to the field is
    added to a field dedicated to holding checkpoint/restart fields.
 #. The field should be initialised with the appropriate mesh and
    function space according to the model configuration.

During the time-step running stage of the model, field collections and
fields set up at initialisation time can be extracted from the
``modeldb`` and passed into the relevant scientific algorithms and
kernels.

.. note::

   Each of the above processes can be written as hard-wired code,
   initialising each field in turn, and adding it to each field
   collection required.

    .. code-block:: fortran

       type(field_collection) :: depository, conv_field, checkpoint_fields
       type(field_type)          :: rain
       type(field_type), pointer :: fld_ptr
       class(pure_abstract_field_type), pointer :: tmp_ptr

       ! Get the field collections from the modeldb
       depository => modeldb%fields%get_field_collection("depository")
       conv_fields => modeldb%fields%get_field_collection("conv_fields")
       checkpoint_fields => modeldb%fields%get_field_collection( &
                                             "checkpoint_fields")

       ! rain is a 2D field
       vector_space => function_space_collection%get_fs( mesh2d, &
                                              element_order, W3 )
       call rain%initialise(vector_space, name='rain')

       ! Add field to depository, then get a pointer to the field
       call depository%add_field(rain)
       call depository%get_field(rain, field_ptr)

       ! Abstract class needed to support all possible field types
       tmp_ptr => field_ptr
       ! Add field pointer to the radiation collection
       call conv_field%add_reference_to_field(tmp_ptr)

       ! If convection scheme is in use,  checkpoint this field
       call checkpoint_fields%add_reference_to_field(tmp_ptr)

   Alternatively, some aspects of the process can be based on sources
   of metadata so as to automate the process and to reduce the risk of
   error. For example, the Momentum\ :sup:`®` atmosphere makes use of
   metadata provided as part of its integration to the XIOS library to
   determine which type of field should be initialised. The use of
   metadata ensures that the type of field created matches up with the
   type of field being written out if the field is output to a
   checkpoint file or as a diagnostic.

Once all the model steps are completed, a checkpoint dump may need to
be written to enable the simulation to be extended in a separate run
of the application. The field_collection object supports an iterator
that allows a function to loop through all the fields in a field
collection. If a checkpoint field collection exists, it can be passed
to a procedure that goes through the collection writing out each field.

Example: The Momentum\ :sup:`®` Model
-------------------------------------

In the Momentum\ :sup:`®` atmosphere model, the set-up of the
prognostics field is based on metadata written in the XIOS
``iodef.xml`` file. All possible prognostic fields are recorded here.
The record includes a string ID for the field, the formal name of the
field (such as the CF name), the units and information about the XIOS
domain used to describe the format of the data in the input or output
file.

In the Momentum\ :sup:`®` atmosphere model, a routine called
``create_physics_prognostics`` includes the code for creating most of
the prognostic fields, including dealing with the steps described in
the preceding section. The ``create_physics_prognostics`` is
responsible for determining whether each potential prognostic field is
required by the model, which scientific field collection it needs to
be added to and whether it should be checkpointed.

A function is then called which both initialises the field and adds it
to the appropriate set of field collections.

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
