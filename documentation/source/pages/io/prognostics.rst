.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section prognostics:

Prognostic Fields
=================

This section describes how to add fields that need to stay in model
state throughout a run and that may need to be written to the
checkpoint dump if the run needs to be stopped and restarted.

Formally, prognostic fields are those fields that capture the physical
quantities within the system being modelled and are therefore core to
the evolution of a numerical model such as wind, density, potential
temperature and moisture. However, the practical difference between
these true prognostic fields other fields that are generated in one
part of a model and need to be passed to another part are small: in
both cases, the fields need to remain in memory through all or part of
the timestep and may need to be stored in the checkpoint dump.

.. note:: Illustrative examples

   #. Fields such as wind, potential temperature, density and moisture
      are physical quantities at the core of the numerical
      simulation. They need to be initialised at the beginning of a
      model run and exist until the run completes. The fields need to
      be written out to any checkpoint file.
   #. Some fields may be read in from file. It is often convenient to read
      them in at the start of a run to keep the IO code separate from
      the scientific code that uses the field data. Hence, like
      prognostics, they need to be initialised at the beginning of the
      run. But they do not need to be written to a checkpoint file as
      the original file can be read in again when the model restarts.
   #. Where a field is passed from one section of science to another
      (such as from radiation to boundary layer) it can be convenient
      to allocate them at the start of the run rather than allocate
      them within the time-step code. If the field is passed across a
      time-step boundary then it will also need to be written to the
      checkpoint file.

This section describes how to add a new prognostic field to an
application. The documentation is targeted at users of the Momentum(R)
atmosphere model which currently uses the XIOS IO library. Commentary
is included for models that use a different IO subsystem.

Overview
--------

In the Momentum(R) atmosphere model, all fields that the model may
read or write, including all prognostic fields, are registered in an
``iodef.xml`` file which is the file recognised by the XIOS IO
library.

Each new prognostic field needs to be added to the ``iodef.xml``
file. The record includes a string ID for the field that is used
within the model, the formal name of the field (such as the CF name),
the units and information about the XIOS domain used to describe the
format of the data in the input or output file.

Prognostic fields should be initialised to the right field type at the
start of a run based on the ``iodef.xml`` metadata. The fields must
exist throughout the model run. Some prognostic fields need to be
written to the checkpoint file during the run or at the end of the
run, to support the ability to restart the model run. Other fields are
used only within the timestep, or can be reinitialised from other data
at the start of a continuation run, and so do not need to be written
to the checkpoint file.

In the Momentum(R) atmosphere model, a routine called
`create_physics_prognostics` includes the code for creating each
field. Logic can be used to create only those fields that are used by
the configuration, and to mark particular fields for writing to the
checkpoint dump.

::

    if (surface == surface_jules .and. albedo_obs) then
      checkpoint_flag = .true.
    else
      checkpoint_flag = .false.
    end if
    call processor%apply(make_spec('albedo_obs_vis', main%radiation,            &
        ckp=checkpoint_flag))
    call processor%apply(make_spec('albedo_obs_nir', main%radiation,            &
        ckp=checkpoint_flag))

The ``processor%apply`` method has two roles. The first is to register
the field with the runtime XIOS system, the second is to query XIOS
for the field metadata held in the ``iodef.xml`` file and then
initialise the field based on that metadata: the metadata will
indicate, for example, whether to create a W3, a Wtheta field, or a 2D
multidata field.
