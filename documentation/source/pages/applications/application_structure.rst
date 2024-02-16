.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section application structure:

Application Structure
=====================

An LFRic application is defined as an executable that runs one, or
possibly more, LFRic models. This section describes the structure of a
simple application that runs just one LFRic model.

Here, we envisage an application that configures, sets up and calls
the model to integrate field data on model meshes over multiple
time-steps, reads input data, outputs diagnostics, and needs to
checkpoint progress.

Furthermore, it is assumed that the model uses the standard features
provided by LFRic core.

TODO: Add a diagram that illustrates a simplified calling tree.

Running the Model
-----------------

The application will call the model via a driver layer. The driver
layer may be standard or bespoke. It comprises an initialise, step
and finalise stage. These stages need to be executed in sequence to
progress the model.

Prior to calling the initialise stage, data structures shared by each
stage of the model need to be set up. An LFRic component provides a
:ref:`modeldb <section modeldb overview>` data structure that aims to store
all the data needed to run the model. Properly encapsulating all the
data allows applications to run two or more models side by side, or to
run ensembles of the same model side by side.

For a model that uses `modeldb`, some aspects of the data structures
held in `modeldb` must be configured prior to calling the
initialisation. Examples include setting up the configuration,
defining the model name and setting the MPI communicator. Configuring
these aspects before calling initialisation can allow multiple
instances of the model to be run alongside each other either
concurrently or sequentially.

Evolution of the model is driven by calling the model driver step the
required number of times, potentially controlled by ticks of the
:ref:`model clock <section model clock>` held in `modeldb`.

Once all steps are executed, the model finalise stage is called after
which processes instantiated by the application prior to
initialisation can be finalised.

The driver layer
----------------

The LFRic infrastructure provides a `driver` component containing code
that can be used when constructing the application driver layer.

Driver Initialise
~~~~~~~~~~~~~~~~~

The driver `initialisation` stage can roughly be divided between
initialising the infrastructure of the model, such as meshes,
coordinates, clocks and calendars, and initialising the initial model
state, including the reading initial data. The model initialisation
will provide procedures for the processes required to complete the
initialisation; separating these processes into multiple procedures
gives applications flexibility in setting up models, for example,
optimising setup where several models use the same or similar meshes.

The model infrastructure typically comprises information about meshes
and coordinates that are fixed throughout the model run, as well as
some fixed data. The model state comprises fields and data that evolve
as the model runs. Data held in the model state needs to be
initialised, either from input data or by computing or choosing
values.

Two similar models running within the same application may share some
constant data.

The driver initialisation may also initialise scientific components
that are used by the model.

Driver Step
~~~~~~~~~~~

The driver `step` stage will execute a single time-step of the model
starting at the input date and lasting for a period defined by a
time-step length.  The driver step is responsible for calling the
model step that integrates the model data forward one time-step, but
will also be responsible for managing infrastructure such as reading
input data, and writing some diagnostic data and writing checkpoint
dumps.

Driver Finalise
~~~~~~~~~~~~~~~

The `finalise` stage will undertake any necessary finalisation
processes, noting that much of the model data may go out of scope as
soon as the driver layer finalise has completed.

The model API
-------------

Mirroring the structure of the driver layer, the model layer will have
initialise, step and finalise stages.

Model Initialise
~~~~~~~~~~~~~~~~

As noted above, the initialise stage may be broken into several
separate procedures to allow for flexibility in application design.

On completion of initialisation, the internal model data structures
should be fully-set up in readiness to run a model timestep.

Model Step
~~~~~~~~~~

The model step will evolve the model prognostics forward by one
timestep.

Model Finalise
~~~~~~~~~~~~~~

The finalise stage will finalise any objects created in the initial
stage.

Data in a model
---------------

Initially, the LFRic field is probably the most important LFRic data
structure to get to know: understanding the role of the field is
critical to understanding LFRic, but the details can be deferred to
the section describing the :ref:`use of PSyclone and the LFRic data
model<section psyclone and the lfric data model>`. For now, it is
sufficient to understand that an LFRic field is a data structure that,
among other things, holds data representing a physical quantity such
as winds, temperature, or density, over the domain of a model mesh.

A complex model such as the Momentum atmosphere requires hundreds of
fields. To simplify the model design, more data structures than just
fields are provided by the LFRIc infrastructure. This section briefly
summarises other data structures that are available, with links to the
more in-depth documentation.

In addition to fields, :ref:`field collections <section field
collections>` can store arbitrarily-large numbers of fields or field
pointers that can be accessed by name. The Momentum atmosphere has
several field collections for each of several major science
components. Use of field collections makes the API of higher-level
science algorithms more manageable by hiding both the large number of
fields and the fact that some fields are not required for all model
configurations.

A :ref:`configuration object <section configuration object>` stores
the model configuration derived from the input namelist, such as input
values for real variables, science options and switches. Settings can
be accessed by a name based on the namelist name and the variable
name.

A :ref:`key-value <section keyvalue pair object>` data structure exist
that stores an arbitrary number of key-value pairs where the value can
be an object of any type. At a basic level, this data structure can
store native fortran types such as real or integer variables and
arrays. More complex abstract or concrete types can also be stored.

The `modeldb` object defined in the `driver` component provides the
ability to store  range of different data structures in addition to
fields. A list of the main data structures declared in `modeldb` is
given here. For more details on how to use these data structures see
the :ref:`modeldb <section modeldb>` documentation.

 - **field** Stores fields and field collections as key-value pairs.
 - **configuration** As described above, stores the model
   configuration: input values, science options, switches and so
   forth.
 - **values** A key-value data structure described above that can
   store any type or class.
 - **mpi** Stores an object that can be used to perform MPI tasks.

While all algorithms in an LFRic model will rely on fields, to retain
a degree of separation between the model and the infrastructure it is
recommended that accesses to `modeldb` do not go too deep into the
code: once an algorithm is sufficiently self-contained, all its inputs
can be extracted from `modeldb` and passed to the algorithm through
the subroutine API.

Operators
~~~~~~~~~

A brief mention of operators is sufficient in this document: an
operator is a data structure that can be used to map a field of one
type onto another type. Its use is relevant to the GungHo mixed finite
element formulation where there is a need to map fields between
different function spaces.

Algorithms and Kernels
----------------------

The work of the initialisation and step processes of the model will be
managed by a set of algorithm and kernel modules. Broadly speaking,
algorithms are higher-level subroutines that deal only with full
fields, and kernels are lower level subroutines that implement the
computation. The structure of algorithms and kernels is governed by
the :ref:`PSyKAl design <section psykal>`.
