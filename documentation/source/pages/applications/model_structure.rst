.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _section model structure:

Model structure
===============

This section focuses on the needs of a typical numerical model, and
the structure of the application that runs the model. Here, we
envisage an application that configures, sets up and calls the model
to integrate field data on model meshes over multiple time-steps,
reads input data, outputs diagnostics, and needs to checkpoint
progress.

Furthermore, we assume that the model is able to use the standard
features provided by LFRic core.

TODO: Add a diagram that illustrates a simplified calling tree.

Running the Model
=================

At the top level of a typical model, there will be a driver layer
comprising an initialise, step and finalise stage. These need to be
executed in sequence to progress the model.

Prior to initialisation, data structures shared by each stage of the
model need to be initialised. An LFRic component provides a `modeldb
<section modeldb overview>` data structure that aims to store all the
data needed to run the model. For a model that uses `modeldb` aspects
of the model data can be configured prior to calling the
initialisation, such as setting the MPI communicator, setting up the
configuration and defining the model name. This can allow multiple
instances of the model to be run alongside each other.

Evolution of the model is driven by executing the model driver step a
required number of times, potentially controlled by ticks of the
`model clock <section model clock>`. The driver step is responsible
for calling the model step that integrates the model data forward one
time-step, but will also be responsible for managing infrastructure
such as reading input data, and writing some diagnostic data and
writing checkpoint dumps.

Once all model steps are executed, the model finalise stage is called
after which processes instantiated by the driver layer can be finalised.

The driver layer
================

Initialise
----------

The `initialisation` stage can roughly be divided between initialising
the infrastructure of the model and initialising the initial model
state.

The model infrastructure can typically comprise information about
meshes and coordinates that are fixed throughout the model run as well
as some fixed data. The model state comprises fields and data that
evolve as the model runs. Commonly, the model state will be
initialised in some way, either from input data or by computing
values.

Model initialisation may also initialise scientific components that
are used by the model.

Step
====

The `step` stage will execute a single time-step of the model starting
at the input date and lasting for a period defined by a time-step
length.

Finalise
========

The `finalise` stage will undertake any necessary finalisation
processes, noting that much of the model data may go out of scope as
soon as the driver layer finalise has completed.

Data in a model
===============

The `modeldb` object provides the ability to store a range of
different data structures. In support of models such as the Momentum
atmosphere model, it provides support for storing many `field
collections <section field collections>` which are arbitrarily large
collections of individual model fields that can be accessed by
name. Field collectons can be populated at run-time meaning the list
of fields can depend on the model configuration.

Using these features, a model such as the Momentum atmosphere can
define a field collection for each of the diverse set of science
components it uses, each containing a number of fields.

Other aspects of `modeldb` enable storing of other data types.

Algorithms and Kernels
======================

The calling structure below the stages above will involve calling
science algorithms. Unless an algorithm requires a very large number
of input fields it is recommended that science algorithm APIs require
lists of individual fields extracted from field collections rather
than the field collections themselves.

At the lowest level, algorithms will call kernels. TODO add reference
to PSyKAl sections.

Use of components
=================
The `driver component <section driver component>` provides code for a
range of functions required at the driver level for setting up data
and methods TODO improve this brief summary.
