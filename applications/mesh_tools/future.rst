..
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.

Ideas For the Future of Mesh Tools
==================================

This document holds thoughts, sketches and notions for potential future work on
the mesh tools.

.. attention::

  The existence of an idea in this document is no guarantee it will ever be
  implemented. Ideas are half-baked and highly speculative. They are captured
  in order that they don't get lost.

More Parallel Partitioning
~~~~~~~~~~~~~~~~~~~~~~~~~~

The use of shared memory OpenMP parallelism has greatly improved the time to
generate and partition (particularly partition) a global mesh. However it is
limited by the I/O bandwidth available on a node. If further speed up were
required we would have to look at spreading the job across multiple nodes in
order to gain more I/O bandwidth.

The obvious approach is to use MPI but partitioning is an embarrassingly
parallel problem and as such there are no messages to pass with the message
passing interface.

In light of this a simple matter of running the mesh generator multiple times,
once per socket, with instructions to generate only a subset of the partitions
would probably be the way to go. It will be a lot simpler than implementing MPI
while not losing any functionality.

Pipeline Fever
~~~~~~~~~~~~~~

The mesh generator has evolved over time as one monolithic slab of code. This is
hard to maintain and harder to extend. What is needed is a more modular
approach, luckily the "pipeline" pattern is right there.

In a pipeline each stage of a process is represented by a stage in the pipeline.
These stages may be re-ordered within sensible constrains. For instance, you
can't partition a global mesh until you have generated a global mesh. Stages can
be included, or not, depending on need. If you do not need the pole rotating you
just omit the pole rotation stage.

From the above it is clear that there are initiating stages which can start a
pipeline. These will be the cubed-sphere and planar mesh generators. There are
also mutation stages which modify a mesh. These may appear in arbitrary order.
Finally there are terminating stages which will end a pipeline, such as writing
out the data to file.

There may be hybrid stages such as partitioning which accepts a global mesh,
like a mutator, but produces a local mesh, like an initiator.

Another mutator/initiator is the flow split which takes a mesh (global or
local), performs a null mutation and passes it on but also takes a copy and
initiates a new pipeline with it. This makes it possible for pipelines to
bifurcate.

With the pipeline structure in place, creating the required meshes becomes a
matter of building an appropriate pipeline and then hitting "play." The correct
files should fall out the end by magic.

..
   This is provided in preparation for migration to new development platform
   and integration with Sphinx, which can support graphviz "dot" syntax.

.. graphviz::
    :caption: Diagram of an example mesh pipeline.
    :alt: Cube-sphere global mesh generator feeds into poll rotation mutator,
          feeds into flow splitter, feeds into both global mesh writer and
          partitioner. Partitioner feeds local meshes into local mesh writer.
    :align: center

    digraph "pipeline-example" {
      rankdir = "LR"
      "cube-sphere" -> "pole rotator" [label="global mesh"]
      "pole rotator" -> "splitter" [label="global mesh"]
      "splitter" -> "global writer" [label="global mesh"]
      "splitter" -> "partitioner" [label="global mesh"]
      "partitioner" -> "partition writer" [label="local mesh"]
    }

Maintenance and extension become easier as each stage is stand-alone and may
be developed and tested as such.

There is one remaining super-power of the pipeline which is that it may be
split between executables at any point.

The entire pipeline may exist in the mesh generator, the result being generated,
modified and partitioned mesh files coming out the end.

Alternatively the pipeline may be split at the point where the global mesh is
complete and written to disc. Reading that mesh and partitioning it can then
exist in the model executable using the same code.

More exotic splits are possible too. If you wanted one executable which
generates a set of different resolution global meshes, a second which mutates
a mesh in particular ways and partitioning happens in the model, that is
possible as well.

It might even be possible to consider various further mesh operations within
the model as part of the pipeline. Applying orography, for instance.

The fact that the pipeline is arbitrarily splittable also points towards easier
migrating from the current situation to a full pipeline. There is no need to
migrate everything all at once but take advantage of the fact that a pipeline
may be split at any point. I will use the example of starting with the
production end but it works just as well starting at the other.

Extract the basic global mesh generation (cube-sphere and planar) from the
monolith and create two potential pipeline initiator stages. Create a pipeline
for them to exist in and logic to instantiate the correct stage in the pipeline.
The pipeline then feeds the resulting global mesh into the existing monolith.

This gives a change set which may be committed to trunk and from the user's
point of view is no different.

The next stage is to extract first mutation from the monolith and implement that
as a pipeline stage. Add logic to instantiate, or not, this stage depending on
configuration. Once again the result feeds into the shrinking monolith and the
user is none the wiser.

Repeat until the monolith has been pared away to nothing and everything is an
optional pipeline stage.
