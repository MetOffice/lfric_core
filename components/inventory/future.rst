..
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.

Ideas For the Future of Inventory Component
===========================================

This document holds thoughts, sketches and notions for potential future work on
the inventory component.

.. attention::

  The existence of an idea in this document is no guarantee it will ever be
  implemented. Ideas are half-baked and highly speculative. They are captured
  in order that they don't get lost.

Now We Know More
----------------

The "inventory" system was put in place to solve an immediate problem without a
full understanding of what the problem was. As such it performs well, however
we may now look at it to better understand the problem and it seems clear
there are issues.

We have an ever expanding pool of global scope data encapsulated in an ever
expanding source file. This is becoming unwieldy and likely to only become
more so.

A chunk of the inventory is for data keyed on the mesh ID. This has a few
problems.

It doesn't meet our aspiration to adopt an object approach with tight
cohesion and loose coupling. To be clear, it could be a lot worse, but it could
be better. In particular we shouldn't be referring to meshes by ID, instead we
should hold a mesh object which is the mesh and does mesh-like things.

Since meshes are singletons shared between users this data should probably be a
property of, or at least associated with, the mesh object. Rather than being
some secondary look-up.

One potential implementation would be a key/value pair store in the mesh object
which could then be populated as needed.

This makes mesh associated data a generic and extensible concept. It also
neatly tackles the problem of Gungho or Atmosphere specific data being held
and described in the inventory component. Instead the code to generate these
values can live in the models which need them. The code for managing this data
becomes an aspect of the mesh.

The down side of this is that it increases the complexity of the mesh object
but it is more conceptually consistent. Properties of the mesh live with, and
are managed by, the mesh. There is also a generic key/value collection
implementation pending.

The same approach could be taken with function-spaces if function-space
specific inventory items exist.

The problem comes if there are any items keyed off multiple objects. In that
case it is not clear who should own them. There may be a requirement for an
"inventory collection" which can "own" the data, then have pointers to it from
the various objects which it refers to. Such data may not exist, in which case
it's not necessary.

Regardless, there is no need to solve all problems at once. Anything removed
from the inventory component is potentially a gain. Even if a residue remains
it may still be worth doing.

Luckily this need not be a "big bang" change. Adding the ability to manage
additional data to the mesh object is a change transparent to the rest of the
system as it extends the API only, it does not alter existing API.

With that in place inventories may be moved over one at a time. This will limit
the size of any one change.
