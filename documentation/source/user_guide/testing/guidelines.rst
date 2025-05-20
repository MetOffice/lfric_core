.. -----------------------------------------------------------------------------
     (c) Crown copyright 2025 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

Functional Testing Guidelines
=============================

There are some things to be aware of when developing functional tests. Some
are relevant regardless of the scale at which you are operating: unit,
integration or system. Others apply more to one scale or another but are still
worth bearing in mind at other scales.

Beware Zero
~~~~~~~~~~~

The problem with an expected output of zero is that any multiplication by
zero will give zero. Therefore there are potentially many wrong ways to come
up with the output zero.

If you can it's best to make your expected output non-zero, this removes the
possibility of an errant multiplication by zero giving a false positive.

Beware Flat Fields
~~~~~~~~~~~~~~~~~~

When testing with array inputs (including ``test_field_type``) it can be
tempting to make them all flat, flooded with a constant value. This makes
calculating the expected output easy since the same value applies everywhere.

The problem with this is that while it will check for the correct value it
can't tell you about order of operations. In particular it can't tell if a
kernel has taken the correct tour around its stencil.

Thus it is best if at least one array going into your test is not flat.

A simple way to achieve this for algorithm tests making use of
``test_helpers_mod`` is to use the ``create_test_field_flood_index`` which
fills each data element with its index in the test data array.

Remember that an unstructured mesh is used so although each cell value is
guaranteed unique it is not guaranteed to appear in any particular order within
the context of the mesh.
