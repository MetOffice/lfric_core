.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _testing:

Testing
=======

Comprehensive and automated testing is an essential aspect of successful
software developement. Comprehensive because we can't find bugs in code we
aren't testing. Automated because tests which aren't run aren't catching bugs.

The LFRic project makes use of functional testing whereby a known set of
input stimuli are presented to a piece of code and the resulting output is
compared to expected results. Any deviation indicates a change in behaviour.

This functional testing is a continuum from fine to coarse grained.

The finest grained end is "unit testing" where code units are tested. This
usually means individual procedures. This testing is very good at isolating
faults to a small piece of code but it can't tell you how these units interact
with each other.

The coarsest grained end is "system testing" where the complete executable
is tested. This is the ultimate test of interaction between units but is poor
at isolating a problem. Our "Rose stem" test suite is an example of system
testing.

Between the two is "integration testing" which considers clusters of units.
This allows a sub-set of interactions to be exercised while still allowing for
reasonable isolation.

.. figure:: /how_to_use_it/images/testing_continuum.svg
   :alt: Diagram showing a continuous ribbon of unit testing to integration
         testing to system testing, with a second, de-emphasised, integration
         testing beyond that.
   :align: center

   The functional testing continuum.

A good testing regime makes use of all these approaches.

Note that the boundaries are not hard drawn. A unit test may well call down to
procedures further down the call tree. Thus the unit test is testing multiple
units. Meanwhile an integration test may be testing a substantial fraction of
the whole system.

Where these boundaries are drawn is a matter of ongoing discussion and debate,
outwith the scope of this document.

Be aware that "integration testing" is also used to refer to the testing of
interaction between systems and so also sits beyond "system testing" as
illustrated above. It is the same basic concept but at different scale. The
systems and their interactions become the units under test.

There is an element of this inter-system integration testing in our "Rose stem"
suite. The mesh generator is built, then used to generate meshes which are
then loaded by models. Thus the interaction between mesh generator and model
is tested.


.. toctree::
    :maxdepth: 1

    unit_testing
    unit_test_canned_data_routines
    integration_testing
    integration_testing_algorithms
    guidelines
