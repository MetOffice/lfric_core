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
at isolating a problem.

Between the two is "integration testing" which considers clusters of units.
This allows a sub-set of interactions to be exercised while still allowing for
reasonable isolation.

A good testing regime makes use of all these approaches.

Note that the boundaries are not hard drawn. A unit under test may well
call down to procedures further down the call tree. Thus the unit test is
testing multiple units. Meanwhile an integration test may be testing a
substantial fraction of the whole system test.

Where these boundaries are drawn is a matter of ongoing discussion and debate
and outwith the scope of this document.

.. toctree::
    :maxdepth: 1

    unit_testing
    unit_test_canned_data_routines
    integration_testing
    integration_testing_algorithms
    guidelines
