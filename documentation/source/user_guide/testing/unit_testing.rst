.. -----------------------------------------------------------------------------
     (c) Crown copyright 2024 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _unitTesting:

Unit Testing
============

The purpose of unit testing is to test the smallest discrete pieces of
functionality that have a clearly specified API. This will typically be a
procedure. It is not uncommon to group all tests pertaining to a single program
unit (for example, a module or type) together.

The unit test code can be found under the ``<project>/unit-test`` directory, in
a tree which mirrors that in ``<project>/source``. Thus the test for a source
file may be easily located by looking in the corresponding location in the unit
test tree.

To aid the writing of unit tests, we use the pFUnit unit testing framework. This
provides a preprocessor to generate test case source and a driver to run the
tests.

Test code can be identified by the ``.pf`` extension and are normal Fortran
files garnished with ``@`` (at symbol) prefixed directives. The processor
substitutes the directives with additional Fortran code and Fortran preprocessor
directives.

It is possible to use a ``.PF`` extension. This causes the Fortran preprocessor
to be run on the file prior to the pFUnit processor. This is necessary for
testing templated source but should otherwise be avoided.

When the resulting program is run, the pFUnit framework marks progress and
reports any failures.

How to Unit Test a Thing
------------------------

A unit test should test only the unit under test. This may seem obvious but it
is easy for a panoply of other things to work their way into a test.

Inputs to the unit under test should never be left to the tender mercies of
undefined or unclear behaviour. Always set them explicitly and ideally to
constants.

The expected result should be constants or calculated from constants. Try to
avoid calculations so as to avoid the possibility that the same (flawed)
algorithm is used in test and unit under test.

Each test procedure should contain a well defined set of tests. The failure of a
test should indicate that a fault lies in a restricted amount of code. A single
test procedure exercising an entire module has some value but if it exercises
only a single procedure under test, your bug hunting effort has been greatly
reduced.

The test should rely on as little external code to set up the test environment
as possible. Ideally none at all. Any external code becomes, by implication,
part of the unit under test. It must work correctly before the test can pass.

Use of LFRic Infrastructure in Unit Tests
-----------------------------------------

Unless you are testing a particular part of the infrastructure the rule is:
Don't use it in a unit test.

This is particularly relevant to testing kernel code. The inputs to kernel
procedures are simple arrays of primitive types. In a model run, these arrays
are derived from infrastructure objects such as fields and function spaces. It
is very easy to rely on these infrastructure objects in the unit tests. However,
this is a poor idea for a number of reasons.

By using the infrastructure you are adding a lot of dependencies into your test.
The effect of using them is that you are not only testing the unit you are
interested in but also all the infrastructure used in the test. This reduces the
locality of any bug discovered.

It also makes the tests much more complicated to implement and to read. This in
itself invites faults.

Canned examples of the arrays used by kernels are provided to set up simple
cases sufficient for most testing. These also add external dependencies which
widen the scope of the test but to a much lesser extent. They are based on
canned data which is merely copied into the appropriate array. Thus they
represent a much lower risk.

Information on how to use the provided helper routines that return the canned
data can be found in the section: :ref:`canned_info`.

What to Test
------------

Obviously you should test the output of the unit against expected good output.
Make sure that when a set of values go in, the correct values come out. When
choosing your test case make sure you don't go for values which can hide
failures. For instance testing that when you get zero out may not be helpful.
There are many ways to get zero out of a calculation.

It may be worth testing a number of good inputs in case by some fluke you happen
to choose the one set of inputs which gives the correct result even though every
other option returns the wrong result.

Where possible, it is also worth testing a few expected failure modes. The
current unit testing framework doesn't handle aborts on error (the whole test
suite will simply abort), but if the error is handled by the code, it can be
tested. For example, if some functionality requires a positive integer, give it
a negative one.

Edge cases are another good area to test. Check for out-by-one errors by passing
values either side of limits.

The idea is to build confidence that the unit not only functions correctly most
of the time but also in stress conditions.

How to Test
-----------

What follows is an introduction to a number of different approaches and
techniques for writing unit tests. In general you should always avoid
implementing a class derived from ``TestCase`` unless you need it to hold
fixture data. Fixtures which do not create data which can be held in this way,
such as namelist feigners, should not cause a ``TestCase`` class to be created.

Even when using fixtures you should prefer standalone ``@before`` and ``@after``
subroutines over ``setUp`` and ``tearDown`` methods.

Simple Test Procedures
----------------------

A minimal example might look like this:

.. code-block::

  module simple_mod_test

    use pFUnit_Mod
    use simple_mod, only: thing, thang

    implicit none

    private
    public test_thing, test_thang

  contains

    @test
    subroutine test_thing

      integer :: result

      result = thing( 1, 2, 3 )
      @assertEqual( 12, result )

    end subroutine test_thing

    @test
    subroutine test_thang

      integer :: result

      call thang( result )
      @assertEqual( 13, result )

    end subroutine test_thang

  end module simple_mod_test

Note the pFUnit directives prefixed with the ``@`` (at) symbol. These are used
to mark out the test cases with ``@test``. They are also used to denote
"assertions". These are the actual business end of a test case. An assertion
must be met in order for the test to pass.

.. note::

   A restriction of pFUnit is that any line that start with the ``@`` symbol
   must all be written on a single line. Continuation lines are not permitted.

In this case ``@assertEqual`` is used. Surprising no one, this requires that the
value from the unit under test (the second argument) must be equal to the
expected result passed in the first argument.

Given the problems inherent in testing for equality between floating point
numbers a fuzzy match may be used. Simply pass tolerance argument, e.g.
``tolerance=0.001``.

A wide variety of other assertions are provided including inequalities such as
``@assertGreaterThan`` and numerical tests such as ``@assertIsNan``.

Also provided is the general purpose ``@assertTrue``. This allows any
unsupported test to be implemented using an expression which yields a boolean
result and testing for its success.

When an assertion fails it will produce a failure message which includes details
about what was expected and what was found. This is usually exactly what is
needed in order to diagnose the problem but in some cases it can be unhelpful.
When this happens the optional "message" argument may be used.

Test Procedures with Fixtures
-----------------------------

The simple test outlined above is fine for situations where there are no
resources to be managed. When there are, "fixtures" are needed.

The most common resource requiring management is memory, allocated space must be
deallocated, but configuration may also be considered a resource in this
context.

Such managed resources are referred to as "fixtures" in testing parlance. They
are things which must exist in order to perform the test but which are not under
test themselves. They are created and initialised before each test and destroyed
after it. This means that each test has a pristine, known, environment in which
to work.

Following is a minimal example of fixtures in use:

.. code-block::

  module simple_mod_test

    use pFUnit_Mod
    use simple_mod, only: thing, thang

    implicit none

    private
    public test_thing, test_thang

    real, allocatable :: data_block(:)

  contains

    @before
    subroutine setUp()

      implicit none

      allocate( data_block(256) )

    end subroutine setUp

    @after
    subroutine tearDown()

      implicit none

      deallocate( data_block )

    end subroutine tearDown

    @test
    subroutine test_thing()

      implicit none

      call thing( 1, 2, 3, data_block )
      @assertEqual( (/12, 13, 14/), data_block )

    end subroutine test_thing

    @test
    subroutine test_thang()

      implicit none

      integer :: result

      data_block = (/-1, -2, -3/)
      result = thang( data_block )
      @assertEqual( 13, result )

    end subroutine test_thang

  end module simple_mod_test

Notice how ``data_block`` is allocated in ``setUp`` and deallocated in
``tearDown``. These procedures are called before and after each test method.
This means that the contents do not carry between tests. Therefore they must be
initialised for each test. This may be done in ``setUp`` if every test needs the
same initial condition or locally to the test if they need different starting
points.

Feigning Configuration
----------------------

If the unit under test makes use of namelist values from configuration modules
they must be suitably initialised. You can not rely on them having been set up
by a previous test as they have not been. You can not rely on them defaulting to
a particular value because they do not.

Failure to do this will lead to the tests failing to pass, in ways that may be
unexpected and unpredictable, as the uninitialised (or initialised elsewhere)
parameters change value.

To perform this initialisation use the "feign" functions provided by
"feign_config_mod" to explicitly set the values needed. These functions work by
creating a temporary namelist, then telling the configuration system to load it.

Bear in mind that you must feign everything the code being tested will use. If
the unit calls down to helper procedures, any configuration they make use of
must also be feigned.

When calling feign functions, use named arguments. This improves the self
documenting nature of the code and protects against new arguments upsetting the
ordering.

For example: 

.. code-block::

  call feign_planet_config( gravity=10.0_r_def,    &
                            radius=6000000_r_def,  &
                            omega=8.0E-5_r_def,    &
                            rd=300.0_r_def,        &
                            cp=1000.0_r_def,       &
                            p_zero=100000.0_r_def, &
                            scaling_factor=1.0_r_def )


MPI Testing
-----------

Our unit test framework supports running tests in an MPI environment. This is
necessary for testing code that will only run in parallel. This type of code
currently only resides in the infrastructure (e.g. the partitioner or local
meshes)

It is easy to use, it just needs the number of processes to use to be specified
for each test.

.. code-block::

  module simple_mod_test

    use pFUnit_Mod
    use simple_mod, only: thing, thang

    implicit none

    private
    public set_up, tear_down, test_thing, test_thang

  contains

    @before
    subroutine set_up( this )

      use feign_config_mod, only : feign_stuff_config

      implicit none

      class(MpiTestMethod), intent(inout) :: this

      !Store the MPI communicator for later use
      call store_comm(this%getMpiCommunicator())

      call feign_stuff_config( name='foo', value=13 )

    end subroutine set_up

    @after
    subroutine tear_down( this )

      use configuration_mod, only : final_configuration

      implicit none

      class(MpiTestMethod), intent(inout) :: this

      call final_configuration()

    end subroutine tear_down

    @test( npes=[1] )
    subroutine test_thing( fixture )

      implicit none

      class(MpiTestMethod), intent(inout) :: fixture

      integer :: data_block(3)

      call thing( 1, 2, 3, data_block )
      @assertEqual( (/12, 13, 14/), data_block )

    end subroutine test_thing

    @test( npes=[1, 2, 4] )
    subroutine test_thang( fixture )

      implicit none

      class(test_simple_type), intent(inout) :: fixture

      integer :: data_block(3)
      integer :: result

      data_block = (/-1, -2, -3/)
      result = thang( fixture%context%getMpiCommunicator(), data_block )
      if (fixture%context%isRootProcess)
        @assertEqual( -13, result )
      else
        @assertEqual( fixture%context%processRank(), result )
      endif

    end subroutine test_thang

  end module simple_mod_test

You can see that a list of numbers of processes is provided for each test. The
test will be run with each number of processes in turn.

The test can discover information about the parallel environment in which it is
running using the ``context`` member of the fixture class. This provides a
number of query functions and some tools allowing you to gather and all-reduce
over the process pool.

Also notice the use of the "feign" procedure to set up configuration for the
unit under test. This is discussed in a previous section.

Test Classes
------------

When the test fixture must hold data for the test cases a test class is used.

.. code-block::

  module simple_mod_test

    use pFUnit_Mod
    use simple_mod, only: thing, thang

    implicit none

    private
    public test_simple_type, test_thing, test_thang

    @TestCase
    type, public, extends(TestCase) :: test_simple_type
      real, allocatable :: data_block(:)
    contains
      procedure SetUp
      procedure tearDown
    end type test_simple_type

  contains

    subroutine setUp( this )

      implicit none

      class(test_simple_type), intent(inout) :: this

      allocate( this%data_block(256) )

    end subroutine setUp

    subroutine tearDown( this )

      implicit none

      class(test_simple_type), intent(inout) :: this

      deallocate( this%data_block )

    end subroutine tearDown

    @test
    subroutine test_thing( context )

      implicit none

      class(test_simple_type), intent(inout) :: context

      call thing( 1, 2, 3, context%data_block )
      @assertEqual( (/12, 13, 14/), context%data_block )

    end subroutine test_thing

    @test
    subroutine test_thang( context )

      implicit none

      class(test_simple_type), intent(inout) :: context

      integer :: result

      this%data_block = (/-1, -2, -3/)
      result = thang( context%data_block )
      @assertEqual( 13, result )

    end subroutine test_thang

  end module simple_mod_test

This is unnecessary for simple tests but can be essential for more complex ones.
It can also be used as a way to organise several sets of tests with different
fixtures in the same source file.

Parameterised Tests
-------------------

This is a more advanced topic but fear not, we have actually come across an
example all ready. MPI tests are parameterised; the test is called once for each
element in the array of number-of-processes. The number of processes is the
parameter.

In the case of MPI tests, the parameterisation is handled for you. However it
can often be useful to implement your own tests in this way. In particular where
you have a series of tests where the logic is identical but the input and output
vary.

Tests like these need a quick and neat way of running the same code with a
number of different conditions. That is what parameterised tests provide.

There is a fair bit of boilerplate but hopefully nothing too off-putting:

.. code-block::

  module simple_mod_test

    use pFUnit_Mod
    use simple_mod, only: thang

    implicit none

    private
    public get_parameters, test_thing, test_thang

    @testParameter
    type, public, extends(MPITestParameter) :: simple_parameter_type
      integer :: input(3)
      integer :: expected
    contains
      procedure :: toString
    end type simple_parameter_type

    @TestCase(npes=[1], testParameters={get_parameters()}, constructor=test_simple_constructor)
    type, public, extends(MPITestCase) :: test_simple_type
      private
      integer :: input(3)
      integer :: expected
      real, allocatable :: data_block(:)
    contains
      procedure setUp
      procedure tearDown
    end type test_simple_type

  contains

    function test_simple_constructor( test_parameter ) result( new_test )

      implicit none

      type(simple_parameter_type), intent( in ) :: test_parameter
      type(test_simple_type) :: new_test

      new_test%input    = test_parameter%input
      new_test%expected = test_parameter%expected

    end function test_simple_constructor

    function toString( this ) result( string )

      implicit none

      class( simple_parameter_type ), intent( in ) :: this
      character(:), allocatable :: string

      character(str_long) :: buffer

      write( buffer, '(3I3)') this%input
      string = trim( buffer )

    end function toString

    function get_parameters() result( parameters )

      implicit none

      type(simple_parameter_type) :: parameters(4)

      parameters = (/simple_parameter_type([0, 0, 0], 0),     &
                     simple_parameter_type([1, 2, 3], 7),     &
                     simple_parameter_type([-1, -1, -1], -1), &
                     simple_parameter_type([3, 2, 1], 14)/)

    end function get_parameters

    subroutine setUp( this )

     implicit none

     class(test_simple_type), intent(inout) :: this

     allocate( this%data_block(256) )

    end subroutine setUp

    subroutine tearDown( this )

      implicit none

      class(test_simple_type), intent(inout) :: this

      deallocate( this%data_block )

    end subroutine tearDown

    @test( npes=[1, 2, 4] )
    subroutine test_thang( fixture )

      implicit none

      class(test_simple_type), intent(inout) :: fixture

      integer :: result

      data_block = fixture%input
      result = thang( fixture%context%getMpiCommunicator(), fixture%data_block )
      @assertEqual( fixture%expected, result )

    end subroutine test_thang

  end module simple_mod_test

As you can see this allows a large number of cases to be tested without a lot of
additional code. It can also work in tandem with MPI testing meaning that each
test case will be run with each number-of-processes. Of course it works just as
well for serial tests if you derive from ``TestCase`` instead.

Despite its verbosity it should be fairly obvious what is going on. The only
thing which really needs explaining is the ``toString`` method on the parameter
type. This is used in the case of a test failure to identify which set of
parameters was in use. It should return a string which uniquely identifies the
case being run.

You may, of course, use the configuration feigning functions in parameterised
tests. They were removed from this example for clarity.

.. _canned_info:

Replacing Infrastructure calls with Canned Information
------------------------------------------------------

When testing kernels, it is often necessary to provide quantities that are quite
difficult to generate - such as dofmaps, basis functions and quadrature
information. When running the full model, these quantities are provided by the
infrastructure, but the infrastructure should not be used to generate them in
unit tests. Canned versions of these quantities are available from helper
routines held in the support directory ``components/science/unit-test/support``

Using Canned Information
^^^^^^^^^^^^^^^^^^^^^^^^

A description of the available canned-data support routines can be found at
:ref:`unit_test_canned_data_routines`.

.. note::

  All arrays returned by the helper routines are allocated inside the
  routines and so, will need to be deallocated in the calling routine when
  they are no longer required.

Adding New Canned Information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to make the use of canned data easier for the unit test writer, it is
important to try to follow a similar API to the already existing canned data
routines. The following rules should, therefore, be followed:

* As noted elsewhere in this page, calls to ``log_event()`` should be avoided.
  We should, therefore, avoid the situation where an error has to be trapped,
  and therefore, logged. For example, different canned dofmaps might be required
  for different orders of function space. A single canned-data support routine
  that takes the required order as an argument could be provided. But this would
  mean that if a non-supported order is passed in, it would have to be caught
  and reported. It would be better to provide separate functions for the
  supported orders. These would then only serve a single purpose and so, need no
  error checking. If anyone tried to call a routine that hadn't been written (a
  currently unsupported order), they would get a linker error.
* The canned-data support routines should have a single mandatory argument
  through which the canned data is returned. Any additional information used to
  control what data is returned should be made through optional arguments. The
  single mandatory argument should be an allocatable that is allocated by the
  support routine. This way, the user knows they always have to deallocate the
  data when they have finished with it. A mixture of some canned data that needs
  deallocating and some that doesn't will cause  confusion and lead to errors.
* Support routines that return the sizes of things (i.e.
  ``get_unit_test_..._sizes_mod.f90``) only return scalar integers and so no
  deallocation of any returned data is required.

Other Considerations
--------------------

There are a few remaining technical issues you might care to know about.

Verbose Output From Testing
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a test suffers an error (rather than a failure) then pFUnit does not always
provide clear details on where it had issues, i.e. a traceback of routines and
line numbers. It does provide the "-v" option to try and mitigate this
shortcoming somewhat.

This option enables "verbose output" which sends a message to the terminal when
it starts and ends each individual tests. This should at least provide a
starting point as to where to look for problems.

The LFRic build system provides the "VERBOSE" switch. When set using ``make test
VERBOSE=1``, it causes the build process to output much more information about
what it is doing. Part of that is to specify the "-v" argument to the unit
tests.

When the unit tests are run as part of the !Rose/Cylc test suites, this option
is specified by default.

Output to the terminal from unit tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Logging using ``log_mod``
"""""""""""""""""""""""""

Using the ``log_event()`` functionality  from ``log_mod`` within a unit test
should be avoided. Calling into ``log_mod`` adds another dependency into your
unit test, so if someone else puts an error in ``log_mod`` your seemingly
unrelated unit test will fail - this makes it harder to debug unit tests.

If you really have to call ``log_event()`` from a unit test (for example in the
unit tests on ``log_mod`` itself!), you must remember to call
``initialise_logging()`` from the unit test code, before you attempt to use
``log_event()``.

Standard Out
""""""""""""

Due to the internal workings of the pFUnit framework writing to standard out
from within tests should also be avoided.

When running in "robust" mode pFUnit actually uses standard out to communicate
between its supervising process and the tests it is running. Therefore,
injecting unexpected data into this stream can trip the framework up.

If you need a quick and easy way to get things to screen for debugging purposes
then there is a solution. If you prepend your debug message with "DEBUG: " then
the framework work will know not to interpret it as internal messaging. Instead
the text will be sent to your terminal.

The message will be inserted into the middle of the progress tally so it should
not be used permanently. However it should be sufficient for temporary debug
purposes.

None of this applies when not in "robust" mode but if you always do it then you
don't have to worry which mode is in use.

Technical Details
-----------------

Some technical details are outlined here, should they interest you.

Robust Mode or Non-robust Mode?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

pFUnit offers two modes of operation: "robust" and "non-robust".

In robust mode each test is run in a subprocess. This is the mode you want to
run in as it means an error (as opposed to a failure) in the test will not cause
the whole test suite to exit. It should also mean that each test is unable to
affect any other test.

Sadly robust mode is rather flawed at the moment. There is leakage between tests
which are supposed to be isolated. This leads to nasty side-effects between
tests.

Also MPI tests, which we use and need, do not work in robust mode.

There is also an odd behaviour whereby errors are reported, then you are told
there were no errors and the suite exits with "Okay".

In light of all this we have chosen to forgo robust mode until it can be made to
work properly and in parallel.

Test Ordering
^^^^^^^^^^^^^

The list of tests to run is automatically generated by the build system.

The way this list is compiled means that tests are run in i-node order. In other
words it is fairly stable on a machine but random between machines. Where there
is an unexpected dependency between tests this can cause a test to fail on one
machine but not on another.
