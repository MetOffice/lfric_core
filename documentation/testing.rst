Testing
=======

A number of different mechanisms and approaches are provided within the Dynamo
build system for testing. Each is discussed here giving details on how to make
use of it and in which circumstances you may wish to do so.

Unit Testing
^^^^^^^^^^^^

Unit testing is intended to exercise a discrete, well contained and defined,
"unit" of code. This will generally be a class, ("type" in Fortran parlance)
module or part thereof.

Tests of this kind are found under ``src/unit-test`` in a tree which mirrors
that in ``src/dynamo``. This is intended to make it easy to find the tests
corresponding to a particular piece of functionality.

We use the pFUnit framework, documentation for which may be found in
``src/pfunit/documentation``.

Tests are identified with the ``.pf`` extension while support code in normal
``.[fF]90`` files may also be included. ``.pf`` files are normal Fortran files
garnished with pFUnit directives. These are preceded with an "at symbol".

A minimal example might look like this::

    module simple_mod_test

    use pFUnit_Mod
    use simple_mod, only: thing, thang

    implicit none

    private
    public test_simple_thing

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

There are a number of things to note here:

  * The module name should follow the form
    ``<module under test>_test[__<subsection>]``
  * Test procedure names should follow the form
    ``test_<procedure name>[__<subsection>]``
  * Inputs to the procedure under test should be provided as constants. This
    prevents them changing unexpectedly. Do not rely on configuration values as
    these can change.
  * If the unit under test makes use of configuration, use the "frig" functions
    provided by "frig_config_mod" to explicitly set the values needed.
  * The expected result given to assertions should be constants or calculated
    from constants. Try to avoid calculations so as to avoid the possibility
    that the same (flawed) algorithm is used in test and unit under test.
  * Each test procedure should contain a single well defined set of tests. This
    is usually for a single procedure under test.

Fully Fledged Test Class
~~~~~~~~~~~~~~~~~~~~~~~~

The simple test outlined above is fine for situations where there are no
resources to be managed. When there are, a full test class is needed.

The most common resource requiring management is memory, allocated space must
be deallocated, but configuration may also be considered a resource in this
case.

These managed resources are referred to as "fixtures" in testing parlance. They
are things which must exist in order to perform the test but which are not
under test themselves. They are created and initialised before each test and
destroyed after it. This means that each test has a pristine, known,
environment in which to work.

Following is a minimal example of a full test class with example fixtures::

    module simple_mod_test

    use pFUnit_Mod
    use simple_mod, only: thing, thang

    implicit none

    private
    public test_simpe_type, test_thing, test_thang

    @TestCase
    type, extends(TestCase) :: test_simple_type
        private
        real :: data_block(:)
    contains
        procedure setUp
        procedure tearDown
        procedure test_thing
        procedure test_thang
    end type test_simple_type

    contains

    subroutine setUp( this )

        use feign_config_mod, only :: feign_stuff_config

        implicit none

        class(test_simple_type), intent(inout) :: this

        call feign_stuff_config( 'foo', 13 )

        allocate( this%data_block(256) )

    end subroutine setUp

    subroutine tearDown( this )

        implicit none

        class(test_simple_type), intent(inout) :: this

        deallocate( this%data_block )

    end subroutine tearDown

    @test
    subroutine test_thing( this )

        implicit none

        class(test_simple_type), intent(inout) :: this

        call thing( 1, 2, 3, data_block )
        @assertEqual( (/12, 13, 14/), data_block )

    end subroutine test_thing

    @test
    subroutine test_thang( this )

        implicit none

        class(test_simple_type), intent(inout) :: this

        integer :: result

        data_block = (/-1, -2, -3/)
        result = thang( data_block )
        @assertEqual( 13, result )

    end subroutine test_thang

    end module simple_mod_test

Notice how ``data_block`` is allocated in ``setUp`` and deallocated in
``tearDown``. These procedures are called before and after each test method.
This means that the contents do not carry between tests. Therefore they must be
initialised for each test. This may be done in ``setUp`` if every test needs
the same initial condition or locally to the test if they need different
starting points.

Also notice the use of the "feign" procedure to set up configuration for the
unit under test. Bear in mind that you must feign everything the code being
tested will use. You can not rely on them having been set up by a previous test
as they have not been. You can not rely on them defaulting to a particular
value because they do not.

If the unit calls down to help procedures any configuration they make use of
must also be feigned.

Problems With pFUnit
~~~~~~~~~~~~~~~~~~~~

Currently there is an issue with pFUnit whereby it will report finding errors,
then tell you no errors were found and exit with "Okay".

This is an issue with running in "robust mode". In this mode each test is run
in a subprocess. This is the mode you want to run in as it means an error (as
opposed to a failure) in the test will not cause the whole test suite to exit.
It should also mean that each test is unable to affect any other test. Sadly
this is not currently the case as there is leakage between tests.

To find out what is failing you can run the tests by hand in non-robust mode.
To do this change to ``src/unit-test`` and run ``../../build/unit-test/test``.

Functional Testing
^^^^^^^^^^^^^^^^^^

Some things can not be easily tested using the unit testing framework. For
instance things which halt execution or interact with command line arguments.
In these cases "functional" tests (which may need a more suitable group name)
may be appropriate.

These tests live in ``src/functional-test``. Each test comprises two parts::

    * A Fortran source file which implements a ``PROGRAM`` which may be
      compiled and linked against Dynamo code
    * A Python script which executes the program and validates the results

A simple framework is provided to aid in writing the test script. Its use is
shown in this example::

    #!/usr/bin/env python2.7

    from __future__ import print_function

    import testframework import Test, TestEngine, TestFailed

    class simple_test( Test ):
    def test( self ):
        question = 'How much wood would a Woodchuck chuck if a Woodchuck could chuck wood'
        expectedTranscript = 'You asked "{}"\nI answer "A Woodchuck would chuck no amount of wood as a Woodchuck can't chuck wood"'.format( question)

        out, err = process.communicate( question )

        if process.returncode != 0:
        raise TestFailed( 'Failed to execute' )

        if err != '':
        raise TestFailed( 'Expected no output on standard error' )

        if out != expectedTranscript:
        raise TestFailed( 'Expected >{}< but found >{}<'.format( expectedTranscript, out ) )

        return 'The question answerer understands the nature of Woodchucks'

    if __name__ == '__main__':
      TestEngine.run( simple_test )

Test Suite
^^^^^^^^^^

`Rose <http://metomi.github.io/rose/doc/rose.html>`_ and
`Cylc <http://cylc.github.io/cylc/>`_ are used to run a suite of tests. There
are two goals for these tests.

  * Ensure the code remains buildable with our supported compilers
  * Ensure the code is tidy
  * Ensure changes to the code do not change the results

The first goal is pursued by automating the process of building with each of
the compilers. The generation of API documentation is also included as this can
sometimes break as well.

Code tidiness is currently pursued only as far as ensuring there are no unused
source files. Other tests may be forthcoming.

The result of runs with known inputs are compared against known good outputs to
check for changes.

The test suite is invoked with ``make test-suite``. This will launch the suite
multiple times, once for each platform listed in the
``DYNAMO_TEST_SUITE_TARGETS`` environment variable. This variable is set up for
you by the LFRic module system.

During development it is often useful to target a single platform. If the test
suite is invoked using ``rose stem`` this will target the default platform.
Out-of-the-box this is the Met Office SPICE server farm. To change the default
target modify ``opts = meto-desktop`` in "rose-stem/rose-suite.conf".

Developers and reviewers should remember not to allow changes to this default
on to trunk.

Test Suite on MONSooN
~~~~~~~~~~~~~~~~~~~~~

The test suite can only be launched from the machine "exvmsrose", within the
MONSooN subnet. Ensure you have suitable, unprotected, keys set up to allow you
to connect from "exvmscylc" back to "exvmsrose" and to "xcml00".

Loading the module "meto-common-environment/dynamo" will set everything up so
you need only type ``make test-suite``.

Nightly Testing with the Test Suite
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are a number of tests which we want to run regularly but which consume
too many resources (especially run-time) to include in the developer test
suite. For this reason the test suite is actually split into a number of
groups.

As of the time of writing available groups are:

- developer
- nightly

Just running ``rose stem``, without specifying a group, causes the developer
tests to be run. If you want to run the nightly tests you must use ``rose
stem --group=nightly``.

These tests will publish various reports to your
``~/public_html/dynamo-<target>`` (where "target" is as described for
``DYNAMO_TEST_SUITE_TARGETS``) directory. If this does not exist you will see
suite failures.

If you want to run the nightly tests as a "Cron" job the following ``crontab``
should give you some pointers::

  MAILTO = joe.bloggs@metoffice.gov.uk
  PATH = /opt/ukmo/utils/bin:/usr/local/sci/bin:/usr/local/bin:/usr/bin:/bin
  DATADIR = /data/users/jbloggs
  TMPDIR = /var/tmp
  ROSE_BUSH_URL = http://fcm1/rose-bush
  # Minute # Hour # Day of Month # Month # Day of Week # Command
  21 03 * * 1-5 $HOME/local/bin/run-nightly-test

Where the ``run-nightly-test`` script looks like this::

  #!/bin/sh

  rm -rf /data/local/jbloggs/dynamo-nightly
  fcm checkout fcm:dynamo.xm-br/dev/joebloggs/r9999_SpecialBranch /data/local/jbloggs/dynamo-nightly
  rose stem --group=nightly --source=/data/local/jbloggs/dynamo-nightly --no-gcontrol

Test Suite Development
~~~~~~~~~~~~~~~~~~~~~~

A Cylc suite consists of an acyclic graph of "tasks". Each task is performed by
an "application" which is a discrete executable unit. The test suite is held
in the "rose-stem" directory.

Adding new science tests
------------------------

It is hoped that adding new science tests should be as simple as you would
hope.

You will need to add a new optional configuration for your test to
``rose-stem/app/dynamo/opt/rose-app-<test name>.conf``. The contents of this
file will override the basic configuration in
``rose-stem/app/dynamo/rose-app.conf``.

Now add the name of this new test to the ``science_configurations`` list at
the top of "suite.rc".

Finally decide which group or groups you would like your test to appear in and
add it to them in the ``groups`` structure just below. The form to use is
``'science(<test name>)'``.

Science tests by default do a kgo comparison so you also need to provide a file
``rose-stem/app/check-dynamo/file/<platform>/checksum.<test name>.<compiler>.kgo.txt``

From that point everything should be taken care of for you.

Adding new tasks
----------------

This requires a little more work. The basic steps to add a task 'new_task' are as follows:

* Add a new app directory ``rose-stem/app/<new_task>``
* Populate the ``rose-stem/app/<new_task>/rose-app.conf`` file to add desired commands 
  (e.g. to run a script)
* If your task runs a script, then put the script in ``rose-stem/app/<new_task>/bin``
* Add ``new_task`` to the appropriate ``groups`` structure near the top of the suite.rc
* Add a Jinja2 macro with the same name as the task, add in any task dependencies::

    {%- macro new_task() %}

       some_other_task => new_task

    {%- endmacro %}

in this example ``new_task`` will run after ``some_other_task`` completes.

* Add Jinja2 code to define the actual task::

    {%- if 'new_task' in scheduledTasks %}

        [[new_task]]
            inherit = None, LOCAL

            script = rose task-run --app-key=<new_task> --command-key=<command>

            [[[environment]]]

            # set up any environment variables here
 
    {%- endif %}

where ``<new_task>`` is the app name and ``<command>`` is the name of a command in
``rose-stem/app/<new_task>/rose-app.conf``. 


suite.rc
--------

The heart of the suite is the "suite.rc" file. This describes the tasks and the
graph which orders them. Each task sets up an execution environment for the
application it invokes. This means multiple tasks may invoke the same
application to gain different effects.

The file itself uses a modified version of the Windows INI file syntax.
Modified to support hierarchical data which is indicated by multiple square
brackets on section headers.

Before it is interpreted by Cylc the suite.rc file is processed using the
Jinja2 template engine. Directives to this process are identified by the use of
curly brackets.

Not all tasks invoke an application. Sometimes called "group tasks" they exist
as a means to consolidate common environment details for multiple "child"
tasks. Be aware, however that this inheritance is only by replacement, not by
extension. This makes sense for environment variables as it means that a child
task's definition of a specific variable will override any in a parent. It is
problematic with command scripting though since it is not possible to have the
parent's script run, then the child.

The task graph is compiled using Jinja scripting. Sadly this rather obfuscates
what is happening but an outline is given below.

The macro ``schedule`` is called to generate the graph. This loops over each
group of tests to be run. For each group it loops over the list of tests and
calls a similarly named macro. These macros generate the set of graph nodes
required to perform that test.

The generating macros do not check to ensure that they are not duplicating
tasks already in the schedule. Thus it is important to deduplicate it and this
is the final stage before the list makes it into the processed "suite.rc".

Having generated the schedule each task must be described. All possible task
definitions are generated by the scripting and are then gated by the schedule
list. If a task isn't actually in the schedule it shouldn't end up in the
processed file.

Suite Style Guide
-----------------

This guide is based on common practice within the Met Office.

Indent your suite.rc file. Each section level should be indented with respect
to the parent section.

Capitalise group task names. Non executive tasks which exist to group other
tasks should have their name in all caps.

Other useful tips
-----------------

It is not always obvious where stuff is located, for example, if you want to copy files
to another location at the end of a task. The ``$CYLC_TASK_WORK_DIR`` environment
variable can be used for the current working directory (e.g. where model output
goes after a run). Other environment variables can be found in any Cylc job output
log.

You can define pre- and post-script commands for tasks to do things like copy
files or create directories. There are many examples of these in the LFRic suite.rc.   
  
