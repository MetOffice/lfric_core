Testing
=======

A number of different mechanisms and approaches are provided within the LFRic
build system for testing. Each is discussed here giving details on how to make
use of it and in which circumstances you may wish to do so.

.. note::
   The canonical version of this document is held as reStructured text in
   the repository at
   `source:LFRic/trunk/infrastructure/documentation/testing.rst`:trac:.
   Any changes in a branch which render this document inaccurate should also
   include updates to the version of this document on that branch. The version
   displayed on the wiki is generated from the head of trunk.

Unit Testing
^^^^^^^^^^^^

Unit testing is intended to exercise a discrete, well contained and defined,
"unit" of code. This will generally be a class, ("type" in Fortran parlance)
module or part thereof.

Tests of this kind are found under ``<project>/unit-test`` in a tree which
mirrors that in ``<project>/source``. This is intended to make it easy to find
the tests corresponding to a particular piece of functionality.

We use the pFUnit framework. Sourceforge holds some `documentation`__.

__ https://sourceforge.net/projects/pfunit/files/Documentation/pFUnit3-ReferenceManual.pdf

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
  ``<module under test>_test[_<subsection>]``
* Test procedure names should follow the form
  ``test_<procedure name>[_<subsection>]``
* Inputs to the procedure under test should be provided as constants. This
  prevents them changing unexpectedly. Do not rely on configuration values as
  these can change.
* If the unit under test makes use of configuration, use the "feign"
  functions provided by "feign_config_mod" to explicitly set the values
  needed.
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

Robust Mode or Non-robust Mode?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pFUnit offers two modes of operation: "robust" and "non-robust".

In robust mode each test is run in a subprocess. This is the mode you want to
run in as it means an error (as opposed to a failure) in the test will not
cause the whole test suite to exit. It should also mean that each test is
unable to affect any other test.

Sadly robust mode is rather flawed at the moment. There is leakage between
tests which are supposed to be isolated. This leads to nasty side-effects
between tests.

Also MPI tests, which we need now and hope to implement soon, do not work in
robust mode.

There is also an odd behaviour where errors are reported, then you are told
there were no errors and the suite exits with "Okay".

In light of all this we have chosen to forgoe robust mode until it can be made
to work properly and in parallel.

Use of pFUnit VERBOSE option
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pFUnit does not always provide clear details on where it had issues, i.e. a 
traceback of routines and line numbers. pFUnit does provide the "-v" option
which sends a message to stdout when it starts and ends individual tests. This
should at least provide a starting point as to where to looking for problems,
especially when the issue is testing order dependant (which can differ on
different platforms).

The VERBOSE option is set to true by default for pFUnit tests on rose
test-suites. If you are running pFUnit tests from your desktop use
``make test VERBOSE=1``.

Integration Testing
^^^^^^^^^^^^^^^^^^^

Some things can not be easily tested using the unit testing framework. For
instance things which halt execution or interact with command line arguments.
Others rely too heavily on additional components to be considered "unit" tests.

In these cases integration tests may be appropriate.

These tests live in ``<project>/integration-test``. Each test comprises two
parts:

* A Fortran source file which implements a ``PROGRAM`` which may be
  compiled and linked against model code
* A Python script which executes the program and validates the results

A simple framework is provided to aid in writing the test script. Its use is
shown in this example::

    #!/usr/bin/env python2.7

    from __future__ import print_function

    from testframework import Test, TestEngine, TestFailed

    class simple_test( Test ):
      def test( self, return_code, out, err ):
        question = 'How much wood would a Woodchuck chuck if a Woodchuck could chuck wood'
        expectedTranscript = 'A Woodchuck would chuck no amount of wood as a Woodchuck can't chuck wood'

        if return_code != 0:
          raise TestFailed( 'Failed to execute' )

        if err != '':
          raise TestFailed( 'Expected no output on standard error' )

        if out != expectedTranscript:
          raise TestFailed( 'Expected >{}< but found >{}<'.format( expectedTranscript, out ) )

        return 'The question answerer understands the nature of Woodchucks'

    if __name__ == '__main__':
      TestEngine.run( simple_test )

Parallel Testing
~~~~~~~~~~~~~~~~

The integration testing framework offers two additional base classes for
parallel testing.

MpiTest
-------

Tests derived from this class will be run using ``mpiexec`` with the
constructor taking the number of processes to use. The class deals with
stripping all the extranious output generated by the MPI system itself so the
test receives only clean output to examine.

Be aware that some systems have faulty MPI installations such that ``mpiexec``
doesn't work. It is preferred as it is the standard but if you are unable to
use it simple set the ``MPIEXEC_BROKEN`` environment variable and ``mpirun``
will be used instead.

LFRicLoggingTest
--------

This class offers further support for tests linked to LFRic Logging. When
run in parallel such tests will tend to send their output to a PET log file
rather than standard out/error. The constructor of this class takes the base
name for these log files and makes their contents available after execution.

Just call ``getLFRicLoggingLog()`` to retrieve the log for a particular process. Be
aware that this class stores the log files in memory so don't let them get
too big. This should not be a problem as tests should not be generating that
much output.

Test Suite
^^^^^^^^^^

`Rose <http://metomi.github.io/rose/doc/rose.html>`_ and
`Cylc <http://cylc.github.io/cylc/>`_ are used to run a suite of tests. There
are a number of goals for these tests.

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

The test suite is invoked with ``make test-suite``. This will launch the test
suite of each sub-project listed in the ``OPERATE_ON`` environment variable.
Each test suite will be launched for each platform listed in the
``TEST_SUITE_TARGETS`` environment variable. Thus the test suite for each
sub-project appearing in ``OPERATE_ON`` must have optional configurations for
each platform appearing in ``TEST_SUITE_TARGETS``.

If ``OPERATE_ON`` is not set a default list will be used. There is not default
for ``TEST_SUITE_TARGETS`` but it will be set up by the Met Office module
system.

`Note:` A clean test-suite can be forced using the ``PURGE_SUITES`` keyword
which will purge any rose suites in the ``~/cycl-run`` directory with the
same name as your suite before launching a new gui
i.e. ``make test-suite PURGE_SUITES=1``.

During development it is often useful to target a single platform.

If the test suite is invoked using ``rose stem --config <project>/rose-stem``
this will target the default platform. Out-of-the-box this is the Met Office
SPICE server farm. To change the default target modify ``opts = meto-desktop``
in ``rose-suite.conf``.

Developers and reviewers should remember not to allow changes to this default
on to trunk.

Most of the time it is more convenient to leave the default as it is and just
specify the desired target on the command line. e.g.
``rose stem --config <project>/rose-stem --opt-conf-key=meto-xc40``

Other target configurations available can be found in ``<project>/rose-stem/opts/``

Test Suite on MONSooN
~~~~~~~~~~~~~~~~~~~~~

Met Office users can launch the test suite on MONSooN from their workstations
but external users must use "exvmsrose", within the MONSooN subnet. Ensure you
have suitable, unprotected, keys set up to allow you to connect from
"exvmscylc" back to "exvmsrose" and to "xcslc0".

Loading the module "meto-common-environment/dynamo" on exvmsrose will set
everything up so you need only type ``make test-suite``.

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
stem --group=nightly``. Alternatively using ``make test-suite
SUITE_GROUP=nightly`` will run the nightly tests on all platforms usually
targeted by ``make test-suite``.

Nightly tests will publish various reports to your
``~/public_html/lfric-<project>-<target>`` (where "target" will be a platform
listed in the ``TEST_SUITE_TARGETS``) directory.

By default nightly tests will be run on the Met Office SPICE server farm.
To run on Cray you need to specify the target explicitly using
``rose stem --group=nightly --opt-conf-key=meto-xc40``.

If the optional configuration key ``meto-mail-list`` is used then the LFRic
mailing list will be sent a message for each failing task in the test suite.
This option is intended for use by timed runs of the suite only. If you want to
be sent e-mails while running the suite interactively be careful. You probably
don't want to send them to the developer mailing list.

If you want to run the nightly tests as a "Cron" job the following ``crontab``
should give you some pointers::

  MAILTO = joe.bloggs@metoffice.gov.uk
  PATH = /opt/ukmo/utils/bin:/usr/local/sci/bin:/usr/local/bin:/usr/bin:/bin
  DATADIR = /data/users/jbloggs
  TMPDIR = /var/tmp
  ROSE_BUSH_URL = http://fcm1/rose-bush
  TROUBLE_EMAIL_ADDRESS=meto-lfric@metoffice.gov.uk
  # Minute # Hour # Day of Month # Month # Day of Week # Command
  21 03 * * 1-5 $HOME/local/bin/run-nightly-test

Where the ``run-nightly-test`` script looks like this::

  #!/bin/sh

  rm -rf /data/local/jbloggs/lfric-nightly
  fcm checkout fcm:lfric.xm-br/dev/joebloggs/r9999_SpecialBranch /data/local/jbloggs/lfric-nightly
  rose stem --group=nightly --source=/data/local/jbloggs/lfric-nightly --no-gcontrol

Testing build with UM Physics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LFRic is currently set to be able to compile with the UM physics codes in order that science 
development can continue with a single code source. This involves extracting the UM code with 
fcm and then building with the native LFRic build system.  A test is provided to ensure this 
build process, but since the analysis of the UM code takes considerably longer than the GungHo
core alone, it is currently separated from the standard developer tests. To run this test use
``make test-umphysics``. 


Test Suite Development
~~~~~~~~~~~~~~~~~~~~~~

A Cylc suite consists of an acyclic graph of "tasks". Each task is performed by
an "application" which is a discrete executable unit. The test suite is held
in the "rose-stem" directory.

Adding new application configuration tests
------------------------------------------

It is hoped that adding new configuration tests should be as simple as you would
hope.

* You will need to add a new optional configuration for your test to
  ``rose-stem/app/<application_name>/opt/rose-app-<configuration_name>.conf``. The
  contents of this file will override the basic configuration in
  ``rose-stem/app/<application_name>/rose-app.conf``.

* Now add the name of this new test configuration to the appropriate application in the
  ``application_configurations`` list at the top of "suite_control.rc". **Note:** if you
  are adding a new application, it should have at least a ``'default'`` configuration
  listed.

* Finally decide which group or groups you would like your test to appear in and
  add it to them in the ``groups`` structure just below. The form to use is
  ``'run_application(<application_name>, <configuration_name>, *keywordArguments)'``.

Configuration tests by default do a known good output (kgo) comparison so you may
also need to provide a kgo file
``rose-stem/app/check_application_kgo/file/<application_name><platform>/checksum.<application_name>_<configuration_name>_<resolution>.<compiler>.kgo.txt``

**Note:** That the kgo filename should include the ``'<resolution>'`` label if a
spatial/temporal resolution has been given for the test configuration (see keyword
arguments for the ``'run_applications'`` macro).

From that point everything should be taken care of for you.

 

Adding new tasks
----------------

This requires a little more work. The basic steps to add a task 'new_task' are
as follows:

* Add a new app directory ``rose-stem/app/<new_task>``
* Populate the ``rose-stem/app/<new_task>/rose-app.conf`` file to add desired
  commands (e.g. to run a script)
* Add a new options directory ``rose-stem/app/<new_task>/opt`` with your predefined
  configurations for your app. These should be named ``rose-app-<configuration_name>.conf``,
  the opt directrory should have at least the default configuration for your application, named
  ``rose-app-default.conf``, these can be an empty file.
* If your task runs a script, then put the script in
  ``rose-stem/app/<new_task>/bin``
* In ``rose-stem/suite_control.rc`` file,  add ``new_task`` to the appropriate ``groups`` structure near the top of
  the file.
* In ``rose-stem/inc/suite_macros.rc``, add a Jinja2 macro with the same name as the task, add in any task
  dependencies::

    {%- macro new_task() %}
      some_other_task => new_task
    {%- endmacro %}

  in this example ``new_task`` will run after ``some_other_task`` completes.
* In ``rose-stem/suite.rc``, add Jinja2 code to define the actual task::

    {%- if 'new_task' in scheduledTasks %}
      [[new_task]]
        inherit = None, LOCAL
        script = rose task-run --app-key=<new_task> --command-key=<command>
        [[[environment]]]
          # set up any environment variables here
    {%- endif %}

  where ``<new_task>`` is the app name and ``<command>`` is the name of a
  command in ``rose-stem/app/<new_task>/rose-app.conf``.


suite.rc
--------

The heart of the suite is the "suite.rc" file. Once processed by Jinja2,
this describes the tasks and the graph which orders them.
Each task sets up an execution environment for the application it invokes.
This means multiple tasks may invoke the same application to gain different
effects.

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

The macro ``schedule`` (found in ``rose-stem/inc/suite_macros.rc``) is called
to generate the graph. This loops over each group of tests to be run. For each
group it loops over the list of tests and calls a similarly named macro.
These macros generate the set of graph nodes required to perform that test.

The generating macros do not check to ensure that they are not duplicating
tasks already in the schedule. Thus it is important to deduplicate it and this
is the final stage before the list makes it into the processed "suite.rc".

Having generated the schedule each task must be described. All possible task
definitions are generated by the scripting and are then gated by the schedule
list. If a task isn't actually in the schedule it shouldn't end up in the
processed file.

Configuring application tests
-----------------------------
Configuration on application tests for the test suite is provided by the 
file ``rose-stem/inc/suite_control.rc``. As described earlier, application
tests are grouped together under a group name for a suite to run,
**e.g.** *developer*. How each application test is configured within a group
is specified using the ``run_application`` macro. The ``run_application``
macro requires at least two arguments, followed by optional
arguments, e.g.::

  {%- set groups = {'<group_1>': [ 'run_application( <application_name>, <configuration_name>, *keywordArguments )',
                                   ... ],
                    '<group_2>': [ 'run_application( <application_name>, <configuration_name>, *keywordArguments )',
                                   ... ],
                    ... } %}

where <**application_name**> is the name of the application in the ``gungho/rose-stem/app/`` directory.
      <**configuration_name**> is one of the configuration options for <**application_name**>
      given in the ``application_configurations`` jinja2 variable at the top of of the ``suite.rc`` file.

accepted keyword arguments are:

* **resolutions**: List of tuples specifying the spatial and time resolutions to run this application configuration at.
* **run_with_compilers**: List of strings giving compilers to run this application configuration at each **resolution** (if applicable). The compiler keys available are given by the ``TARGET_COMPILER_SETUP`` jinja2 variable. If this keyword is not provided, the task will be run with the compiler given by the jinja2 variable ``TARGET_SCIENCE_COMPILER``
* **support meshes**: Meshes than must be generated before the application configuration can be run, i.e. multimesh configurations.
* **env**: Dictionary specifying environment variables for the task
* **checkkgo**: Whether the output of the task requires a check against *Known Good Output*, default=True
* **plotstr**: String input which is a command line call to a require plotting routine.

The keyword arguments, **resolutions** and **env** follow the CSAR functionality for the rose stem suite, details of how to specify these keyword arguments are given `here <https://code.metoffice.gov.uk/trac/lfric/wiki/GungHo/RoseStemCsarFunctionality>`_ 

Suite Style Guide
-----------------

This guide is based on common practice within the Met Office.

Indent your suite.rc file. Each section level should be indented with respect
to the parent section. Following LFRic standard we use 2 spaces per indent
level.

Capitalise group task names. Non executive tasks which exist to group other
tasks should have their name in all caps.

Other useful tips
-----------------

It is not always obvious where stuff is located, for example, if you want to
copy files to another location at the end of a task. The
``$CYLC_TASK_WORK_DIR`` environment variable can be used for the current
working directory (e.g. where model output goes after a run). Other environment
variables can be found in any Cylc job output log.

You can define pre- and post-script commands for tasks to do things like copy
files or create directories. There are many examples of these in the LFRic
suite.rc.
