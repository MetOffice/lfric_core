.. -----------------------------------------------------------------------------
     (c) Crown copyright 2025 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

Testing Algorithms
==================

Science kernels are an ideal candidate for unit testing, being small and self
contained with well understood input and output. On the other hand algorithms
are too unwieldy requiring substantial supporting set-up and covering too many
interacting sub-components.

This is the sort of situation integration testing is well suited.

Support for Tests
~~~~~~~~~~~~~~~~~

One of the things which makes testing algorithms difficult is the fact that
they operate on field objects. This implies a substantial amount of LFRic
infrastructure in order to create these fields.

In order to make life easier some support code is provided in
``gungho/integration-test/support``.

Tiny World
----------

The file ``support/test_tiny_world_mod.f90`` provides
``intialise_tiny_world()``. This sets up everything needed to create fields on
a C3 mesh. ``finalise_tiny_world()`` is provided to clean up afterwards.

Field Data Access
-----------------

Once some fields have been created it will be necessary to interact with them.
This runs counter to the O.O. principles of data hiding and separation of
concerns. As such ``support/test_field_mod.f90`` provides a clean interface
allowing data to be written to and read from fields without breaking
encapsulation. It provides ``test_field_xx_type`` where `xx` is 32 or 64
denoting the precision of the field.

These test fields map an LFRic field to an array of floating-point values which
may be accessed using ``test_field_type%get_data()``. Once the array has been
set up the way you like you may call ``test_field_type%copy_to_lfric()`` which
will copy the content into the field.

Once the algorithm under test has been called any output fields may have their
``test_field_type%copy_from_lfric()`` method called to copy the field data into
the associated array. This may then be examined for correct behaviour.

Helping Hands
-------------

When setting up the input fields for an algorithm there are a few initial
patterns which are very common. The module in ``support/test_helpers_mod.f90``
provides convenience procedures to support these common patterns.

Creating Fields
...............

The subroutine ``create_test_field()`` takes a function space identifier and
name. It will return an LFRic field constructed with those parameters and a
test field associated with it.

If no other arguments are given the field will be initialised with a simple
cell count. i.e. the first cell holds 1, the second 2 and so on. There is no
gurantee over which cell is associated with which number but each number will
appear only once.

If an initialisation value is passed to the subroutine the field will be
flooded with this value.

Printing Fields
...............

In order to give the Python script visibility of field content a
``print_test_field`` subroutine is provided. This prints out the contents of a
field in a fixed format. Obviously care is needed with this as it is printing
a full 3D field so even a C3 mesh generates a lot of output.

Command Line Arguments
----------------------

If your test requires interaction with the command-line then
``support/test_cli.f90`` provides ``get_cli_argument()`` which takes an index
(default the first argument) and returns a string holding that argument.

Implementing Tests
~~~~~~~~~~~~~~~~~~

As with any other integration test there will be a Fortran program which
creates any necessary fields, calls the algorithm, then prints suitible
diagnostic quantities to standard out. This will be received by the Python
script which checks for correctness and reports success or failure.

Use of the support utilities should make writing your test considerably easier.

Example
-------

To illustrate this, consider the following example from
``gungho/integration-test/algorithm/transport/ffsl/ffsl_advective_updates_alg_mod.f90``:

.. code-block:: fortran

  program test_swift_outer_update_tracer

    use, intrinsic :: iso_fortran_env,  only : output_unit

    use constants_mod,                  only : r_second
    use ffsl_advective_updates_alg_mod, only : swift_outer_update_tracer
    use fs_continuity_mod,              only : W3
    use function_space_mod,             only : function_space_type
    use function_space_collection_mod,  only : function_space_collection
    use test_helpers_mod,               only : create_test_field, print_test_field
    use test_tiny_world_mod,            only : initialise_tiny_world, &
                                               finalise_tiny_world,   &
                                               mesh

    implicit none

    character(:), allocatable :: option

    call initialise_tiny_world()
    call run_test_64()
    call finalise_tiny_world()

Here we see the use of "tiny world" to avoid a lot of boiler-plate.

.. code-block:: fortran

  contains

    subroutine run_test_64()

      use, intrinsic :: iso_fortran_env, only : real64

      use field_real64_mod, only : field_real64_type
      use test_field_mod,   only : test_field_64_type

      implicit none

      type(function_space_type), pointer :: w3_fs

      type(field_real64_type),  allocatable, target :: post_x, post_y
      type(test_field_64_type), allocatable         :: post_x_test, post_y_test
      type(field_real64_type),  allocatable, target :: tracer_field
      type(test_field_64_type), allocatable         :: tracer_test_field
      type(field_real64_type),  allocatable, target :: dry_mass_next
      type(test_field_64_type), allocatable         :: dry_mass_next_test
      type(field_real64_type),  allocatable, target :: dry_mass_x, dry_mass_y
      type(test_field_64_type), allocatable         :: dry_mass_x_test, dry_mass_y_test
      type(field_real64_type),  allocatable, target :: inc_x, inc_y
      type(test_field_64_type), allocatable         :: inc_x_test, inc_y_test

These are the fields needed by the algorithm to be tested and their associated
``test_field_type`` objects.

.. code-block:: fortran

      w3_fs => function_space_collection%get_fs(mesh, 0, 0, W3)

      call create_test_field( w3_fs, "tracer", &
                              tracer_field, tracer_test_field, 1.0_real64 )

      call create_test_field( w3_fs, "post x", post_x, post_x_test, 0.2_real64 )
      call create_test_field( w3_fs, "post y", post_y, post_y_test, 0.4_real64 )

These fields are initialised by flooding them with the specified value.

.. code-block:: fortran

      call create_test_field( w3_fs, "dry mass next", &
                              dry_mass_next, dry_mass_next_test )

This field is initialised with cell count numbers.

.. code-block:: fortran

      call create_test_field( w3_fs, "dry mass x", &
                              dry_mass_x, dry_mass_x_test, 0.3_real64 )
      call create_test_field( w3_fs, "dry mass y", &
                              dry_mass_y, dry_mass_y_test, 0.5_real64 )

      call create_test_field( w3_fs, "incrament x", &
                              inc_x, inc_x_test, 0.7_real64 )
      call create_test_field( w3_fs, "incrament y", &
                              inc_y, inc_y_test, 0.9_real64 )

      call swift_outer_update_tracer( tracer_field,           &
                                      post_x, post_y,         &
                                      dry_mass_next,          &
                                      dry_mass_x, dry_mass_y, &
                                      inc_x, inc_y,           &
                                      dt=0.5_r_second )

      call print_test_field( tracer_test_field )

Finally the algorithm is called and its output field printed.

.. code-block:: fortran

    end subroutine run_test_64

  end program test_swift_outer_update_tracer

Meanwhile, in
``gungho/integration-test/algorithm/transport/ffsl/ffsl_advective_updates_alg_mod.py``:

.. code-block:: python

  from re import compile as re_compile
  from sys import argv as sys_argv
  from testframework import MpiTest, TestEngine, TestFailed
  from typing import Dict

  class SwiftOuterUpdateTracerTest(MpiTest):

      def __init__(self):
          super().__init__([sys_argv[1]], processes=6)
          self.__pattern = re_compile(r'tracer \((\d+)\) = (.*)')

This test will be run in parallel to better simulate its use int the model. It
will use 6 processes.

The regular expression is used to locate field data in standard output.

.. code-block:: python

      def test(self, return_code: int, out: str, err: str) -> str:
          expected = {
              0: [-0.270, -0.135, -0.090, -0.068, -0.054, -0.045, -0.039, -0.034,
                  -0.030, -0.027, -0.025, -0.023, -0.021, -0.019, -0.018, -0.017,
                  -0.016, -0.015, -0.014, -0.014, -0.013, -0.012, -0.012, -0.011,
                  -0.011, -0.010, -0.010],
              1: [-0.135, -0.090, -0.068, -0.054, -0.045, -0.039, -0.034, -0.030,
                  -0.027, -0.025, -0.023, -0.021, -0.019, -0.018, -0.017, -0.016,
                  -0.015, -0.014, -0.014, -0.013, -0.012, -0.012, -0.011, -0.011,
                  -0.010, -0.010, -0.010],
              2: [-0.090, -0.068, -0.054, -0.045, -0.039, -0.034, -0.030, -0.027,
                  -0.025, -0.023, -0.021, -0.019, -0.018, -0.017, -0.016, -0.015,
                  -0.014, -0.014, -0.013, -0.012, -0.012, -0.011, -0.011, -0.010,
                  -0.010, -0.010, -0.009],
              3: [-0.068, -0.054, -0.045, -0.039, -0.034, -0.030, -0.027, -0.025,
                  -0.023, -0.021, -0.019, -0.018, -0.017, -0.016, -0.015, -0.014,
                  -0.014, -0.013, -0.012, -0.012, -0.011, -0.011, -0.010, -0.010,
                  -0.010, -0.009, -0.009],
              4: [-0.054, -0.045, -0.039, -0.034, -0.030, -0.027, -0.025, -0.023,
                  -0.021, -0.019, -0.018, -0.017, -0.016, -0.015, -0.014, -0.014,
                  -0.013, -0.012, -0.012, -0.011, -0.011, -0.010, -0.010, -0.010,
                  -0.009, -0.009, -0.009],
              5: [-0.045, -0.039, -0.034, -0.030, -0.027, -0.025, -0.023, -0.021,
                  -0.019, -0.018, -0.017, -0.016, -0.015, -0.014, -0.014, -0.013,
                  -0.012, -0.012, -0.011, -0.011, -0.010, -0.010, -0.010, -0.009,
                  -0.009, -0.009, -0.008]
          }

Remember we are making use of the C3 "tiny world" mesh. Even this small mesh
means 3×3×3×6 cells. The data is presented here in layers from the ground
up.

.. code-block:: python

          result: Dict[float] = {}
          for line in out.splitlines():
              match = self.__pattern.match(line)
              if match:
                  result[int(match.group(1))] = [float(value) for value in match.group(2).split()]

Isolate and extract the field data from standard output.

.. code-block:: python

          for rank in range(0, 5):
              if result[rank] != expected[rank]:
                  message = "Tracer field does not match expectation for rank " + rank
                  raise TestFailed(message, stdout=out, stderr=err)

Check the field data against the expected result and report an error if they
don't match.

.. code-block:: python

          return "Swift outer update tracer"


  if __name__ == '__main__':
      TestEngine.run(SwiftOuterUpdateTracerTest())
