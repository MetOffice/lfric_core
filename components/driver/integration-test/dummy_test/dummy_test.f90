!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise.

!> @brief Placeholder integration test for the driver component
!> @details
!>   TODO: Remove this test and replace with a proper test of the field caching
!>         system.
!>
!>   This test currently creates a simple constant-valued test field on the W3
!>   function space and prints it to verify basic infrastructure functionality.

program dummy_test

  use constants_mod,                  only : r_tran, r_second
  use fs_continuity_mod,              only : W3
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use test_helpers_mod,               only : create_test_field_flood_constant, &
                                             print_test_field
  use test_tiny_world_mod,            only : initialise_tiny_world, &
                                             finalise_tiny_world,   &
                                             mesh

  implicit none

  character(:), allocatable :: option

  ! Initialize the tiny world test infrastructure
  call initialise_tiny_world()
  call run_test_64()
  ! Clean up test infrastructure
  call finalise_tiny_world()

contains

  !> @brief Runs the test case for 64-bit real fields
  !> @details
  !>   Creates a constant-valued test field (all values = 1.0) on the W3
  !>   function space and prints the field values for verification across
  !>   all MPI processes.
  subroutine run_test_64()

    use, intrinsic :: iso_fortran_env, only : real64

    use field_real64_mod, only : field_real64_type
    use test_field_mod,   only : test_field_64_type

    implicit none

    type(function_space_type), pointer :: w3_fs

    type(field_real64_type),  allocatable, target :: dummy_field
    type(test_field_64_type), allocatable         :: dummy_test_field

    ! Retrieve the W3 function space from the mesh
    w3_fs => function_space_collection%get_fs(mesh, 0, 0, W3)

    ! Create a test field filled with constant value 1.0_real64
    call create_test_field_flood_constant( w3_fs, "dummy_field",                 &
                                           dummy_field, dummy_test_field, &
                                           1.0_real64 )

    ! Print the field values to stdout for verification
    call print_test_field( dummy_test_field )

  end subroutine run_test_64

end program dummy_test