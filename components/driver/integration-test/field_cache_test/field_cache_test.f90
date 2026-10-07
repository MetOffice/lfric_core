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

program field_cache_test

  use constants_mod,                  only : r_tran, r_second, i_def
  use field_cache_alg_mod,            only : field_cache_type
  use field_collection_mod,           only : field_collection_type
  use fs_continuity_mod,              only : W3
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use test_field_mod,                 only : test_field_64_type
  use test_helpers_mod,               only : print_test_field
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

    type(field_real64_type) :: field_cache_field
    type(field_real64_type), pointer :: field_cache_field_ptr
    type(test_field_64_type), allocatable :: field_cache_test_field
    real(real64), pointer :: field_values(:)

    type(field_collection_type) :: depository

    type(field_cache_type) :: field_cache

    call depository%initialise(name = "depository", table_len = 100_i_def)

    ! Retrieve the W3 function space from the mesh
    w3_fs => function_space_collection%get_fs(mesh, 0, 0, W3)
    call field_cache_field%initialise(w3_fs, name = "field_cache_field")
    call depository%add_field(field_cache_field)
    call depository%get_field("field_cache_field", field_cache_field_ptr)

    ! Create a test field wrapper for the LFRic field
    allocate(field_cache_test_field, source=test_field_64_type( field_cache_field_ptr ))
    field_values => field_cache_test_field%get_data()
    field_values = 1.1_real64
    call field_cache_test_field%copy_to_lfric()

    ! Initialize the field cache
    call field_cache%initialise(depository, (/"field_cache_field"/))

    ! Print the initial field values to stdout for verification
    call print_test_field( field_cache_test_field, "field_before_caching" )

    call field_cache%cache_fields()

    field_values = 2.2_real64
    call field_cache_test_field%copy_to_lfric()

    ! Print the field values to stdout for verification
    call print_test_field( field_cache_test_field, "field_after_caching" )

    ! Reset the field cache, which should restore the field values to those cached
    call field_cache%reset_fields()

    call print_test_field( field_cache_test_field, "field_after_reset" )

  end subroutine run_test_64

end program field_cache_test