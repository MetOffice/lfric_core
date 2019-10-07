!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module init_ancils_mod

  use constants_mod,                   only: i_def
  use log_mod,                         only: log_event,         &
                                             log_scratch_space, &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_TRACE
  use finite_element_config_mod,       only: element_order
  use init_tstar_analytic_alg_mod,     only: init_tstar_analytic_alg

  use field_mod,                       only: field_type,     &
                                             read_interface
  use io_mod,                          only: xios_read_field_single_face
  use field_collection_mod,            only: field_collection_type
  use function_space_collection_mod,   only: function_space_collection
  use fs_continuity_mod,               only: W3
  
  implicit none

  private
  public :: init_analytic_ancils, &
            init_aquaplanet_ancils

contains


  !> @details Initialises ancillary fields analytically
  !> @param[in,out] twod_fields the 2D field collection
  subroutine init_analytic_ancils(twod_fields)

    implicit none

    type( field_collection_type ), intent( inout ) :: twod_fields

    type( field_type ), pointer ::    tstar_ptr  => null()

    tstar_ptr => twod_fields%get_field('tstar')

    call init_tstar_analytic_alg(tstar_ptr)

    nullify(tstar_ptr)

  end subroutine init_analytic_ancils

  !> @details Initialises ancillary fields from
  !>          an aquaplanet dump
  !> @param[in,out] twod_fields the 2D field collection
  subroutine init_aquaplanet_ancils(twod_fields)


    implicit none

    type( field_collection_type ), intent( inout ) :: twod_fields

    ! local variables
    ! Pointer to the 2D tstar in the twod field collection
    type( field_type ), pointer ::    tstar_ptr  => null()

    procedure(read_interface), pointer  :: tmp_read_ptr => null()

    call log_event("Reading tstar from dump", LOG_LEVEL_INFO)

    tstar_ptr => twod_fields%get_field('tstar')


    ! Need to set the I/O handler for read. Any ancils here
    ! are currently read from a UM2LFRic dump

    tmp_read_ptr => xios_read_field_single_face
    call tstar_ptr%set_read_behaviour(tmp_read_ptr)
    call tstar_ptr%read_field("read_tstar")

    nullify( tstar_ptr, tmp_read_ptr )

  end subroutine init_aquaplanet_ancils


end module init_ancils_mod
