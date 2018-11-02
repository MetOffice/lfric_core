!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Initialises functionality for transport.

!> @details Handles initialisation of wind, density and chi fields.
!> transport_init_fields_alg() is used to initialise density.

module init_transport_mod

  use constants_mod,                  only: i_def
  use field_mod,                      only: field_type,               &
                                            write_interface,          &
                                            checkpoint_interface,     &
                                            restart_interface
  use finite_element_config_mod,      only: element_order
  use fs_continuity_mod,              only: W2, W3
  use restart_control_mod,            only: restart_type
  use runtime_constants_mod,          only: create_runtime_constants
  use function_space_mod,             only: function_space_type
  use function_space_collection_mod,  only: function_space_collection
  use output_config_mod,              only: write_xios_output
  use io_mod,                         only: xios_write_field_face,    &
                                            checkpoint_xios,          &
                                            checkpoint_netcdf,        &
                                            restart_netcdf,           &
                                            restart_xios
  use log_mod,                        only: log_event,                &
                                            LOG_LEVEL_INFO
  use transport_init_fields_alg_mod,  only: transport_init_fields_alg

  implicit none

  contains

  !> @param[in] mesh_id        Mesh-id
  !> @param[in,out] chi        Coordinate field
  !> @param[in,out] wind       Wind field
  !> @param[in,out] density    Density field
  !> @param[in,out] dep_pts_x  Departure points in the x-direction
  !> @param[in,out] dep_pts_y  Departure points in the y-direction
  !> @param[in,out] dep_pts_z  Departure points in the z-direction
  subroutine init_transport( mesh_id, chi, wind, density, dep_pts_x, dep_pts_y, dep_pts_z )

    integer(i_def), intent(in)        :: mesh_id
    type(field_type), intent(inout)   :: chi(:)
    type(field_type), intent(inout)   :: wind
    type(field_type), intent(inout)   :: density
    type(field_type), intent(inout)   :: dep_pts_x
    type(field_type), intent(inout)   :: dep_pts_y
    type(field_type), intent(inout)   :: dep_pts_z

    type(function_space_type), pointer       :: function_space => null()
    procedure(write_interface), pointer      :: tmp_write_ptr => null()
    procedure(restart_interface), pointer    :: tmp_restart_ptr => null()

    wind    = field_type( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, W2 ), &
                          output_space = W3 )
    density = field_type( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, W3 ) )

    dep_pts_x  = field_type( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, W2 ), &
                          output_space = W3 )
    dep_pts_y  = field_type( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, W2 ), &
                          output_space = W3 )
    dep_pts_z  = field_type( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, W2 ), &
                          output_space = W3 )


    ! Create runtime_constants object.
    call create_runtime_constants( mesh_id, chi )

    ! Initialise density field
    call transport_init_fields_alg( density )

    ! Set I/O behaviours for diagnostic output

    if ( write_xios_output ) then
       ! Fields that are output on the XIOS face domain
       tmp_write_ptr => xios_write_field_face
       call wind%set_write_field_behaviour( tmp_write_ptr )
       call density%set_write_field_behaviour( tmp_write_ptr )
    end if

    call density%set_restart_behaviour( tmp_restart_ptr )
    call wind%set_restart_behaviour( tmp_restart_ptr )

    nullify( function_space, tmp_write_ptr, tmp_restart_ptr )

  end subroutine init_transport

end module init_transport_mod
