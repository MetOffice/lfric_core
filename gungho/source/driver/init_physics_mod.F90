!-------------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Init functionality for physics
!> @details Handles initialization of prognostic fields
module init_physics_mod

  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type, &
                                             write_interface
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use field_collection_mod,           only : field_collection_type
  use fs_continuity_mod,              only : W0, W1, W2, W3, Wtheta
  use function_space_mod,             only : function_space_type
  use init_prognostic_fields_alg_mod, only : init_prognostic_fields_alg
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO,         &
                                             LOG_LEVEL_ERROR
  use formulation_config_mod,         only : transport_only
  use transport_config_mod,           only : scheme, &
                                             operators, &
                                             transport_scheme_method_of_lines, &
                                             transport_operators_fv
  use mr_indices_mod,                 only : nummr
  use runtime_constants_mod,          only : create_runtime_constants

  implicit none


contains
  !>@brief Routine to initialise the field objects required by the physics
  !> @param[in] mesh_id Identifier of the mesh
  !> @param[in] twod_mesh_id Identifier of the 2D (surface) mesh
  !> @param[inout] u Wind field
  !> @param[inout] exner Exner pressure field
  !> @param[inout] rho Density field
  !> @param[inout] theta Potential temperature field
  !> @param[out]   derived_fields Collection of FD fields derived from FE fields 
  !> @param[out]   cloud_fields Collection of FD cloud fields
  !> @param[out]   twod_fields Collection of two fields
  subroutine init_physics(mesh_id, twod_mesh_id,                      &
                          u, exner, rho, theta,                       &
                          derived_fields, cloud_fields, twod_fields)

    integer(i_def), intent(in)               :: mesh_id
    integer(i_def), intent(in)               :: twod_mesh_id
    ! Prognostic fields
    type( field_type ), intent(inout)        :: u, exner, rho, theta

    ! Collections of fields
    type(field_collection_type), intent(out) :: cloud_fields
    type(field_collection_type), intent(out) :: twod_fields
    type(field_collection_type), intent(out) :: derived_fields

    ! pointers to vector spaces
    type(function_space_type), pointer  :: vector_space => null()

    integer(i_def) :: theta_space

    call log_event( 'Physics: initialisation...', LOG_LEVEL_INFO )
    
    theta_space=theta%which_function_space()
    
    if (theta_space /= Wtheta)then
      call log_event( 'Physics: requires theta variable to be in Wtheta', LOG_LEVEL_ERROR )
    end if

    if (element_order > 0)then
      call log_event( 'Physics: requires lowest order elements', LOG_LEVEL_ERROR )
    end if

    !========================================================================
    ! Here we create some field collections
    !========================================================================
    derived_fields  =  field_collection_type()
   
    ! Wtheta fields
    vector_space=>function_space_collection%get_fs(mesh_id, 0, Wtheta)

    call add_physics_field(derived_fields, 'w_physics', vector_space)
    call add_physics_field(derived_fields, 'rho_in_wth', vector_space)
    call add_physics_field(derived_fields, 'exner_in_wth', vector_space)
    call add_physics_field(derived_fields, 'w_physics_star', vector_space)
    call add_physics_field(derived_fields, 'theta_star', vector_space)

    ! W3 fields
    vector_space=> function_space_collection%get_fs(mesh_id, 0, W3) 

    call add_physics_field(derived_fields, 'u1_in_w3', vector_space)
    call add_physics_field(derived_fields, 'u2_in_w3', vector_space)
    call add_physics_field(derived_fields, 'u3_in_w3', vector_space)
    call add_physics_field(derived_fields, 'theta_in_w3', vector_space)
    call add_physics_field(derived_fields, 'u1_in_w3_star', vector_space)
    call add_physics_field(derived_fields, 'u2_in_w3_star', vector_space)
    call add_physics_field(derived_fields, 'u3_in_w3_star', vector_space)

    ! W2 fields
    vector_space=>function_space_collection%get_fs(mesh_id, 0, W2)

    call add_physics_field(derived_fields, 'u_physics', vector_space)
    call add_physics_field(derived_fields, 'u_physics_star', vector_space)
    call add_physics_field(derived_fields, 'u_star', vector_space)


    !========================================================================
    ! Here we create some 2d fields for the UM physics
    !========================================================================
    twod_fields = field_collection_type()
    vector_space=> function_space_collection%get_fs(twod_mesh_id, 0, W3) 

    call add_physics_field(twod_fields, 'tstar', vector_space)
    call add_physics_field(twod_fields, 'zh', vector_space)
    call add_physics_field(twod_fields, 'z0msea', vector_space)


    !========================================================================
    ! ...and some fields for use with clouds
    !========================================================================
    cloud_fields = field_collection_type()
    vector_space=>function_space_collection%get_fs(mesh_id, 0, Wtheta)

    call add_physics_field(cloud_fields, 'area_fraction', vector_space)
    call add_physics_field(cloud_fields, 'ice_fraction', vector_space)
    call add_physics_field(cloud_fields, 'liquid_fraction', vector_space)
    call add_physics_field(cloud_fields, 'bulk_fraction', vector_space)

    call log_event( 'Physics initialised', LOG_LEVEL_INFO )

  end subroutine init_physics

  subroutine add_physics_field(field_collection, name, vector_space)

    use output_config_mod,  only : write_xios_output
    use io_mod,             only : xios_write_field_face

    implicit none
    
    type(field_collection_type), intent(inout) :: field_collection
    character(*), intent(in) :: name
    type(function_space_type), intent(in)  :: vector_space

    !Local variables
    type(field_type) :: new_field

    ! pointers for xios write interface
    procedure(write_interface), pointer :: write_behaviour => null()

    ! All physics fields currently require output on faces...
    write_behaviour => xios_write_field_face

    new_field = field_type( vector_space, name=trim(name) )

    if (write_xios_output) then
      call new_field%set_write_field_behaviour(write_behaviour)
    end if

    call field_collection%add_field(new_field)

  end subroutine add_physics_field

end module init_physics_mod
