!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief init functionality for dynamo

!> @details Handles init of prognostic and coordinate fields

module init_dynamo_mod

  use assign_coordinate_field_mod,    only : assign_coordinate_field
  use base_mesh_config_mod,          only: geometry, &
                                         base_mesh_geometry_spherical
  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type
  use finite_element_config_mod,      only : element_order, coordinate_order, wtheta_on
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W0, W1, W2, W3, Wtheta, Wchi
  use init_prognostic_fields_alg_mod, only : init_prognostic_fields_alg
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO
  use restart_control_mod,            only : restart_type
  use formulation_config_mod,         only : transport_only
  use transport_config_mod,           only : scheme, &
                                             operators, &
                                             transport_scheme_method_of_lines, &
                                             transport_operators_fv

  implicit none


  contains

  subroutine init_dynamo(mesh_id, chi, u, rho, theta, xi, restart)

    integer(i_def), intent(in)               :: mesh_id
    ! coordinate field
    type( field_type ), intent(inout)        :: chi(3)
    ! prognostic fields
    type( field_type ), intent(inout)        :: u, rho, theta, xi
    type(restart_type), intent(in)           :: restart

    integer(i_def)                           :: coord


    call log_event( 'Dynamo: initialisation...', LOG_LEVEL_INFO )


    ! Calculate coordinates
    if ( geometry == base_mesh_geometry_spherical ) then
      do coord = 1,3
          chi(coord) = field_type (vector_space =                                    &
                         function_space_collection%get_fs(mesh_id,element_order,W0) )
      end do
      ! Assign coordinate field
      call log_event( "Dynamo: Computing W0 coordinate fields", LOG_LEVEL_INFO )
    else
      do coord = 1,3
          chi(coord) = field_type (vector_space =                                    &
                         function_space_collection%get_fs(mesh_id,coordinate_order,Wchi) )
      end do
      ! Assign coordinate field
      call log_event( "Dynamo: Computing Wchi coordinate fields", LOG_LEVEL_INFO )

    end if

    call assign_coordinate_field(chi, mesh_id)


    ! Create prognostic fields
    if ( (transport_only .and. &
         scheme == transport_scheme_method_of_lines .and. &
         operators == transport_operators_fv) .or. &
         wtheta_on  ) then
      ! Only use Wtheta for fv method of lines transport or if wtheta_on
      theta = field_type( vector_space = &
                             function_space_collection%get_fs(mesh_id, element_order, Wtheta) )
    else
      theta = field_type( vector_space = &
                         function_space_collection%get_fs(mesh_id, element_order, W0) )
    end if
    xi    = field_type( vector_space = &
                       function_space_collection%get_fs(mesh_id, element_order, W1) )
    u     = field_type( vector_space = &
                       function_space_collection%get_fs(mesh_id, element_order, W2) )
    rho   = field_type( vector_space = &
                       function_space_collection%get_fs(mesh_id, element_order, W3) )

    ! Initialise prognostic fields
    call init_prognostic_fields_alg( mesh_id, chi, u, rho, theta, xi, restart)


    call log_event( 'Dynamo initialised', LOG_LEVEL_INFO )

  end subroutine init_dynamo

end module init_dynamo_mod
