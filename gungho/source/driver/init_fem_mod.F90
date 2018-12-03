!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief init functionality for fem specific choices for model

!> @details Creates a collection of function spaces and the coordinate field (chi)

module init_fem_mod

  use constants_mod,                  only : i_def, i_native
  use finite_element_config_mod,      only : element_order, coordinate_order
  use field_mod,                      only : field_type
  use fs_continuity_mod,              only : W0, W1, W2, W3, Wtheta, Wchi
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection_type, &
                                             function_space_collection
  use assign_coordinate_field_mod,    only : assign_coordinate_field
  use assign_orography_field_mod,     only : assign_orography_field
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO


  implicit none


  contains

  !>@brief Initialises the coordinate field (chi) and (optionally) the vertically
  !>        shifted coordinate field.
  !> @param[in] mesh_id         Identifier of the mesh
  !> @param[in,out] chi         Coordinate field
  !> @param[in] shifted_mesh_id Identifier of vertically shifted mesh with an extra level
  !> @param[in,out] shifted_chi Spatial coordinates of vertically shifted mesh
  subroutine init_fem( mesh_id, chi, shifted_mesh_id, shifted_chi )

    integer(i_def), intent(in)                  :: mesh_id
    ! Coordinate field
    type( field_type ), intent(inout)           :: chi(:)
    integer(i_def), intent(in), optional        :: shifted_mesh_id
    type( field_type ), intent(inout), optional :: shifted_chi(:)

    integer(i_native), parameter :: fs_list(5) = [W0, W1, W2, W3, Wtheta]

    type(function_space_type), pointer :: fs => null()
    type(function_space_type), pointer :: shifted_fs => null()
    integer(i_native)                  :: fs_index
    integer(i_def)                     :: chi_space
    integer(i_def)                     :: coord


    call log_event( 'FEM specifics: creating function spaces...', LOG_LEVEL_INFO )

    allocate( function_space_collection, &
              source = function_space_collection_type() )

    ! Create function spaces from W0 to Wtheta
    do fs_index = 1, size(fs_list)
      fs => function_space_collection%get_fs( mesh_id,       &
                                              element_order, &
                                              fs_list(fs_index) )
    end do


    ! Compute coordinates
    if ( coordinate_order == 0 ) then
      chi_space = W0
      call log_event( "FEM specifics: Computing W0 coordinate fields", LOG_LEVEL_INFO )
    else
      chi_space = Wchi
      call log_event( "FEM specifics: Computing Wchi coordinate fields", LOG_LEVEL_INFO )
    end if

    fs => function_space_collection%get_fs(mesh_id, coordinate_order, chi_space)

    do coord = 1, size(chi)
      chi(coord) = field_type (vector_space = fs )
    end do

    call assign_coordinate_field(chi, mesh_id)
    call assign_orography_field(chi, mesh_id)


    ! Create shifted vertical mesh extrusion.
    if (present(shifted_mesh_id)) then
      shifted_fs => function_space_collection%get_fs(shifted_mesh_id, coordinate_order, chi_space)

      do coord = 1, size(chi)
        shifted_chi(coord) = field_type(vector_space = shifted_fs)
      end do

      call assign_coordinate_field(shifted_chi, shifted_mesh_id)
      call assign_orography_field(shifted_chi, shifted_mesh_id)

      nullify( shifted_fs )
    end if

    nullify( fs )
    call log_event( 'FEM specifics created', LOG_LEVEL_INFO )

  end subroutine init_fem

end module init_fem_mod
