!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief init functionality for fem specific choices for model

!> @details Creates a collection of function spaces and the coordinate field (chi)

module init_fem_mod

  use constants_mod,                  only : i_def
  use finite_element_config_mod,      only : element_order, coordinate_order
  use field_mod,                      only : field_type
  use fs_continuity_mod,              only : W0, Wtheta, Wchi
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection_type, &
                                             function_space_collection
  use assign_coordinate_field_mod,    only : assign_coordinate_field
  use assign_orography_field_mod,     only : assign_orography_field
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO


  implicit none


  contains

  subroutine init_fem( mesh_id, chi )

    integer(i_def), intent(in)                :: mesh_id
    ! Coordinate field
    type( field_type ), target, intent(inout) :: chi(:)

    type(function_space_type), pointer        :: fs => null()
    integer(i_def)                            :: fs_id
    integer(i_def)                            :: chi_space
    integer(i_def)                            :: coord


    call log_event( 'FEM specifics: creating function spaces...', LOG_LEVEL_INFO )

    allocate( function_space_collection, &
              source = function_space_collection_type() )

    ! Create function spaces from W0 to Wtheta

    do fs_id = W0, Wtheta,1

      fs => function_space_collection%get_fs(mesh_id, element_order, fs_id)

    end do


    ! Compute coordinates
    if ( coordinate_order == 0 ) then
      chi_space = W0
      call log_event( "FEM specifics: Computing W0 coordinate fields", LOG_LEVEL_INFO )
    else
      chi_space = Wchi
      call log_event( "FEM specifics: Computing Wchi coordinate fields", LOG_LEVEL_INFO )
    end if

    do coord = 1,3
      chi(coord) = field_type (vector_space = &
                   function_space_collection%get_fs(mesh_id, coordinate_order, chi_space) )
    end do

    call assign_coordinate_field(chi, mesh_id)
    call assign_orography_field(chi, mesh_id)

    call log_event( 'FEM specifics created', LOG_LEVEL_INFO )

  end subroutine init_fem

end module init_fem_mod
