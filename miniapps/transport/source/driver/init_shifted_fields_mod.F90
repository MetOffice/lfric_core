!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Initialises functionality for transport.
!> @details Initialises vertically shifted wind and density fields.

module init_shifted_fields_mod

  use constants_mod,                  only: i_def
  use field_mod,                      only: field_type
  use finite_element_config_mod,      only: element_order
  use fs_continuity_mod,              only: W2, W3
  use restart_control_mod,            only: restart_type
  use function_space_mod,             only: function_space_type
  use function_space_collection_mod,  only: function_space_collection

  implicit none

  contains

  !> @param[in] mesh_id      Mesh-id
  !> @param[in,out] wind     Wind field
  !> @param[in,out] density  Density field
  subroutine init_shifted_fields( mesh_id, wind, density )

    integer(i_def), intent(in)        :: mesh_id
    type(field_type), intent(inout)   :: wind
    type(field_type), intent(inout)   :: density

    wind    = field_type( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, W2 ) )

    density = field_type( vector_space = &
                          function_space_collection%get_fs( mesh_id, element_order, W3 ) )

  end subroutine init_shifted_fields

end module init_shifted_fields_mod
