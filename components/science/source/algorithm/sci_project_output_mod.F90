!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Creates a new field of one function space and projects an existing field on
!> a different function space to it.
!>
module sci_project_output_mod

  implicit none

  private
  public :: project_output

contains

  !> @brief Project a field to a new field.
  !>
  !> @details This procedure uses the galerkin projection and a precomputed
  !>          mass matrix to project a field
  !>
  !> @param[in]    config           Application configuration object
  !> @param[in]    field            To be projected.
  !> @param[inout] projected_field  Receives projection.
  !> @param[in]    chi              Field entity co-ordinates.
  !> @param[in]    panel_id         Cell orientation map.
  !> @param[in]    output_fs        Desired output function space.
  !>
  subroutine project_output( config, field, projected_field, &
                             chi, panel_id,          &
                             output_fs )

    use config_mod,                only: config_type
    use constants_mod,             only: r_def, str_max_filename, i_def
    use field_mod,                 only: field_type
    use field_parent_mod,          only: write_interface
    use operator_mod,              only: operator_type
    use function_space_collection_mod,  only: function_space_collection
    use fs_continuity_mod,         only: W0, W1, W2, W3
    use quadrature_xyoz_mod,               only: quadrature_xyoz_type
    use quadrature_rule_gaussian_mod,      only: quadrature_rule_gaussian_type
    use sci_galerkin_projection_alg_mod, only: galerkin_projection_algorithm

    implicit none

    type(config_type), intent(in) :: config

    ! Input field to project from
    type(field_type),         intent(in)    :: field

    ! Output field to project to
    type(field_type),         intent(inout) :: projected_field(:)

    ! Co-ordinate system
    type(field_type),         intent(in)    :: chi(:)
    type(field_type),         intent(in)    :: panel_id

    ! Output function space
    integer(i_def),           intent(in)    :: output_fs

    type( quadrature_xyoz_type )          :: qr
    type( quadrature_rule_gaussian_type ) :: quadrature_rule
    integer(i_def)                        :: idx, fs_handle
    integer(i_def)                        :: element_order_h
    integer(i_def)                        :: element_order_v
    procedure(write_interface), pointer   :: tmp_write_ptr => null()

    element_order_h = field%get_element_order_h()
    element_order_v = field%get_element_order_v()

    qr = quadrature_xyoz_type( element_order_h+3, element_order_h+3, &
                               element_order_v+3, quadrature_rule )

    ! Determine the input function space
    fs_handle = field%which_function_space()

    ! Create the output field
    call field%get_write_behaviour( tmp_write_ptr )
    do idx = 1, size(projected_field)
      call projected_field(idx)%initialise( &
        vector_space = function_space_collection%get_fs( &
          field%get_mesh(), element_order_h, element_order_v, output_fs &
        ) &
      )
      !
      ! set the write field behaviour based upon what is set in the original
      ! field
      !
      call projected_field(idx)%set_write_behaviour(tmp_write_ptr)
    end do

    ! do the projection
    call galerkin_projection_algorithm(config, projected_field, field, &
                                       chi, panel_id, qr)

  end subroutine project_output

end module sci_project_output_mod
