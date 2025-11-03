!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module create_field_set_mod

  use constants_mod,           only: i_def
  use field_collection_mod,    only: field_collection_type
  use field_mod,               only: field_type
  use integer_field_mod,       only: integer_field_type
  use namelist_collection_mod, only: namelist_collection_type
  use namelist_mod,            only: namelist_type
  use mesh_collection_mod,     only: mesh_collection
  use mesh_mod,                only: mesh_type
  use fs_continuity_mod,       only: W0, W2H, W2V, W3, Wtheta
  use extrusion_mod,           only: twod
  use create_field_mod,        only: create_field
  use log_mod,                 only: log_event, log_level_info

  use lbc_demo_config_mod, only: field_type_real, field_type_integer

  implicit none

  private

  public :: create_field_set

contains

!> @brief Instantiates field set for lbc_demo application
!! @param[in, out]  fld_collection   Field collection to add field set
!! @param[in]       mesh             Mesh to use for field set
!! @param[in]       configuration    Configuration namelist
subroutine create_field_set( fld_collection, mesh, configuration )

  implicit none


  type(field_collection_type),    pointer, intent(inout) :: fld_collection
  type(mesh_type),                pointer, intent(in)    :: mesh
  type(namelist_collection_type),          intent(in)    :: configuration

  type(namelist_type), pointer :: finite_element_nml
  type(namelist_type), pointer :: lbc_demo_nml

  type(field_type)         :: fld
  type(integer_field_type) :: int_fld

  integer(i_def) :: ndata

  type(mesh_type), pointer :: mesh_2d

  integer(i_def) :: order_h
  integer(i_def) :: order_v

  ! Enumerations
  integer :: test_field_type

  lbc_demo_nml       => configuration%get_namelist('lbc_demo')
  finite_element_nml => configuration%get_namelist('finite_element')

  call lbc_demo_nml%get_value( 'field_type', test_field_type )
  call finite_element_nml%get_value( 'element_order_h', order_h )
  call finite_element_nml%get_value( 'element_order_v', order_v )

  mesh_2d => mesh_collection%get_mesh(mesh, twod)

  !----------------------------------------------------------------------------
  ! Create core fields to send/recieve data from file and set I/O behaviours
  !----------------------------------------------------------------------------
  select case ( test_field_type )

  case (field_type_real )

    ndata = 1

    ! W0 (node) field
    call create_field( fld_collection, fld, "w0_field", &
                       mesh, W0, order_h, order_v, ndata )

    ! W2 (edge) fields
    call create_field( fld_collection, fld, "w2h_field", &
                       mesh, W2H, order_h, order_v, ndata )
    call create_field( fld_collection, fld, "w2v_field", &
                       mesh, W2V, order_h, order_v, ndata )

    ! W3 (face) fields
    call create_field( fld_collection, fld, "w3_field", &
                       mesh, W3, order_h, order_v, ndata )
    call create_field( fld_collection, fld, "w3_2d_field", &
                       mesh_2d, W3, order_h, order_v, ndata )

    ! Wtheta (face) fields
    call create_field( fld_collection, fld, "wtheta_field", &
                       mesh, wtheta, order_h, order_v, ndata )

    ! Note: This needs to tally with the 'n_glo' in the
    !       iodef.xml file which will size the field.
    ndata = 2

    call create_field( fld_collection, fld, "multi_data_field", &
                       mesh_2d, W3, order_h, order_v, ndata )


  case ( field_type_integer )

    ndata = 1

    ! W0 (node) field
    call create_field( fld_collection, int_fld, "w0_field", &
                       mesh, W0, order_h, order_v, ndata )

    ! W2 (edge) fields
    call create_field( fld_collection, int_fld, "w2h_field", &
                         mesh, W2H, order_h, order_v, ndata )
    call create_field( fld_collection, int_fld, "w2v_field", &
                       mesh, W2V, order_h, order_v, ndata )

    ! W3 (face) fields
    call create_field( fld_collection, int_fld, "w3_field", &
                       mesh, W3, order_h, order_v, ndata )
    call create_field( fld_collection, int_fld, "w3_2d_field", &
                       mesh_2d, W3, order_h, order_v, ndata )

    ! Wtheta (face) fields
    call create_field( fld_collection, int_fld, "wtheta_field", &
                       mesh, wtheta, order_h, order_v, ndata )

    ! Note: This needs to tally with the 'n_glo' in the
    !       iodef.xml file which will size the field.
    ndata = 2

    call create_field( fld_collection, int_fld, "multi_data_field", &
                       mesh_2d, W3, order_h, order_v, ndata )

    end select

  end subroutine create_field_set

end module create_field_set_mod
