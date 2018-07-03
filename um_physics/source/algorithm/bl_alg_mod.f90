!-------------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief Interface to the UM Boundary Layer scheme

module bl_alg_mod

  use field_mod,            only: field_type
  use field_collection_mod, only: field_collection_type
  use mr_indices_mod,       only: nummr, imr_v, imr_c

  implicit none

  public bl_alg_step

contains

  !>@brief Run the UM Boundary Layer scheme
  !>@details The UM Boundary Layer scheme does:
  !>             vertical mixing of heat and momentum (no moisture yet),
  !>             as documented in UMDP24
  !>         NB This version uses winds in w3 space (i.e. A-grid)
  !>@param[inout] dtheta_bl Theta increment 
  !>@param[in]    theta Theta in its native space
  !>@param[in]    rho    Dry density in its native space
  !>@param[in]    mr Water species mixing ratios
  !>@param[in]    exner  Exner Pressure in the w3 space (i.e. colocated with density)
  !>@param[in]    height_w3  Height in w3
  !>@param[in]    height_wth Height in wth
  !>@param[inout] twod_fields Group of 2D fields
  !>@param[inout] derived_fields Group of derived fields
  subroutine bl_alg_step(dtheta_bl, theta, rho, mr, exner, &
                         height_w3, height_wth,            &
                         derived_fields, twod_fields)

    use psykal_lite_phys_mod, only: invoke_bl_kernel

    implicit none

    type( field_type ), intent( inout )     :: dtheta_bl, mr(nummr)
    type( field_type ), intent( in )        :: theta, rho, exner
    type( field_type ), intent( in )        :: height_w3, height_wth
    type( field_collection_type ), intent(inout) :: derived_fields
    type( field_collection_type ), intent(inout) :: twod_fields

    ! Temporary fields unpacked fields from collections
    type( field_type ), pointer :: exner_in_wth  => null()
    type( field_type ), pointer :: rho_in_wth => null()
    type( field_type ), pointer :: u1_in_w3 => null()
    type( field_type ), pointer :: u2_in_w3  => null()
    type( field_type ), pointer :: tstar_2d  => null()
    type( field_type ), pointer :: zh_2d => null()
    type( field_type ), pointer :: z0msea_2d => null()

    !Unpack fields
    exner_in_wth => derived_fields%get_field('exner_in_wth')
    rho_in_wth => derived_fields%get_field('rho_in_wth')
    u1_in_w3 => derived_fields%get_field('u1_in_w3')
    u2_in_w3 => derived_fields%get_field('u2_in_w3')
    tstar_2d => twod_fields%get_field('tstar')
    zh_2d => twod_fields%get_field('zh')
    z0msea_2d => twod_fields%get_field('z0msea')
    
    call invoke_bl_kernel( theta, rho, rho_in_wth,                        &
                           exner, exner_in_wth, u1_in_w3, u2_in_w3,       &
                           height_w3, height_wth, tstar_2d, zh_2d,        &
                           z0msea_2d, dtheta_bl, mr(imr_v), mr(imr_c) )

  end subroutine bl_alg_step

end module bl_alg_mod
