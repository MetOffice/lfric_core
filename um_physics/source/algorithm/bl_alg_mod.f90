!-------------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief Interface to the UM Boundary Layer scheme

module bl_alg_mod

  use constants_mod,        only: i_def
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
  !>@param[inout] mr        Mixing ratios
  !>@param[in]    outer     outer loop counter
  !>@param[in]    theta     Theta in its native space
  !>@param[in]    rho       Dry density in its native space
  !>@param[in]    mr_n      mixing ratios at time level n
  !>@param[in]    exner     Exner Pressure in the w3 space (i.e. colocated with density)
  !>@param[in]    height_w3  Height in w3
  !>@param[in]    height_wth Height in wth
  !>@param[in]    derived_fields Group of derived fields
  !>@param[inout] twod_fields Group of 2D fields
  subroutine bl_alg_step(dtheta_bl, mr, outer, theta, rho, mr_n, exner, &
                         height_w3, height_wth, derived_fields, twod_fields)

    use psykal_lite_phys_mod, only: invoke_bl_kernel

    implicit none

    type( field_type ), intent( inout )     :: dtheta_bl, mr(nummr)
    type( field_type ), intent( in )        :: theta, rho, exner, mr_n(nummr)
    type( field_type ), intent( in )        :: height_w3, height_wth
    type( field_collection_type ), intent(inout) :: derived_fields
    type( field_collection_type ), intent(inout) :: twod_fields
    integer(kind=i_def), intent(in)         :: outer

    ! Temporary fields unpacked fields from collections
    type( field_type ), pointer :: exner_in_wth  => null()
    type( field_type ), pointer :: rho_in_wth => null()
    type( field_type ), pointer :: u1_in_w3 => null()
    type( field_type ), pointer :: u2_in_w3  => null()
    type( field_type ), pointer :: w_physics => null()
    type( field_type ), pointer :: theta_star => null()
    type( field_type ), pointer :: u1_in_w3_star => null()
    type( field_type ), pointer :: u2_in_w3_star => null()
    type( field_type ), pointer :: w_physics_star => null()
    type( field_type ), pointer :: tstar_2d  => null()
    type( field_type ), pointer :: zh_2d => null()
    type( field_type ), pointer :: z0msea_2d => null()

    !Unpack fields
    exner_in_wth => derived_fields%get_field('exner_in_wth')
    rho_in_wth => derived_fields%get_field('rho_in_wth')
    u1_in_w3 => derived_fields%get_field('u1_in_w3')
    u2_in_w3 => derived_fields%get_field('u2_in_w3')
    w_physics => derived_fields%get_field('w_physics')
    theta_star => derived_fields%get_field('theta_star')
    u1_in_w3_star => derived_fields%get_field('u1_in_w3_star')
    u2_in_w3_star => derived_fields%get_field('u2_in_w3_star')
    w_physics_star => derived_fields%get_field('w_physics_star')

    tstar_2d => twod_fields%get_field('tstar')
    zh_2d => twod_fields%get_field('zh')
    z0msea_2d => twod_fields%get_field('z0msea')
    
    call invoke_bl_kernel( outer, theta, rho, rho_in_wth,                 &
                           exner, exner_in_wth, u1_in_w3, u2_in_w3,       &
                           w_physics, mr_n(imr_v), mr_n(imr_c),           &
                           theta_star, u1_in_w3_star, u2_in_w3_star,      &
                           w_physics_star,                                &
                           height_w3, height_wth, tstar_2d, zh_2d,        &
                           z0msea_2d, dtheta_bl, mr(imr_v), mr(imr_c) )

  end subroutine bl_alg_step

end module bl_alg_mod
