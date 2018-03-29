!-------------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief Interface to the UM Boundary Layer scheme

module bl_alg_mod

  use field_mod, only: field_type

  implicit none

  public bl_alg_step

contains

  !>@brief Run the UM Boundary Layer scheme
  !>@details The UM Boundary Layer scheme does:
  !>             vertical mixing of heat and momentum (no moisture yet),
  !>             as documented in UMDP24
  !>         NB This version uses winds in w3 space (i.e. A-grid)
  !>@param[in]    theta_in_wth Theta in its native space
  !>@param[in]    rho_in_w3  Density in its native space
  !>@param[in]    rho_in_wth Density in the theta space
  !>@param[in]    exner_in_w3 Exner Pressure in the w3 space (i.e. colocated with density)
  !>@param[in]    exner_in_wth Exner Pressure in the theta space
  !>@param[in]    u1_in_w3   1st horizontal component of wind averaged to w3 space
  !>@param[in]    u2_in_w3   2st horizontal component of wind averaged to w3 space
  !>@param[in]    height_w3  Height in w3
  !>@param[in]    height_wth Height in wth
  !>@param[inout] tstar_2d   Surface tempature
  !>@param[inout] zh_2d      Boundary layer depth
  !>@param[inout] z0msea_2d  Roughness length
  !>@param[inout] theta_inc  BL theta increment
  !>@param[inout] m_v        Vapour mixing ratio
  !>@param[inout] m_cl       Cloud liquid mixing ratio
  subroutine bl_alg_step(theta_in_wth, rho_in_w3, rho_in_wth, exner_in_w3,  &
                         exner_in_wth, u1_in_w3, u2_in_w3, height_w3,       &
                         height_wth, tstar_2d, zh_2d, z0msea_2d, theta_inc, &
                         m_v, m_cl)

    use psykal_lite_phys_mod, only: invoke_bl_kernel

    implicit none

    type( field_type ), intent( in ) :: theta_in_wth, rho_in_w3, rho_in_wth, &
         exner_in_w3, exner_in_wth, u1_in_w3, u2_in_w3, height_w3, height_wth
    type( field_type ), intent( inout ) :: theta_inc, tstar_2d,         &
         zh_2d, z0msea_2d, m_v, m_cl

    call invoke_bl_kernel( theta_in_wth, rho_in_w3, rho_in_wth,           &
                           exner_in_w3, exner_in_wth, u1_in_w3, u2_in_w3, &
                           height_w3, height_wth, tstar_2d, zh_2d,        &
                           z0msea_2d, theta_inc, m_v, m_cl )

  end subroutine bl_alg_step

end module bl_alg_mod
