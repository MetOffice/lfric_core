!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to boundary layer scheme.
!>
module bl_kernel_mod

  use argument_mod,       only : arg_type,                    &
                                 GH_FIELD, GH_READ, GH_WRITE, &
                                 CELLS
  use constants_mod,      only : r_def, i_def, r_double
  use formulation_config_mod, only: use_moisture
  use fs_continuity_mod,  only : W3, Wtheta
  use kernel_mod,         only : kernel_type
  use physics_config_mod, only: l_flux_bc, fixed_flux_e, fixed_flux_h

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: bl_kernel_type
    private
    type(arg_type) :: meta_args(15) = (/         &
        arg_type(GH_FIELD,   GH_READ,   WTHETA), &
        arg_type(GH_FIELD,   GH_READ,   W3),     &
        arg_type(GH_FIELD,   GH_READ,   WTHETA), &
        arg_type(GH_FIELD,   GH_READ,   W3),     &
        arg_type(GH_FIELD,   GH_READ,   WTHETA), &
        arg_type(GH_FIELD,   GH_READ,   W3),     &
        arg_type(GH_FIELD,   GH_READ,   W3),     &
        arg_type(GH_FIELD,   GH_READ,   W3),     &
        arg_type(GH_FIELD,   GH_READ,   WTHETA), &
        arg_type(GH_FIELD,   GH_WRITE,  W3),     &
        arg_type(GH_FIELD,   GH_WRITE,  W3),     &
        arg_type(GH_FIELD,   GH_WRITE,  W3),     &
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA), &
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA), &
        arg_type(GH_FIELD,   GH_WRITE,  WTHETA)  &
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass ::bl_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! overload the default structure constructor for function space
  interface bl_kernel_type
    module procedure bl_kernel_constructor
  end interface

  public bl_code

contains

type(bl_kernel_type) function bl_kernel_constructor() result(self)
  return
end function bl_kernel_constructor

!> @brief Interface to the BL scheme
!! @param[in] nlayers Number of layers
!! @param[in] nlayers_2d Number of layers for 2d field (i.e. 1)
!! @param[in] theta_in_wth Potential temperature field
!! @param[in] rho_in_w3 density field in potential density space
!! @param[in] rho_in_wth density field in potential temperature space
!! @param[in] exner_in_w3 Exner pressure field in density space
!! @param[in] exner_in_wth Exner pressure field in potential temperature space
!! @param[in] u1_in_w3 'zonal' wind in density space
!! @param[in] u2_in_w3 'meridional' wind in density space
!! @param[in] heigh_w3 Height of density space levels above surface
!! @param[in] heigh_wth Height of temperature space levels above surface
!! @param[inout] tstar_2d   Surface tempature
!! @param[inout] zh_2d      Boundary layer depth
!! @param[inout] z0msea_2d  Roughness length
!! @param[inout] theta_inc  BL theta increment
!! @param[inout] m_v        Vapour mixing ratio
!! @param[inout] m_cl       Cloud liquid mixing ratio
!! @param[in] ndf_wth Number of degrees of freedom per cell for potential temperature space
!! @param[in] undf_wth Number unique of degrees of freedom  for potential temperature space
!! @param[in] map_wth Dofmap for the cell at the base of the column for potential temperature space
!! @param[in] ndf_w3 Number of degrees of freedom per cell for density space
!! @param[in] undf_w3 Number unique of degrees of freedom  for density space
!! @param[in] map_w3 Dofmap for the cell at the base of the column for density space
!! @param[in] ndf_2d Number of degrees of freedom per cell for 2D fields
!! @param[in] undf_2d Number unique of degrees of freedom  for 2D fields
!! @param[in] map_2d Dofmap for the cell at the base of the column for 2D fields
subroutine bl_code(nlayers, & 
                   nlayers_2d, & 
                   theta_in_wth, &
                   rho_in_w3, &
                   rho_in_wth, &
                   exner_in_w3, &
                   exner_in_wth, &
                   u1_in_w3, &
                   u2_in_w3, &
                   height_w3, &
                   height_wth, &
                   tstar_2d, &
                   zh_2d, &
                   z0msea_2d, &
                   theta_inc, &
                   m_v, &
                   m_cl, &
                   ndf_wth, &
                   undf_wth, &
                   map_wth, &
                   ndf_w3, &
                   undf_w3, &
                   map_w3, &
                   ndf_2d, &
                   undf_2d, &
                   map_2d)

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    ! structures holding diagnostic arrays - not used
    use bl_diags_mod, only: BL_diag
    use sf_diags_mod, only: sf_diag
    ! other modules containing stuff passed to BL
    use atm_fields_bounds_mod, only: tdims, udims, vdims, udims_s, vdims_s, &
         pdims
    use atmos_physics2_alloc_mod !everything
    use bdy_expl3_mod, only: bdy_expl3
    use conv_diag_6a_mod, only: conv_diag_6a
    use jules_sea_seaice_mod, only: nice_use
    use jules_surface_types_mod, only: npft, ntype
    use level_heights_mod, only: r_theta_levels, r_rho_levels, eta_theta_levels
    use ni_bl_ctl_mod, only: ni_bl_ctl
    use ni_imp_ctl_mod, only: ni_imp_ctl
    use nlsizes_namelist_mod, only: row_length, rows, model_levels, land_field,&
         sm_levels, ntiles, bl_levels, tr_vars
    use planet_constants_mod, only: p_zero, kappa, planet_radius
    use timestep_mod, only: timestep

    implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, nlayers_2d
  integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3
  integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
  integer(kind=i_def), intent(in) :: undf_wth, undf_w3

  integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth
  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_2d),  intent(in) :: map_2d

  real(kind=r_def), dimension(undf_wth), intent(inout):: theta_inc, m_v, m_cl
  real(kind=r_def), dimension(undf_w3),  intent(in)   :: rho_in_w3,          &
                                                         exner_in_w3,        &
                                                         u1_in_w3, u2_in_w3, &
                                                         height_w3
  real(kind=r_def), dimension(undf_wth), intent(in)   :: rho_in_wth,         &
                                                         exner_in_wth,       &
                                                         height_wth,         &
                                                         theta_in_wth
  real(kind=r_def), dimension(undf_2d), intent(inout) :: tstar_2d, zh_2d,    &
                                                         z0msea_2d

  ! Local variables for the kernel
    integer :: k

    ! switches and model parameters/dimensions/time etc
    integer :: cycleno, error_code
    integer :: val_year, val_day_number, val_hour, val_minute, val_second
    integer :: co2_dim_len, co2_dim_row, rhc_row_length, rhc_rows,          &
         cloud_levels, n_cca_levels
    integer :: asteps_since_triffid
    integer :: nscmdpkgs
    logical :: l_scrn, l_aero_classic, l_spec_z0, l_plsp,        &
         l_mixing_ratio, l_extra_call, l_calc_at_p
    logical :: l_scmdiags(0)

    ! profile fields from level 1 upwards
    real(r_double), dimension(row_length,rows,model_levels) ::              &
         p, rho_wet_rsq, rho_wet, rho_dry, z_rho, z_theta,                  &
         bulk_cloud_fraction, bl_w_var, rhcpt, t_latest, q_latest,          &
         qcl_latest, qcf_latest, cf_latest, cfl_latest, cff_latest, cca,    &
         cca0, ccw0, area_cloud_fraction, cloud_fraction_liquid,            &
         cloud_fraction_frozen
    real(r_double), dimension(row_length,rows,bl_levels) ::                 &
         e_trb, tsq_trb, qsq_trb, cov_trb, tau_fd_x, tau_fd_y, rhogamu,     &
         rhogamv, f_ngstress, fqw, ftl, rhokh
    real(r_double), dimension(row_length,rows,model_levels) :: rho_wet_tq
                             !IB removed -1 from model_levels here
    real(r_double), dimension(row_length,rows,model_levels) :: exner_rho_levels
                             !IB removed +1 from model_levels here
    ! profile fields on u/v points
    real(r_double), dimension(row_length,rows,model_levels) :: u, v, r_u, r_v
    real(r_double), dimension(row_length,rows,bl_levels) ::                 &
         taux_fd_u, tauy_fd_v, rhogamu_u, rhogamv_v, f_ngstress_u,          &
         f_ngstress_v, taux, tauy
    ! profile fields from level 0 upwards
    real(r_double), dimension(row_length,rows,0:model_levels) ::            &
         p_layer_centres, w_copy, etadot_copy, R_w, p_layer_boundaries, w
    ! profile fields with 0 level which isn't used
    real(r_double), dimension(row_length,rows,0:model_levels) ::            &
         q, qcl, qcf, theta, exner_theta_levels, aerosol, dust_div1,        &
         dust_div2, dust_div3, dust_div4, dust_div5, dust_div6, so2, dms,   &
         so4_aitken, so4_accu, so4_diss, nh3, soot_new, soot_aged, soot_cld,&
         bmass_new, bmass_aged, bmass_cld, ocff_new, ocff_aged, ocff_cld,   &
         nitr_acc, nitr_diss, ozone_tracer, conv_prog_precip
    ! multi tracer variables
    real(r_double), dimension(row_length,rows,0:model_levels,tr_vars) ::    &
         free_tracers
    ! profile fields with a hard-wired 2
    real(r_double), dimension(row_length,rows,2,bl_levels) ::               &
         rad_hr, micro_tends
    ! single level fields
    real(r_double), dimension(row_length,rows) ::                           &
         p_star, lw_down, cos_zenith_angle, ti_gb, tstar, zh_prev, ddmfx,   &
         zlcl, zhpar, z0h_scm, z0m_scm, flux_e, flux_h, z0msea,             &
         photosynth_act_rad, soil_clay, soil_sand, dust_mrel1, dust_mrel2,  &
         dust_mrel3, dust_mrel4, dust_mrel5, dust_mrel6, tstar_sea, zh, dzh,&
         zhpar_shcu, flandfac, fseafac, rhokm_land, rhokm_ssi, cdr10m,      &
         tstar_land, tstar_ssi, dtstar_sea, t1_sd, q1_sd, wstar, wthvs,     &
         xx_cos_theta_latitude, ice_fract, ls_rain, ls_snow, conv_rain,     &
         conv_snow, cca0_2d, qcl_inv_top, co2_emits, co2flux, tscrndcl_ssi, &
         tstbtrans, sum_eng_fluxes, sum_moist_flux, drydep2, olr,           &
         surf_ht_flux_land, zlcl_mixed, theta_star_surf, qv_star_surf,      &
         snowmelt, tstar_sice, u_0_p, v_0_p, w_max, deep_flag, past_precip, &
         past_conv_ht, zlcl_uv, ql_ad, cin_undilute, cape_undilute,         &
         entrain_coef, qsat_lcl, delthvu, dtstar_sice
    ! single level fields on u/v points
    real(r_double), dimension(row_length,rows) :: u_0, v_0, rhokm_u_land,   &
         rhokm_u_ssi, rhokm_v_land, rhokm_v_ssi, flandfac_u, flandfac_v,    &
         fseafac_u, fseafac_v, taux_land, tauy_land, taux_ssi, tauy_ssi
    ! single level fields
    integer, dimension(row_length,rows) :: nlcl, ntml, ntpar, lcbase, ccb0, &
         cct0, conv_type, nbdsc, ntdsc
    ! single level fields
    logical, dimension(row_length,rows) :: land_sea_mask, cumulus, l_shallow,&
         l_pc2_diag_sh_pts, no_cumulus, l_congestus, l_congestus2
    ! fields on ice categories
    real(r_double), dimension(row_length,rows,nice_use) ::                  &
         ice_fract_cat_use, k_sice, co2, ti, tstar_sice_cat, radnet_sice,   &
         radnet_sea, fqw_ice, ftl_ice, di_ncat, ice_fract_ncat
    ! field on land points and soil levels
    real(r_double), dimension(land_field,sm_levels) :: soil_layer_moisture, &
         smvccl_levs, smvcwt_levs, smvcst_levs, sthf, sthu, ext
    ! fields on land points
    real(r_double), dimension(land_field) :: hcon, sil_orog_land,           &
         ho2r2_orog, sd_orog, z0m_soil, albsoil, gs
    integer, dimension(land_field) :: land_index
    ! fields on land points and tiles
    real(r_double), dimension(land_field,ntiles) :: canopy, catch,          &
         catch_snow, snow_tile, z0_tile, z0h_tile_bare, sw_tile, tstar_tile,&
         cs, frac, canht_ft, lai_ft, t_soil, g_leaf_acc, npp_ft_acc,        &
         resp_w_ft_acc, resp_s_acc, ftl_tile, fqw_tile, epot_tile,          &
         dtstar_tile, tsurf_elev_surft, tscrndcl_tile, ei_tile, ecan_tile,  &
         melt_tile, surf_htf_tile
    ! stashwork arrays
    real(r_double) :: stashwork3(1), stashwork9(1)    

    !-----------------------------------------------------------------------
    ! Initialisation of variables and arrays
    ! These will need to be configured/passed in by LFRic, but for
    ! the initial implementation these will be set here...
    !-----------------------------------------------------------------------
    ! endgame loop counter
    cycleno=1
    ! model dimensions
    co2_dim_len=1
    co2_dim_row=1
    rhc_row_length=1
    rhc_rows=1
    cloud_levels=model_levels
    n_cca_levels=model_levels
    ! surface ancils
    land_sea_mask=.false.
    flandg=0.0
    ice_fract=0.0
    ice_fract_cat_use=0.0
    ice_fract_ncat=0.0
    ! diagnostic flags
    l_scrn=.false.
    nscmdpkgs=0
    error_code=0
    l_plsp=.false.
    ! other logicals
    l_aero_classic=.false.
    l_mixing_ratio=.false.
    l_extra_call=.false.
    ! surface forcing
    if ( l_flux_bc ) then
      flux_e(:,:)=fixed_flux_e
      flux_h(:,:)=fixed_flux_h
    end if
    l_spec_z0=.false.

    if ( .not. allocated(uhalo) ) then
      call atmos_physics2_alloc(land_field,ntiles,1,1,npft,ntype,       &
                                1,1,sm_levels,bl_levels,nice_use)
    end if

    !-----------------------------------------------------------------------
    ! Main model prognostics - passed in from Gungho
    ! For the initial implementation we pass each individual column
    ! of data to an array sized (1,1,k) to match the UMs (i,j,k) data 
    ! layout.  
    !-----------------------------------------------------------------------

    ! assuming map_wth(1) points to level 0
    ! and map_w3(1) points to level 1

    do k = 1, model_levels
      ! potential temperature on theta levels
      theta(1,1,k) = theta_in_wth(map_wth(1) + k)
      ! wet density on theta levels
      rho_wet_tq(1,1,k) = rho_in_wth(map_wth(1) + k)
      ! wet density on rho levels
      rho_wet(1,1,k) = rho_in_w3(map_w3(1) + k-1)
      ! pressure on rho levels
      p(1,1,k) = p_zero*(exner_in_w3(map_w3(1) + k-1))**(1.0_r_def/kappa)
      ! pressure on theta levels
      p_layer_centres(1,1,k) = p_zero*(exner_in_wth(map_wth(1) + k))**(1.0_r_def/kappa)
      ! u wind on rho levels
      u_p(1,1,k) = u1_in_w3(map_w3(1) + k-1)
      ! v wind on rho levels
      v_p(1,1,k) = u2_in_w3(map_w3(1) + k-1)
      ! height of rho levels from centre of planet
      r_rho_levels(1,1,k) = height_w3(map_w3(1) + k-1)                  &
                          + planet_radius
      ! height of theta levels from centre of planet
      r_theta_levels(1,1,k) = height_wth(map_wth(1) + k)                &
                            + planet_radius
    end do

    ! surface pressure
    p_layer_centres(1,1,0) = p_zero*(exner_in_wth(map_wth(1) + 0))**(1.0_r_def/kappa)
    p_star(1,1) = p_layer_centres(1,1,0)
    ! exner pressure on rho and theta levels
    exner_rho_levels = (p/p_zero)**kappa
    exner_theta_levels = (p_layer_centres/p_zero)**kappa
    ! near surface potential temperature
    theta(1,1,0) = theta_in_wth(map_wth(1) + 0)
    ! wet density multiplied by planet radius squared on rho levs
    rho_wet_rsq = rho_wet * r_rho_levels**2
    ! extended halo u and v winds
    u_px = u_p
    v_px = v_p
    ! u and v winds on their native grid
    u = u_p
    v = v_p
    ! surface currents
    u_0_px = 0.0
    v_0_px = 0.0
    u_0 = 0.0
    v_0 = 0.0
    u_0_p = 0.0
    v_0_p = 0.0
    ! moisture fields
    if (use_moisture) then
      do k = 1, model_levels
        ! water vapour mixing ratio
        q(1,1,k) = m_v(map_wth(1) + k)
        ! cloud liquid mixing ratio
        qcl(1,1,k) = m_cl(map_wth(1) + k)
      end do
      q(1,1,0) = m_v(map_wth(1) + 0)
      qcl(1,1,0) = m_cl(map_wth(1) + 0)
    else
      q = 0.0
      qcl = 0.0
    end if
    qcf = 0.0
    bulk_cloud_fraction = 0.0
    ! surface height
    r_theta_levels(1,1,0) = height_wth(map_wth(1) + 0) + planet_radius
    ! eta space, 0-1 scaled height
    eta_theta_levels(:) = (r_theta_levels(1,1,:)-r_theta_levels(1,1,0))  &
                        /(r_theta_levels(1,1,model_levels)-planet_radius)
    ! height of levels above surface
    z_rho = r_rho_levels-r_theta_levels(1,1,0)
    z_theta(1,1,:) = r_theta_levels(1,1,1:model_levels)-r_theta_levels(1,1,0)
    ! vertical velocity - needs setting------------------------------
    w = 0.0
    w_copy = w
    etadot_copy = w_copy / r_theta_levels(1,1,model_levels)

    !-----------------------------------------------------------------------
    ! Things saved from one timestep to the next
    !-----------------------------------------------------------------------
    ! sea surface temperature
    tstar(1,1) = tstar_2d(map_2d(1))
    tstar_sea = tstar
    tstar_ssi = tstar
    ! not needed but given values to avoid qsat crash in Jules
    tstar_land = 300.0
    tstar_tile = 300.0
    tstar_sice = 300.0
    tstar_sice_cat = 300.0
    ! previous BL height
    zh(1,1) = zh_2d(map_2d(1))
    zh_prev = zh
    ! surface roughness
    z0msea(1,1) = z0msea_2d(map_2d(1))
    ! downdraft at cloud base
    ddmfx = 0.0

    !-----------------------------------------------------------------------
    ! Things saved from other parametrization schemes on this timestep
    !-----------------------------------------------------------------------
    rad_hr = 0.0        ! From radiation
    micro_tends = 0.0   ! From microphysics

    !-----------------------------------------------------------------------
    ! code below here should mimic the call from the UMs atmos_physics2
    !-----------------------------------------------------------------------

    CALL conv_diag_6a(                                                  &

    !     IN Parallel variables
            row_length, rows                                            &

    !     IN model dimensions.
          , bl_levels                                                   &
          , p, p_layer_centres(1,1,1),exner_rho_levels                  &
          , rho_wet, rho_wet_tq, z_theta, z_rho                         &

    !     IN Model switches
          , l_mixing_ratio                                              &
          , l_extra_call                                                &
          , no_cumulus                                                  &

    !     IN cloud data
          , qcf(1:row_length,1:rows,1:tdims%k_end)                      &
          , qcl(1:row_length,1:rows,1:tdims%k_end), bulk_cloud_fraction &

    !     IN everything not covered so far :

          , p_star, q(1:row_length,1:rows,1:tdims%k_end)                &
          , theta(tdims%i_start:tdims%i_end,                            &
              tdims%j_start:tdims%j_end,1:tdims%k_end)                  &
          , exner_theta_levels(tdims%i_start:tdims%i_end,               &
              tdims%j_start:tdims%j_end, 1:tdims%k_end)                 &
          , u_p, v_p, u_0_p, v_0_p                                      &
          , tstar_land, tstar_sea, tstar_sice, z0msea                   &
          , L_flux_bc, flux_e, flux_h, L_spec_z0, z0m_scm, z0h_scm      &
          , tstar, land_sea_mask, flandg, ice_fract                     &
          , w_copy, w_max, deep_flag, past_precip, past_conv_ht         &
          , conv_prog_precip                                            &
    !
    !     SCM Diagnostics (dummy values in full UM)
          , nSCMDpkgs,L_SCMDiags                                        &

    !     OUT data required elsewhere in UM system :
          , zh,zhpar,dzh,qcl_inv_top,zlcl,zlcl_uv,delthvu,ql_ad         &
          , ntml,ntpar,nlcl                                             &
          , cumulus,l_shallow,l_congestus,l_congestus2                  &
          , conv_type                                                   &
          , CIN_undilute,CAPE_undilute, wstar, wthvs                    &
          , entrain_coef, qsat_lcl                                      &
          , Error_code                                                  &
            )

    CALL NI_bl_ctl (                                                    &
    !     IN parameters for SISL scheme
         CycleNo,                                                       &
    !     IN time stepping information
         val_year, val_day_number, val_hour, val_minute, val_second,    &
    !     IN model dimensions.
         land_field, ntiles, bl_levels,                                 &
    !     IN switches
         L_scrn, L_aero_classic,                                        &
    !     IN data fields.
         p, p_layer_centres, rho_wet_rsq,rho_wet,rho_dry, u_p, v_p,     &
         u_px, v_px, u_0_px, v_0_px,                                    &
         land_sea_mask, q, qcl, qcf, p_star, theta, exner_theta_levels, rad_hr,&
         micro_tends, soil_layer_moisture, rho_wet_tq, z_rho, z_theta,  &
    !     IN ancillary fields and fields needed to be kept from tstep to tstep
         hcon,smvccl_levs, smvcwt_levs, smvcst_levs, sthf,sthu,sil_orog_land,&
    !-------------------------------------------------------------------------
         ho2r2_orog, sd_orog, ice_fract_cat_use, k_sice(:,:,1:nice_use),&
         land_index, photosynth_act_rad,                                &
         soil_clay,soil_sand,dust_mrel1,dust_mrel2,                     &
         dust_mrel3,dust_mrel4,dust_mrel5,dust_mrel6,                   &
    !     IN additional variables for JULES
         canopy, catch, catch_snow, snow_tile, z0_tile, z0h_tile_bare,  &
         z0m_soil, lw_down, sw_tile, tstar_tile, tsurf_elev_surft,      &
         co2(1:co2_dim_len,1:co2_dim_row,1),                            &
         asteps_since_triffid,                                          &
         cs,frac,canht_ft,lai_ft,fland,flandg,albsoil,cos_zenith_angle, &
    !     IN everything not covered so far
             t_soil,ti_gb,                                              &
             ti,tstar,zh_prev,ddmfx,bulk_cloud_fraction,zhpar,zlcl,     &
    !     IN SCM namelist data
         L_spec_z0, z0m_scm, z0h_scm, flux_e, flux_h, L_flux_bc,        &
    !     SCM diagnostics and STASH
         nSCMDpkgs, L_SCMDiags, BL_diag, sf_diag,                       &
    !     INOUT data
         gs,z0msea,w_copy,etadot_copy,tstar_sea,tstar_sice_cat,zh,dzh,cumulus,&
         ntml,ntpar,l_shallow,error_code,                               &
    !     INOUT additional variables for JULES
         g_leaf_acc,npp_ft_acc,resp_w_ft_acc,resp_s_acc,                &
    !     INOUT variables for TKE based turbulence schemes
         e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                  &
    !     OUT variables for message passing
         flandfac, fseafac,rhokm_land, rhokm_ssi,                       &
         cdr10m, cdr10m_n, cd10m_n, tau_fd_x, tau_fd_y,                 &
         rhogamu, rhogamv, f_ngstress,                                  &
    !     OUT variables required in IMP_SOLVER
         alpha1_sea, alpha1_sice, ashtf_sea, ashtf, bq_gb, bt_gb,       &
         dtrdz_charney_grid, rdz_charney_grid,                          &
         dtrdz_u, dtrdz_v, rdz_u, rdz_v,                                &
         k_blend_tq, k_blend_uv,uStarGBM,                               &
    !     OUT diagnostics (done after implicit solver)
         fqw, ftl, rib_gb, vshr, zht, zhnl, shallowc, cu_over_orog,     &
         bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6,   &
         bl_type_7, z0m_eff_gb, z0h_eff_gb, bl_w_var,                   &
    !     OUT diagnostics required for soil moisture nudging scheme :
         wt_ext,                                                        &
    !     OUT data required for tracer mixing :
         rho_aresist,aresist,resist_b,                                  &
    !     OUT variables required for mineral dust scheme
         r_b_dust,dust_flux,dust_emiss_frac,                            &
         u_s_t_tile,u_s_t_dry_tile,u_s_std_tile, kent, we_lim, t_frac, zrzi,&
         kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,              &
    !     OUT additional variables for JULES
         ftl_tile,radnet_sea,radnet_sice,rib_tile,rho_aresist_tile,     &
         aresist_tile,resist_b_tile,alpha1,ashtf_tile,fqw_tile,epot_tile,&
         fqw_ice,ftl_ice,fraca,resfs,resft,rhokh_tile,rhokh_sice,rhokh_sea,&
         z0hssi,z0h_tile,z0m_gb,z0mssi,z0m_tile,chr1p5m,chr1p5m_sice,smc,&
         gpp,npp,resp_p,g_leaf,gpp_ft,npp_ft,resp_p_ft,resp_s,resp_s_tot,&
         resp_w_ft,gc,canhc_tile,wt_ext_tile,flake,tile_index,tile_pts, &
         tile_frac,fsmc,rib_ssi, vshr_land,vshr_ssi,tstar_land,tstar_ssi,&
         dtstar_tile,dtstar_sea,dtstar_sice,hcons,emis_tile,emis_soil,  &
    !     OUT fields
         t1_sd,q1_sd,nbdsc,ntdsc,wstar,wthvs,uw0,vw0,taux_p,tauy_p,     &
         rhokm,rhokh,rhcpt                                              &
     )

    !-----------------------------------------------------------------------
    ! re-grid stuff onto u/v points
    ! The UM requires this regridding.  This will be dealt with 
    ! differently in LFRic (see #1059)
    !-----------------------------------------------------------------------
    flandg_u=flandg
    flandg_v=flandg
    rhokm_u=rhokm
    rhokm_v=rhokm
    rhokm_u_ssi=rhokm_ssi
    rhokm_v_ssi=rhokm_ssi
    rhokm_u_land=rhokm_land
    rhokm_v_land=rhokm_land
    f_ngstress_u=f_ngstress
    f_ngstress_v=f_ngstress
    taux_fd_u=tau_fd_x
    tauy_fd_v=tau_fd_y

    l_calc_at_p = .FALSE.  ! Flag for calling this routine on p-grid
      CALL bdy_expl3 (                                                 &
      ! IN grid related variables
           bl_levels, l_calc_at_p,                                     &
           udims, vdims, udims_s, vdims_s, pdims,                      &
           udims, vdims,                                               &
      ! IN SCM diags
           nSCMDpkgs,L_SCMDiags,                                       &
      ! IN variables used in flux calculations
           u, v, u_0, v_0, rhokm_u_land, rhokm_v_land, flandfac_u,     &
           flandfac_v, rhokm_u_ssi, rhokm_v_ssi, fseafac_u, fseafac_v, &
           flandg_u, flandg_v, zhnl, rdz_u, rdz_v, rhokm_u, rhokm_v,   &
           taux_fd_u, tauy_fd_v,                                       &
           rhogamu_u, rhogamv_v, f_ngstress_u, f_ngstress_v,           &
      ! OUT explicit momentum fluxes
           taux_land, tauy_land, taux_ssi, tauy_ssi, taux, tauy        &
                   )

    !-----------------------------------------------------------------------
    ! increments / fields with increments added
    ! Currently we're only modifying temperature.
    !-----------------------------------------------------------------------
    t_latest(1,1,:) = theta(1,1,1:model_levels)                         &
                    * exner_theta_levels(1,1,1:model_levels)
    q_latest(1,1,:)   = q(1,1,1:model_levels)
    qcl_latest(1,1,:) = qcl(1,1,1:model_levels)
    qcf_latest = 0.0
    cf_latest = 0.0
    cfl_latest = 0.0
    cff_latest = 0.0
    r_u = 0.0
    r_v = 0.0
    r_w = 0.0

    CALL NI_imp_ctl (                                                   &
    ! IN model dimensions.
            rhc_row_length, rhc_rows, land_field                        &
          , ntiles, bl_levels                                           &
          , cloud_levels, n_cca_levels                                  &
    ! IN Model switches
          , CycleNo                                                     &
    ! IN model Parameters
          , tr_vars                                                     &
    ! IN trig arrays
          , xx_cos_theta_latitude                                       &
    ! IN data fields.
          , p_layer_centres, p_layer_boundaries, rho_wet_rsq, rho_wet_tq&
          , u, v, w                                                     &
          , land_sea_mask, q, qcl, qcf, p_star, theta, exner_theta_levels&
    ! IN ancillary fields and fields needed to be kept from tstep to tstep
          , sil_orog_land, ho2r2_orog                                   &
          , ice_fract, di_ncat, ice_fract_ncat, k_sice, u_0, v_0, land_index&
          , cca, lcbase, ccb0, cct0                                     &
          , ls_rain, ls_snow, conv_rain, conv_snow, L_scrn, L_plsp      &
    ! IN variables required from BDY_LAYR
          , alpha1_sea, alpha1_sice, ashtf_sea, ashtf, bq_gb, bt_gb     &
          , dtrdz_charney_grid, rdz_charney_grid, dtrdz_u, dtrdz_v      &
          , rdz_u, rdz_v, cdr10m_u, cdr10m_v, cdr10m_n_u, cdr10m_n_v    &
          , cd10m_n_u, cd10m_n_v, z_theta                               &
          , k_blend_tq, k_blend_uv, uStarGBM, rhokm_u, rhokm_v          &
    ! IN diagnostics (started or from) BDY_LAYR
          , rib_gb,zlcl,zht,zhnl,dzh,qcl_inv_top,zh                     &
          , bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6 &
          , bl_type_7, z0m_gb, z0m_eff_gb, z0h_eff_gb                   &
          , ntml, cumulus, l_pc2_diag_sh_pts                            &
    ! IN logical for scm surface forcing
          , L_flux_bc                                                   &
    ! IN data required for tracer mixing :
          , rho_aresist,aresist,r_b_dust                                &
          , kent, we_lim, t_frac, zrzi                                  &
          , kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                  &
          , zhsc,z_rho,dust_flux,dust_emiss_frac                        &
          , u_s_t_tile,u_s_t_dry_tile,u_s_std_tile                      &
     ! IN additional variables for JULES. Now includes lai_ft, canht_ft.
          , tile_pts,tile_index,tile_frac,canopy                        &
          , alpha1,fraca,rhokh_tile,smc,chr1p5m,resfs,z0hssi,z0mssi     &
          , canhc_tile,flake,wt_ext_tile,lw_down,lai_ft,canht_ft        &
          , sw_tile,ashtf_tile,gc,aresist_tile                          &
          , resft,rhokh_sice,rhokh_sea,z0h_tile,z0m_tile                &
          , chr1p5m_sice                                                &
          , fland, flandg, flandg_u,flandg_v,vshr_land,vshr_ssi         &
          , emis_tile, t_soil, snow_tile, rib_ssi                       &
    ! IN JULES variables for STASH
          , gs,gpp,npp,resp_p,gpp_ft,npp_ft,resp_p_ft,resp_s            &
          , resp_s_tot,cs                                               &
          , rib_tile,fsmc,catch,g_leaf                                  &
          , co2_emits, co2flux                                          &
    !IN additional variables for soil moisture nudging scheme
          , wt_ext,                                                     &
    ! INOUT diagnostic info
           STASHwork3, STASHwork9                                       &
    ! SCM Diagnostics (dummy in full UM) & bl diags
          , nSCMDpkgs, L_SCMDiags, bl_diag, sf_diag                     &
    ! INOUT (Note ti and ti_gb are IN only if l_sice_multilayers=T)
          , TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans                      &
          , cca0,ccw0,cca0_2d,fqw,ftl,taux,tauy, rhokh                  &
          , fqw_ice,ftl_ice,dtstar_tile,dtstar_sea,dtstar_sice,ti       &
          , area_cloud_fraction, bulk_cloud_fraction                    &
          , T_latest, q_latest, qcl_latest, qcf_latest                  &
          , cf_latest, cfl_latest, cff_latest                           &
          , R_u, R_v, R_w, cloud_fraction_liquid, cloud_fraction_frozen &
          , sum_eng_fluxes,sum_moist_flux, rhcpt                        &
    ! INOUT tracer fields
          , aerosol, free_tracers,  resist_b,  resist_b_tile            &
          , dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6 &
          , drydep2, so2, dms, so4_aitken, so4_accu, so4_diss, nh3      &
          , soot_new, soot_aged, soot_cld, bmass_new, bmass_aged        &
          , bmass_cld, ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss&
          , co2, ozone_tracer                                           &
    ! INOUT additional variables for JULES
          , tstar_tile,fqw_tile,epot_tile,ftl_tile                      &
          , radnet_sea,radnet_sice,olr,tstar_sice_cat,tstar_ssi         &
          , tstar_sea,taux_land,taux_ssi,tauy_land,tauy_ssi,Error_code  &
    ! OUT fields
          , surf_ht_flux_land, zlcl_mixed                               &
          , theta_star_surf, qv_star_surf                               &
    ! OUT additional variables for JULES
          , tstar, ti_gb, ext, snowmelt,tstar_land,tstar_sice, ei_tile  &
          , ecan_tile,melt_tile, surf_htf_tile                          &
            )

    deallocate(bl_diag%q_incr)
    deallocate(bl_diag%qcl_incr)
    deallocate(bl_diag%qcf_incr)

    !-----------------------------------------------------------------------
    ! update main model prognostics
    ! Currently we're only modifying temperature.
    !-----------------------------------------------------------------------

    do k = 1, model_levels
      ! potential temperature increment on theta levels
      theta_inc(map_wth(1) + k) = t_latest(1,1,k)                       &
                                   / exner_theta_levels(1,1,k)          &
                                - theta_in_wth(map_wth(1)+k)
      ! water vapour on theta levels
      m_v(map_wth(1) + k)  = q_latest(1,1,k)
      ! cloud liquid water on theta levels
      m_cl(map_wth(1) + k) = qcl_latest(1,1,k)
    end do
    ! copy down lowest level to surface as done in UM
    theta_inc(map_wth(1) + 0) = theta_inc(map_wth(1) + 1)
    m_v(map_wth(1) + 0)  = m_v(map_wth(1) + 1)
    m_cl(map_wth(1) + 0) = m_cl(map_wth(1) + 1)

    ! update BL prognostics
    ! There is no need to loop over nlayers_2d, since here
    ! nlayers_2d=1
    tstar_2d(map_2d(1))=tstar(1,1)
    zh_2d(map_2d(1))=zh(1,1)
    z0msea_2d(map_2d(1))=z0msea(1,1)

end subroutine bl_code

end module bl_kernel_mod
