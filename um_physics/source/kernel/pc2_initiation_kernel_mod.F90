!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to cloud scheme.

module pc2_initiation_kernel_mod

use argument_mod, only: arg_type,                                &
                        GH_FIELD, GH_READ, GH_WRITE,             &
                        GH_READWRITE, CELLS, ANY_SPACE_1
use fs_continuity_mod, only: WTHETA, W3

use kernel_mod,   only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: pc2_initiation_kernel_type
  private
  type(arg_type) :: meta_args(23) = (/                &
       arg_type(GH_FIELD,   GH_READ,    WTHETA),      & ! mv_wth
       arg_type(GH_FIELD,   GH_READ,    WTHETA),      & ! ml_wth
       arg_type(GH_FIELD,   GH_READ,    WTHETA),      & ! mi_wth
       arg_type(GH_FIELD,   GH_READ,    WTHETA),      & ! cfl_wth
       arg_type(GH_FIELD,   GH_READ,    WTHETA),      & ! cff_wth
       arg_type(GH_FIELD,   GH_READ,    WTHETA),      & ! bcf_wth
       arg_type(GH_FIELD,   GH_READ,    WTHETA),      & ! theta_wth
       arg_type(GH_FIELD,   GH_READ,    WTHETA),      & ! exner_wth
       arg_type(GH_FIELD,   GH_READ,    W3),          & ! exner_w3
       arg_type(GH_FIELD,   GH_READ,    WTHETA),      & ! mv_n_wth
       arg_type(GH_FIELD,   GH_READ,    WTHETA),      & ! ml_n_wth
       arg_type(GH_FIELD,   GH_READ,    WTHETA),      & ! theta_n_wth
       arg_type(GH_FIELD,   GH_READ,    WTHETA),      & ! exner_n_wth
       arg_type(GH_FIELD,   GH_READ,    ANY_SPACE_1), & ! zlcl_mixed
       arg_type(GH_FIELD,   GH_READ,    ANY_SPACE_1), & ! r_cumulus
       arg_type(GH_FIELD,   GH_WRITE,   WTHETA),      & ! dtheta_inc
       arg_type(GH_FIELD,   GH_WRITE,   WTHETA),      & ! dqv_inc_wth
       arg_type(GH_FIELD,   GH_WRITE,   WTHETA),      & ! dqcl_inc_wth
       arg_type(GH_FIELD,   GH_WRITE,   WTHETA),      & ! dqcf_inc_wth
       arg_type(GH_FIELD,   GH_WRITE,   WTHETA),      & ! dcfl_inc_wth
       arg_type(GH_FIELD,   GH_WRITE,   WTHETA),      & ! dcff_inc_wth
       arg_type(GH_FIELD,   GH_WRITE,   WTHETA),      & ! dbcf_inc_wth
       arg_type(GH_FIELD,   GH_READ,    WTHETA)       & ! rh_crit_wth
       /)
   integer :: iterates_over = CELLS
contains
  procedure, nopass :: pc2_initiation_code
end type

public pc2_initiation_code
contains

!> @brief Interface to pc2 initiation 
!> @details Calculates whether cloud should be created from clear-sky conditions
!>          or whether overcast conditions should be broken up.
!>          More info is in UMDP 30.
!> @param[in] nlayers         Number of layers
!> @param[in] mv_wth          Vapour mass mixing ratio
!> @param[in] ml_wth          Liquid cloud mass mixing ratio
!> @param[in] mi_wth          Liquid cloud mass mixing ratio
!> @param[in] cfl_wth         Liquid cloud fraction
!> @param[in] cff_wth         Ice cloud fraction
!> @param[in] bcf_wth         Bulk cloud fraction
!> @param[in] theta_wth       Potential temperature field
!> @param[in] exner_wth       Exner pressure in theta space
!> @param[in] exner_w3        Exner pressure in w3 space
!> @param[in] mv_n_wth        Start of timestep vapour mass mixing ratio
!> @param[in] ml_n_wth        Start of timestep liquid cloud mass mixing ratio
!> @param[in] theta_n_wth     Start of timestep theta in theta space
!> @param[in] exner_n_wth     Start of timestep exner in theta space
!> @param[in] zlcl_mixed      The height of the lifting condensation level in the mixed layer
!> @param[in] r_cumulus       A real number representing the logical cumulus flag
!> @param[out] dtheta_inc_wth Increment to theta in theta space
!> @param[out] dqv_inc_wth    Increment to water vapour in theta space
!> @param[out] dqcl_inc_wth   Increment to liquid water content in theta space
!> @param[out] dqcf_inc_wth   Increment to ice water content in theta space
!> @param[out] dcfl_inc_wth   Increment to liquid cloud fraction in theta space
!> @param[out] dcff_inc_wth   Increment to ice cloud fraction in theta space
!> @param[out] dbcf_inc_wth   Increment to bulk cloud fraction in theta space
!> @param[in] rh_crit_wth     Critical relative humidity in theta space
!> @param[in] ndf_wth         Number of degrees of freedom per cell for theta space
!> @param[in] undf_wth        Number unique of degrees of freedom  for theta space
!> @param[in] map_wth         Dofmap for the cell at the base of the column for theta space
!> @param[in] ndf_w3          Number of degrees of freedom per cell for density space
!> @param[in] undf_w3         Number unique of degrees of freedom  for density space
!> @param[in] map_w3          Dofmap for the cell at the base of the column for density space
!> @param[in] ndf_2d          Number of degrees of freedom per cell for density space
!> @param[in] undf_2d         Number of unique degrees of freedom for density space
!> @param[in] map_2d          Dofmap for the cell at the base of the column for density space

subroutine pc2_initiation_code( nlayers,                           &
                                mv_wth,                            &
                                ml_wth,                            &
                                mi_wth,                            &
                                cfl_wth,                           &
                                cff_wth,                           &
                                bcf_wth,                           &
                                theta_wth,                         &
                                exner_wth,                         &
                                exner_w3,                          &
                                ! Start of timestep values for RHt
                                mv_n_wth,                          &
                                ml_n_wth,                          &
                                theta_n_wth,                       &
                                exner_n_wth,                       &
                                zlcl_mixed,                        &
                                r_cumulus,                         &
                                ! Responses
                                dtheta_inc_wth,                    &
                                dqv_inc_wth,                       &
                                dqcl_inc_wth,                      &
                                dqcf_inc_wth,                      &
                                dcfl_inc_wth,                      &
                                dcff_inc_wth,                      &
                                dbcf_inc_wth,                      &
                                rh_crit_wth,                       &
                                ! Other
                                ndf_wth, undf_wth, map_wth,        &
                                ndf_w3,  undf_w3,  map_w3,         &
                                ndf_2d,  undf_2d,  map_2d   )

    use constants_mod, only: r_def, i_def, r_um, i_um

    !---------------------------------------
    ! UM modules
    !---------------------------------------

    use nlsizes_namelist_mod,       only: row_length, rows, model_levels
    use atm_step_local,             only: rhc_row_length, rhc_rows
    use pc2_initiation_ctl_mod,     only: pc2_initiation_ctl
    use planet_constants_mod,       only: p_zero, kappa, lcrcp
    use gen_phys_inputs_mod,        only: l_mr_physics

    ! Redirect routine names to avoid clash with existing qsat routines
    use qsat_mod, only: qsat_wat_new     => qsat_wat,                       &
                        qsat_wat_mix_new => qsat_wat_mix

    implicit none

    ! Arguments

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth , ndf_w3
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d

    real(kind=r_def), intent(in), dimension(undf_wth) :: mv_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: ml_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: mi_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: bcf_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: cfl_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: cff_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: theta_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: exner_wth
    real(kind=r_def), intent(in), dimension(undf_w3)  :: exner_w3

    real(kind=r_def), intent(in), dimension(undf_2d) :: zlcl_mixed
    real(kind=r_def), intent(in), dimension(undf_2d) :: r_cumulus

    logical, dimension(row_length,rows) :: l_cumulus

    ! Start of timestep values
    real(kind=r_def), intent(in),     dimension(undf_wth) :: mv_n_wth
    real(kind=r_def), intent(in),     dimension(undf_wth) :: ml_n_wth
    real(kind=r_def), intent(in),     dimension(undf_wth) :: theta_n_wth
    real(kind=r_def), intent(in),     dimension(undf_wth) :: exner_n_wth

    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_w3)  :: map_w3
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    ! The changes to the fields as a result
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dtheta_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dqv_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dqcl_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dqcf_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcfl_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcff_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dbcf_inc_wth

    real(kind=r_def), intent(in), dimension(undf_wth) :: rh_crit_wth

    real(r_um), dimension(row_length,rows,model_levels) ::   &
                  qv_work, qcl_work, qcf_work,               &
                  cfl_work, cff_work, bcf_work,              &
                  t_work, theta_work, rhts,                  &
                  t_incr, qv_incr, qcl_incr, qcf_incr,       &
                  cfl_incr, cff_incr, bcf_incr,              &
                  rhcpt, zeros

    real(r_um), dimension(row_length,rows,model_levels) ::   &
                  tlts, qtts, ptts, qsl_tl

    real(r_um), dimension(row_length,rows,model_levels) :: &
                  p_theta_levels

    real(r_um), dimension(row_length,rows,model_levels+1) ::   &
                  p_rho_levels

    real(r_um), dimension(row_length,rows) :: t_n, p_star

    integer(i_um) :: k

    ! Hardwired things for PC2
    !
    integer(i_um), dimension(row_length,rows) :: i_dummy
    integer(i_um), parameter :: nSCMDpkgs=15
    logical,       parameter :: l_scmdiags(nscmdpkgs) = .false.
    logical,       parameter :: calculate_increments  = .false.

    ! Convert cumulus flag from real to logical
    if (r_cumulus(map_2d(1))>0.5_r_def) then
      l_cumulus(1,1) = .true.
    else
      l_cumulus(1,1) = .false.
    endif

    i_dummy=1
    zeros=0.0_r_um
    !-----------------------------------------------------------------------
    ! Initialisation of prognostic variables and arrays
    !-----------------------------------------------------------------------

    p_star(1,1) = p_zero*(exner_wth(map_wth(1) + 0))     &
                                   **(1.0_r_def/kappa)

    p_rho_levels(1,1,model_levels+1) = 0.0_r_um

    do k = 1, model_levels

      ! Calculate temperature, this array will be updated.
      t_work(1,1,k)  = theta_wth(map_wth(1) + k) *                  &
                       exner_wth(map_wth(1) + k)

      ! Pressure at centre of theta levels
      p_theta_levels(1,1,k)    = p_zero*(exner_wth(map_wth(1) + k)) &
                                           **(1.0_r_def/kappa)

      ! pressure at layer boundaries
      p_rho_levels(1,1,k) = p_zero*( exner_w3(map_w3(1) + k-1))     &
                                           **(1.0_r_def/kappa)

      ! Moist prognostics
      qv_work(1,1,k)   = mv_wth(map_wth(1) + k)
      qcl_work(1,1,k)  = ml_wth(map_wth(1) + k)
      qcf_work(1,1,k)  = mi_wth(map_wth(1) + k)

      ! Critical relative humidity
      rhcpt(1,1,k)     = rh_crit_wth(map_wth(1) + k)

      ! Some start of timestep fields:
      ! Pressure on theta levels at start of time-step
      ptts(1,1,k) = p_zero*(exner_n_wth(map_wth(1) + k))     &
                                           **(1.0_r_def/kappa)

      ! Temperature at start of timestep
      t_n(1,1)    = theta_n_wth(map_wth(1) + k) *                &
                    exner_n_wth(map_wth(1) + k)

      ! Q total at start of time-step
      qtts(1,1,k) = mv_n_wth(map_wth(1) + k) + ml_n_wth(map_wth(1) + k)

      ! Liquid temperature
      tlts(1,1,k) = t_n(1,1) - ( lcrcp * ml_n_wth(map_wth(1) + k) )

    end do     ! k

    ! Calculate qsat(TL) with respect to liquid water, operate on whole column.
    if ( l_mr_physics ) then
      call qsat_wat_mix_new(qsl_tl(1,1,:), tlts(1,1,:), ptts(1,1,:), model_levels )
    else
      call qsat_wat_new(qsl_tl(1,1,:), tlts(1,1,:), ptts(1,1,:), model_levels )
    end if

    do k = 1, model_levels
      ! Total relative humidity (using start of time-step liquid temperature).
      rhts(1,1,k) = qtts(1,1,k) / qsl_tl(1,1,k)
      ! Recast LFRic cloud fractions onto cloud fraction work arrays.
      bcf_work(1,1,k) = bcf_wth(map_wth(1) + k)
      cfl_work(1,1,k) = cfl_wth(map_wth(1) + k)
      cff_work(1,1,k) = cff_wth(map_wth(1) + k)
    end do

    ! Initialize
    t_incr   = 0.0_r_um
    qv_incr  = 0.0_r_um
    qcl_incr = 0.0_r_um
    qcf_incr = 0.0_r_um
    bcf_incr = 0.0_r_um
    cfl_incr = 0.0_r_um
    cff_incr = 0.0_r_um

    CALL pc2_initiation_ctl(                               &
                            ! Dimensions of Rh crit array
                            rhc_row_length,                &
                            rhc_rows,                      &
                            ! Pass in zlcl_mixed=zero for NOW
                            zlcl_mixed(map_2d(1)),         &
                            ! Model switches
                            l_mr_physics,                  &
                            ! SCM diagnostics switches
                            nSCMDpkgs,                     &
                            L_SCMDiags,                    &
                            ! Primary fields passed IN/OUT
                            t_work,                        &
                            qv_work,                       &
                            qcl_work,                      &
                            qcf_work,                      &
                            bcf_work,                      &
                            cfl_work,                      &
                            cff_work,                      &
                            rhts,                          &
                            tlts,                          &
                            qtts,                          &
                            ptts,                          &
                            zeros,                         &
                            ! Primary fields passed IN
                            p_rho_levels,                  &
                            p_star,                        &
                            p_theta_levels,                &
                            i_dummy,                       &
                            l_cumulus,                     &
                            rhcpt,                         &
                            calculate_increments,          &
                            ! Output increments
                            t_incr,                        &
                            qv_incr,                       &
                            qcl_incr,                      &
                            qcf_incr,                      &
                            bcf_incr,                      &
                            cfl_incr,                      &
                            cff_incr                       )

    ! Recast back to LFRic space
    do k = 1, model_levels
      ! *_work arrays have been updated

      ! New theta found from new temperature.
      theta_work(1,1,k) = t_work(1,1,k) /                  &
                          exner_wth(map_wth(1) + k)
      ! All increments found from difference between the updated *_work values
      ! and the values that were intent in.
      dtheta_inc_wth(map_wth(1) + k) = theta_work(1,1,k)   &
                                     - theta_wth(map_wth(1) + k)
      !
      dqv_inc_wth (map_wth(1)+k) = qv_work(1,1,k)  - mv_wth(map_wth(1) + k)
      dqcl_inc_wth(map_wth(1)+k) = qcl_work(1,1,k) - ml_wth(map_wth(1) + k)
      dqcf_inc_wth(map_wth(1)+k) = qcf_work(1,1,k) - mi_wth(map_wth(1) + k)
      dcfl_inc_wth(map_wth(1)+k) = cfl_work(1,1,k) - cfl_wth(map_wth(1) + k)
      dcff_inc_wth(map_wth(1)+k) = cff_work(1,1,k) - cff_wth(map_wth(1) + k)
      dbcf_inc_wth(map_wth(1)+k) = bcf_work(1,1,k) - bcf_wth(map_wth(1) + k)

    end do

    dtheta_inc_wth(map_wth(1)+0) = dtheta_inc_wth(map_wth(1)+1)
    dqv_inc_wth   (map_wth(1)+0) = dqv_inc_wth   (map_wth(1)+1)
    dqcl_inc_wth  (map_wth(1)+0) = dqcl_inc_wth  (map_wth(1)+1)
    dqcf_inc_wth  (map_wth(1)+0) = dqcf_inc_wth  (map_wth(1)+1)
    dcfl_inc_wth  (map_wth(1)+0) = dcfl_inc_wth  (map_wth(1)+1)
    dcff_inc_wth  (map_wth(1)+0) = dcff_inc_wth  (map_wth(1)+1)
    dbcf_inc_wth  (map_wth(1)+0) = dbcf_inc_wth  (map_wth(1)+1)

end subroutine pc2_initiation_code

end module pc2_initiation_kernel_mod
