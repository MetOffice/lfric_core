! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_stash_to_lfric_map_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64
! lfric modules
USE constants_mod, ONLY: str_def, imdi
USE log_mod,  ONLY : log_event, log_scratch_space, LOG_LEVEL_ERROR

IMPLICIT NONE

PRIVATE

PUBLIC :: lfricinp_init_stash_to_lfric_map, get_field_name,                    &
          get_lfric_field_kind, w2h_field, w3_field, w3_field_2d,              &
          w3_soil_field, wtheta_field

INTEGER(KIND=int64), PARAMETER :: max_lfric_field_names = 500
CHARACTER(LEN=str_def) :: field_name(max_lfric_field_names)
INTEGER(KIND=int64), SAVE :: field_counter = 0

! The get_index array is dimensioned to the maximum number of stashcodes
! possible in the UM. The array index will correspond to the
! stashcode. First two digits are the section code, remaining 3 digits
! are the items code. E.g. Section code 0, item code 3 will be at
! get_index(3), section code 2, item code 21 would be get_index(2021)
! The contents of the array will an index to the character array of lfric
! field names. This intermediate step is used to avoid having a hugely
! oversized character array
INTEGER(KIND=int64) :: get_index(99999) = int(imdi, int64)

INTEGER(KIND=int64), PARAMETER :: w2h_field     = 1
INTEGER(KIND=int64), PARAMETER :: w3_field      = 2
INTEGER(KIND=int64), PARAMETER :: w3_field_2d   = 3
INTEGER(KIND=int64), PARAMETER :: w3_soil_field = 4
INTEGER(KIND=int64), PARAMETER :: wtheta_field  = 5

CONTAINS

SUBROUTINE lfricinp_init_stash_to_lfric_map()

! Description:
!  Map stashcode to lfric field name

USE lfricinp_stashmaster_mod, ONLY: &
    stashcode_u, stashcode_v, stashcode_theta, stashcode_soil_moist,          &
    stashcode_soil_temp, stashcode_tstar, stashcode_bl_depth, stashcode_orog, &
    stashcode_ozone, stashcode_w, stashcode_can_conduct,                      &
    stashcode_unfrozen_soil, stashcode_frozen_soil, stashcode_snow_soot_tile, &
    stashcode_can_water_tile, stashcode_rgrain, stashcode_tstar_tile,         &
    stashcode_snow_tile, stashcode_snow_grnd, stashcode_area_cf,              &
    stashcode_bulk_cf, stashcode_liquid_cf, stashcode_frozen_cf, stashcode_zw,&
    stashcode_fsat, stashcode_fwetl, stashcode_sthzw,                         &
    stashcode_snowdep_grd_tile, stashcode_snowpack_bk_dens,                   &
    stashcode_nsnow_layrs_tiles, stashcode_snow_laythk_tiles,                 &
    stashcode_snow_ice_tile, stashcode_snow_liq_tile, stashcode_snow_T_tile,  &
    stashcode_snow_grnsiz_tiles, stashcode_dry_rho, stashcode_mv,             &
    stashcode_mcl, stashcode_mcf, stashcode_mr, stashcode_ddmfx,              &
    stashcode_tstar_sea, stashcode_tstar_sice, stashcode_sea_ice_temp,        &
    stashcode_z0, stashcode_q, stashcode_qcf, stashcode_qcl, stashcode_qrain, &
    stashcode_rhor2, stashcode_lsm, stashcode_icefrac, stashcode_icethick,    &
    stashcode_total_aero, stashcode_z0h_tile, stashcode_dust1_mmr,            &
    stashcode_dust2_mmr, stashcode_ls_snow_rate, stashcode_conv_rain_rate,    &
    stashcode_qt , stashcode_exner
USE lfricinp_regrid_options_mod, ONLY: winds_on_w3

IMPLICIT NONE

field_name(:) = TRIM('unset')
field_counter = 0

! PLEASE KEEP THIS LIST OF SUPPORTED STASHCODES IN NUMERICAL ORDER
IF (winds_on_w3) THEN
  CALL map_field_name(stashcode_u, 'ew_wind')                        ! stash 2
  CALL map_field_name(stashcode_v, 'ns_wind')                        ! stash 3
ELSE
  CALL map_field_name(stashcode_u, 'h_wind')                         ! stash 2&3
                                                 !combined into single W2H field
  CALL map_field_name(stashcode_v, 'h_wind')
ENDIF
CALL map_field_name(stashcode_theta, 'theta')                        ! stash 4
CALL map_field_name(stashcode_soil_moist, 'soil_moisture')           ! stash 9
CALL map_field_name(stashcode_q, 'q')                                ! stash 10
CALL map_field_name(stashcode_qcf, 'qcf')                            ! stash 12
CALL map_field_name(stashcode_soil_temp, 'soil_temperature')         ! stash 20
CALL map_field_name(stashcode_tstar, 'tstar')                        ! stash 24
CALL map_field_name(stashcode_bl_depth, 'zh')                        ! stash 25
CALL map_field_name(stashcode_lsm, 'land_mask')                      ! stash 30
CALL map_field_name(stashcode_icefrac, 'icefrac')                    ! stash 31
CALL map_field_name(stashcode_icethick, 'icethick')                  ! stash 32
CALL map_field_name(stashcode_orog, 'surface_altitude')              ! stash 33
CALL map_field_name(stashcode_ozone, 'ozone')                        ! stash 60
CALL map_field_name(stashcode_z0, 'z0msea')                          ! stash 26
CALL map_field_name(stashcode_sea_ice_temp, 'sea_ice_temperature')   ! stash 49
CALL map_field_name(stashcode_total_aero, 'total_aero')              ! stash 90
CALL map_field_name(stashcode_w, 'upward_wind')                      ! stash 150
CALL map_field_name(stashcode_ls_snow_rate, 'ls_snow_rate')          ! stash 187
CALL map_field_name(stashcode_conv_rain_rate, 'conv_rain_rate')      ! stash 188
CALL map_field_name(stashcode_can_conduct, 'surface_conductance')    ! stash 213
CALL map_field_name(stashcode_unfrozen_soil, 'unfrozen_soil_moisture')
                                                                     ! stash 214
CALL map_field_name(stashcode_frozen_soil, 'frozen_soil_moisture')   ! stash 215
CALL map_field_name(stashcode_snow_soot_tile, 'snow_soot')           ! stash 221
CALL map_field_name(stashcode_can_water_tile, 'tile_canopy_water')   ! stash 229
CALL map_field_name(stashcode_rgrain, 'tile_snow_rgrain')            ! stash 231
CALL map_field_name(stashcode_tstar_tile, 'tile_temperature')        ! stash 233
CALL map_field_name(stashcode_snow_tile , 'tile_snow_mass')          ! stash 240
CALL map_field_name(stashcode_snow_grnd, 'tile_snow_under_canopy')   ! stash 242
CALL map_field_name(stashcode_z0h_tile, 'z0h_tile')                  ! stash 246
CALL map_field_name(stashcode_rhor2, 'rho_r2')                       ! stash 253
CALL map_field_name(stashcode_qcl, 'qcl')                            ! stash 254
CALL map_field_name(stashcode_exner, 'exner')                            ! stash 255
CALL map_field_name(stashcode_area_cf, 'area_fraction')              ! stash 265
CALL map_field_name(stashcode_bulk_cf, 'bulk_fraction')              ! stash 266
CALL map_field_name(stashcode_liquid_cf, 'liquid_fraction')          ! stash 267
CALL map_field_name(stashcode_frozen_cf, 'frozen_fraction')          ! stash 268
CALL map_field_name(stashcode_qrain, 'qrain')                        ! stash 272
CALL map_field_name(stashcode_fsat, 'soil_sat_frac')                 ! stash 279
CALL map_field_name(stashcode_fwetl, 'soil_wet_frac')                ! stash 280
CALL map_field_name(stashcode_zw, 'water_table')                     ! stash 278
CALL map_field_name(stashcode_sthzw, 'wetness_under_soil')           ! stash 281
CALL map_field_name(stashcode_snowdep_grd_tile, 'tile_snow_depth')   ! stash 376
CALL map_field_name(stashcode_snowpack_bk_dens, 'tile_snowpack_density')
                                                                     ! stash 377
CALL map_field_name(stashcode_nsnow_layrs_tiles, 'tile_n_snow_layers')
                                                                     ! stash 380
CALL map_field_name(stashcode_snow_laythk_tiles,                    &! stash 381
     'tile_snow_layer_thickness')
CALL map_field_name(stashcode_snow_liq_tile, 'tile_snow_layer_liq_mass')
                                                                     ! stash 383
CALL map_field_name(stashcode_snow_ice_tile, 'tile_snow_layer_ice_mass')
                                                                     ! stash 382
CALL map_field_name(stashcode_snow_T_tile, 'tile_snow_layer_temp')   ! stash 384
CALL map_field_name(stashcode_snow_grnsiz_tiles,                    &! stash 386
     'tile_snow_layer_rgrain')
CALL map_field_name(stashcode_dry_rho, 'rho')                        ! stash 389
CALL map_field_name(stashcode_mv, 'm_v')                             ! stash 391
CALL map_field_name(stashcode_mcl, 'm_cl')                           ! stash 392
CALL map_field_name(stashcode_mcf, 'm_cf')                           ! stash 393
CALL map_field_name(stashcode_mr, 'm_r')                             ! stash 394
CALL map_field_name(stashcode_dust1_mmr, 'dust1_mmr')                ! stash 431
CALL map_field_name(stashcode_dust2_mmr, 'dust2_mmr')                ! stash 432
CALL map_field_name(stashcode_ddmfx , 'dd_mf_cb')                    ! stash 493
CALL map_field_name(stashcode_tstar_sea, 'tstar_sea')                ! stash 507
CALL map_field_name(stashcode_tstar_sice, 'tstar_sea_ice')           ! stash 508
CALL map_field_name(stashcode_qt, 'qt')                              ! stash 16207
! PLEASE KEEP THIS LIST OF SUPPORTED STASHCODES IN NUMERICAL ORDER

END SUBROUTINE lfricinp_init_stash_to_lfric_map


SUBROUTINE map_field_name(stashcode, lfric_field_name)

IMPLICIT NONE

INTEGER(KIND=int64), INTENT(IN) :: stashcode
CHARACTER(LEN=*), INTENT(IN) :: lfric_field_name

field_counter = field_counter + 1

IF (field_counter >  max_lfric_field_names) THEN
  WRITE(log_scratch_space, '(A)') "field_counter is greater than" // &
       " max_lfric_field_names. Recompile to increase."
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
ELSE
  field_name(field_counter) = lfric_field_name
  get_index(stashcode) = field_counter
END IF

END SUBROUTINE map_field_name


FUNCTION get_field_name(stashcode) RESULT(name)

IMPLICIT NONE

INTEGER(KIND=int64), INTENT(IN) :: stashcode
CHARACTER(LEN=str_def) :: name ! Result

IF (get_index(stashcode) == imdi) THEN
  WRITE(log_scratch_space, '(A,I0,A)') "Stashcode ", stashcode, &
       " has not been mapped to lfric field name"
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
ELSE
  name = TRIM(field_name(get_index(stashcode)))
END IF

END FUNCTION get_field_name


FUNCTION get_lfric_field_kind(stashcode) RESULT(lfric_field_kind)

USE lfricinp_regrid_options_mod, ONLY: winds_on_w3
USE lfricinp_stashmaster_mod,    ONLY: get_stashmaster_item, levelt,       &
                                       rho_levels, theta_levels,           &
                                       single_level, soil_levels,          &
                                       stashcode_u, stashcode_v
IMPLICIT NONE

INTEGER(KIND=int64), INTENT(IN) :: stashcode
INTEGER(KIND=int64) :: lfric_field_kind, level_code

level_code = get_stashmaster_item(stashcode, levelt)

IF ((stashcode == stashcode_u .OR. stashcode == stashcode_v) .AND.     &
    (.NOT. winds_on_w3)) THEN
  lfric_field_kind = w2h_field

ELSE IF (level_code == rho_levels) THEN
  lfric_field_kind = w3_field

ELSE IF (level_code == theta_levels) THEN
  lfric_field_kind = wtheta_field

ELSE IF (level_code == single_level) THEN
  lfric_field_kind = w3_field_2d

ELSE IF (level_code == soil_levels) THEN
  lfric_field_kind = w3_soil_field

ELSE

  WRITE(log_scratch_space, '(A,I0,A)') "Stashcode ", stashcode, &
       " is not mapped to a lfric field type"
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)

END IF

END FUNCTION get_lfric_field_kind

END MODULE lfricinp_stash_to_lfric_map_mod
