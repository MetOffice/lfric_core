!!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!
!> @brief ENUM containing different model levels
!> @details This enumerator contains all the possible model levels
!> for an LFRic field. It is used by vertical_dimensions_mod.f90 and
!> field meta data definition files for specifying model levels

module levels_enum_mod

  implicit none
  !> THE NUMBERS IN THIS FILE ARE ARBITARY.
  !> DO NOT ATTEMPT TO USE THEM IN ANY CODE
  enum, bind(c)
    enumerator :: BOTTOM_ATMOSPHERIC_LEVEL = 1001,                          &
                  TOP_ATMOSPHERIC_LEVEL = 1002,                             &
                  TOP_WET_LEVEL = 1003,                                     &
                  TOP_ATMOSPHERIC_LEVEL_MINUS_1 = 1004,                     &
                  BOTTOM_LEVEL_IN_BOUNDARY_LAYER = 1005,                    &
                  TOP_LEVEL_IN_BOUNDARY_LAYER = 1006,                       &
                  BOTTOM_LEVEL_ABOVE_BOUNDARY_LAYER = 1007,           &
                  BOTTOM_SOIL_LEVEL = 1008,                                 &
                  TOP_SOIL_LEVEL = 1009,                                    &
                  BOTTOM_TRACER_LEVEL = 1010,                               &
                  TOP_TRACER_LEVEL = 1011,                                  &
                  TOP_GRAVITY_WAVE_DRAG_LEVEL_PLUS_ONE = 1012,              &
                  BOTTOM_GRAVITY_WAVE_DRAG_LEVEL = 1013,                    &
                  TOP_GRAVITY_WAVE_DRAG_LEVEL = 1014,                       &
                  BOTTOM_LEVEL_IN_VERTICAL_DIFFUSION_ROUTINE = 1015,        &
                  TOP_LEVEL_IN_VERTICAL_DIFFUSION_ROUTINE_MINUS_ONE = 1016, &
                  TOP_LEVEL_IN_VERTICAL_DIFFUSION_ROUTINE = 1017,           &
                  TOP_LEVEL_IN_BOUNDARY_LAYER_MINUS_ONE = 1018,             &
                  TOP_LEVEL_OF_ATMOSPHERE_PLUS_ONE = 1019,                  &
                  BOTTOM_SOIL_LEVEL_PLUS_ONE = 1020,                        &
                  TOP_OZONE_LEVEL = 1021,                                   &
                  NUMBER_OF_ATMOS_LEVELS_MULTIPLY_SW_BANDS = 1022,          &
                  NUMBER_OF_ATMOS_LEVELS_PLUS_ONE_MULTIPLY_SW_BANDS = 1023, &
                  NUMBER_OF_WET_LEVELS_MULTIPLY_SW_BANDS = 1024,            &
                  NUMBER_OF_ATMOS_LEVELS_MULTIPLY_LW_BANDS = 1025,          &
                  NUMBER_OF_ATMOS_LEVELS_PLUS_ONE_MULTIPLY_LW_BANDS = 1026, &
                  NUMBER_OF_WET_LEVELS_MULTIPLY_LW_BANDS = 1027,            &
                  NUMBER_OF_SW_RADIATION_BANDS = 1028,                      &
                  NUMBER_OF_LW_RADIATION_BANDS = 1029,                      &
                  NUMBER_OF_CLOUDY_LEVELS = 1030

  end enum
end module levels_enum_mod