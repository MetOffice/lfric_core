!!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!
!> @brief ENUM containing different time steps
!> @details This enumarator contains all the timesteps available for an LFRic
!> field. It is used by vertical_dimensions_mod.f90 and field meta data
!> definition files for specifying the timestep that the field is output on
module time_step_enum_mod

  implicit none
  !> THE NUMBERS IN THIS FILE ARE ARBITARY. DO NOT ATTEMPT TO USE THEM IN ANY CODE
  enum, bind(c)
    enumerator :: STANDARD_TIMESTEP = 3001, &
                  RADIATION_TIMESTEP = 3002
  end enum
end module time_step_enum_mod