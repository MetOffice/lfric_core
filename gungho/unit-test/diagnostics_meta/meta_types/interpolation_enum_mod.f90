!!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!
!> @brief ENUM containing different interpolation methods
!> @details This enumerator contains all the possible interpolation
!> methods for an LFRic field. It is used by meta_data_mod.f90 for
!> specifying the recommended interpolation type

module interpolation_enum_mod

  implicit none
  !> THE NUMBERS IN THIS FILE ARE ARBITARY.
  !> DO NOT ATTEMPT TO USE THEM IN ANY CODE
  enum, bind(c)

    enumerator :: BILINEAR = 4002

  end enum

end module interpolation_enum_mod