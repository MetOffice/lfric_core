!!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!
!> @brief ENUM containing different values of positive
!> @details This enumerator contains all the values for the attribute positive
!> of an LFRic field. It is used by vertical_dimensions_mod.f90 and
!> field meta data definition files for specifying the positive direction

module positive_enum_mod

  implicit none
  !> THE NUMBERS IN THIS FILE ARE ARBITARY.
  !> DO NOT ATTEMPT TO USE THEM IN ANY CODE
  enum, bind(c)
    enumerator :: POSITIVE_UP = 2001, &
                  POSITIVE_DOWN = 2002
  end enum
end module positive_enum_mod