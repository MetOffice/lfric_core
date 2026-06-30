!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Shifts a mass field from W3 to shifted W3.
!> @details Map a W3 mass field to the shifted W3 function space.
!!          This is only designed to worked for the lowest-order elements.
module sci_shift_mass_w3_kernel_mod

  use argument_mod,            only : arg_type,                  &
                                      GH_FIELD, GH_REAL,         &
                                      GH_READ, GH_WRITE,         &
                                      ANY_DISCONTINUOUS_SPACE_3, &
                                      CELL_COLUMN
  use constants_mod,           only : i_def
  use, intrinsic :: iso_fortran_env, only: real32, real64
  use fs_continuity_mod,       only : W3

  use kernel_mod,              only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: shift_mass_w3_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                      &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3),   & ! mass_shifted
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3)                           & ! mass_prime
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: shift_mass_w3_code

  ! Generic interface for real32 and real64 types
  interface shift_mass_w3_code
    module procedure  &
      shift_mass_w3_code_real32, &
      shift_mass_w3_code_real64
  end interface

contains

!> @brief Maps a W3 mass field from the prime mesh to the shifted mesh
!> @param[in]     nlayers_sh    Number of layers in the vertically-shifted mesh.
!> @param[in,out] mass_shifted  Mass field on the shifted mesh
!> @param[in]     mass_prime    Mass field on the prime mesh
!> @param[in]     ndf_w3_sh     Num DoFs per cell for W3 shifted
!> @param[in]     undf_w3_sh    Num DoFs for this partition for W3 shifted
!> @param[in]     map_w3_sh     Map of DoFs in the bottom layer for W3 shifted
!> @param[in]     ndf_w3        Num DoFs per cell for W3 prime
!> @param[in]     undf_w3       Num DoFs for this partition for W3 prime
!> @param[in]     map_w3        Map of DoFs in the bottom layer for W3 prime

! REAL32 PRECISION
! ==================
subroutine shift_mass_w3_code_real32(                  &
                                        nlayers_sh,    &
                                        mass_shifted,  &
                                        mass_prime,    &
                                        ndf_w3_sh,     &
                                        undf_w3_sh,    &
                                        map_w3_sh,     &
                                        ndf_w3,        &
                                        undf_w3,       &
                                        map_w3         &
                                      )

  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in) :: nlayers_sh
  integer(kind=i_def),                           intent(in) :: ndf_w3_sh, ndf_w3
  integer(kind=i_def),                           intent(in) :: undf_w3_sh, undf_w3
  integer(kind=i_def), dimension(ndf_w3_sh),     intent(in) :: map_w3_sh
  integer(kind=i_def), dimension(ndf_w3),        intent(in) :: map_w3

  real(kind=real32),      dimension(undf_w3_sh), intent(inout) :: mass_shifted
  real(kind=real32),      dimension(undf_w3),       intent(in) :: mass_prime

  ! Internal variables
  integer(kind=i_def) :: k

  ! Assume lowest order so only a single DoF per cell
  ! Bottom boundary value
  mass_shifted(map_w3_sh(1)) = 0.5_real32 * mass_prime(map_w3(1))

  ! All interior levels
  do k = 1, nlayers_sh - 2
    mass_shifted(map_w3_sh(1)+k) = &
      0.5_real32 * (mass_prime(map_w3(1)+k) + mass_prime(map_w3(1)+k-1))
  end do

  ! Top boundary value
  k = nlayers_sh - 1
  mass_shifted(map_w3_sh(1)+k) = 0.5_real32 * mass_prime(map_w3(1)+k-1)

end subroutine shift_mass_w3_code_real32

! REAL64 PRECISION
! ==================
subroutine shift_mass_w3_code_real64(                  &
                                        nlayers_sh,    &
                                        mass_shifted,  &
                                        mass_prime,    &
                                        ndf_w3_sh,     &
                                        undf_w3_sh,    &
                                        map_w3_sh,     &
                                        ndf_w3,        &
                                        undf_w3,       &
                                        map_w3         &
                                      )

  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in) :: nlayers_sh
  integer(kind=i_def),                           intent(in) :: ndf_w3_sh, ndf_w3
  integer(kind=i_def),                           intent(in) :: undf_w3_sh, undf_w3
  integer(kind=i_def), dimension(ndf_w3_sh),     intent(in) :: map_w3_sh
  integer(kind=i_def), dimension(ndf_w3),        intent(in) :: map_w3

  real(kind=real64),      dimension(undf_w3_sh), intent(inout) :: mass_shifted
  real(kind=real64),      dimension(undf_w3),       intent(in) :: mass_prime

  ! Internal variables
  integer(kind=i_def) :: k

  ! Assume lowest order so only a single DoF per cell
  ! Bottom boundary value
  mass_shifted(map_w3_sh(1)) = 0.5_real64 * mass_prime(map_w3(1))

  ! All interior levels
  do k = 1, nlayers_sh - 2
    mass_shifted(map_w3_sh(1)+k) = &
      0.5_real64 * (mass_prime(map_w3(1)+k) + mass_prime(map_w3(1)+k-1))
  end do

  ! Top boundary value
  k = nlayers_sh - 1
  mass_shifted(map_w3_sh(1)+k) = 0.5_real64 * mass_prime(map_w3(1)+k-1)

end subroutine shift_mass_w3_code_real64

end module sci_shift_mass_w3_kernel_mod
