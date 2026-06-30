!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Applies an operator to take a field from W3 to shifted W3.
!> @details Applies a pre-computed operator to map a W3 field to the shifted W3
!!          function space. The operator is stored in two W3 fields, which
!!          represent contributions from the lower and upper halves of cells.
!!          This is only designed to worked for the lowest-order elements.
module sci_apply_w3_to_sh_w3_kernel_mod

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
  type, public, extends(kernel_type) :: apply_w3_to_sh_w3_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                      &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3), & ! rhs_w3_sh
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! field_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! detj_shifted
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3)                         & ! detj_prime
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: apply_w3_to_sh_w3_code

  ! Generic interface for real32 and real64 types
  interface apply_w3_to_sh_w3_code
    module procedure  &
      apply_w3_to_sh_w3_code_real32, &
      apply_w3_to_sh_w3_code_real64
  end interface

contains

!> @brief Applies an operator to take a field from W3 to shifted W3.
!> @param[in] nlayers_sh  Number of layers in the vertically-shifted mesh.
!> @param[in,out] rhs_w3_sh Right hand side vector to compute.
!> It is a field in w3 shifted.
!> @param[in] field_w3 The original wind in w3 from which to compute the RHS.
!> @param[in] detj_shifted The volume of shifted mesh cells stored in W3 shifted.
!> @param[in] detj_prime The volume of prime mesh cells stored in W3.
!> @param[in] ndf_w3_sh Number of degrees of freedom per cell for w3 shifted
!> @param[in] undf_w3_sh Number of (local) unique degrees of freedom for w3 shifted
!> @param[in] map_w3_sh Dofmap for the cell at the base of the column for w3 shifted
!> @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!> @param[in] undf_w3 Number of (local) unique degrees of freedom for w3
!> @param[in] map_w3 Dofmap for the cell at the base of the column for w3

! REAL32 PRECISION
! ==================
subroutine apply_w3_to_sh_w3_code_real32(                  &
                                            nlayers_sh,    &
                                            rhs_w3_sh,     &
                                            field_w3,      &
                                            detj_shifted,  &
                                            detj_prime,    &
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

  real(kind=real32),      dimension(undf_w3_sh), intent(inout) :: rhs_w3_sh
  real(kind=real32),      dimension(undf_w3),       intent(in) :: field_w3
  real(kind=real32),      dimension(undf_w3_sh),    intent(in) :: detj_shifted
  real(kind=real32),      dimension(undf_w3),       intent(in) :: detj_prime

  ! Internal variables
  integer(kind=i_def) :: k

  ! Assume lowest order so only a single DoF per cell
  ! Bottom boundary value
  rhs_w3_sh(map_w3_sh(1)) = 0.5_real32 * field_w3(map_w3(1)) &
                    * detj_prime(map_w3(1)) / detj_shifted(map_w3_sh(1))

  ! All interior levels
  do k = 1, nlayers_sh - 2
    rhs_w3_sh(map_w3_sh(1)+k) = &
      (0.5_real32 * detj_prime(map_w3(1)+k) * field_w3(map_w3(1)+k) +    &
       0.5_real32 * detj_prime(map_w3(1)+k-1) * field_w3(map_w3(1)+k-1)) &
       / detj_shifted(map_w3_sh(1)+k)
  end do

  ! Top boundary value
  k = nlayers_sh - 1
  rhs_w3_sh(map_w3_sh(1)+k) = 0.5_real32 * field_w3(map_w3(1) + k-1) &
                    * detj_prime(map_w3(1)+k-1) / detj_shifted(map_w3_sh(1)+k)

end subroutine apply_w3_to_sh_w3_code_real32

! REAL64 PRECISION
! ==================
subroutine apply_w3_to_sh_w3_code_real64(                  &
                                            nlayers_sh,    &
                                            rhs_w3_sh,     &
                                            field_w3,      &
                                            detj_shifted,  &
                                            detj_prime,    &
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

  real(kind=real64),      dimension(undf_w3_sh), intent(inout) :: rhs_w3_sh
  real(kind=real64),      dimension(undf_w3),       intent(in) :: field_w3
  real(kind=real64),      dimension(undf_w3_sh),    intent(in) :: detj_shifted
  real(kind=real64),      dimension(undf_w3),       intent(in) :: detj_prime

  ! Internal variables
  integer(kind=i_def) :: k

  ! Assume lowest order so only a single DoF per cell
  ! Bottom boundary value
  rhs_w3_sh(map_w3_sh(1)) = 0.5_real64 * field_w3(map_w3(1)) &
                    * detj_prime(map_w3(1)) / detj_shifted(map_w3_sh(1))

  ! All interior levels
  do k = 1, nlayers_sh - 2
    rhs_w3_sh(map_w3_sh(1)+k) = &
      (0.5_real64 * detj_prime(map_w3(1)+k) * field_w3(map_w3(1)+k) +    &
       0.5_real64 * detj_prime(map_w3(1)+k-1) * field_w3(map_w3(1)+k-1)) &
       / detj_shifted(map_w3_sh(1)+k)
  end do

  ! Top boundary value
  k = nlayers_sh - 1
  rhs_w3_sh(map_w3_sh(1)+k) = 0.5_real64 * field_w3(map_w3(1) + k-1) &
                    * detj_prime(map_w3(1)+k-1) / detj_shifted(map_w3_sh(1)+k)

end subroutine apply_w3_to_sh_w3_code_real64

end module sci_apply_w3_to_sh_w3_kernel_mod
