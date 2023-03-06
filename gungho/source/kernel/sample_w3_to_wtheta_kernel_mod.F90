!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Obtains a Wtheta field from a W3 field by sampling at Wtheta DoFs.
!>
!> @details Takes a W3 field and returns a Wtheta field, whose values are
!>          obtained by interpolating the Wtheta values to the Wtheta DoFs.
!>          The values at the top and bottom boundaries are obtained by linear
!>          extrapolation.
!>
!>          Only designed to work with the lowest order finite elements.
!>
module sample_w3_to_wtheta_kernel_mod

  use argument_mod,            only : arg_type,          &
                                      GH_FIELD, GH_REAL, &
                                      GH_WRITE, GH_READ, &
                                      CELL_COLUMN
  use constants_mod,           only : r_def, i_def, r_double, r_single
  use fs_continuity_mod,       only : Wtheta, W3
  use kernel_mod,              only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: sample_w3_to_wtheta_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                 &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),     &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3)      &
         /)
    integer :: operates_on = CELL_COLUMN
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: sample_w3_to_wtheta_code

  ! Generic interface for real32 and real64 types
  interface sample_w3_to_wtheta_code
    module procedure  &
      sample_w3_to_wtheta_code_r_single, &
      sample_w3_to_wtheta_code_r_double
  end interface

contains

!> @brief Obtains a Wtheta field from a W3 field by sampling at W3 DoFs.
!! @param[in]     nlayers    Number of layers
!! @param[in,out] field_wt   Output Wtheta field
!! @param[in]     field_w3   Input W3 field
!! @param[in]     height_wt  Height above surface of Wtheta DoFs
!! @param[in]     height_w3  Height above surface of W3 DoFs
!! @param[in]     ndf_wt     Number of degrees of freedom per cell for Wtheta
!! @param[in]     undf_wt    Total number of degrees of freedom for Wtheta
!! @param[in]     map_wt     Dofmap for the cell at the base of the column for Wtheta
!! @param[in]     ndf_w3     Number of degrees of freedom per cell for W3
!! @param[in]     undf_w3    Total number of degrees of freedom for W3
!! @param[in]     map_w3     Dofmap for the cell at the base of the column for W3

! R_SINGLE PRECISION
! ==================
! height fields are always r_def
subroutine sample_w3_to_wtheta_code_r_single( nlayers,   &
                                              field_wt,  &
                                              field_w3,  &
                                              height_wt, &
                                              height_w3, &
                                              ndf_wt,    &
                                              undf_wt,   &
                                              map_wt,    &
                                              ndf_w3,    &
                                              undf_w3,   &
                                              map_w3     )

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in)    :: nlayers
  integer(kind=i_def),                     intent(in)    :: undf_wt, ndf_wt
  integer(kind=i_def),                     intent(in)    :: undf_w3, ndf_w3

  integer(kind=i_def), dimension(ndf_w3),  intent(in)    :: map_w3
  integer(kind=i_def), dimension(ndf_wt),  intent(in)    :: map_wt

  real(kind=r_single), dimension(undf_wt), intent(inout) :: field_wt
  real(kind=r_single), dimension(undf_w3), intent(in)    :: field_w3
  real(kind=r_def),    dimension(undf_wt), intent(in)    :: height_wt
  real(kind=r_def),    dimension(undf_w3), intent(in)    :: height_w3

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_single) :: weight_lower, weight_upper, weight_denom
  real(kind=r_single) :: top_value, log_top_value

  ! At top and bottom do linear extrapolation to get values on boundaries
  ! Bottom first
  weight_denom = real(height_w3(map_w3(1)+1) - height_w3(map_w3(1)), r_single)
  weight_upper = real(height_wt(map_wt(1)) - height_w3(map_w3(1)), r_single)
  weight_lower = real(height_w3(map_w3(1)+1) - height_wt(map_wt(1)), r_single)

  field_wt(map_wt(1)) = weight_upper / weight_denom * field_w3(map_w3(1)+1) &
                      + weight_lower / weight_denom * field_w3(map_w3(1))

  do k = 1, nlayers - 1

    weight_denom = real(height_w3(map_w3(1)+k) - height_w3(map_w3(1)+k-1), r_single)
    weight_upper = real(height_wt(map_wt(1)+k) - height_w3(map_w3(1)+k-1), r_single)
    weight_lower = real(height_w3(map_w3(1)+k) - height_wt(map_wt(1)+k), r_single)

    field_wt(map_wt(1)+k) = weight_upper / weight_denom * field_w3(map_w3(1)+k) &
                        + weight_lower / weight_denom * field_w3(map_w3(1)+k-1)

  end do

  ! Now top
  k = nlayers
  weight_denom = real(height_w3(map_w3(1)+k-1) - height_w3(map_w3(1)+k-2), r_single)
  weight_upper = real(height_wt(map_wt(1)+k) - height_w3(map_w3(1)+k-2), r_single)
  weight_lower = real(height_w3(map_w3(1)+k-1) - height_wt(map_wt(1)+k), r_single)

  if ( field_w3(map_w3(1)+k-1) > 0.0_r_single .and. &
       field_w3(map_w3(1)+k-2) > 0.0_r_single ) then
    ! Danger of getting a negative value from extrapolation, so perform the
    ! calculation for the top value in log space (which cannot be negative)
    log_top_value = weight_upper / weight_denom * log(field_w3(map_w3(1)+k-1)) &
                    + weight_lower / weight_denom * log(field_w3(map_w3(1)+k-2))
    top_value = exp(log_top_value)
  else
    ! If the field was originally negative, allow extrapolation to be negative
    top_value = weight_upper / weight_denom * field_w3(map_w3(1)+k-1) &
                + weight_lower / weight_denom * field_w3(map_w3(1)+k-2)
  end if

  field_wt(map_wt(1)+k) = top_value

end subroutine sample_w3_to_wtheta_code_r_single

! R_DOUBLE PRECISION
! ==================
! height fields are always r_def
subroutine sample_w3_to_wtheta_code_r_double( nlayers,   &
                                              field_wt,  &
                                              field_w3,  &
                                              height_wt, &
                                              height_w3, &
                                              ndf_wt,    &
                                              undf_wt,   &
                                              map_wt,    &
                                              ndf_w3,    &
                                              undf_w3,   &
                                              map_w3     )

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in)    :: nlayers
  integer(kind=i_def),                     intent(in)    :: undf_wt, ndf_wt
  integer(kind=i_def),                     intent(in)    :: undf_w3, ndf_w3

  integer(kind=i_def), dimension(ndf_w3),  intent(in)    :: map_w3
  integer(kind=i_def), dimension(ndf_wt),  intent(in)    :: map_wt

  real(kind=r_double), dimension(undf_wt), intent(inout) :: field_wt
  real(kind=r_double), dimension(undf_w3), intent(in)    :: field_w3
  real(kind=r_def),    dimension(undf_wt), intent(in)    :: height_wt
  real(kind=r_def),    dimension(undf_w3), intent(in)    :: height_w3

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_double) :: weight_lower, weight_upper, weight_denom
  real(kind=r_double) :: top_value, log_top_value

  ! At top and bottom do linear extrapolation to get values on boundaries
  ! Bottom first
  weight_denom = real(height_w3(map_w3(1)+1) - height_w3(map_w3(1)), r_double)
  weight_upper = real(height_wt(map_wt(1)) - height_w3(map_w3(1)), r_double)
  weight_lower = real(height_w3(map_w3(1)+1) - height_wt(map_wt(1)), r_double)

  field_wt(map_wt(1)) = weight_upper / weight_denom * field_w3(map_w3(1)+1) &
                      + weight_lower / weight_denom * field_w3(map_w3(1))

  do k = 1, nlayers - 1
    weight_denom = real(height_w3(map_w3(1)+k) - height_w3(map_w3(1)+k-1), r_double)
    weight_upper = real(height_wt(map_wt(1)+k) - height_w3(map_w3(1)+k-1), r_double)
    weight_lower = real(height_w3(map_w3(1)+k) - height_wt(map_wt(1)+k), r_double)

    field_wt(map_wt(1)+k) = weight_upper / weight_denom * field_w3(map_w3(1)+k) &
                        + weight_lower / weight_denom * field_w3(map_w3(1)+k-1)

  end do

  ! Now top
  k = nlayers
  weight_denom = real(height_w3(map_w3(1)+k-1) - height_w3(map_w3(1)+k-2), r_double)
  weight_upper = real(height_wt(map_wt(1)+k) - height_w3(map_w3(1)+k-2), r_double)
  weight_lower = real(height_w3(map_w3(1)+k-1) - height_wt(map_wt(1)+k), r_double)

  if ( field_w3(map_w3(1)+k-1) > 0.0_r_double .and. &
       field_w3(map_w3(1)+k-2) > 0.0_r_double ) then
    ! Danger of getting a negative value from extrapolation, so perform the
    ! calculation for the top value in log space (which cannot be negative)
    log_top_value = weight_upper / weight_denom * log(field_w3(map_w3(1)+k-1)) &
                    + weight_lower / weight_denom * log(field_w3(map_w3(1)+k-2))
    top_value = exp(log_top_value)
  else
    ! If the field was originally negative, allow extrapolation to be negative
    top_value = weight_upper / weight_denom * field_w3(map_w3(1)+k-1) &
                + weight_lower / weight_denom * field_w3(map_w3(1)+k-2)
  end if

  field_wt(map_wt(1)+k) = top_value

end subroutine sample_w3_to_wtheta_code_r_double

end module sample_w3_to_wtheta_kernel_mod
