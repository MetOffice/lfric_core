!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Performs a simple insertion sort on a Wtheta field
!! @details Performs a sort on a Wtheta field in a column to remove static
!!          instability, leaving a field that is increasing in the vertical
!!          direction. This is done only above a particular height (and below
!!          that height the field is not sorted). This has been implemented in
!!          single and double precision.
!>          Only written for the lowest order elements.
module sci_sort_column_above_kernel_mod

  use argument_mod,      only: arg_type,             &
                               GH_FIELD, GH_SCALAR,  &
                               GH_REAL, GH_READ,     &
                               GH_READWRITE, CELL_COLUMN
  use constants_mod,     only: i_def, r_def
  use, intrinsic :: iso_fortran_env, only: real32, real64
  use fs_continuity_mod, only: Wtheta
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: sort_column_above_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                      &
         arg_type(GH_FIELD,  GH_REAL, GH_READWRITE, Wtheta), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,      Wtheta), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)               &
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: sort_column_above_code

  ! Generic interface for real32 and real64 types
  interface sort_column_above_code
    module procedure  &
      sort_column_above_code_real32, &
      sort_column_above_code_real64
  end interface
contains

!> @brief Performs a simple sort on a Wtheta field. There are single and double
!!        precision implementations.
!> @param[in]     nlayers                Number of layers
!> @param[in,out] theta                  Theta field
!> @param[in]     height_wth             Field containing heights of Wtheta DoFs
!> @param[in]     height_to_sort_above   The height to sort above
!> @param[in]     ndf_wth                Num of DoFs per cell for Wtheta
!> @param[in]     undf_wth               Num of DoFs per partition for Wtheta
!> @param[in]     map_wth                DoFmap for Wtheta

! REAL32 PRECISION
! ==================
subroutine sort_column_above_code_real32(nlayers,                   &
                                           theta,                     &
                                           height_wth,                &
                                           height_to_sort_above,      &
                                           ndf_wth, undf_wth, map_wth &
                                          )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
  real(kind=r_def),    intent(in) :: height_to_sort_above

  real(kind=real32), dimension(undf_wth), intent(inout) :: theta
  real(kind=r_def),    dimension(undf_wth), intent(in)    :: height_wth
  integer(kind=i_def), dimension(ndf_wth),  intent(in)    :: map_wth

  ! Internal variables
  integer(kind=i_def) :: k, kcnt, k_low
  real(kind=real32) :: theta_k

  ! Work out level to sort about
  k_low = nlayers + 1
  do k = 1, nlayers
    if (height_wth(map_wth(1)+k) >= height_to_sort_above) then
      k_low = k + 1
      exit
    end if
  end do

  do k = k_low, nlayers

    theta_k = theta(map_wth(1) + k)
    kcnt = k

    do while (theta(map_wth(1) + kcnt -1) > theta_k)
      theta(map_wth(1) + kcnt) = theta(map_wth(1) + kcnt -1)
      kcnt = kcnt - 1
      if (kcnt < k_low) exit

    end do
    theta(map_wth(1) + kcnt) = theta_k
  end do

end subroutine sort_column_above_code_real32

! REAL64 PRECISION
! ==================
subroutine sort_column_above_code_real64(nlayers,                   &
                                           theta,                     &
                                           height_wth,                &
                                           height_to_sort_above,      &
                                           ndf_wth, undf_wth, map_wth &
                                          )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
  real(kind=r_def),    intent(in) :: height_to_sort_above

  real(kind=real64), dimension(undf_wth), intent(inout) :: theta
  real(kind=r_def),    dimension(undf_wth), intent(in)    :: height_wth
  integer(kind=i_def), dimension(ndf_wth),  intent(in)    :: map_wth

  ! Internal variables
  integer(kind=i_def) :: k, kcnt, k_low
  real(kind=real64) :: theta_k

  ! Work out level to sort about
  k_low = nlayers + 1
  do k = 1, nlayers
    if (height_wth(map_wth(1)+k) >= height_to_sort_above) then
      k_low = k + 1
      exit
    end if
  end do

  do k = k_low, nlayers

    theta_k = theta(map_wth(1) + k)
    kcnt = k

    do while (theta(map_wth(1) + kcnt -1) > theta_k)
      theta(map_wth(1) + kcnt) = theta(map_wth(1) + kcnt -1)
      kcnt = kcnt - 1
      if (kcnt < k_low) exit

    end do
    theta(map_wth(1) + kcnt) = theta_k
  end do

end subroutine sort_column_above_code_real64

end module sci_sort_column_above_kernel_mod
