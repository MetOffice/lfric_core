!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Performs a simple insert sort on reference theta field to remove any
!>        static instability.
!>        Only written for the lowest order elements.
module sci_sort_ref_kernel_mod

  use argument_mod,      only: arg_type,          &
                               GH_FIELD, GH_REAL, &
                               GH_READWRITE, CELL_COLUMN
  use constants_mod,     only: i_def
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
  type, public, extends(kernel_type) :: sort_ref_kernel_type
    private
    type(arg_type) :: meta_args(1) = (/                    &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta) &
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: sort_ref_code

  ! Generic interface for real32 and real64 types
  interface sort_ref_code
    module procedure  &
      sort_ref_code_real32, &
      sort_ref_code_real64
  end interface
contains

!> @brief The subroutine which is called directly by the psy layer
!! @param[in] nlayers Number of layers
!! @param[in,out] theta_ref Theta reference state
!! @param[in] ndf_wth Number of degrees of freedom per cell for Wtheta
!! @param[in] undf_wth The number of unique degrees of freedom for Wtheta
!! @param[in] map_wth Dofmap for the cell at the base of the column for Wtheta

! REAL32 PRECISION
! ==================
subroutine sort_ref_code_real32(  nlayers,                   &
                                  theta_ref,                 &
                                  ndf_wth, undf_wth, map_wth &
                                 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth

  real(kind=real32),   dimension(undf_wth), intent(inout) :: theta_ref
  integer(kind=i_def), dimension(ndf_wth),  intent(in)    :: map_wth

  ! Internal variables
  integer(kind=i_def) :: k, kcnt
  real(kind=real32)   :: theta_k

  do k = 1, nlayers

    theta_k=theta_ref(map_wth(1) + k)
    kcnt=k

    do while (theta_ref(map_wth(1) + kcnt -1) > theta_k)
      theta_ref(map_wth(1) + kcnt) = theta_ref(map_wth(1) + kcnt -1)
      kcnt = kcnt - 1
      if (kcnt == 0) exit

    end do
    theta_ref(map_wth(1) + kcnt) = theta_k
  end do

end subroutine sort_ref_code_real32

! REAL64 PRECISION
! ==================
subroutine sort_ref_code_real64(  nlayers,                   &
                                  theta_ref,                 &
                                  ndf_wth, undf_wth, map_wth &
                                 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth

  real(kind=real64),   dimension(undf_wth), intent(inout) :: theta_ref
  integer(kind=i_def), dimension(ndf_wth),  intent(in)    :: map_wth

  ! Internal variables
  integer(kind=i_def) :: k, kcnt
  real(kind=real64)  :: theta_k

  do k = 1, nlayers

    theta_k=theta_ref(map_wth(1) + k)
    kcnt=k

    do while (theta_ref(map_wth(1) + kcnt -1) > theta_k)
      theta_ref(map_wth(1) + kcnt) = theta_ref(map_wth(1) + kcnt -1)
      kcnt = kcnt - 1
      if (kcnt == 0) exit

    end do
    theta_ref(map_wth(1) + kcnt) = theta_k
  end do

end subroutine sort_ref_code_real64

end module sci_sort_ref_kernel_mod
