!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief The restriction operation from a fine W2V field to a coarse W2V field
!> @details Restrict the W2V field from a fine mesh into a W2V field on a coarse
!!          mesh. The fields are extensive -- i.e. this works on flux values and
!!          not pointwise values.
!!          This method is only designed for the lowest order W2V space.
module sci_restrict_w2v_kernel_mod

use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_WRITE,         &
                                   GH_COARSE, GH_FINE,        &
                                   ANY_DISCONTINUOUS_SPACE_2, &
                                   CELL_COLUMN
use constants_mod,           only: i_def, r_def
use, intrinsic :: iso_fortran_env, only: real32, real64
use fs_continuity_mod,       only: W2V
use kernel_mod,              only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!> Psy layer.
!>

type, public, extends(kernel_type) :: restrict_w2v_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                          &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, W2V, mesh_arg=GH_COARSE),         &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2,        &
                                                          mesh_arg=GH_FINE)    &
       /)
  integer :: operates_on = CELL_COLUMN
end type restrict_w2v_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

public :: restrict_w2v_code

  ! Generic interface for real32 and real64 types
  interface restrict_w2v_code
    module procedure  &
      restrict_w2v_code_real32, &
      restrict_w2v_code_real64
  end interface

contains

  !> @brief Restrict a fine W2V field to a coarse mesh
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] coarse_field             Coarse grid W2V field to compute
  !> @param[in]     fine_field               Fine grid  W2V field to restrict
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid

  ! REAL32 PRECISION
  ! ==================
  subroutine restrict_w2v_code_real32(nlayers,                 &
                                        cell_map,                &
                                        ncell_fine_per_coarse_x, &
                                        ncell_fine_per_coarse_y, &
                                        ncell_fine,              &
                                        coarse_field,            &
                                        fine_field,              &
                                        undf_coarse,             &
                                        map_coarse,              &
                                        ndf,                     &
                                        undf_fine,               &
                                        map_fine                 )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in) :: ncell_fine
    integer(kind=i_def), intent(in) :: ndf
    integer(kind=i_def), intent(in) :: undf_fine, undf_coarse

    ! Fields
    real(kind=real32), intent(inout) :: coarse_field(undf_coarse)
    real(kind=real32), intent(in)    :: fine_field(undf_fine)

    ! Maps
    integer(kind=i_def), intent(in) :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in) :: map_coarse(ndf)
    integer(kind=i_def), intent(in) :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)

    ! Internal variables
    integer(kind=i_def) :: df, k, x_idx, y_idx
    real(kind=real32) :: new_coarse(nlayers+1)

    !---------------------------------------------------------------------------
    ! Vertical components
    !---------------------------------------------------------------------------

    ! Only do bottom value of cell
    ! Loop over an extra layer to get the very top
    df = 1
    new_coarse(:) = 0.0_real32

    ! Build up 1D array of new coarse values for this column
    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do k = 0, nlayers
          new_coarse(k+1) = new_coarse(k+1) + fine_field(map_fine(df,cell_map(x_idx,y_idx))+k)
        end do
      end do
    end do

    ! Copy over values into coarse field
    do k = 0, nlayers
      coarse_field(map_coarse(df)+k) = new_coarse(k+1)
    end do

  end subroutine restrict_w2v_code_real32

  ! REAL64 PRECISION
  ! ==================
  subroutine restrict_w2v_code_real64(nlayers,                 &
                                        cell_map,                &
                                        ncell_fine_per_coarse_x, &
                                        ncell_fine_per_coarse_y, &
                                        ncell_fine,              &
                                        coarse_field,            &
                                        fine_field,              &
                                        undf_coarse,             &
                                        map_coarse,              &
                                        ndf,                     &
                                        undf_fine,               &
                                        map_fine                 )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in) :: ncell_fine
    integer(kind=i_def), intent(in) :: ndf
    integer(kind=i_def), intent(in) :: undf_fine, undf_coarse

    ! Fields
    real(kind=real64), intent(inout) :: coarse_field(undf_coarse)
    real(kind=real64), intent(in)    :: fine_field(undf_fine)

    ! Maps
    integer(kind=i_def), intent(in) :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in) :: map_coarse(ndf)
    integer(kind=i_def), intent(in) :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)

    ! Internal variables
    integer(kind=i_def) :: df, k, x_idx, y_idx
    real(kind=real64) :: new_coarse(nlayers+1)

    !---------------------------------------------------------------------------
    ! Vertical components
    !---------------------------------------------------------------------------

    ! Only do bottom value of cell
    ! Loop over an extra layer to get the very top
    df = 1
    new_coarse(:) = 0.0_real64

    ! Build up 1D array of new coarse values for this column
    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do k = 0, nlayers
          new_coarse(k+1) = new_coarse(k+1) + fine_field(map_fine(df,cell_map(x_idx,y_idx))+k)
        end do
      end do
    end do

    ! Copy over values into coarse field
    do k = 0, nlayers
      coarse_field(map_coarse(df)+k) = new_coarse(k+1)
    end do

  end subroutine restrict_w2v_code_real64

end module sci_restrict_w2v_kernel_mod
