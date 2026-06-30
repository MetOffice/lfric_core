!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief The restriction operation for mass fields.
!> @details Restrict a mass field on a fine grid to a coarse grid field
!!          This kernel only works for the lowest-order W3 and Wtheta spaces

module sci_restrict_mass_w3_kernel_mod

use constants_mod,           only: i_def
use, intrinsic :: iso_fortran_env, only: real32, real64
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_WRITE,         &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_2, &
                                   GH_COARSE, GH_FINE, CELL_COLUMN

implicit none

private

type, public, extends(kernel_type) :: restrict_mass_w3_kernel_type
   private
   type(arg_type) :: meta_args(2) = (/                                   &
        arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1, &
                                              mesh_arg=GH_COARSE),       &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2, &
                                              mesh_arg=GH_FINE  )        &
        /)
  integer :: operates_on = CELL_COLUMN
end type restrict_mass_w3_kernel_type

public :: restrict_mass_w3_kernel_code

  ! Generic interface for real32 and real64 types
  interface restrict_mass_w3_kernel_code
    module procedure  &
      restrict_mass_w3_code_real32, &
      restrict_mass_w3_code_real64
  end interface

contains

  !> @brief The restriction operation for mass fields.
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_y  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] coarse_field             Coarse grid mass field to write to
  !> @param[in]     fine_field               Fine grid mass field to restrict
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid

  ! REAL32 PRECISION
  ! ==================
  subroutine restrict_mass_w3_code_real32(                            &
                                             nlayers,                 &
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
                                             map_fine)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in)    :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)
    integer(kind=i_def), intent(in)    :: ncell_fine
    integer(kind=i_def), intent(in)    :: ndf
    integer(kind=i_def), intent(in)    :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in)    :: map_coarse(ndf)
    integer(kind=i_def), intent(in)    :: undf_fine, undf_coarse
    real(kind=real32),   intent(inout) :: coarse_field(undf_coarse)
    real(kind=real32),   intent(in)    :: fine_field(undf_fine)

    integer(kind=i_def) :: df, k, x_idx, y_idx, top_df
    real(kind=real32)   :: coarse_value(nlayers-1+ndf)

    ! Assume lowest order W3 or Wtheta space
    df = 1
    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_df = nlayers - 2 + ndf
    coarse_value(:) = 0.0_real32

    ! Build up 1D array of new coarse values for this column
    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do k = 0, top_df
          coarse_value(k+1) = coarse_value(k+1) + &
            fine_field(map_fine(df,cell_map(x_idx,y_idx))+k)
        end do
      end do
    end do

    ! Copy over values into coarse field
    do k = 0, top_df
      coarse_field(map_coarse(df) + k) = coarse_value(k+1)
    end do

  end subroutine restrict_mass_w3_code_real32

  ! REAL64 PRECISION
  ! ==================
  subroutine restrict_mass_w3_code_real64(                            &
                                             nlayers,                 &
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
                                             map_fine)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in)    :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)
    integer(kind=i_def), intent(in)    :: ncell_fine
    integer(kind=i_def), intent(in)    :: ndf
    integer(kind=i_def), intent(in)    :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in)    :: map_coarse(ndf)
    integer(kind=i_def), intent(in)    :: undf_fine, undf_coarse
    real(kind=real64),   intent(inout) :: coarse_field(undf_coarse)
    real(kind=real64),   intent(in)    :: fine_field(undf_fine)

    integer(kind=i_def) :: df, k, x_idx, y_idx, top_df
    real(kind=real64)   :: coarse_value(nlayers-1+ndf)

    ! Assume lowest order W3 or Wtheta space
    df = 1
    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_df = nlayers - 2 + ndf
    coarse_value(:) = 0.0_real64

    ! Build up 1D array of new coarse values for this column
    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do k = 0, top_df
          coarse_value(k+1) =  coarse_value(k+1) + &
            fine_field(map_fine(df,cell_map(x_idx,y_idx))+k)
        end do
      end do
    end do

    ! Copy over values into coarse field
    do k = 0, top_df
      coarse_field(map_coarse(df) + k) = coarse_value(k+1)
    end do

  end subroutine restrict_mass_w3_code_real64


end module sci_restrict_mass_w3_kernel_mod
