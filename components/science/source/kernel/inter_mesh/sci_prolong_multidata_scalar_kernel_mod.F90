!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Perform the unweighted injection-prolongation operation of a coarse
!!        grid scalar multidata field to a fine grid multidata field
!> @details Prolong the coarse grid correction into all cells on the fine grid
!!          that are contained in that coarse grid cell. This is an injection,
!!          so values for the fine field are simple assigned from the coarse
!!          field, and no weighting is used.
!!          This kernel only works for the lowest-order W3 and Wtheta spaces
!!          and for column first multidata fields.


module sci_prolong_multidata_scalar_kernel_mod

use constants_mod,           only: i_def
use, intrinsic :: iso_fortran_env, only: real32, real64
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_READWRITE,     &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_2, &
                                   GH_COARSE, GH_FINE, CELL_COLUMN

implicit none

private

type, public, extends(kernel_type) :: prolong_multidata_scalar_kernel_type
   private
   type(arg_type) :: meta_args(2) = (/                                       &
        arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1, &
                                                  mesh_arg=GH_FINE),         &
        arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2, &
                                                  mesh_arg=GH_COARSE )       &
        /)
  integer :: operates_on = CELL_COLUMN
end type prolong_multidata_scalar_kernel_type

public :: prolong_multidata_scalar_kernel_code

  ! Generic interface for real32 and real64 types
  interface prolong_multidata_scalar_kernel_code
    module procedure  &
      prolong_multidata_scalar_code_real32, &
      prolong_multidata_scalar_code_real64
  end interface

contains

  !> @brief Performs the unweighted injection-prolongation for scalar multidata fields
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] fine_field               The fine grid field to write to
  !> @param[in]     coarse_field             Coarse grid field to prolong
  !> @param[in]     ndata                    Number fields in multidata field
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid

  ! REAL32 PRECISION
  ! ==================
  subroutine prolong_multidata_scalar_code_real32(               &
                                          nlayers,                 &
                                          cell_map,                &
                                          ncell_fine_per_coarse_x, &
                                          ncell_fine_per_coarse_y, &
                                          ncell_fine,              &
                                          fine_field,              &
                                          coarse_field,            &
                                          ndata,                   &
                                          ndf,                     &
                                          undf_fine,               &
                                          map_fine,                &
                                          undf_coarse,             &
                                          map_coarse)

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
    real(kind=real32), intent(in)    :: coarse_field(undf_coarse)
    real(kind=real32), intent(inout) :: fine_field(undf_fine)
    integer(kind=i_def), intent(in)    :: ndata

    integer(kind=i_def) :: k, x_idx, y_idx, k_start, i, top_df
    integer(kind=i_def), parameter :: df = 1 ! Lowest order function space

    ! Assume lowest order W3 or Wtheta space
    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_df = nlayers - 2 + ndf

    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do i = 1, ndata
          k_start = (i-1)*(top_df+1)+1
          do k = 0, top_df
            fine_field(map_fine(df,cell_map(x_idx,y_idx))+k+k_start-1) =   &
              fine_field(map_fine(df,cell_map(x_idx,y_idx))+k+k_start-1) + &
              coarse_field(map_coarse(df)+k+k_start-1)
          end do
        end do
      end do
    end do

  end subroutine prolong_multidata_scalar_code_real32

  ! REAL64 PRECISION
  ! ==================
  subroutine prolong_multidata_scalar_code_real64(               &
                                          nlayers,                 &
                                          cell_map,                &
                                          ncell_fine_per_coarse_x, &
                                          ncell_fine_per_coarse_y, &
                                          ncell_fine,              &
                                          fine_field,              &
                                          coarse_field,            &
                                          ndata,                   &
                                          ndf,                     &
                                          undf_fine,               &
                                          map_fine,                &
                                          undf_coarse,             &
                                          map_coarse)

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
    real(kind=real64), intent(in)    :: coarse_field(undf_coarse)
    real(kind=real64), intent(inout) :: fine_field(undf_fine)
    integer(kind=i_def), intent(in)    :: ndata

    integer(kind=i_def) :: k, x_idx, y_idx, k_start, i, top_df
    integer(kind=i_def), parameter :: df = 1 ! Lowest order function space

    ! Assume lowest order W3 or Wtheta space
    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_df = nlayers - 2 + ndf

    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do i = 1, ndata
          k_start = (i-1)*(top_df+1)+1
          do k = 0, top_df
            fine_field(map_fine(df,cell_map(x_idx,y_idx))+k+k_start-1) =   &
              fine_field(map_fine(df,cell_map(x_idx,y_idx))+k+k_start-1) + &
              coarse_field(map_coarse(df)+k+k_start-1)
          end do
        end do
      end do
    end do

  end subroutine prolong_multidata_scalar_code_real64


end module sci_prolong_multidata_scalar_kernel_mod
