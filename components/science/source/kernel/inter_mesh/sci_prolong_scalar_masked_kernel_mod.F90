!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Perform the prolongation operation from a coarse grid field to a fine
!!        grid field, using a mask, for an intensive field.
!> @details Prolong the coarse grid correction into a subset of corresponding
!!          sub-cells on the fine grid, as defined by the masking field on the
!!          fine mesh. e.g. using a limited area mask field.
!!          Can only be used with finite elements at lowest-order.
!!          For a coarse mesh cell j, that contains Nf fine mesh cells
!!          fine(i) = mask(i) * coarse(j), for i=1,Nf
module sci_prolong_scalar_masked_kernel_mod

use constants_mod, only: i_def
use, intrinsic :: iso_fortran_env, only: real32, real64
use kernel_mod,    only: kernel_type
use argument_mod,  only: arg_type,                  &
                         GH_FIELD, GH_REAL,         &
                         GH_READ, GH_READWRITE,     &
                         ANY_DISCONTINUOUS_SPACE_1, &
                         ANY_DISCONTINUOUS_SPACE_2, &
                         GH_COARSE, GH_FINE, CELL_COLUMN

implicit none

private

type, public, extends(kernel_type) :: prolong_scalar_masked_kernel_type
   private
   type(arg_type) :: meta_args(3) = (/                                       &
        arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1, & ! fine_field
                                                  mesh_arg=GH_FINE),         &
        arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_2, & ! coarse_field
                                                  mesh_arg=GH_COARSE ),      &
        arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1, & ! mask_fine
                                                  mesh_arg=GH_FINE )         &
        /)
  integer :: operates_on = CELL_COLUMN
end type prolong_scalar_masked_kernel_type

public :: prolong_scalar_masked_kernel_code

  ! Generic interface for real32 and real64 types
  interface prolong_scalar_masked_kernel_code
    module procedure  &
      prolong_scalar_masked_code_real32, &
      prolong_scalar_masked_code_real64
  end interface

contains

  !> @brief Performs the masked injection-prolongation for scalar fields
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
  !> @param[in]     mask_fine                Mask field on fine grid. Only the
  !!                                         bottom layer of cells are used.
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid

  ! REAL32 PRECISION
  ! ==================
  subroutine prolong_scalar_masked_code_real32(   nlayers,                 &
                                                  cell_map,                &
                                                  ncell_fine_per_coarse_x, &
                                                  ncell_fine_per_coarse_y, &
                                                  ncell_fine,              &
                                                  fine_field,              &
                                                  coarse_field,            &
                                                  mask_fine,               &
                                                  ndf,                     &
                                                  undf_fine,               &
                                                  map_fine,                &
                                                  undf_coarse,             &
                                                  map_coarse )

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
    real(kind=real32),   intent(in)    :: coarse_field(undf_coarse)
    real(kind=real32),   intent(inout) :: fine_field(undf_fine)
    real(kind=real32),   intent(in)    :: mask_fine(undf_fine)

    integer(kind=i_def) :: df, k, x_idx, y_idx, top_df

    ! Loop over vertical layers from bottom to top
    ! Wtheta field with ndf=2: update the lower df and loops from k=0 to nlayers
    ! This uses the fact that for Wtheta, map(upper df)+ nlayers-1 = map(lower df) + nlayers
    ! W3 field with ndf=1: update the only (cell centre) df and loop from k=0 to nlayers-1
    ! Use the mask_fine from the bottom-layer (at present mask_fine is stored
    ! as a 3D field - but it may be stored as a 2D field in the future).
    df = 1
    top_df = nlayers - 2 + ndf

    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do k = 0, top_df
          fine_field( map_fine(df,cell_map(x_idx,y_idx)) + k ) =     &
              fine_field( map_fine(df,cell_map(x_idx,y_idx)) + k ) + &
              mask_fine( map_fine(df,cell_map(x_idx,y_idx)) ) *      &
              coarse_field( map_coarse(df)+k )
        end do
      end do
    end do

  end subroutine prolong_scalar_masked_code_real32

  ! REAL64 PRECISION
  ! ==================
  subroutine prolong_scalar_masked_code_real64(   nlayers,                 &
                                                  cell_map,                &
                                                  ncell_fine_per_coarse_x, &
                                                  ncell_fine_per_coarse_y, &
                                                  ncell_fine,              &
                                                  fine_field,              &
                                                  coarse_field,            &
                                                  mask_fine,               &
                                                  ndf,                     &
                                                  undf_fine,               &
                                                  map_fine,                &
                                                  undf_coarse,             &
                                                  map_coarse )

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
    real(kind=real64),   intent(in)    :: coarse_field(undf_coarse)
    real(kind=real64),   intent(inout) :: fine_field(undf_fine)
    real(kind=real64),   intent(in)    :: mask_fine(undf_fine)

    integer(kind=i_def) :: df, k, x_idx, y_idx, top_df

    ! Loop over vertical layers from bottom to top
    ! Wtheta field with ndf=2: update the lower df and loops from k=0 to nlayers
    ! This uses the fact that for Wtheta, map(upper df)+ nlayers-1 = map(lower df) + nlayers
    ! W3 field with ndf=1: update the only (cell centre) df and loop from k=0 to nlayers-1
    ! Use the mask_fine from the bottom-layer (at present mask_fine is stored
    ! as a 3D field - but it may be stored as a 2D field in the future).
    df = 1
    top_df = nlayers - 2 + ndf

    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do k = 0, top_df
          fine_field( map_fine(df,cell_map(x_idx,y_idx)) + k ) =     &
              fine_field( map_fine(df,cell_map(x_idx,y_idx)) + k ) + &
              mask_fine( map_fine(df,cell_map(x_idx,y_idx)) ) *      &
              coarse_field( map_coarse(df)+k )
        end do
      end do
    end do

  end subroutine prolong_scalar_masked_code_real64

end module sci_prolong_scalar_masked_kernel_mod
