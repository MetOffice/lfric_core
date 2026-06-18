!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   The weighted mapping from a higher to a lower-order W2 space.
!> @details Take a W2 space which is higher-order in the horizontal but
!!          shares the same number of levels, and map it to a lowest-order
!!          W2 space on a finer mesh.

module sci_map_w2_fe_to_fv_kernel_mod

use argument_mod,          only: arg_type,                  &
                                 GH_FIELD, GH_REAL,         &
                                 GH_INTEGER,                &
                                 GH_READ, GH_READWRITE,     &
                                 ANY_DISCONTINUOUS_SPACE_1, &
                                 ANY_DISCONTINUOUS_SPACE_2, &
                                 ANY_DISCONTINUOUS_SPACE_3, &
                                 ANY_DISCONTINUOUS_SPACE_4, &
                                 GH_COARSE, GH_FINE, CELL_COLUMN
use constants_mod,         only: i_def, r_def, IMDI
use kernel_mod,            only: kernel_type
use reference_element_mod, only: E, N, T

implicit none

private

type, public, extends(kernel_type) :: map_w2_fe_to_fv_kernel_type
   private
   type(arg_type) :: meta_args(5) = (/                                          &
        arg_type(GH_FIELD, GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1, &
                 mesh_arg=GH_FINE),                                             &
        arg_type(GH_FIELD, GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2, &
                 mesh_arg=GH_COARSE ),                                          &
        arg_type(GH_FIELD, GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_3, &
                 mesh_arg=GH_COARSE ),                                          &
        arg_type(GH_FIELD, GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_4, &
                 mesh_arg=GH_COARSE ),                                          &
        arg_type(GH_FIELD, GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_4, &
                 mesh_arg=GH_COARSE )                                           &
        /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: map_w2_fe_to_fv_code
end type map_w2_fe_to_fv_kernel_type

public :: map_w2_fe_to_fv_code

contains

  !> @brief Performs the mapping between W2 fields of different orders
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_y  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] fine_field               Lowest-order fine grid field to
  !!                                         map to
  !> @param[in]     coarse_field             Higher-order coarse field to map
  !!                                         from
  !> @param[in]     weights                  Weights for the mapping from
  !!                                         coarse to fine grid
  !> @param[in]     ndf_fine                 Num of DoFs per cell on the fine
  !!                                         grid
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  !> @param[in]     ndf_coarse               Num of DoFs per cell on the coarse
  !!                                         grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  !> @param[in]     undf_weights             Num of DoFs on the weights field
  !> @param[in]     map_weights              DoFmap of cells on the weights
  !!                                         field

  subroutine map_w2_fe_to_fv_code(nlayers,                 &
                                  cell_map,                &
                                  ncell_fine_per_coarse_x, &
                                  ncell_fine_per_coarse_y, &
                                  ncell_fine,              &
                                  fine_field,              &
                                  coarse_field,            &
                                  weights,                 &
                                  face_selector_ew,        &
                                  face_selector_ns,        &
                                  ndf_fine,                &
                                  undf_fine,               &
                                  map_fine,                &
                                  ndf_coarse,              &
                                  undf_coarse,             &
                                  map_coarse,              &
                                  undf_weights,            &
                                  map_weights,             &
                                  undf_w3_2d,              &
                                  map_w3_2d)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in)    :: ncell_fine
    integer(kind=i_def), intent(in)    :: ndf_fine, ndf_coarse
    integer(kind=i_def), intent(in)    :: undf_fine
    integer(kind=i_def), intent(in)    :: undf_coarse
    integer(kind=i_def), intent(in)    :: undf_weights
    integer(kind=i_def), intent(in)    :: undf_w3_2d

    ! Fields
    real(kind=r_def),    intent(in)    :: coarse_field(undf_coarse)
    real(kind=r_def),    intent(inout) :: fine_field(undf_fine)
    real(kind=r_def),    intent(in)    :: weights(undf_weights)
    integer(kind=i_def), intent(in)    :: face_selector_ew(undf_w3_2d)
    integer(kind=i_def), intent(in)    :: face_selector_ns(undf_w3_2d)

    ! Maps
    integer(kind=i_def), intent(in) :: cell_map(ncell_fine_per_coarse_x, &
                                                ncell_fine_per_coarse_y)
    integer(kind=i_def), intent(in) :: map_fine(ndf_fine, ncell_fine)
    integer(kind=i_def), intent(in) :: map_coarse(ndf_coarse)
    integer(kind=i_def), intent(in) :: map_weights(ndf_coarse)
    integer(kind=i_def), intent(in) :: map_w3_2d(undf_w3_2d)

    ! Internal variables
    real(kind=r_def) :: new_fine

    integer(kind=i_def) :: k, x_idx, y_idx, df_c, j, df_f

    integer(kind=i_def) :: map_idx(ndf_fine, ncell_fine_per_coarse_x, &
                                   ncell_fine_per_coarse_y)

    !---------------------------------------------------------------------------
    ! Redefine optimised mapping used in weights kernel
    !---------------------------------------------------------------------------
    j = 1
    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do df_f = 1, ndf_fine
          if (df_f == E .and. x_idx /= ncell_fine_per_coarse_x) then
            map_idx(df_f, x_idx, y_idx) = IMDI
          else if (df_f == N .and. y_idx /= 1) then
            map_idx(df_f, x_idx, y_idx) = IMDI
          else
            map_idx(df_f, x_idx, y_idx) = j
            j = j + 1
          end if
        end do
      end do
    end do

    !---------------------------------------------------------------------------
    ! Calculate the new fine field
    !---------------------------------------------------------------------------
    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do df_f = 1, ndf_fine

          ! Only evaluate East and North dofs at the East and North edges
          ! of the cell to avoid overcounting
          if (df_f == E .and. x_idx /= ncell_fine_per_coarse_x) cycle
          if (df_f == N .and. y_idx /= 1) cycle

          ! Only set dofs prescribed by the face selectors
          if (df_f == E .and. (face_selector_ew(map_w3_2d(1)) == 0             &
                               .or. face_selector_ew(map_w3_2d(1)) == 1)) cycle
          if (df_f == N .and. (face_selector_ns(map_w3_2d(1)) == 0             &
                               .or. face_selector_ns(map_w3_2d(1)) == 1)) cycle

          do k = 0, nlayers-1
            new_fine = 0.0_r_def
            ! Only evaluate the top dofs at the top of the column
            if (df_f == T .and. k /= nlayers-1) cycle

            j = map_idx(df_f, x_idx, y_idx) - 1

            ! Calculate the new value
            do df_c = 1, ndf_coarse
              new_fine = new_fine                                              &
                + coarse_field(map_coarse(df_c)+k)*weights(map_weights(df_c) + j)
            end do

            fine_field(map_fine(df_f,cell_map(x_idx,y_idx))+k) = new_fine
          end do
        end do
      end do
    end do

  end subroutine map_w2_fe_to_fv_code

end module sci_map_w2_fe_to_fv_kernel_mod
