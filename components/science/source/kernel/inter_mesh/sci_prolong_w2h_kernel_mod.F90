!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Perform the prolongation operation from a coarse grid W2h field to a
!!        fine grid W2h field
!> @details Prolong the W2h coarse grid correction into all cells in a W2h field
!!          on a fine grid. The fine grid cells are contained in that coarse
!!          grid cell. Values are found by multiplying by weights. DoFs internal
!!          to the coarse cell are treated differently to those DoFs on the
!!          edges of the coarse cell:
!!          - internal DoFs: values interpolated from coarse cell values
!!          - edge DoFs: divide coarse cell flux by fine cell area fraction

module sci_prolong_w2h_kernel_mod

use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_WRITE,         &
                                   GH_COARSE, GH_FINE,        &
                                   ANY_SPACE_2, CELL_COLUMN,  &
                                   ANY_DISCONTINUOUS_SPACE_2
use constants_mod,           only: i_def, r_def, r_single, r_double
use kernel_mod,              only: kernel_type
use reference_element_mod,   only: W, S, E, N

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!> Psy layer.
!>

type, public, extends(kernel_type) :: prolong_w2h_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                           &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2,         &
                                                          mesh_arg=GH_FINE ),   &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_2, mesh_arg=GH_COARSE ), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2,         &
                                                          mesh_arg=GH_FINE )    &
       /)
  integer :: operates_on = CELL_COLUMN
end type prolong_w2h_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

public :: prolong_w2h_kernel_code

  ! Generic interface for real32 and real64 types
  interface prolong_w2h_kernel_code
    module procedure  &
      prolong_w2h_kernel_code_r_single, &
      prolong_w2h_kernel_code_r_double
  end interface

contains

  !> @brief Subroutine to perform the W2h prolongation operation
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] fine_field               Fine grid  W2h field to restrict
  !> @param[in]     coarse_field             Coarse grid W2h field to compute
  !> @param[in]     weights                  A fine grid W2h field containing
  !!                                         weights for the mapping process
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid

  ! R_SINGLE PRECISION
  ! ==================
  subroutine prolong_w2h_kernel_code_r_single(nlayers,                 &
                                              cell_map,                &
                                              ncell_fine_per_coarse_x, &
                                              ncell_fine_per_coarse_y, &
                                              ncell_fine,              &
                                              fine_field,              &
                                              coarse_field,            &
                                              weights,                 &
                                              ndf,                     &
                                              undf_fine,               &
                                              map_fine,                &
                                              undf_coarse,             &
                                              map_coarse               )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in) :: ncell_fine
    integer(kind=i_def), intent(in) :: ndf
    integer(kind=i_def), intent(in) :: undf_fine, undf_coarse

    ! Fields
    real(kind=r_single), intent(inout) :: fine_field(undf_fine)
    real(kind=r_single), intent(in)    :: coarse_field(undf_coarse)
    real(kind=r_def),    intent(in)    :: weights(undf_fine)

    ! Maps
    integer(kind=i_def), intent(in) :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in) :: map_coarse(ndf)
    integer(kind=i_def), intent(in) :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)

    ! Internal arguments
    integer(kind=i_def) :: df, opp_df, k, x_idx, y_idx
    real(kind=r_single)    :: edge_weight, internal_weight

    !===========================================================================
    ! Horizontal components
    !===========================================================================

    do k = 0, nlayers-1

      ! The rows and columns forming the cell map match the arrangment
      ! of fine cells within the coarse cell

      ! These are aligned as follows with the LFRic directions:
      !         N
      !   |--------------|
      !   |    row 1     |
      !   |c            c|
      !   |o            o|
      ! W |l            l| E
      !   |              |
      !   |1           nx|
      !   |    row ny    |
      !   |--------------|
      !          S

      !-------------------------------------------------------------------------
      ! Edges of coarse cell
      !-------------------------------------------------------------------------

      ! The weights here describe the flux through the edge of the fine cell
      ! This should be the fractional area of the coarse cell edge taken up
      ! by the fine cell. If fine cells have equal area then this will be the
      ! reciprocal of the number of fine cells along that edge

      do x_idx = 1, ncell_fine_per_coarse_x
        ! N edge is first row of cell map
        y_idx = 1
        df = N
        fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) =                                        &
                                    real(weights(map_fine(df,cell_map(x_idx,y_idx))+k), r_single) &
                                    * coarse_field(map_coarse(df)+k)

        ! S edge is last row of cell map
        y_idx = ncell_fine_per_coarse_y
        df = S
        fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) =                                        &
                                    real(weights(map_fine(df,cell_map(x_idx,y_idx))+k), r_single) &
                                    * coarse_field(map_coarse(df)+k)
      end do

      do y_idx = 1, ncell_fine_per_coarse_y
        ! W edge is first column of cell map
        x_idx = 1
        df = W
        fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) =                                        &
                                    real(weights(map_fine(df,cell_map(x_idx,y_idx))+k), r_single) &
                                    * coarse_field(map_coarse(df)+k)

        ! E edge is last column of cell map
        x_idx = ncell_fine_per_coarse_x
        df = E
        fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) =                                        &
                                    real(weights(map_fine(df,cell_map(x_idx,y_idx))+k), r_single) &
                                    * coarse_field(map_coarse(df)+k)
      end do

      !-------------------------------------------------------------------------
      ! DoFs internal to coarse cell
      !-------------------------------------------------------------------------

      ! For internal cells, we interpolate. The weights at these DoFs give the
      ! relative weighting of the WEST and NORTH components. The contributions
      ! from the east and south components are given by 1 minus those weights.
      ! All must be combined with the area fraction weights to give the true
      ! values. As the fine cells are internal to a coarse cell, they will all
      ! have the same orientation as one another.

      ! Do West/East DoFs
      ! These loops will only execute if there are internal DoFs
      df = W
      opp_df = E
      do y_idx = 1, ncell_fine_per_coarse_y
        edge_weight = real(weights(map_fine(df,cell_map(1,y_idx))+k), r_single)

        do x_idx = 2, ncell_fine_per_coarse_x
          internal_weight = real(weights(map_fine(df,cell_map(x_idx,y_idx))+k), r_single)

          fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) =                      &
               edge_weight * ( internal_weight * coarse_field(map_coarse(df)+k) + &
               (1.0_r_single - internal_weight) * coarse_field(map_coarse(opp_df)+k) )
        end do
      end do

      ! Do North/South DoFs
      ! These loops will only execute if there are internal DoFs
      df = N
      opp_df = S
      do x_idx = 1, ncell_fine_per_coarse_x
        edge_weight = real(weights(map_fine(df,cell_map(x_idx,1))+k), r_single)

        do y_idx = 2, ncell_fine_per_coarse_y
          internal_weight = real(weights(map_fine(df,cell_map(x_idx,y_idx))+k), r_single)

          fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) =                      &
               edge_weight * ( internal_weight * coarse_field(map_coarse(df)+k) + &
               (1.0_r_single - internal_weight) * coarse_field(map_coarse(opp_df)+k) )

        end do
      end do

    end do

  end subroutine prolong_w2h_kernel_code_r_single

  ! R_DOUBLE PRECISION
  ! ==================
  subroutine prolong_w2h_kernel_code_r_double(nlayers,                 &
                                              cell_map,                &
                                              ncell_fine_per_coarse_x, &
                                              ncell_fine_per_coarse_y, &
                                              ncell_fine,              &
                                              fine_field,              &
                                              coarse_field,            &
                                              weights,                 &
                                              ndf,                     &
                                              undf_fine,               &
                                              map_fine,                &
                                              undf_coarse,             &
                                              map_coarse               )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in) :: ncell_fine
    integer(kind=i_def), intent(in) :: ndf
    integer(kind=i_def), intent(in) :: undf_fine, undf_coarse

    ! Fields
    real(kind=r_double), intent(inout) :: fine_field(undf_fine)
    real(kind=r_double), intent(in)    :: coarse_field(undf_coarse)
    real(kind=r_def),    intent(in)    :: weights(undf_fine)

    ! Maps
    integer(kind=i_def), intent(in) :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in) :: map_coarse(ndf)
    integer(kind=i_def), intent(in) :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)

    ! Internal arguments
    integer(kind=i_def) :: df, opp_df, k, x_idx, y_idx
    real(kind=r_double)    :: edge_weight, internal_weight

    !===========================================================================
    ! Horizontal components
    !===========================================================================

    do k = 0, nlayers-1

      ! The rows and columns forming the cell map match the arrangment
      ! of fine cells within the coarse cell

      ! These are aligned as follows with the LFRic directions:
      !         N
      !   |--------------|
      !   |    row 1     |
      !   |c            c|
      !   |o            o|
      ! W |l            l| E
      !   |              |
      !   |1           nx|
      !   |    row ny    |
      !   |--------------|
      !          S

      !-------------------------------------------------------------------------
      ! Edges of coarse cell
      !-------------------------------------------------------------------------

      ! The weights here describe the flux through the edge of the fine cell
      ! This should be the fractional area of the coarse cell edge taken up
      ! by the fine cell. If fine cells have equal area then this will be the
      ! reciprocal of the number of fine cells along that edge

      do x_idx = 1, ncell_fine_per_coarse_x
        ! N edge is first row of cell map
        y_idx = 1
        df = N
        fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) =                                        &
                                    real(weights(map_fine(df,cell_map(x_idx,y_idx))+k), r_double) &
                                    * coarse_field(map_coarse(df)+k)

        ! S edge is last row of cell map
        y_idx = ncell_fine_per_coarse_y
        df = S
        fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) =                                        &
                                    real(weights(map_fine(df,cell_map(x_idx,y_idx))+k), r_double) &
                                    * coarse_field(map_coarse(df)+k)
      end do

      do y_idx = 1, ncell_fine_per_coarse_y
        ! W edge is first column of cell map
        x_idx = 1
        df = W
        fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) =                                        &
                                    real(weights(map_fine(df,cell_map(x_idx,y_idx))+k), r_double) &
                                    * coarse_field(map_coarse(df)+k)

        ! E edge is last column of cell map
        x_idx = ncell_fine_per_coarse_x
        df = E
        fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) =                                        &
                                    real(weights(map_fine(df,cell_map(x_idx,y_idx))+k), r_double) &
                                    * coarse_field(map_coarse(df)+k)
      end do

      !-------------------------------------------------------------------------
      ! DoFs internal to coarse cell
      !-------------------------------------------------------------------------

      ! For internal cells, we interpolate. The weights at these DoFs give the
      ! relative weighting of the WEST and NORTH components. The contributions
      ! from the east and south components are given by 1 minus those weights.
      ! All must be combined with the area fraction weights to give the true
      ! values. As the fine cells are internal to a coarse cell, they will all
      ! have the same orientation as one another.

      ! Do West/East DoFs
      ! These loops will only execute if there are internal DoFs
      df = W
      opp_df = E
      do y_idx = 1, ncell_fine_per_coarse_y
        edge_weight = real(weights(map_fine(df,cell_map(1,y_idx))+k), r_double)

        do x_idx = 2, ncell_fine_per_coarse_x
          internal_weight = real(weights(map_fine(df,cell_map(x_idx,y_idx))+k), r_double)

          fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) =                      &
               edge_weight * ( internal_weight * coarse_field(map_coarse(df)+k) + &
               (1.0_r_double - internal_weight) * coarse_field(map_coarse(opp_df)+k) )
        end do
      end do

      ! Do North/South DoFs
      ! These loops will only execute if there are internal DoFs
      df = N
      opp_df = S
      do x_idx = 1, ncell_fine_per_coarse_x
        edge_weight = real(weights(map_fine(df,cell_map(x_idx,1))+k), r_double)

        do y_idx = 2, ncell_fine_per_coarse_y
          internal_weight = real(weights(map_fine(df,cell_map(x_idx,y_idx))+k), r_double)

          fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) =                      &
               edge_weight * ( internal_weight * coarse_field(map_coarse(df)+k) + &
               (1.0_r_double - internal_weight) * coarse_field(map_coarse(opp_df)+k) )

        end do
      end do

    end do

  end subroutine prolong_w2h_kernel_code_r_double


end module sci_prolong_w2h_kernel_mod
