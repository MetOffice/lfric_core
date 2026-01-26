!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Computes weights for the prologation operation on a W2h field.
!> @details Computes weights for performing the prolongation operation to take a
!!          W2h field from a coarse mesh to a fine mesh. The weights are stored
!!          in a W2h field on the fine mesh. For horizontal W2h components, the
!!          fine cell DoFs on the edges of coarse cells have weight given by
!!          the reciprocal of the number of fine cells touching that edge of
!!          the coarse cell. At fine cell DoFs on the interior of the coarse
!!          cell, weights correspond to interpolation between the edges of the
!!          coarse cell. 
!!          This is only designed to work with the lowest order W2h fields and
!!          for fine cells that are exactly nested within a coarse cell.

module sci_weights_prolong_w2h_kernel_mod

use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_WRITE,         &
                                   GH_COARSE, GH_FINE,        &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_2, &
                                   CELL_COLUMN
use constants_mod,           only: i_def, r_def
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

type, public, extends(kernel_type) :: weights_prolong_w2h_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                           &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1, mesh_arg=GH_FINE   ), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2, mesh_arg=GH_COARSE )  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: weights_prolong_w2h_kernel_code
end type weights_prolong_w2h_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

public :: weights_prolong_w2h_kernel_code

contains

  !> @brief Computes weights for the prologation operation on a W2h field.
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] weights_fine             The fine grid W2h field that contains
  !!                                         the weights to be computed here for
  !!                                         prolongation
  !> @param[in]     dummy_coarse             An unused coarse W2h field supplied
  !!                                         to ensure that appropriate coarse
  !!                                         field properties are provided.
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  subroutine weights_prolong_w2h_kernel_code(nlayers,                 &
                                             cell_map,                &
                                             ncell_fine_per_coarse_x, &
                                             ncell_fine_per_coarse_y, &
                                             ncell_fine,              &
                                             weights_fine,            &
                                             dummy_coarse,            &
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
    real(kind=r_def),    intent(inout) :: weights_fine(undf_fine)
    real(kind=r_def),    intent(in)    :: dummy_coarse(undf_coarse)

    ! Maps
    integer(kind=i_def), intent(in) :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in) :: map_coarse(ndf)
    integer(kind=i_def), intent(in) :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)

    ! Internal arguments
    integer(kind=i_def) :: df, k, x_idx, y_idx

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
        weights_fine(map_fine(df,cell_map(x_idx,y_idx))+k) =                      &
                                        1.0_r_def / real(ncell_fine_per_coarse_x, r_def)

        ! S edge is last row of cell map
        y_idx = ncell_fine_per_coarse_y
        df = S
        weights_fine(map_fine(df,cell_map(x_idx,y_idx))+k) =                      &
                                        1.0_r_def / real(ncell_fine_per_coarse_x, r_def)
      end do

      do y_idx = 1, ncell_fine_per_coarse_y
        ! W edge is first column of cell map
        x_idx = 1
        df = W
        weights_fine(map_fine(df,cell_map(x_idx,y_idx))+k) =                      &
                                        1.0_r_def / real(ncell_fine_per_coarse_y, r_def)

        ! E edge is last column of cell map
        x_idx = ncell_fine_per_coarse_x
        df = E
        weights_fine(map_fine(df,cell_map(x_idx,y_idx))+k) =                      &
                                        1.0_r_def / real(ncell_fine_per_coarse_y, r_def)
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
      do y_idx = 1, ncell_fine_per_coarse_y
        do x_idx = 2, ncell_fine_per_coarse_x

          weights_fine(map_fine(df,cell_map(x_idx,y_idx))+k) =                    &
              real(ncell_fine_per_coarse_x + 1 - x_idx, r_def) / real(ncell_fine_per_coarse_x, r_def)
        end do
      end do

      ! Do North/South DoFs
      ! These loops will only execute if there are internal DoFs
      df = N
      do x_idx = 1, ncell_fine_per_coarse_x
        do y_idx = 2, ncell_fine_per_coarse_y

          weights_fine(map_fine(df,cell_map(x_idx,y_idx))+k) =                    &
              real(ncell_fine_per_coarse_y + 1 - y_idx, r_def) / real(ncell_fine_per_coarse_y, r_def)

        end do
      end do

    end do

  end subroutine weights_prolong_w2h_kernel_code

end module sci_weights_prolong_w2h_kernel_mod
