!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Applies a correction to the process of averaging from W3 to W2
!> @details The averaging operation of a scalar W3 field to W2 points is not
!!          accurate at the edges of cubed sphere panels. This kernel applies a
!!          correction to the naive W2 field, based on a pre-computed
!!          displacement field. The correction depends on the gradient of the
!!          original field parallel to the panel boundary.
!!          This kernel only works for the lowest-order elements.
module sci_w3_to_w2_correction_kernel_mod

  use argument_mod,      only : arg_type,                    &
                                GH_FIELD, GH_REAL,           &
                                GH_READ, GH_INC,             &
                                STENCIL, CROSS, CELL_COLUMN, &
                                ANY_DISCONTINUOUS_SPACE_9,   &
                                ANY_SPACE_2
  use constants_mod,     only : r_def, i_def, l_def
  use fs_continuity_mod, only : W2, W3
  use kernel_mod,        only : kernel_type

  implicit none
  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: w3_to_w2_correction_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                        &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                           &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W3, STENCIL(CROSS)),           &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_SPACE_2),                  &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_9,     &
                                                            STENCIL(CROSS))    &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: w3_to_w2_correction_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: w3_to_w2_correction_code

contains

!> @brief Applies a correction to the process of averaging from W3 to W2
!> @param[in]     nlayers          Number of layers in the mesh
!> @param[in,out] field_w2         The field already naively averaged to W2, to
!!                                 be corrected
!> @param[in]     field_w3         The original scalar field in W3
!> @param[in]     stencil_size_w3  The size of the stencil for the W3 field
!> @param[in]     stencil_map_w3   The stencil DoFmap for W3
!> @param[in]     displacement     2D W2H field containing the displacements
!!                                 corresponding to the averaging error. This is
!!                                 dimensionless, divided by the cell width
!> @param[in]     panel_id         ID for panels of the underlying mesh
!> @param[in]     stencil_size_pid The size of the stencil for the panel ID
!> @param[in]     stencil_map_pid  The stencil DoFmap for the panel ID
!> @param[in]     ndf_w2           Number of DoFs for W2 per cell
!> @param[in]     undf_w2          Number of unique DoFs for W2 per partition
!> @param[in]     map_w2           The DoF map for bottom layer cells for W2
!> @param[in]     ndf_w3           Number of DoFs for W3 per cell
!> @param[in]     undf_w3          Number of unique DoFs for W3 per partition
!> @param[in]     map_w3           The DoF map for bottom layer cells for W3
!> @param[in]     ndf_w2h_2d       Number of DoFs for W2H per cell
!> @param[in]     undf_w2h_2d      Num of unique DoFs for 2D W2H per partition
!> @param[in]     map_w2h_2d       The DoF map for bottom layer cells for 2D W2H
!> @param[in]     ndf_pid          Number of DoFs for panel id per cell
!> @param[in]     undf_pid         Num of unique DoFs for panel id per partition
!> @param[in]     map_pid          DoF map for bottom layer cells for panel ID
subroutine w3_to_w2_correction_code(                                        &
                                     nlayers,                               &
                                     field_w2,                              &
                                     field_w3,                              &
                                     stencil_size_w3, stencil_map_w3,       &
                                     displacement,                          &
                                     panel_id,                              &
                                     stencil_size_pid, stencil_map_pid,     &
                                     ndf_w2, undf_w2, map_w2,               &
                                     ndf_w3, undf_w3, map_w3,               &
                                     ndf_w2h_2d, undf_w2h_2d, map_w2h_2d,   &
                                     ndf_pid, undf_pid, map_pid             &
                                   )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in)    :: ndf_w2h_2d, undf_w2h_2d
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), intent(in)    :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in)    :: stencil_size_w3, stencil_size_pid
  integer(kind=i_def), intent(in)    :: map_w2h_2d(ndf_w2h_2d)
  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: stencil_map_w3(ndf_w3, stencil_size_w3)
  integer(kind=i_def), intent(in)    :: stencil_map_pid(ndf_pid, stencil_size_pid)
  real(kind=r_def),    intent(inout) :: field_w2(undf_w2)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  real(kind=r_def),    intent(in)    :: field_w3(undf_w3)
  real(kind=r_def),    intent(in)    :: displacement(undf_w2h_2d)

  ! Internal variables
  integer(kind=i_def) :: i, k, face
  integer(kind=i_def) :: cell_panel, stencil_panel, perp_panel
  integer(kind=i_def) :: perp_cell, stencil_cell, same_panel_idx, sign_recon
  logical(kind=l_def) :: panel_corner
  real(kind=r_def)    :: gradient

  ! Indices of cell through each face in the cross stencil
  integer(kind=i_def) :: stencil_cells_by_face(4)

  ! Indices of cells in cross stencil that are perpendicular to the cell of interest
  integer(kind=i_def) :: perp_cells(2,4)

  ! Assumed direction for derivatives in this kernel is:
  !  y
  !  ^
  !  |-> x
  !

  ! Layout of dofs for the stencil map
  ! dimensions of map are (ndf, ncell)
  ! Horizontally:
  !
  !   -- 4 --
  !   |     |
  !   1     3
  !   |     |
  !   -- 2 --
  !
  ! df = 5 is in the centre on the bottom face
  ! df = 6 is in the centre on the top face

  ! The layout of the cells in the stencil is:
  !
  !          -----
  !          |   |
  !          | 5 |
  !     ---------------
  !     |    |   |    |
  !     |  2 | 1 |  4 |
  !     ---------------
  !          |   |
  !          | 3 |
  !          -----

  stencil_cells_by_face(:) = (/ 2, 3, 4, 5 /)
  perp_cells(:,:) = reshape([5, 3, 4, 2, 5, 3, 4, 2], [2,4])

  ! If the full stencil isn't available, we must be at the domain edge.
  ! Don't need to increment the field here so just exit the kernel.
  if (stencil_size_w3 < 5_i_def) then
    return
  end if

  cell_panel = int(panel_id(map_pid(1)), i_def)

  do face = 1, 4

    stencil_cell = stencil_cells_by_face(face)
    stencil_panel = int(panel_id(stencil_map_pid(1, stencil_cell)), i_def)

    if (stencil_panel /= cell_panel) then

      ! Check panels of other stencil cells -- are we at the corner of a panel?
      panel_corner = .false.
      do i = 1, 2
        perp_cell = perp_cells(i,face)
        perp_panel = int(panel_id(stencil_map_pid(1, perp_cell)), i_def)
        if (perp_panel /= cell_panel) then
          panel_corner = .true.
          same_panel_idx = MOD(i, 2) + 1  ! Gives 1 if i == 2, or 2 if i == 1
          sign_recon = (-1)**(same_panel_idx-1)  ! 1 if same_panel_idx == 1 or -1 if same_panel_idx == 2
        end if
      end do

      if (panel_corner) then
        ! If panel corner then reconstruction is based on nearest cell
        do k = 0, nlayers - 1
          gradient = sign_recon *                                              &
            (field_w3(stencil_map_w3(1,perp_cells(same_panel_idx,face)) + k)   &
            - field_w3(stencil_map_w3(1,1) + k))

            ! Add contribution to the averaged W2 field. Factor of 0.5 as there
            ! is a contribution from each face
            field_w2(map_w2(face)+k) = field_w2(map_w2(face)+k) -              &
              0.5_r_def * gradient * displacement(map_w2h_2d(face))
        end do

      else
        ! If not panel corner, then fit quadratic through neighbouring cells
        ! Assume cells are of equal width, so gradient is linear coefficient
        do k = 0, nlayers - 1
          gradient = (field_w3(stencil_map_w3(1,perp_cells(1,face)) + k)       &
            - field_w3(stencil_map_w3(1,perp_cells(2,face)) + k)) * 0.5_r_def

          ! Add contribution to the averaged W2 field. Factor of 0.5 as there
          ! is a contribution from each face
          field_w2(map_w2(face)+k) = field_w2(map_w2(face)+k) -                &
            0.5_r_def * gradient * displacement(map_w2h_2d(face))
        end do
      end if
    end if
  end do

end subroutine w3_to_w2_correction_code

end module sci_w3_to_w2_correction_kernel_mod
