!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the orientations of halo cells in a partition.
!>
!> The kernel computes the orientation of core and halo cells in a partition.
!> The kernel assumes that the number of mpi partitions is a multiple of 6
!> which implies that no core cells of a partition overlap a cubed-sphere
!> panel edge. Thus it is assumed that all core cells have the same
!> orientation which we assign orientation 1. The kernel computes the
!> orientation of cells using the W2 dof values and the orientation can take
!> the value 1, 2, 3 or 4.
!>
module calc_cell_orientation_kernel_mod

  use argument_mod,      only : arg_type, func_type, &
                                GH_FIELD, GH_WRITE,  &
                                GH_BASIS, CELLS
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W3
  use kernel_mod,        only : kernel_type

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: calc_cell_orientation_kernel_type
    private
    type(arg_type) :: meta_args(1) = (/  &
        arg_type(GH_FIELD, GH_WRITE, W3) &
        /)
    type(func_type) :: meta_funcs(1) = (/ &
        func_type(W3, GH_BASIS)           &
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass ::calc_cell_orientation_code
  end type

  public calc_cell_orientation_code

contains

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
!>  @brief  Function which returns the index depending on orientation of cell
!!
!!  @param[in]      branch                    branch of the cross-stencil
!!  @param[in]      orientation               orientation of cell
!!  @param[return]  index_of_interest_map     returns index
!--------------------------------------------------------------------------------
function index_of_interest_map(branch, orientation)

  implicit none

  integer(kind=i_def) :: index_of_interest_map, branch, orientation

  index_of_interest_map = mod(orientation + 2 + branch, 4) + 1

end function index_of_interest_map

!--------------------------------------------------------------------------------
!>  @brief  Function which returns the cells orientation
!!
!!  @param[in]      ii        relates to the W2 dof
!!  @param[in]      branch    relates to the branch of the cross-stencil
!!  @param[return]  orientation_of_cell    orientation of the cell
!--------------------------------------------------------------------------------
function orientation_of_cell(branch,ii)

  implicit none

  integer(kind=i_def) :: branch, ii
  real(kind=r_def)    :: orientation_of_cell

  orientation_of_cell = real(mod(ii + 6 - branch, 4) + 1,r_def)

end function orientation_of_cell

!--------------------------------------------------------------------------------
!>  @brief  Subroutine which calculates the orientation of cells
!!
!> @details The kernel computes the orientation of cells which is dependent on
!>          the numbering of the local W2 dof values.
!>          The kernel iterates over all core cells in a partition and uses a
!>          cross-stencil inorder to calculate the orientation of cells in the halo.
!>          The orientation of cells in the cross-stencil is calculated by
!>          iterating over all the cells in a branch of the stencil and then for
!>          each of the subsequent branches.
!!
!! @param[in] nlayers              Number of layers
!! @param[inout] orientation       W3 field containing orientation of cells
!! @param[in] undf_w3              Number of unique degrees of freedom
!! @param[in] ndf_w3               Number of degrees of freedom per cell
!! @param[in] map_w3               Dofmap for the cell at the base of the column
!! @param[in] ndf_w2               Number of degrees of freedom per cell
!! @param[in] size_of_stencil      Size of cross-stencil
!! @param[in] stencil_map_w2_cross Dofmap for the stencil
!! @param[in] stencil_map_w3_cross Dofmap for the stencil
!--------------------------------------------------------------------------------
subroutine calc_cell_orientation_code(  nlayers,                       &
                                        orientation,                   &
                                        undf_w3,                       &
                                        ndf_w3,                        &
                                        map_w3,                        &
                                        ndf_w2,                        &
                                        size_of_stencil,               &
                                        stencil_map_w2_cross,          &
                                        stencil_map_w3_cross )

  use log_mod, only : log_event, LOG_LEVEL_ERROR

  implicit none

  integer(kind=i_def), intent(in)                     :: nlayers
  integer(kind=i_def), intent(in)                     :: undf_w3
  integer(kind=i_def), intent(in)                     :: ndf_w3
  real(kind=r_def), intent(inout)                     :: orientation(1:undf_w3)
  integer(kind=i_def), dimension(ndf_w3), intent(in)  :: map_w3
  integer(kind=i_def), intent(in)                     :: ndf_w2
  integer(kind=i_def), intent(in)                     :: size_of_stencil
  integer(kind=i_def), intent(in)                     :: stencil_map_w2_cross(1:ndf_w2,1:size_of_stencil)
  integer(kind=i_def), intent(in)                     :: stencil_map_w3_cross(1:ndf_w3,1:size_of_stencil)

  integer(kind=i_def) :: k, ii, branch, stencil ! iteration counters
  integer(kind=i_def) :: cell_inner, cell_outer
  integer(kind=i_def) :: cell_inner_index_of_interest
  integer(kind=i_def) :: shared_dof_of_interest
  integer(kind=i_def) :: halo_depth, cross_stencil_branches

  halo_depth = (size_of_stencil-1)/4
  cross_stencil_branches = 4

  do k=0, nlayers-1

    ! Kernel loops over core cells only, all core cells are assumed to have orientation 1
    orientation(stencil_map_w3_cross(1,1)+k) = 1_r_def
    do branch = 1,cross_stencil_branches
      cell_inner = 1
      do stencil = 1,halo_depth
        cell_outer = 1 + (halo_depth * (branch - 1)) + stencil

        if ( orientation(stencil_map_w3_cross(1,cell_inner)+k) < 0.5_r_def ) then
          call log_event( "Orientation unassigned ", LOG_LEVEL_ERROR )
        end if

        cell_inner_index_of_interest = index_of_interest_map(branch,int(orientation(stencil_map_w3_cross(1,cell_inner)+k)))
        shared_dof_of_interest = int(stencil_map_w2_cross(cell_inner_index_of_interest,cell_inner), i_def)

        do ii=1,ndf_w2
          if ( shared_dof_of_interest == int(stencil_map_w2_cross(ii,cell_outer),i_def) ) then
            orientation(stencil_map_w3_cross(1,cell_outer)+k) = orientation_of_cell(branch,ii)
          end if

        end do
        cell_inner = cell_outer
      end do
    end do
  end do

end subroutine calc_cell_orientation_code

end module calc_cell_orientation_kernel_mod
