
!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-----------------------------------------------------------------------------
!> @details Creates "vertical-edge selector" fields to avoid calculation
!!          duplication on vertical edges when iterating over columns
module sci_vert_edge_selector_halo_kernel_mod

  use argument_mod,                       only: arg_type,                      &
                                                GH_FIELD, GH_WRITE,            &
                                                GH_INTEGER,                    &
                                                ANY_SPACE_1,                   &
                                                HALO_CELL_COLUMN
  use constants_mod,                      only: i_def
  use fs_continuity_mod,                  only: W3
  use kernel_mod,                         only: kernel_type
  use sci_vert_edge_selector_support_mod, only: compute_vert_edge_selector

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  !>
  type, public, extends(kernel_type) :: vert_edge_selector_halo_kernel_type
    private
    ! TODO: replace ANY_SPACE_1 with W1V once PSyclone's LFRic API recognises
    ! W1V as a valid function space for kernel metadata
    type(arg_type) :: meta_args(2) = (/                                        &
        arg_type(GH_FIELD, GH_INTEGER, GH_WRITE, W3),                          &
        arg_type(GH_FIELD, GH_INTEGER, GH_WRITE, ANY_SPACE_1)                  &
    /)
    integer :: operates_on = HALO_CELL_COLUMN
  contains
    procedure, nopass :: vert_edge_selector_halo_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: vert_edge_selector_halo_code

contains

!> @details Creates "vertical-edge selector" fields to avoid calculation
!!          duplication on vertical edges when iterating over columns
!> @param[in]     nlayers             The number of layers
!> @param[in,out] vert_edge_selector  The vertical-edge selector. It is a W3
!!                                    integer field, encoding which vertical
!!                                    edges of a column this cell should
!!                                    process.
!> @param[in,out] vert_edge_counter   An integer W1V field, counting the
!!                                    number of times each vertical edge has
!!                                    been claimed
!> @param[in]     ndf_w3              Num of DoFs for W3 per cell
!> @param[in]     undf_w3             Num of DoFs for this partition for W3
!> @param[in]     map_w3              DoF-map for W3 in base cells
!> @param[in]     ndf_w1v             Num of DoFs for W1v per cell
!> @param[in]     undf_w1v            Num of DoFs for this partition for W1v
!> @param[in]     map_w1v             DoF-map for W1v in base cells
subroutine vert_edge_selector_halo_code( nlayers,                   &
                                         vert_edge_selector,        &
                                         vert_edge_counter,         &
                                         ndf_w3, undf_w3, map_w3,   &
                                         ndf_w1v, undf_w1v, map_w1v &
                                       )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: ndf_w1v, undf_w1v

  integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in) :: map_w1v(ndf_w1v)

  integer(kind=i_def), intent(inout) :: vert_edge_selector(undf_w3)
  integer(kind=i_def), intent(inout) :: vert_edge_counter(undf_w1v)

  call compute_vert_edge_selector(                                             &
      vert_edge_selector,                                                      &
      vert_edge_counter,                                                       &
      ndf_w3, undf_w3, map_w3,                                                 &
      ndf_w1v, undf_w1v, map_w1v                                               &
  )

end subroutine vert_edge_selector_halo_code

end module sci_vert_edge_selector_halo_kernel_mod
