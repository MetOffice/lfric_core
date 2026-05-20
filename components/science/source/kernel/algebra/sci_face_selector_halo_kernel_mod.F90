
!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @details Creates "face selector" fields to avoid calculation duplication on
!!          faces when iterating over columns
module sci_face_selector_halo_kernel_mod

  use argument_mod,                  only: arg_type,                           &
                                           GH_FIELD, GH_WRITE,                 &
                                           GH_REAL, GH_INTEGER,                &
                                           HALO_CELL_COLUMN
  use constants_mod,                 only: i_def
  use fs_continuity_mod,             only: W3, W2H
  use kernel_mod,                    only: kernel_type
  use sci_face_selector_support_mod, only: compute_face_selector

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  !>
  type, public, extends(kernel_type) :: face_selector_halo_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
        arg_type(GH_FIELD, GH_INTEGER, GH_WRITE, W3),                          &
        arg_type(GH_FIELD, GH_INTEGER, GH_WRITE, W3),                          &
        arg_type(GH_FIELD, GH_INTEGER, GH_WRITE, W2H)                          &
    /)
    integer :: operates_on = HALO_CELL_COLUMN
  contains
    procedure, nopass :: face_selector_halo_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: face_selector_halo_code

contains

!> @details Creates "face selector" fields to avoid calculation duplication on
!!          faces when iterating over columns
!> @param[in]     nlayers           The number of layers
!> @param[in,out] face_selector_ew  The East-West face selector. It is a W3
!!                                  integer field, which contains 2 when
!!                                  iterating over the East and West faces of a
!!                                  cell, or 1 when just the West face.
!> @param[in,out] face_selector_ns  The North-South face selector. It is a W3
!!                                  integer field, which contains 2 when
!!                                  iterating over the North and South faces of
!!                                  a cell, or 1 when just the South face.
!> @param[in,out] face_counter      An integer W2H field, counting the number of
!!                                  times that each face has been iterated over
!> @param[in]     ndf_w3            Num of DoFs for W3 per cell
!> @param[in]     undf_w3           Num of DoFs for this partition for W3
!> @param[in]     map_w3            DoF-map for W3 in base cells
!> @param[in]     ndf_w2h           Num of DoFs for W2h per cell
!> @param[in]     undf_w2h          Num of DoFs for this partition for W2h
!> @param[in]     map_w2h           DoF-map for W2h in base cells
subroutine face_selector_halo_code( nlayers,                   &
                                    face_selector_ew,          &
                                    face_selector_ns,          &
                                    face_counter,              &
                                    ndf_w3, undf_w3, map_w3,   &
                                    ndf_w2h, undf_w2h, map_w2h &
                                  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: ndf_w2h, undf_w2h

  integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)

  integer(kind=i_def), intent(inout) :: face_selector_ew(undf_w3)
  integer(kind=i_def), intent(inout) :: face_selector_ns(undf_w3)
  integer(kind=i_def), intent(inout) :: face_counter(undf_w2h)

  call compute_face_selector(                                                  &
      face_selector_ew,                                                        &
      face_selector_ns,                                                        &
      face_counter,                                                            &
      ndf_w3, undf_w3, map_w3,                                                 &
      ndf_w2h, undf_w2h, map_w2h                                               &
  )

end subroutine face_selector_halo_code

end module sci_face_selector_halo_kernel_mod
