!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-----------------------------------------------------------------------------
!> @brief Split a W1 field into the component W1h and W1v fields
!> @details This kernel is only valid for the lowest-order W1 space
!>          (ndf_w1 == 12). It relies on a dof-ordering assumption that
!>          does not hold for higher-order W1 spaces.
module sci_split_w1_field_kernel_mod

  use argument_mod,                       only : arg_type,                     &
                                                 GH_FIELD, GH_REAL,            &
                                                 GH_WRITE, GH_INTEGER,         &
                                                 GH_READ,                      &
                                                 ANY_DISCONTINUOUS_SPACE_1,    &
                                                 ANY_DISCONTINUOUS_SPACE_2,    &
                                                 ANY_DISCONTINUOUS_SPACE_3,    &
                                                 CELL_COLUMN
  use constants_mod,                      only : r_double, r_single, i_def
  use fs_continuity_mod,                  only : W1
  use kernel_mod,                         only : kernel_type
  use reference_element_mod,              only : SW
  use sci_face_selector_support_mod,      only : face_from_face_selector
  use sci_vert_edge_selector_support_mod, only : num_vert_edges_to_do,         &
                                                 decode_vert_edge_selector

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: split_w1_field_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                        &
         ! NB: This is to be used to write to a continuous W1h field, but
         ! using a discontinuous data pattern, so use discontinuous metadata.
         ! TODO: replace ANY_DISCONTINUOUS_SPACE_2 with W1h once PSyclone's
         ! LFRic API recognises W1h as a valid function space for kernel
         ! metadata
        arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2),   &
         ! NB: as above, but for a continuous W1v field.
         ! TODO: replace ANY_DISCONTINUOUS_SPACE_1 with W1v once PSyclone's
         ! LFRic API recognises W1v as a valid function space for kernel
         ! metadata
        arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   &
        arg_type(GH_FIELD, GH_REAL,    GH_READ,  W1),                          &
        arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),   &
        arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),   &
        arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)    &
    /)
    integer :: operates_on = CELL_COLUMN
  end type
  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: split_w1_field_code

  ! Generic interface for real32 and real64 types
  interface split_w1_field_code
    module procedure  &
      split_w1_field_code_r_single, &
      split_w1_field_code_r_double
  end interface

contains

!> @brief Split a W1 field into the component W1h and W1v fields
!> @param[in,out] uv Horizontal-edge component
!> @param[in,out] w Vertical-edge component
!> @param[in]     uvw Full W1 field
!> @param[in]     face_selector_ew  2D field indicating which W/E horizontal
!!                                  edges to loop over in this column
!> @param[in]     face_selector_ns  2D field indicating which N/S horizontal
!!                                  edges to loop over in this column
!> @param[in]     vert_edge_selector 2D field indicating which vertical edges
!!                                   to loop over in this column
!> @param[in]     ndf_w1h    Number of degrees of freedom per cell for W1h
!> @param[in]     undf_w1h   Number of unique degrees of freedom for W1h
!> @param[in]     map_w1h    Dofmap for the base cell of the column for W1h
!> @param[in]     ndf_w1v    Number of degrees of freedom per cell for W1v
!> @param[in]     undf_w1v   Number of unique degrees of freedom for W1v
!> @param[in]     map_w1v    Dofmap for the base cell of the column for W1v
!> @param[in]     ndf_w1     Number of degrees of freedom per cell for W1
!> @param[in]     undf_w1    Number of unique degrees of freedom for W1
!> @param[in]     map_w1     Dofmap for the base cell of the column for W1
!> @param[in]     ndf_w3_2d  Num of DoFs for 2D W3 per cell
!> @param[in]     undf_w3_2d Num of DoFs for this partition for 2D W3
!> @param[in]     map_w3_2d  Map for 2D W3

! R_DOUBLE PRECISION
! ==================
subroutine split_w1_field_code_r_double(nlayers,                               &
                                        uv, w, uvw,                            &
                                        face_selector_ew,                      &
                                        face_selector_ns,                      &
                                        vert_edge_selector,                    &
                                        ndf_w1h, undf_w1h, map_w1h,            &
                                        ndf_w1v, undf_w1v, map_w1v,            &
                                        ndf_w1, undf_w1, map_w1,               &
                                        ndf_w3_2d, undf_w3_2d, map_w3_2d)
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w1, ndf_w1h, ndf_w1v
  integer(kind=i_def), intent(in) :: undf_w1, undf_w1h, undf_w1v
  integer(kind=i_def), intent(in) :: ndf_w3_2d, undf_w3_2d
  integer(kind=i_def), intent(in) :: map_w1(ndf_w1)
  integer(kind=i_def), intent(in) :: map_w1h(ndf_w1h)
  integer(kind=i_def), intent(in) :: map_w1v(ndf_w1v)
  integer(kind=i_def), intent(in) :: map_w3_2d(ndf_w3_2d)
  real(kind=r_double), intent(inout) :: uv(undf_w1h)
  real(kind=r_double), intent(inout) :: w(undf_w1v)
  real(kind=r_double), intent(in) :: uvw(undf_w1)
  integer(kind=i_def), intent(in) :: face_selector_ew(undf_w3_2d)
  integer(kind=i_def), intent(in) :: face_selector_ns(undf_w3_2d)
  integer(kind=i_def), intent(in) :: vert_edge_selector(undf_w3_2d)

  ! Internal variables
  integer(kind=i_def) :: df, k, j, edge, base_df
  integer(kind=i_def) :: hori_dofs_to_do, vert_dofs_to_do

  ! Loop over horizontal W1 DoFs (WB/SB/EB/NB and WT/ST/ET/NT). These are
  ! shared horizontally between (up to) two neighbouring columns, so use
  ! the face selector to avoid duplicating work. Only the "bottom" local
  ! dof of each direction is needed, as the vertical stacking from k=0 to
  ! nlayers reaches the "top" edge of the column automatically
  hori_dofs_to_do = ABS(face_selector_ew(map_w3_2d(1))) +                      &
                    ABS(face_selector_ns(map_w3_2d(1)))

  do j = 1, hori_dofs_to_do
    df = face_from_face_selector(j, face_selector_ew(map_w3_2d(1)),            &
                                    face_selector_ns(map_w3_2d(1)))

    do k = 0, nlayers
      uv(map_w1h(df) + k) = uvw(map_w1(df) + k)
    end do
  end do

  ! Loop over vertical W1 DoFs (SW/SE/NE/NW). These are shared horizontally
  ! between (up to) four neighbouring columns, so use the vertical-edge
  ! selector to avoid duplicating work. Each vertical edge only spans a
  ! single layer, so is looped over once per layer
  vert_dofs_to_do = num_vert_edges_to_do(vert_edge_selector(map_w3_2d(1)))

  do j = 1, vert_dofs_to_do
    edge = decode_vert_edge_selector(j, vert_edge_selector(map_w3_2d(1)))
    df = edge - SW + 1
    base_df = edge

    do k = 0, nlayers-1
      w(map_w1v(df) + k) = uvw(map_w1(base_df) + k)
    end do
  end do

end subroutine split_w1_field_code_r_double

! R_SINGLE PRECISION
! ==================
subroutine split_w1_field_code_r_single(nlayers,                               &
                                        uv, w, uvw,                            &
                                        face_selector_ew,                      &
                                        face_selector_ns,                      &
                                        vert_edge_selector,                    &
                                        ndf_w1h, undf_w1h, map_w1h,            &
                                        ndf_w1v, undf_w1v, map_w1v,            &
                                        ndf_w1, undf_w1, map_w1,               &
                                        ndf_w3_2d, undf_w3_2d, map_w3_2d)
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w1, ndf_w1h, ndf_w1v
  integer(kind=i_def), intent(in) :: undf_w1, undf_w1h, undf_w1v
  integer(kind=i_def), intent(in) :: ndf_w3_2d, undf_w3_2d
  integer(kind=i_def), intent(in) :: map_w1(ndf_w1)
  integer(kind=i_def), intent(in) :: map_w1h(ndf_w1h)
  integer(kind=i_def), intent(in) :: map_w1v(ndf_w1v)
  integer(kind=i_def), intent(in) :: map_w3_2d(ndf_w3_2d)
  real(kind=r_single), intent(inout) :: uv(undf_w1h)
  real(kind=r_single), intent(inout) :: w(undf_w1v)
  real(kind=r_single), intent(in) :: uvw(undf_w1)
  integer(kind=i_def), intent(in) :: face_selector_ew(undf_w3_2d)
  integer(kind=i_def), intent(in) :: face_selector_ns(undf_w3_2d)
  integer(kind=i_def), intent(in) :: vert_edge_selector(undf_w3_2d)

  ! Internal variables
  integer(kind=i_def) :: df, k, j, edge, base_df
  integer(kind=i_def) :: hori_dofs_to_do, vert_dofs_to_do

  ! Loop over horizontal W1 DoFs (WB/SB/EB/NB and WT/ST/ET/NT). These are
  ! shared horizontally between (up to) two neighbouring columns, so use
  ! the face selector to avoid duplicating work. Only the "bottom" local
  ! dof of each direction is needed, as the vertical stacking from k=0 to
  ! nlayers reaches the "top" edge of the column automatically
  hori_dofs_to_do = ABS(face_selector_ew(map_w3_2d(1))) +                      &
                    ABS(face_selector_ns(map_w3_2d(1)))

  do j = 1, hori_dofs_to_do
    df = face_from_face_selector(j, face_selector_ew(map_w3_2d(1)),            &
                                    face_selector_ns(map_w3_2d(1)))

    do k = 0, nlayers
      uv(map_w1h(df) + k) = uvw(map_w1(df) + k)
    end do
  end do

  ! Loop over vertical W1 DoFs (SW/SE/NE/NW). These are shared horizontally
  ! between (up to) four neighbouring columns, so use the vertical-edge
  ! selector to avoid duplicating work. Each vertical edge only spans a
  ! single layer, so is looped over once per layer
  vert_dofs_to_do = num_vert_edges_to_do(vert_edge_selector(map_w3_2d(1)))

  do j = 1, vert_dofs_to_do
    edge = decode_vert_edge_selector(j, vert_edge_selector(map_w3_2d(1)))
    df = edge - SW + 1
    base_df = edge

    do k = 0, nlayers-1
      w(map_w1v(df) + k) = uvw(map_w1(base_df) + k)
    end do
  end do

end subroutine split_w1_field_code_r_single

end module sci_split_w1_field_kernel_mod
