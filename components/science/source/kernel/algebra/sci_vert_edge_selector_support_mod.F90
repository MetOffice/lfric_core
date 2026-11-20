
!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-----------------------------------------------------------------------------
!> @brief Support routines for computing and using the vertical-edge selector
!!        field
module sci_vert_edge_selector_support_mod

use constants_mod,         only: i_def
use reference_element_mod, only: SW, SE, NE, NW

implicit none

private

public :: compute_vert_edge_selector, num_vert_edges_to_do, &
          decode_vert_edge_selector

! Bit values used to encode, in the vert_edge_selector field, which vertical
! edges of a column this cell should process. The order matches the local
! dof ordering of the W1V function space at order (0,0), i.e. SW, SE, NE, NW.
integer(kind=i_def), parameter :: sw_bit = 1
integer(kind=i_def), parameter :: se_bit = 2
integer(kind=i_def), parameter :: ne_bit = 4
integer(kind=i_def), parameter :: nw_bit = 8

contains

!> @details Creates a "vertical-edge selector" field to avoid calculation
!!          duplication on vertical edges when iterating over columns
!> @param[in,out] vert_edge_selector  The vertical-edge selector. It is a W3
!!                                    integer field. Bits sw_bit, se_bit,
!!                                    ne_bit and nw_bit are set if this cell
!!                                    should process the corresponding
!!                                    vertical edge of the column.
!> @param[in,out] vert_edge_counter   An integer W1V field, counting the
!!                                    number of times each vertical edge has
!!                                    been claimed
!> @param[in]     ndf_w3              Num of DoFs for W3 per cell
!> @param[in]     undf_w3             Num of DoFs for this partition for W3
!> @param[in]     map_w3              DoF-map for W3 in base cells
!> @param[in]     ndf_w1v             Num of DoFs for W1v per cell
!> @param[in]     undf_w1v            Num of DoFs for this partition for W1v
!> @param[in]     map_w1v             DoF-map for W1v in base cells
subroutine compute_vert_edge_selector(                                         &
    vert_edge_selector,                                                        &
    vert_edge_counter,                                                         &
    ndf_w3, undf_w3, map_w3,                                                   &
    ndf_w1v, undf_w1v, map_w1v                                                 &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: ndf_w1v, undf_w1v

  integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in) :: map_w1v(ndf_w1v)

  integer(kind=i_def), intent(inout) :: vert_edge_selector(undf_w3)
  integer(kind=i_def), intent(inout) :: vert_edge_counter(undf_w1v)

  ! Internal variables
  integer(kind=i_def) :: idx, s

  idx = map_w3(1)

  if (vert_edge_selector(idx) == 0) then

    s = 0

    if (vert_edge_counter(map_w1v(1)) == 0) then
      s = s + sw_bit
      vert_edge_counter(map_w1v(1)) = 1
    end if

    if (vert_edge_counter(map_w1v(2)) == 0) then
      s = s + se_bit
      vert_edge_counter(map_w1v(2)) = 1
    end if

    if (vert_edge_counter(map_w1v(3)) == 0) then
      s = s + ne_bit
      vert_edge_counter(map_w1v(3)) = 1
    end if

    if (vert_edge_counter(map_w1v(4)) == 0) then
      s = s + nw_bit
      vert_edge_counter(map_w1v(4)) = 1
    end if

    vert_edge_selector(idx) = s

  end if

end subroutine compute_vert_edge_selector

!> @details Returns the number of vertical edges that a kernel should
!!          iterate over for a given value of the vertical-edge selector.
!!          Intended to be called before a loop over vertical edges.
!> @param[in] vert_edge_selector  The value of the vertical-edge selector in
!!                                this column
function num_vert_edges_to_do(vert_edge_selector) result(n)

  implicit none

  integer(kind=i_def), intent(in) :: vert_edge_selector
  integer(kind=i_def)             :: n

  n = 0
  if (iand(vert_edge_selector, sw_bit) /= 0) n = n + 1
  if (iand(vert_edge_selector, se_bit) /= 0) n = n + 1
  if (iand(vert_edge_selector, ne_bit) /= 0) n = n + 1
  if (iand(vert_edge_selector, nw_bit) /= 0) n = n + 1

end function num_vert_edges_to_do

!> @details Determines which vertical edge to perform calculations on, from
!!          the vertical-edge selector field and the index of the loop over
!!          vertical edges
!> @param[in] j                   The index of the loop over vertical edges,
!!                                between 1 and num_vert_edges_to_do
!> @param[in] vert_edge_selector  The value of the vertical-edge selector in
!!                                this column
function decode_vert_edge_selector(j, vert_edge_selector) result(edge)

  implicit none

  integer(kind=i_def), intent(in) :: j
  integer(kind=i_def), intent(in) :: vert_edge_selector
  integer(kind=i_def)             :: edge

  ! Internal variables
  integer(kind=i_def) :: count

  count = 0

  if (iand(vert_edge_selector, sw_bit) /= 0) then
    count = count + 1
    if (count == j) then
      edge = SW
      return
    end if
  end if

  if (iand(vert_edge_selector, se_bit) /= 0) then
    count = count + 1
    if (count == j) then
      edge = SE
      return
    end if
  end if

  if (iand(vert_edge_selector, ne_bit) /= 0) then
    count = count + 1
    if (count == j) then
      edge = NE
      return
    end if
  end if

  if (iand(vert_edge_selector, nw_bit) /= 0) then
    count = count + 1
    if (count == j) then
      edge = NW
      return
    end if
  end if

end function decode_vert_edge_selector

end module sci_vert_edge_selector_support_mod
