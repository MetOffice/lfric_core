
!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Support routines for computing and using the face selector fields
module sci_face_selector_support_mod

use constants_mod,         only: i_def
use reference_element_mod, only: W, S, E, N

implicit none

private

public :: compute_face_selector, face_from_face_selector

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
subroutine compute_face_selector(                                              &
    face_selector_ew,                                                          &
    face_selector_ns,                                                          &
    face_counter,                                                              &
    ndf_w3, undf_w3, map_w3,                                                   &
    ndf_w2h, undf_w2h, map_w2h                                                 &
)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: ndf_w2h, undf_w2h

  integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)

  integer(kind=i_def), intent(inout) :: face_selector_ew(undf_w3)
  integer(kind=i_def), intent(inout) :: face_selector_ns(undf_w3)
  integer(kind=i_def), intent(inout) :: face_counter(undf_w2h)

  ! Internal variables
  integer(kind=i_def) :: idx, w_val, s_val, e_val, n_val

  idx = map_w3(1)
  w_val = face_counter(map_w2h(W))
  s_val = face_counter(map_w2h(S))
  e_val = face_counter(map_w2h(E))
  n_val = face_counter(map_w2h(N))

  if (face_selector_ew(idx) == 0) then
    if (w_val == 0 .and. e_val == 0) then
      face_selector_ew(idx) = 2
      face_counter(map_w2h(W)) = 1
      face_counter(map_w2h(E)) = 1
    else if (w_val == 0 .and. e_val == 1) then
      face_selector_ew(idx) = 1
      face_counter(map_w2h(W)) = 1
    else if (w_val == 1 .and. e_val == 0) then
      face_selector_ew(idx) = -1
      face_counter(map_w2h(E)) = 1
    end if
  end if

  if (face_selector_ns(idx) == 0) then
    if (s_val == 0 .and. n_val == 0) then
      face_selector_ns(idx) = 2
      face_counter(map_w2h(S)) = 1
      face_counter(map_w2h(N)) = 1
    else if (s_val == 0 .and. n_val == 1) then
      face_selector_ns(idx) = 1
      face_counter(map_w2h(S)) = 1
    else if (s_val == 1 .and. n_val == 0) then
      face_selector_ns(idx) = -1
      face_counter(map_w2h(N)) = 1
    end if
  end if

end subroutine compute_face_selector

!> @details Determines the face to perform calculations on from the face
!!          selector fields, and the index of the loop.
!> @param[in] idx               The index of the loop over faces
!> @param[in] face_selector_ew  The value of the East-West face selector
!> @param[in] face_selector_ns  The value of the North-South face selector
function face_from_face_selector(                                              &
    idx,                                                                       &
    face_selector_ew,                                                          &
    face_selector_ns                                                           &
) result(face)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: idx
  integer(kind=i_def), intent(in) :: face_selector_ew
  integer(kind=i_def), intent(in) :: face_selector_ns

  ! Internal variables
  integer(kind=i_def) :: ew_or_ns, face

  ! Switch to determine whether to iterate over East-West or North-South faces
  ! 1: for E/W faces, 0: for N/S faces
  ew_or_ns = (1 + SIGN(1, ABS(face_selector_ew) - idx)) / 2

  ! Faces are looped over from 1 to ABS(face_selector_ew)+ABS(face_selector_ns)
  ! This equation forms a unique mapping from the loop index and face selector
  ! values to the faces (W,S,E,N) respectively, where the face selectors take
  ! the values:
  ! 1: W face / S face
  ! -1: E face / N face
  ! 0: no faces
  ! 2: W and E faces / S and N faces
  ! Faces are always iterated over in the order W,E,S,N
  face = (                                                                     &
      -1                                                                       &
      + ew_or_ns * 2*(idx - MIN(0, face_selector_ew))                          &
      + (1 - ew_or_ns) * (                                                     &
          1 + 2*(idx - ABS(face_selector_ew) - MIN(0, face_selector_ns))       &
      )                                                                        &
  )

end function face_from_face_selector

end module sci_face_selector_support_mod
