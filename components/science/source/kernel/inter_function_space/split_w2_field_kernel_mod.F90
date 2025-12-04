!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Split a W2 field into the component W2v and W2h fields
module split_w2_field_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    GH_WRITE, GH_INTEGER,      &
                                    GH_READ, ANY_SPACE_1,      &
                                    ANY_DISCONTINUOUS_SPACE_2, &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    CELL_COLUMN
  use constants_mod,         only : r_double, r_single, i_def, l_def
  use fs_continuity_mod,     only : W2, W2h, W2v
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : N

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: split_w2_field_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                       &
         ! NB: This is to be used to write to a continuous W2 field, but using
         ! a discontinuous data pattern, so use discontinuous metadata
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, W2v),                       &
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  W2),                        &
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)  &
         /)
    integer :: operates_on = CELL_COLUMN
  end type
  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: split_w2_field_code

  ! Generic interface for real32 and real64 types
  interface split_w2_field_code
    module procedure  &
      split_w2_field_code_r_single, &
      split_w2_field_code_r_double
  end interface

contains

!> @brief Split a W2 field into the component W2v and W2h fields
!> @param[in,out] uv Horizontal wind
!> @param[in,out] w Vertical wind
!> @param[in]     uvw 3D wind
!> @param[in]     face_selector_ew  2D field indicating which W/E faces
!!                                  to loop over in this column
!> @param[in]     face_selector_ns  2D field indicating which N/S faces
!!                                  to loop over in this column
!> @param[in]     ndf_w2h    Number of degrees of freedom per cell for W2h
!> @param[in]     undf_w2h   Number of unique degrees of freedom for W2h
!> @param[in]     map_w2h    Dofmap for the cell at the base of the column for W2h
!> @param[in]     ndf_w2v    Number of degrees of freedom per cell for W2v
!> @param[in]     undf_w2v   Number of unique degrees of freedom for W2v
!> @param[in]     map_w2v    Dofmap for the cell at the base of the column for W2v
!> @param[in]     ndf_w2     Number of degrees of freedom per cell for W2
!> @param[in]     undf_w2    Number of unique degrees of freedom for W2
!> @param[in]     map_w2     Dofmap for the cell at the base of the column for W2
!> @param[in]     ndf_w3_2d  Num of DoFs for 2D W3 per cell
!> @param[in]     undf_w3_2d Num of DoFs for this partition for 2D W3
!> @param[in]     map_w3_2d  Map for 2D W3

! R_DOUBLE PRECISION
! ==================
subroutine split_w2_field_code_r_double(nlayers,                         &
                                        uv, w, uvw,                      &
                                        face_selector_ew,                &
                                        face_selector_ns,                &
                                        ndf_w2h,   undf_w2h,   map_w2h,  &
                                        ndf_w2v,   undf_w2v,   map_w2v,  &
                                        ndf_w2,    undf_w2,    map_w2,   &
                                        ndf_w3_2d, undf_w3_2d, map_w3_2d )
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, ndf_w2h, ndf_w2v
  integer(kind=i_def), intent(in) :: undf_w2, undf_w2h, undf_w2v
  integer(kind=i_def), intent(in) :: ndf_w3_2d, undf_w3_2d
  integer(kind=i_def), dimension(ndf_w2),     intent(in)    :: map_w2
  integer(kind=i_def), dimension(ndf_w2h),    intent(in)    :: map_w2h
  integer(kind=i_def), dimension(ndf_w2v),    intent(in)    :: map_w2v
  integer(kind=i_def), dimension(ndf_w3_2d),  intent(in)    :: map_w3_2d
  real(kind=r_double), dimension(undf_w2h),   intent(inout) :: uv
  real(kind=r_double), dimension(undf_w2v),   intent(inout) :: w
  real(kind=r_double), dimension(undf_w2),    intent(in)    :: uvw
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ew
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ns

  ! Internal variables
  integer(kind=i_def) :: df, k, j
  integer(kind=i_def) :: hori_dofs_to_do
  logical(kind=l_def) :: lowest_order
  logical(kind=l_def) :: dof3_is_N

  if (ndf_w2 == 6) then
    lowest_order = .true.
  else
    lowest_order = .false.
  end if

  if (lowest_order) then
    hori_dofs_to_do = face_selector_ew(map_w3_2d(1)) + face_selector_ns(map_w3_2d(1))
    if (face_selector_ns(map_w3_2d(1)) == 2 .and. face_selector_ew(map_w3_2d(1)) == 1) then
      dof3_is_N = .true.
    else
      dof3_is_N = .false.
    end if
  else
    hori_dofs_to_do = ndf_w2h
    dof3_is_N = .false.
  end if

  ! Loop over horizontal W2 DoFs
  do j = 1, hori_dofs_to_do
    df = j
    if (j == 3 .and. dof3_is_N) df = N

    ! Loop over layers
    do k = 0, nlayers-1
      uv(map_w2h(df) + k) = uvw(map_w2(df) + k)
    end do
  end do

  ! Loop over vertical W2 DoFs
  if (ndf_w2v == 2) then
    ! Optimal looping for lowest-order space
    do k = 0, nlayers
      w(map_w2v(1) + k) = uvw(map_w2(5) + k)
    end do

  else
    ! General looping for not lowest-order
    do df = 1, ndf_w2v
      do k = 0, nlayers-1
        w(map_w2v(df) + k) = uvw(map_w2(ndf_w2h+df) + k)
      end do
    end do
  end if

end subroutine split_w2_field_code_r_double

! R_SINGLE PRECISION
! ==================
subroutine split_w2_field_code_r_single(nlayers,                         &
                                        uv, w, uvw,                      &
                                        face_selector_ew,                &
                                        face_selector_ns,                &
                                        ndf_w2h,   undf_w2h,   map_w2h,  &
                                        ndf_w2v,   undf_w2v,   map_w2v,  &
                                        ndf_w2,    undf_w2,    map_w2,   &
                                        ndf_w3_2d, undf_w3_2d, map_w3_2d )
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, ndf_w2h, ndf_w2v
  integer(kind=i_def), intent(in) :: undf_w2, undf_w2h, undf_w2v
  integer(kind=i_def), intent(in) :: ndf_w3_2d, undf_w3_2d
  integer(kind=i_def), dimension(ndf_w2),     intent(in)    :: map_w2
  integer(kind=i_def), dimension(ndf_w2h),    intent(in)    :: map_w2h
  integer(kind=i_def), dimension(ndf_w2v),    intent(in)    :: map_w2v
  integer(kind=i_def), dimension(ndf_w3_2d),  intent(in)    :: map_w3_2d
  real(kind=r_single), dimension(undf_w2h),   intent(inout) :: uv
  real(kind=r_single), dimension(undf_w2v),   intent(inout) :: w
  real(kind=r_single), dimension(undf_w2),    intent(in)    :: uvw
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ew
  integer(kind=i_def), dimension(undf_w3_2d), intent(in)    :: face_selector_ns

  ! Internal variables
  integer(kind=i_def) :: df, k, j
  integer(kind=i_def) :: hori_dofs_to_do
  logical(kind=l_def) :: lowest_order
  logical(kind=l_def) :: dof3_is_N

  if (ndf_w2 == 6) then
    lowest_order = .true.
  else
    lowest_order = .false.
  end if

  if (lowest_order) then
    hori_dofs_to_do = face_selector_ew(map_w3_2d(1)) + face_selector_ns(map_w3_2d(1))
    if (face_selector_ns(map_w3_2d(1)) == 2 .and. face_selector_ew(map_w3_2d(1)) == 1) then
      dof3_is_N = .true.
    else
      dof3_is_N = .false.
    end if
  else
    hori_dofs_to_do = ndf_w2h
    dof3_is_N = .false.
  end if

  ! Loop over horizontal W2 DoFs
  do j = 1, hori_dofs_to_do
    df = j
    if (j == 3 .and. dof3_is_N) df = N

    ! Loop over layers
    do k = 0, nlayers-1
      uv(map_w2h(df) + k) = uvw(map_w2(df) + k)
    end do
  end do

  ! Loop over vertical W2 DoFs
  if (ndf_w2v == 2) then
    ! Optimal looping for lowest-order space
    do k = 0, nlayers
      w(map_w2v(1) + k) = uvw(map_w2(5) + k)
    end do

  else
    ! General looping for not lowest-order
    do df = 1, ndf_w2v
      do k = 0, nlayers-1
        w(map_w2v(df) + k) = uvw(map_w2(ndf_w2h+df) + k)
      end do
    end do
  end if

end subroutine split_w2_field_code_r_single

end module split_w2_field_kernel_mod
