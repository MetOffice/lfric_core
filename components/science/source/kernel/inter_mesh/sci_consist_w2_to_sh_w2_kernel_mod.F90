!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Maps a field from W2 to shifted W2.
!> @details Calculates a shifted W2 field from a W2 field. The new horizontal
!!          and vertical fluxes take the average from each of the half-levels
!!          from the original field. In conjunction with the mapping from W3 to
!!          shifted W3 this provides a consistent mapping of fluxes and
!!          densities to the shifted mesh.
!!          This kernel only works for the lowest-order elements.
module sci_consist_w2_to_sh_w2_kernel_mod

  use argument_mod,                  only : arg_type, GH_INTEGER,      &
                                            GH_FIELD, GH_REAL,         &
                                            GH_READ, GH_WRITE,         &
                                            ANY_DISCONTINUOUS_SPACE_2, &
                                            ANY_DISCONTINUOUS_SPACE_3, &
                                            CELL_COLUMN
  use constants_mod,                 only : r_def, i_def, r_single, r_double
  use fs_continuity_mod,             only : W2
  use kernel_mod,                    only : kernel_type
  use reference_element_mod,         only : T, B
  use sci_face_selector_support_mod, only : face_from_face_selector

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !! Psy layer.
  !!
  type, public, extends(kernel_type) :: consist_w2_to_sh_w2_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                       &
         ! NB: This is to be used to write to a continuous W2 field, but using
         ! a discontinuous data pattern, so use discontinuous metadata
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  W2),                        &
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)  &
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: consist_w2_to_sh_w2_code

  ! Generic interface for real32 and real64 types
  interface consist_w2_to_sh_w2_code
    module procedure  &
      consist_w2_to_sh_w2_code_single, &
      consist_w2_to_sh_w2_code_double
  end interface

contains

!> @brief Maps a field from W2 to the W2 shifted space.
!> @param[in] nlayers_sh Number of layers in the shifted mesh
!> @param[in,out] field_w2_sh Field in the shifted W2 space to be returned.
!> @param[in] field_w2 Original field in W2 to be used.
!> @param[in] face_selector_ew  2D field indicating which W/E faces
!!                              to loop over in this column
!> @param[in] face_selector_ns  2D field indicating which N/S faces
!!                              to loop over in this column
!> @param[in] ndf_w2_sh Number of degrees of freedom per cell for W2 shifted
!> @param[in] undf_w2_sh Number of (local) unique degrees of freedom for W2 shifted
!> @param[in] map_w2_sh Dofmap for the cell at the base of the column for W2 shifted
!> @param[in] ndf_w2 Number of degrees of freedom per cell for W2
!> @param[in] undf_w2 Number of (local) unique degrees of freedom for W2
!> @param[in] map_w2 Dofmap for the cell at the base of the column for W2
!> @param[in] ndf_w3_2d  Num of DoFs for 2D W3 per cell
!> @param[in] undf_w3_2d Num of DoFs for this partition for 2D W3
!> @param[in] map_w3_2d  Map for 2D W3

! R_SINGLE PRECISION
! ==================
subroutine consist_w2_to_sh_w2_code_single(  nlayers_sh,       &
                                             field_w2_sh,      &
                                             field_w2,         &
                                             face_selector_ew, &
                                             face_selector_ns, &
                                             ndf_w2_sh,        &
                                             undf_w2_sh,       &
                                             map_w2_sh,        &
                                             ndf_w2,           &
                                             undf_w2,          &
                                             map_w2,           &
                                             ndf_w3_2d,        &
                                             undf_w3_2d,       &
                                             map_w3_2d         &
                                             )

  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in) :: nlayers_sh
  integer(kind=i_def),                           intent(in) :: ndf_w2_sh, ndf_w2
  integer(kind=i_def),                           intent(in) :: undf_w2_sh, undf_w2
  integer(kind=i_def),                           intent(in) :: ndf_w3_2d, undf_w3_2d
  integer(kind=i_def), dimension(ndf_w2_sh),     intent(in) :: map_w2_sh
  integer(kind=i_def), dimension(ndf_w2),        intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w3_2d),     intent(in) :: map_w3_2d

  real(kind=r_single), dimension(undf_w2_sh), intent(inout) :: field_w2_sh
  real(kind=r_single), dimension(undf_w2),       intent(in) :: field_w2
  integer(kind=i_def), dimension(undf_w3_2d),    intent(in) :: face_selector_ew
  integer(kind=i_def), dimension(undf_w3_2d),    intent(in) :: face_selector_ns

  ! Internal variables
  integer(kind=i_def) :: df, k, j

  ! Loop over horizontal W2 DoFs
  do j = 1, ABS(face_selector_ew(map_w3_2d(1))) + ABS(face_selector_ns(map_w3_2d(1)))
    df = face_from_face_selector(j, face_selector_ew(map_w3_2d(1)), face_selector_ns(map_w3_2d(1)))

    ! Bottom boundary value
    field_w2_sh(map_w2_sh(df)) = 0.5_r_single * field_w2(map_w2(df))

    ! Loop over all interior layers of shifted mesh
    do k = 1, nlayers_sh - 2
      field_w2_sh(map_w2_sh(df)+k) =                                           &
          0.5_r_single * field_w2(map_w2(df)+k-1)                              &
          + 0.5_r_single * field_w2(map_w2(df)+k)
    end do

    ! Top boundary value
    k = nlayers_sh - 1
    field_w2_sh(map_w2_sh(df)+k) = 0.5_r_single * field_w2(map_w2(df)+k-1)
  end do

  do k = 1, nlayers_sh - 1
    ! Loop over vertical W2 DoFs. Only need to do bottom DoF of each cell.
    ! Values are the average from the overlapping cells on the original mesh.
    df = B
    field_w2_sh(map_w2_sh(df)+k) = &
      0.5_r_single * (field_w2(map_w2(df)+k-1) + field_w2(map_w2(df)+k) )
  end do

  ! Top and bottom values are the same as the original space
  field_w2_sh(map_w2_sh(B)) = field_w2(map_w2(B))
  field_w2_sh(map_w2_sh(T)+nlayers_sh-1) = field_w2(map_w2(T)+nlayers_sh-2)

end subroutine consist_w2_to_sh_w2_code_single

! R_DOUBLE PRECISION
! ==================
subroutine consist_w2_to_sh_w2_code_double(  nlayers_sh,       &
                                             field_w2_sh,      &
                                             field_w2,         &
                                             face_selector_ew, &
                                             face_selector_ns, &
                                             ndf_w2_sh,        &
                                             undf_w2_sh,       &
                                             map_w2_sh,        &
                                             ndf_w2,           &
                                             undf_w2,          &
                                             map_w2,           &
                                             ndf_w3_2d,        &
                                             undf_w3_2d,       &
                                             map_w3_2d         &
                                             )

  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in) :: nlayers_sh
  integer(kind=i_def),                           intent(in) :: ndf_w2_sh, ndf_w2
  integer(kind=i_def),                           intent(in) :: undf_w2_sh, undf_w2
  integer(kind=i_def),                           intent(in) :: ndf_w3_2d, undf_w3_2d
  integer(kind=i_def), dimension(ndf_w2_sh),     intent(in) :: map_w2_sh
  integer(kind=i_def), dimension(ndf_w2),        intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w3_2d),     intent(in) :: map_w3_2d

  real(kind=r_double), dimension(undf_w2_sh),    intent(inout) :: field_w2_sh
  real(kind=r_double), dimension(undf_w2),       intent(in) :: field_w2
  integer(kind=i_def), dimension(undf_w3_2d),    intent(in) :: face_selector_ew
  integer(kind=i_def), dimension(undf_w3_2d),    intent(in) :: face_selector_ns

  ! Internal variables
  integer(kind=i_def) :: df, k, j


  ! Loop over horizontal W2 DoFs
  do j = 1, ABS(face_selector_ew(map_w3_2d(1))) + ABS(face_selector_ns(map_w3_2d(1)))
    df = face_from_face_selector(j, face_selector_ew(map_w3_2d(1)), face_selector_ns(map_w3_2d(1)))

    ! Bottom boundary value
    field_w2_sh(map_w2_sh(df)) = 0.5_r_double * field_w2(map_w2(df))

    ! Loop over all interior layers of shifted mesh
    do k = 1, nlayers_sh - 2
      field_w2_sh(map_w2_sh(df)+k) =                                           &
          0.5_r_double * field_w2(map_w2(df)+k-1)                              &
          + 0.5_r_double * field_w2(map_w2(df)+k)
    end do

    ! Top boundary value
    k = nlayers_sh - 1
    field_w2_sh(map_w2_sh(df)+k) = 0.5_r_double * field_w2(map_w2(df)+k-1)
  end do

  do k = 1, nlayers_sh - 1
    ! Loop over vertical W2 DoFs. Only need to do bottom DoF of each cell.
    ! Values are the average from the overlapping cells on the original mesh.
    df = B
    field_w2_sh(map_w2_sh(df)+k) = &
      0.5_r_double * (field_w2(map_w2(df)+k-1) + field_w2(map_w2(df)+k) )
  end do

  ! Top and bottom values are the same as the original space
  field_w2_sh(map_w2_sh(B)) = field_w2(map_w2(B))
  field_w2_sh(map_w2_sh(T)+nlayers_sh-1) = field_w2(map_w2(T)+nlayers_sh-2)

end subroutine consist_w2_to_sh_w2_code_double

end module sci_consist_w2_to_sh_w2_kernel_mod
