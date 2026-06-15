!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Samples a field from W2 to shifted W2.
!> @details Calculates a shifted W2 field from a W2 by sampling the velocity
!!          and using area weightings. This preserves a horizontal velocity that
!!          varies linearly in the vertical.
!!          This kernel only works with the lowest-order finite element spaces
!!          on quadrilateral cells.
module sci_sample_w2_to_sh_w2_kernel_mod

  use argument_mod,                  only : arg_type, GH_INTEGER,              &
                                            GH_FIELD, GH_REAL,                 &
                                            GH_READ, GH_WRITE,                 &
                                            ANY_DISCONTINUOUS_SPACE_2,         &
                                            ANY_DISCONTINUOUS_SPACE_3,         &
                                            CELL_COLUMN
  use constants_mod,                 only : r_def, i_def
  use fs_continuity_mod,             only : W2
  use kernel_mod,                    only : kernel_type
  use reference_element_mod,         only : N, E, S, W, T, B
  use sci_face_selector_support_mod, only : face_from_face_selector

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: sample_w2_to_sh_w2_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                       &
         ! NB: This is to be used to write to a continuous W2 field, but using
         ! a discontinuous data pattern, so use discontinuous metadata
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! field_w2_sh
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  W2),                        & ! field_w2
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! area_w2_sh
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  W2),                        & ! area_w2
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face_selector_ew
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)  & ! face_selector_ns
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: sample_w2_to_sh_w2_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: sample_w2_to_sh_w2_code

contains

!> @brief Samples a W2 field in the W2 shifted space.
!>
!> @param[in] nlayers_sh Number of layers in the shifted mesh
!> @param[in,out] field_w2_sh Field in the shifted W2 space to be returned.
!> @param[in] field_w2 Original field in W2 to be used.
!> @param[in] area_w2_sh The areas of cell faces of the shifted mesh. In W2 shifted.
!> @param[in] area_w2 The areas of cell faces of the primal mesh. In W2.
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
subroutine sample_w2_to_sh_w2_code( nlayers_sh,       &
                                    field_w2_sh,      &
                                    field_w2,         &
                                    area_w2_sh,       &
                                    area_w2,          &
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

  real(kind=r_def),    dimension(undf_w2_sh), intent(inout) :: field_w2_sh
  real(kind=r_def),    dimension(undf_w2),       intent(in) :: field_w2
  real(kind=r_def),    dimension(undf_w2_sh),    intent(in) :: area_w2_sh
  real(kind=r_def),    dimension(undf_w2),       intent(in) :: area_w2
  integer(kind=i_def), dimension(undf_w3_2d),    intent(in) :: face_selector_ew
  integer(kind=i_def), dimension(undf_w3_2d),    intent(in) :: face_selector_ns

  ! Internal variables
  integer(kind=i_def) :: df, k, j


  ! Loop over horizontal W2 DoFs
  do j = 1, ABS(face_selector_ew(map_w3_2d(1))) + ABS(face_selector_ns(map_w3_2d(1)))
    df = face_from_face_selector(j, face_selector_ew(map_w3_2d(1)), face_selector_ns(map_w3_2d(1)))

    ! Bottom layer
    field_w2_sh(map_w2_sh(df)) =                                               &
      (area_w2_sh(map_w2_sh(df)) / area_w2(map_w2(df)))                        &
      * field_w2(map_w2(df))

    ! Loop over inter layers (but not top and bottom)
    do k = 1, nlayers_sh - 2

      field_w2_sh(map_w2_sh(df)+k) =                                           &
        0.5_r_def * area_w2_sh(map_w2_sh(df)+k)                                &
        * ( field_w2(map_w2(df)+k) / area_w2(map_w2(df)+k)                     &
            + field_w2(map_w2(df)+k-1) / area_w2(map_w2(df)+k-1) )
    end do

    ! Top layer
    k = nlayers_sh-1
    field_w2_sh(map_w2_sh(df)+k) =                                             &
      (area_w2_sh(map_w2_sh(df)+k) / area_w2(map_w2(df)+k-1))                  &
      * field_w2(map_w2(df)+k-1)
  end do

  ! Loop over vertical W2 DoFs. Only need to do bottom DoF of each cell.
  ! Top and bottom values are the same as the original space
  field_w2_sh(map_w2_sh(B)) = field_w2(map_w2(B))

  do k = 1, nlayers_sh - 1
    ! Values are obtained from linear interpolation.
    df = B
    field_w2_sh(map_w2_sh(df)+k) = 0.5_r_def * area_w2_sh(map_w2_sh(df)+k)     &
      * ( field_w2(map_w2(df)+k-1) / area_w2(map_w2(df)+k-1)                   &
          + field_w2(map_w2(df)+k) / area_w2(map_w2(df)+k) )
  end do

  field_w2_sh(map_w2_sh(T)+nlayers_sh-1) = field_w2(map_w2(T)+nlayers_sh-2)

end subroutine sample_w2_to_sh_w2_code

end module sci_sample_w2_to_sh_w2_kernel_mod
