!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Maps a field from W2v to shifted W2v.
!> @details Calculates a shifted W2v field from a W2v field. The new vertical
!!          fluxes take the average from each of the half-levels from the
!!          original field. In conjunction with the mapping from W3 to
!!          shifted W3 this provides a consistent mapping of fluxes and
!!          densities to the shifted mesh.
!!          This kernel only works for the lowest-order elements.
module sci_consist_w2v_to_sh_w2v_kernel_mod

  use, intrinsic :: iso_fortran_env, only: real32, real64

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    GH_READ, GH_WRITE,         &
                                    ANY_DISCONTINUOUS_SPACE_2, &
                                    CELL_COLUMN
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : W2v
  use kernel_mod,            only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !! Psy layer.
  !!
  type, public, extends(kernel_type) :: consist_w2v_to_sh_w2v_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W2v)                        &
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: consist_w2v_to_sh_w2v_code

  ! Generic interface for real32 and real64 types
  interface consist_w2v_to_sh_w2v_code
    module procedure  &
      consist_w2v_to_sh_w2v_code_real32, &
      consist_w2v_to_sh_w2v_code_real64
  end interface

contains

!> @brief Maps a field from W2v to the W2v shifted space.
!> @param[in] nlayers_sh Number of layers in the shifted mesh
!> @param[in,out] field_w2v_sh Field in the shifted W2v space to be returned.
!> @param[in] field_w2v Original field in W2v to be used.
!> @param[in] ndf_w2v_sh Number of degrees of freedom per cell for W2v shifted
!> @param[in] undf_w2v_sh Number of (local) unique degrees of freedom for W2v shifted
!> @param[in] map_w2v_sh Dofmap for the cell at the base of the column for W2v shifted
!> @param[in] ndf_w2v Number of degrees of freedom per cell for W2v
!> @param[in] undf_w2v Number of (local) unique degrees of freedom for W2v
!> @param[in] map_w2v Dofmap for the cell at the base of the column for W2v

! REAL32 PRECISION
! ==================
subroutine consist_w2v_to_sh_w2v_code_real32(    nlayers_sh,        &
                                                 field_w2v_sh,      &
                                                 field_w2v,         &
                                                 ndf_w2v_sh,        &
                                                 undf_w2v_sh,       &
                                                 map_w2v_sh,        &
                                                 ndf_w2v,           &
                                                 undf_w2v,          &
                                                 map_w2v            &
                                                 )

  implicit none

  ! Arguments
  integer(kind=i_def),                            intent(in) :: nlayers_sh
  integer(kind=i_def),                            intent(in) :: ndf_w2v_sh, ndf_w2v
  integer(kind=i_def),                            intent(in) :: undf_w2v_sh, undf_w2v
  integer(kind=i_def), dimension(ndf_w2v_sh),     intent(in) :: map_w2v_sh
  integer(kind=i_def), dimension(ndf_w2v),        intent(in) :: map_w2v

  real(kind=real32),   dimension(undf_w2v_sh), intent(inout) :: field_w2v_sh
  real(kind=real32),   dimension(undf_w2v),       intent(in) :: field_w2v

  ! Internal variables
  integer(kind=i_def) :: bottom_df, top_df, k

  bottom_df = 1
  top_df = 2

  do k = 1, nlayers_sh - 1
    ! Loop over vertical W2 DoFs. Only need to do bottom DoF of each cell.
    ! Values are the average from the overlapping cells on the original mesh.
    field_w2v_sh(map_w2v_sh(bottom_df)+k) = &
      0.5_real32 * (field_w2v(map_w2v(bottom_df)+k-1) + field_w2v(map_w2v(bottom_df)+k) )
  end do

  ! Top and bottom values are the same as the original space
  field_w2v_sh(map_w2v_sh(bottom_df)) = field_w2v(map_w2v(bottom_df))
  field_w2v_sh(map_w2v_sh(top_df)+nlayers_sh-1) = field_w2v(map_w2v(top_df)+nlayers_sh-2)

end subroutine consist_w2v_to_sh_w2v_code_real32


! REAL64 PRECISION
! ==================
subroutine consist_w2v_to_sh_w2v_code_real64(    nlayers_sh,        &
                                                 field_w2v_sh,      &
                                                 field_w2v,         &
                                                 ndf_w2v_sh,        &
                                                 undf_w2v_sh,       &
                                                 map_w2v_sh,        &
                                                 ndf_w2v,           &
                                                 undf_w2v,          &
                                                 map_w2v            &
                                                 )

  implicit none

  ! Arguments
  integer(kind=i_def),                            intent(in) :: nlayers_sh
  integer(kind=i_def),                            intent(in) :: ndf_w2v_sh, ndf_w2v
  integer(kind=i_def),                            intent(in) :: undf_w2v_sh, undf_w2v
  integer(kind=i_def), dimension(ndf_w2v_sh),     intent(in) :: map_w2v_sh
  integer(kind=i_def), dimension(ndf_w2v),        intent(in) :: map_w2v

  real(kind=real64),   dimension(undf_w2v_sh), intent(inout) :: field_w2v_sh
  real(kind=real64),   dimension(undf_w2v),       intent(in) :: field_w2v

  ! Internal variables
  integer(kind=i_def) :: bottom_df, top_df, k

  bottom_df = 1
  top_df = 2

  do k = 1, nlayers_sh - 1
    ! Loop over vertical W2 DoFs. Only need to do bottom DoF of each cell.
    ! Values are the average from the overlapping cells on the original mesh.
    field_w2v_sh(map_w2v_sh(bottom_df)+k) = &
      0.5_real64 * (field_w2v(map_w2v(bottom_df)+k-1) + field_w2v(map_w2v(bottom_df)+k) )
  end do

  ! Top and bottom values are the same as the original space
  field_w2v_sh(map_w2v_sh(bottom_df)) = field_w2v(map_w2v(bottom_df))
  field_w2v_sh(map_w2v_sh(top_df)+nlayers_sh-1) = field_w2v(map_w2v(top_df)+nlayers_sh-2)

end subroutine consist_w2v_to_sh_w2v_code_real64

end module sci_consist_w2v_to_sh_w2v_kernel_mod
