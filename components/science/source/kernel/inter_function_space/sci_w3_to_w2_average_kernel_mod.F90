!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes a W3 scalar field at W2 points by averaging.
!> @details Kernel to average a W3 lowest-order field to W2 or W2H points.

module sci_w3_to_w2_average_kernel_mod

  use argument_mod,          only : arg_type,          &
                                    GH_FIELD, GH_REAL, &
                                    GH_INC, GH_READ,   &
                                    CELL_COLUMN
  use constants_mod,         only : i_def, r_def
  use, intrinsic :: iso_fortran_env, only: real32, real64
  use fs_continuity_mod,     only : W2, W3
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : B

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: w3_to_w2_average_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/            &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W2), &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3), &
         arg_type(GH_FIELD, GH_REAL, GH_READ, W2)  &
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  ! -----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: w3_to_w2_average_kernel_code

  ! Generic interface for real32 and real64 types
  interface w3_to_w2_average_kernel_code
    module procedure  &
      w3_to_w2_average_kernel_code_real32, &
      w3_to_w2_average_kernel_code_real64
  end interface


  contains

  !> @brief Computes a W3 scalar field at W2 points by averaging.
  !> @param[in]     nlayers       Number of layers in the mesh
  !> @param[in,out] field_w2      Output field in W2 or W2H space
  !> @param[in]     field_w3      Input field in W3 space
  !> @param[in]     rmultiplicity Reciprocal multiplicity of W2 DoFs
  !> @param[in]     ndf_w2        Number of degrees of freedom per cell for W2
  !> @param[in]     undf_w2       Number of DoFs per partition for W2
  !> @param[in]     map_w2        Map of bottom-layer DoFs for W2
  !> @param[in]     ndf_w3        Number of degrees of freedom per cell for W3
  !> @param[in]     undf_w3       Number of DoFs per partition for W3
  !> @param[in]     map_w3        Map of bottom-layer DoFs for W3
  subroutine w3_to_w2_average_kernel_code_real32( nlayers,                 &
                                                    field_w2,                &
                                                    field_w3,                &
                                                    rmultiplicity,           &
                                                    ndf_w2, undf_w2, map_w2, &
                                                    ndf_w3, undf_w3, map_w3  )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ndf_w2, undf_w2
    integer(kind=i_def), intent(in)    :: ndf_w3, undf_w3

    integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)
    integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)

    real(kind=real32), intent(inout) :: field_w2(undf_w2)
    real(kind=r_def),    intent(in)    :: rmultiplicity(undf_w2)
    real(kind=real32), intent(in)    :: field_w3(undf_w3)

    ! Internal variables
    integer(kind=i_def) :: df, k

    ! Loop over horizontal faces
    do df = 1, 4
      do k = 0, nlayers-1
        ! Multiplicity gives averaging factor. Can just use multplicity from
        ! the bottom layer as this should be the same for all layers
        field_w2(map_w2(df) + k) = field_w2(map_w2(df) + k) +                  &
            field_w3(map_w3(1) + k) * real(rmultiplicity(map_w2(df)), real32)
      end do
    end do

    ! Loop over vertical DoFs, if the output field is in W2 rather than W2H
    if (ndf_w2 == 6) then
      ! At bottom boundary, take the value from the bottom cell
      field_w2(map_w2(B)) = field_w3(map_w3(1))
      ! Loop through internal layers
      do k = 1, nlayers - 1
        field_w2(map_w2(B)+k) =                                                &
              0.5_real32*(field_w3(map_w3(1)+k) + field_w3(map_w3(1)+k-1))
      end do
      ! At bottom boundary, take the value from the top cell
      k = nlayers
      field_w2(map_w2(B)+k) = field_w3(map_w3(1)+k-1)
    end if

  end subroutine w3_to_w2_average_kernel_code_real32

  !> @brief Computes a W3 scalar field at W2 points by averaging.
  !> @param[in]     nlayers       Number of layers in the mesh
  !> @param[in,out] field_w2      Output field in W2 or W2H space
  !> @param[in]     field_w3      Input field in W3 space
  !> @param[in]     rmultiplicity Reciprocal multiplicity of W2 DoFs
  !> @param[in]     ndf_w2        Number of degrees of freedom per cell for W2
  !> @param[in]     undf_w2       Number of DoFs per partition for W2
  !> @param[in]     map_w2        Map of bottom-layer DoFs for W2
  !> @param[in]     ndf_w3        Number of degrees of freedom per cell for W3
  !> @param[in]     undf_w3       Number of DoFs per partition for W3
  !> @param[in]     map_w3        Map of bottom-layer DoFs for W3
  subroutine w3_to_w2_average_kernel_code_real64( nlayers,                 &
                                                    field_w2,                &
                                                    field_w3,                &
                                                    rmultiplicity,           &
                                                    ndf_w2, undf_w2, map_w2, &
                                                    ndf_w3, undf_w3, map_w3  )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ndf_w2, undf_w2
    integer(kind=i_def), intent(in)    :: ndf_w3, undf_w3

    integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)
    integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)

    real(kind=real64), intent(inout) :: field_w2(undf_w2)
    real(kind=r_def),    intent(in)    :: rmultiplicity(undf_w2)
    real(kind=real64), intent(in)    :: field_w3(undf_w3)

    ! Internal variables
    integer(kind=i_def) :: df, k

    ! Loop over horizontal faces
    do df = 1, 4
      do k = 0, nlayers-1
        ! Multiplicity gives averaging factor. Can just use multplicity from
        ! the bottom layer as this should be the same for all layers
        field_w2(map_w2(df) + k) = field_w2(map_w2(df) + k) +                  &
            field_w3(map_w3(1) + k) * real(rmultiplicity(map_w2(df)), real64)
      end do
    end do

    ! Loop over vertical DoFs, if the output field is in W2 rather than W2H
    if (ndf_w2 == 6) then
      ! At bottom boundary, take the value from the bottom cell
      field_w2(map_w2(B)) = field_w3(map_w3(1))
      ! Loop through internal layers
      do k = 1, nlayers - 1
        field_w2(map_w2(B)+k) =                                                &
              0.5_real64*(field_w3(map_w3(1)+k) + field_w3(map_w3(1)+k-1))
      end do
      ! At bottom boundary, take the value from the top cell
      k = nlayers
      field_w2(map_w2(B)+k) = field_w3(map_w3(1)+k-1)
    end if

  end subroutine w3_to_w2_average_kernel_code_real64

end module sci_w3_to_w2_average_kernel_mod
