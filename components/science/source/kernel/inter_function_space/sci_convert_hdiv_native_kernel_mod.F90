!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes three physical Cartesian component fields from a HDiv field
!> @details Applies the div-conforming Piola transform to a computational HDiv
!!          vector field, and returns the 3 Cartesian components of the
!!          corresponding physical field as separate W2 fields.
!!
!!          This version differs from "dg_convert_hdiv_native", in that the
!!          target field must be in W2, and so horizontally continuous. It also
!!          differs from the "convert_hdiv_field" kernel as the coordinate
!!          fields must be the mesh's "native" coordinates, which allows some
!!          optimisation.

module sci_convert_hdiv_native_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,       &
                                    GH_FIELD, GH_REAL, GH_INC, &
                                    GH_READ, ANY_SPACE_9,      &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    GH_DIFF_BASIS, GH_BASIS,   &
                                    CELL_COLUMN, GH_EVALUATOR
use constants_mod,           only : r_def, i_def
use fs_continuity_mod,       only : W2

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: convert_hdiv_native_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                          &
      arg_type(GH_FIELD*3, GH_REAL, GH_INC,  W2),                              &
      arg_type(GH_FIELD,   GH_REAL, GH_READ, W2),                              &
      arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),                     &
      arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3)        &
  /)
  type(func_type) :: meta_funcs(2) = (/                                        &
      func_type(W2, GH_BASIS),                                                 &
      func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                          &
  /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: convert_hdiv_native_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: convert_hdiv_native_code
contains

!> @brief Computes three physical Cartesian component fields from a HDiv field
!> @param[in]     nlayers            Number of layers in the mesh
!> @param[in,out] physical_field_1   X component of the vector field
!> @param[in,out] physical_field_2   Y component of the vector field
!> @param[in,out] physical_field_3   Z component of the vector field
!> @param[in]     hdiv_field         Input HDiv (W2) field
!> @param[in]     chi_1              Native coordinates in the first direction
!> @param[in]     chi_2              Native coordinates in the second direction
!> @param[in]     chi_3              Native coordinates in the third direction
!> @param[in]     panel_id           Field storing mesh panel ID for each column
!> @param[in]     ndf_w2             Num DoFs per cell for W2
!> @param[in]     undf_w2            Num DoFs in this partition for W2
!> @param[in]     map_w2             Map of lowest-level DoFs for W2
!> @param[in]     basis_w2           W2 basis functions evaluated at W2 DoFs
!> @param[in]     ndf_chi            Num DoFs per cell for Wchi
!> @param[in]     undf_chi           Num DoFs in this partition for Wchi
!> @param[in]     map_chi            Map of lowest-level DoFs for Wchi
!> @param[in]     basis_chi          Wchi basis functions evaluated at W2 DoFs
!> @param[in]     diff_basis_chi     Derivatives of Wchi basis functions,
!!                                   evaluated at W2 DoFs
!> @param[in]     ndf_pid            Num DoFs per cell for panel ID
!> @param[in]     undf_pid           Num DoFs in this partition for panel ID
!> @param[in]     map_pid            Map of lowest-level DoFs for panel ID
subroutine convert_hdiv_native_code(nlayers,                                   &
                                    physical_field_1,                          &
                                    physical_field_2,                          &
                                    physical_field_3,                          &
                                    hdiv_field,                                &
                                    chi_1, chi_2, chi_3,                       &
                                    panel_id,                                  &
                                    ndf_w2, undf_w2, map_w2,                   &
                                    basis_w2,                                  &
                                    ndf_chi, undf_chi, map_chi,                &
                                    basis_chi, diff_basis_chi,                 &
                                    ndf_pid, undf_pid, map_pid                 &
                                    )

  use sci_native_jacobian_mod, only: native_jacobian

  use base_mesh_config_mod,      only: geometry, topology
  use finite_element_config_mod, only: coord_system
  use planet_config_mod,         only: scaled_radius

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in)    :: ndf_chi
  integer(kind=i_def), intent(in)    :: undf_chi
  integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid

  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)

  real(kind=r_def),    intent(in)    :: hdiv_field(undf_w2)
  real(kind=r_def),    intent(in)    :: chi_1(undf_chi)
  real(kind=r_def),    intent(in)    :: chi_2(undf_chi)
  real(kind=r_def),    intent(in)    :: chi_3(undf_chi)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  real(kind=r_def),    intent(inout) :: physical_field_1(undf_w2)
  real(kind=r_def),    intent(inout) :: physical_field_2(undf_w2)
  real(kind=r_def),    intent(inout) :: physical_field_3(undf_w2)
  real(kind=r_def),    intent(in)    :: basis_chi(1,ndf_chi,ndf_w2)
  real(kind=r_def),    intent(in)    :: diff_basis_chi(3,ndf_chi,ndf_w2)
  real(kind=r_def),    intent(in)    :: basis_w2(3,ndf_w2,ndf_w2)

  ! Internal variables
  integer(kind=i_def) :: i, j
  integer(kind=i_def) :: df_w2, df_w2_in, df_chi
  integer(kind=i_def) :: w2_idx, chi_idx
  real(kind=r_def) :: jacobian(nlayers,3,3), dj(nlayers)
  real(kind=r_def) :: vector_in(nlayers,3), vector_out(nlayers,3)
  real(kind=r_def) :: chi_1_e(ndf_chi), chi_2_e(ndf_chi), chi_3_e(nlayers,ndf_chi)

  integer(kind=i_def) :: ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Fill arrays of coordinate values
  do df_chi = 1, ndf_chi
    chi_idx = map_chi(df_chi)
    chi_1_e(df_chi) = chi_1(chi_idx)
    chi_2_e(df_chi) = chi_2(chi_idx)
    chi_3_e(:,df_chi) = chi_3(chi_idx : chi_idx+nlayers-1)
  end do

  ! Loop through W2 DoFs
  do df_w2 = 1, ndf_w2

    ! Compute Jacobian for whole column at this DoF
    call native_jacobian(                                        &
            coord_system, geometry, topology, scaled_radius,     &
            ndf_chi, nlayers, chi_1_e, chi_2_e, chi_3_e, ipanel, &
            basis_chi(:,:,df_w2),  diff_basis_chi(:,:,df_w2),    &
            jacobian, dj )

    ! Create vector of HDiv values at this point
    vector_in(:,:) = 0.0_r_def
    do df_w2_in = 1, ndf_w2
      w2_idx = map_w2(df_w2_in)
      do i = 1, 3
        vector_in(:,i) = vector_in(:,i)                                        &
            + hdiv_field(w2_idx : w2_idx+nlayers-1)*basis_w2(i,df_w2_in,df_w2)
      end do
    end do

    ! Calculate the contribution to components at df_w2
    vector_out(:,:) = 0.0_r_def
    do i = 1, 3
      do j = 1, 3
        vector_out(:,i) = vector_out(:,i) + jacobian(:,i,j)*vector_in(:,j)
      end do
      vector_out(:,i) = vector_out(:,i) / dj(:)
    end do

    w2_idx = map_w2(df_w2)
    physical_field_1(w2_idx : w2_idx+nlayers-1) = (                            &
        physical_field_1(w2_idx : w2_idx+nlayers-1) + vector_out(:,1)          &
    )
    physical_field_2(w2_idx : w2_idx+nlayers-1) = (                            &
        physical_field_2(w2_idx : w2_idx+nlayers-1) + vector_out(:,2)          &
    )
    physical_field_3(w2_idx : w2_idx+nlayers-1) = (                            &
        physical_field_3(w2_idx : w2_idx+nlayers-1) + vector_out(:,3)          &
    )
  end do

end subroutine convert_hdiv_native_code

end module sci_convert_hdiv_native_kernel_mod
