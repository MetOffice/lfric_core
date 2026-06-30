!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Makes tridiagonal matrix to get mixing ratio from Wtheta to density in
!> shifted W3
!> @details To transform from a mixing ratio in Wtheta to a density in shifted W3
!>  we need a  tridiagonal transform matrix whose elements are the integrals
!> \f[ \int( \phi_i * \Phi_j \Psi_k dV \f]
!> where \phi_i are the basis functions of the W3 shifted space, \Phi_j are
!> the basis functions of Wtheta and \Psi_k are the basis functions for W3.
!> In this kernel we get the triadiagonal elements of the matrix by multiplying
!> the integrals by the coefficients of the dry density.
!> Only valid for the lowest order elements.
!>
module sci_proj_mr_to_sh_rho_rhs_update_kernel_mod

  use argument_mod,      only : arg_type,                  &
                                GH_FIELD, GH_REAL,         &
                                GH_WRITE, GH_READ,         &
                                ANY_DISCONTINUOUS_SPACE_3, &
                                CELL_COLUMN
  use constants_mod,     only : i_def
  use, intrinsic :: iso_fortran_env, only: real32, real64
  use fs_continuity_mod, only : W3
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  !>
  type, public, extends(kernel_type) :: proj_mr_to_sh_rho_rhs_update_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                      &
         arg_type(GH_FIELD*3, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3), & ! tri_below/diag/above
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                        & ! rho_d
         arg_type(GH_FIELD*4, GH_REAL, GH_READ,  W3)                         & ! I_lower/upper
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: proj_mr_to_sh_rho_rhs_update_code

  ! Generic interface for real32 and real64 types
  interface proj_mr_to_sh_rho_rhs_update_code
    module procedure  &
      proj_mr_to_sh_rho_rhs_update_code_real32, &
      proj_mr_to_sh_rho_rhs_update_code_real64
  end interface
contains

!> @brief Compute the terms of the tridiagonal matrix for transforming from
!> mixing ratio in Wtheta to density in shifted W3.
!! @param[in] nlayers_shifted Number of layers in the shifted mesh
!! @param[in,out] tri_below The below-diagonal elements of the tridiagonal matrix to
!! be returned. It is a field in shifted W3.
!! @param[in,out] tri_diag The central diagonal elements of the tridiagonal matrix to
!! be returned. It is a field in shifted W3.
!! @param[in,out] tri_above The above-diagonal elements of the tridiagonal matrix to
!! be returned. It is a field in shifted W3.
!! @param[in] rho_d The dry density in W3 on the original mesh.
!! @param[in] I_lower_i_ip1 The integral of the (i+1)-th Wtheta basis function on
!! the lower half of the original mesh. Is a W3 field.
!! @param[in] I_lower_i_i The integral of the i-th Wtheta basis function on
!! the lower half of the original mesh. Is a W3 field.
!! @param[in] I_upper_i_i The integral of the i-th Wtheta basis function on
!! the upper half of the original mesh. Is a W3 field.
!! @param[in] I_upper_i_im1 The integral of the (i-1)-th Wtheta basis function on
!! the upper half of the original mesh. Is a W3 field.
!! @param[in] ndf_sh_w3 The number of degrees of freedom per cell for shifted w3
!! @param[in] undf_sh_w3 The number of unique degrees of freedom for shifted w3
!! @param[in] map_sh_w3 Dofmap for the cell at the base of the column for shifted w3
!! @param[in] ndf_w3 The number of degrees of freedom per cell for w3
!! @param[in] undf_w3 The number of unique degrees of freedom for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3

! REAL32 PRECISION
! ==================
subroutine proj_mr_to_sh_rho_rhs_update_code_real32(                             &
                                              nlayers_shifted,                   &
                                              tri_below,                         &
                                              tri_diag,                          &
                                              tri_above,                         &
                                              rho_d,                             &
                                              I_lower_i_ip1,                     &
                                              I_lower_i_i,                       &
                                              I_upper_i_i,                       &
                                              I_upper_i_im1,                     &
                                              ndf_sh_w3, undf_sh_w3, map_sh_w3,  &
                                              ndf_w3, undf_w3, map_w3            &
                                            )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers_shifted
  integer(kind=i_def), intent(in) :: ndf_w3, ndf_sh_w3
  integer(kind=i_def), intent(in) :: undf_w3, undf_sh_w3
  integer(kind=i_def), dimension(ndf_w3),    intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_sh_w3), intent(in) :: map_sh_w3

  real(kind=real32),   dimension(undf_sh_w3),  intent(inout) :: tri_below
  real(kind=real32),   dimension(undf_sh_w3),  intent(inout) :: tri_diag
  real(kind=real32),   dimension(undf_sh_w3),  intent(inout) :: tri_above
  real(kind=real32),   dimension(undf_w3),     intent(in)    :: rho_d
  real(kind=real32),   dimension(undf_w3),     intent(in)    :: I_lower_i_ip1
  real(kind=real32),   dimension(undf_w3),     intent(in)    :: I_upper_i_i
  real(kind=real32),   dimension(undf_w3),     intent(in)    :: I_lower_i_i
  real(kind=real32),   dimension(undf_w3),     intent(in)    :: I_upper_i_im1

  ! Internal variables
  integer(kind=i_def) :: df, k

  ! Calculation for top and bottom layers (bottom is k=0)
  k = nlayers_shifted - 1
  do df = 1, ndf_sh_w3
    tri_below(map_sh_w3(df)) = 0.0_real32
    tri_diag(map_sh_w3(df)) = rho_d(map_w3(df)) * I_lower_i_i(map_w3(df))
    tri_above(map_sh_w3(df)) = rho_d(map_w3(df)) * I_lower_i_ip1(map_w3(df))

    tri_below(map_sh_w3(df)+k) = rho_d(map_w3(df)+k-1) * I_upper_i_im1(map_w3(df)+k-1)
    tri_diag(map_sh_w3(df)+k) = rho_d(map_w3(df)+k-1) * I_upper_i_i(map_w3(df)+k-1)
    tri_above(map_sh_w3(df)+k) = 0.0_real32
  end do

  ! Calculation for generic internal layers
  do k = 1, nlayers_shifted-2
    do df = 1, ndf_sh_w3
      tri_below(map_sh_w3(df)+k) = rho_d(map_w3(df)+k-1) * I_upper_i_im1(map_w3(df)+k-1)
      tri_diag(map_sh_w3(df)+k) = (rho_d(map_w3(df)+k) * I_lower_i_i(map_w3(df)+k) &
                                  + rho_d(map_w3(df)+k-1) * I_upper_i_i(map_w3(df)+k-1))
      tri_above(map_sh_w3(df)+k) = rho_d(map_w3(df)+k) * I_lower_i_ip1(map_w3(df)+k)
    end do
  end do

end subroutine proj_mr_to_sh_rho_rhs_update_code_real32

! REAL64 PRECISION
! ==================
subroutine proj_mr_to_sh_rho_rhs_update_code_real64(                             &
                                              nlayers_shifted,                   &
                                              tri_below,                         &
                                              tri_diag,                          &
                                              tri_above,                         &
                                              rho_d,                             &
                                              I_lower_i_ip1,                     &
                                              I_lower_i_i,                       &
                                              I_upper_i_i,                       &
                                              I_upper_i_im1,                     &
                                              ndf_sh_w3, undf_sh_w3, map_sh_w3,  &
                                              ndf_w3, undf_w3, map_w3            &
                                            )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers_shifted
  integer(kind=i_def), intent(in) :: ndf_w3, ndf_sh_w3
  integer(kind=i_def), intent(in) :: undf_w3, undf_sh_w3
  integer(kind=i_def), dimension(ndf_w3),    intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_sh_w3), intent(in) :: map_sh_w3

  real(kind=real64),   dimension(undf_sh_w3),  intent(inout) :: tri_below
  real(kind=real64),   dimension(undf_sh_w3),  intent(inout) :: tri_diag
  real(kind=real64),   dimension(undf_sh_w3),  intent(inout) :: tri_above
  real(kind=real64),   dimension(undf_w3),     intent(in)    :: rho_d
  real(kind=real64),   dimension(undf_w3),     intent(in)    :: I_lower_i_ip1
  real(kind=real64),   dimension(undf_w3),     intent(in)    :: I_upper_i_i
  real(kind=real64),   dimension(undf_w3),     intent(in)    :: I_lower_i_i
  real(kind=real64),   dimension(undf_w3),     intent(in)    :: I_upper_i_im1

  ! Internal variables
  integer(kind=i_def) :: df, k

  ! Calculation for top and bottom layers (bottom is k=0)
  k = nlayers_shifted - 1
  do df = 1, ndf_sh_w3
    tri_below(map_sh_w3(df)) = 0.0_real64
    tri_diag(map_sh_w3(df)) = rho_d(map_w3(df)) * I_lower_i_i(map_w3(df))
    tri_above(map_sh_w3(df)) = rho_d(map_w3(df)) * I_lower_i_ip1(map_w3(df))

    tri_below(map_sh_w3(df)+k) = rho_d(map_w3(df)+k-1) * I_upper_i_im1(map_w3(df)+k-1)
    tri_diag(map_sh_w3(df)+k) = rho_d(map_w3(df)+k-1) * I_upper_i_i(map_w3(df)+k-1)
    tri_above(map_sh_w3(df)+k) = 0.0_real64
  end do

  ! Calculation for generic internal layers
  do k = 1, nlayers_shifted-2
    do df = 1, ndf_sh_w3
      tri_below(map_sh_w3(df)+k) = rho_d(map_w3(df)+k-1) * I_upper_i_im1(map_w3(df)+k-1)
      tri_diag(map_sh_w3(df)+k) = (rho_d(map_w3(df)+k) * I_lower_i_i(map_w3(df)+k) &
                                  + rho_d(map_w3(df)+k-1) * I_upper_i_i(map_w3(df)+k-1))
      tri_above(map_sh_w3(df)+k) = rho_d(map_w3(df)+k) * I_lower_i_ip1(map_w3(df)+k)
    end do
  end do

end subroutine proj_mr_to_sh_rho_rhs_update_code_real64

end module sci_proj_mr_to_sh_rho_rhs_update_kernel_mod
