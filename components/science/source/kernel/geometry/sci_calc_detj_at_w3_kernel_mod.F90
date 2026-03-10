!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the values of detj at w3 locations.
!>
module sci_calc_detj_at_w3_kernel_mod

  use argument_mod,      only : arg_type, func_type,         &
                                GH_FIELD, GH_REAL, GH_WRITE, &
                                GH_READ, ANY_SPACE_1,        &
                                GH_DIFF_BASIS, GH_BASIS,     &
                                ANY_DISCONTINUOUS_SPACE_3,   &
                                OWNED_AND_HALO_CELL_COLUMN,  &
                                GH_EVALUATOR,                &
                                GH_SCALAR, GH_INTEGER

  use constants_mod,     only : r_def, r_single, r_double, i_def
  use fs_continuity_mod, only : W3
  use kernel_mod,        only : kernel_type

  use base_mesh_config_mod,      only: geometry, topology
  use finite_element_config_mod, only: coord_system
  use planet_config_mod,         only: scaled_radius

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: calc_detj_at_w3_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
         arg_type(GH_FIELD,   GH_REAL,    GH_WRITE,  W3),                      &
         arg_type(GH_FIELD*3, GH_REAL,    GH_READ, ANY_SPACE_1),               &
         arg_type(GH_FIELD,   GH_REAL,    GH_READ, ANY_DISCONTINUOUS_SPACE_3)  &
         /)
    type(func_type) :: meta_funcs(1) = (/                                      &
         func_type(ANY_SPACE_1, GH_BASIS, GH_DIFF_BASIS)                       &
         /)
    integer :: operates_on = OWNED_AND_HALO_CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: calc_detj_at_w3_code

  ! Generic interface for real32 and real64 types
  interface calc_detj_at_w3_code
    module procedure  &
      calc_detj_at_w3_code_r_single, &
      calc_detj_at_w3_code_r_double
  end interface

contains

!> @param[in]  nlayers        Integer the number of layers
!> @param[in,out] detj_w3     The output field containing the detj values at w3 locations
!> @param[in]  chi1           1st coordinate field in Wchi
!> @param[in]  chi2           2nd coordinate field in Wchi
!> @param[in]  chi3           3rd coordinate field in Wchi
!> @param[in]  panel_id       Field giving the ID for mesh panels.
!> @param[in]  ndf_w3         The number of degrees of freedom per cell for the output field
!> @param[in]  undf_w3        The number of unique degrees of freedom for the output field
!> @param[in]  map_w3         Integer array holding the dofmap for the cell at the base of the column for the output field
!> @param[in]  nodal_basis_w3 The nodal basis functions evaluated at the nodal points for the w3 field
!> @param[in]  ndf_chi        The number of degrees of freedom per cell for the coordinate field
!> @param[in]  undf_chi       The number of unique degrees of freedom for the coordinate field
!> @param[in]  map_chi        Integer array holding the dofmap for the cell at the base of the column for the coordinate field
!> @param[in]  diff_basis_chi Basis functions of the coordinate space evaluated at the nodal points
!> @param[in]  diff_basis_chi The diff basis functions of the coordinate space evaluated at the nodal points
!> @param[in]  ndf_pid        Number of degrees of freedom per cell for panel_id
!> @param[in]  undf_pid       Number of unique degrees of freedom for panel_id
!> @param[in]  map_pid        Dofmap for the cell at the base of the column for panel_id

! R_SINGLE PRECISION
! ==================
subroutine calc_detj_at_w3_code_r_single( nlayers,                             &
                                          detj_w3,                             &
                                          chi1, chi2, chi3,                    &
                                          panel_id,                            &
                                          ndf_w3, undf_w3, map_w3,             &
                                          ndf_chi, undf_chi, map_chi,          &
                                          basis_chi, diff_basis_chi,           &
                                          ndf_pid, undf_pid, map_pid           &
                                        )

  use sci_coordinate_jacobian_mod, only: pointwise_coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def),                            intent(in)    :: nlayers
  integer(kind=i_def),                            intent(in)    :: ndf_w3
  integer(kind=i_def),                            intent(in)    :: undf_w3
  integer(kind=i_def),                            intent(in)    :: ndf_chi
  integer(kind=i_def),                            intent(in)    :: undf_chi
  integer(kind=i_def),                            intent(in)    :: ndf_pid
  integer(kind=i_def),                            intent(in)    :: undf_pid
  real(kind=r_single), dimension(undf_w3),        intent(inout) :: detj_w3
  real(kind=r_def), dimension(undf_chi),          intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_pid),          intent(in)    :: panel_id
  integer(kind=i_def), dimension(ndf_w3),         intent(in)    :: map_w3
  integer(kind=i_def), dimension(ndf_chi),        intent(in)    :: map_chi
  integer(kind=i_def), dimension(ndf_pid),        intent(in)    :: map_pid
  real(kind=r_def), dimension(3,ndf_chi,ndf_w3),  intent(in)    :: diff_basis_chi
  real(kind=r_def), dimension(1,ndf_chi,ndf_w3),  intent(in)    :: basis_chi

  ! Internal variables
  integer(kind=i_def)                  :: df, k
  integer(kind=i_def)                  :: ipanel
  real(kind=r_def), dimension(ndf_chi) :: chi1_e, chi2_e, chi3_e
  real(kind=r_def), dimension(3,3)     :: jacobian
  real(kind=r_def)                     :: detj

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1

    do df = 1,ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do

    do df = 1,ndf_w3
      call pointwise_coordinate_jacobian(coord_system, geometry,          &
                                         topology, scaled_radius,         &
                                         ndf_chi, chi1_e, chi2_e, chi3_e, &
                                         ipanel, basis_chi(:,:,df),       &
                                         diff_basis_chi(:,:,df),          &
                                         jacobian, detj)
      detj_w3(map_w3(df)+k) = real(detj, r_single)
    end do

  end do

end subroutine calc_detj_at_w3_code_r_single

! R_DOUBLE PRECISION
! ==================
subroutine calc_detj_at_w3_code_r_double( nlayers,                             &
                                          detj_w3,                             &
                                          chi1, chi2, chi3,                    &
                                          panel_id,                            &
                                          ndf_w3, undf_w3, map_w3,             &
                                          ndf_chi, undf_chi, map_chi,          &
                                          basis_chi, diff_basis_chi,           &
                                          ndf_pid, undf_pid, map_pid           &
                                        )

  use sci_coordinate_jacobian_mod, only: pointwise_coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def),                            intent(in)    :: nlayers
  integer(kind=i_def),                            intent(in)    :: ndf_w3
  integer(kind=i_def),                            intent(in)    :: undf_w3
  integer(kind=i_def),                            intent(in)    :: ndf_chi
  integer(kind=i_def),                            intent(in)    :: undf_chi
  integer(kind=i_def),                            intent(in)    :: ndf_pid
  integer(kind=i_def),                            intent(in)    :: undf_pid
  real(kind=r_double), dimension(undf_w3),        intent(inout) :: detj_w3
  real(kind=r_def), dimension(undf_chi),          intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_pid),          intent(in)    :: panel_id
  integer(kind=i_def), dimension(ndf_w3),         intent(in)    :: map_w3
  integer(kind=i_def), dimension(ndf_chi),        intent(in)    :: map_chi
  integer(kind=i_def), dimension(ndf_pid),        intent(in)    :: map_pid
  real(kind=r_def), dimension(3,ndf_chi,ndf_w3),  intent(in)    :: diff_basis_chi
  real(kind=r_def), dimension(1,ndf_chi,ndf_w3),  intent(in)    :: basis_chi

  ! Internal variables
  integer(kind=i_def)                  :: df, k
  integer(kind=i_def)                  :: ipanel
  real(kind=r_def), dimension(ndf_chi) :: chi1_e, chi2_e, chi3_e
  real(kind=r_def), dimension(3,3)     :: jacobian
  real(kind=r_def)                     :: detj

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1

    do df = 1,ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do

    do df = 1,ndf_w3
      call pointwise_coordinate_jacobian(coord_system, geometry,          &
                                         topology, scaled_radius,         &
                                         ndf_chi, chi1_e, chi2_e, chi3_e, &
                                         ipanel, basis_chi(:,:,df),       &
                                         diff_basis_chi(:,:,df),          &
                                         jacobian, detj)
      detj_w3(map_w3(df)+k) = real(detj, r_double)
    end do

  end do

end subroutine calc_detj_at_w3_code_r_double

end module sci_calc_detj_at_w3_kernel_mod
