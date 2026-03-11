!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes a HDiv wind field from co-located physical components
!>
!> @details Converts a wind field from physical components, located at W2 DoFs,
!!          into a computational wind field in W2.
module sci_convert_phys_to_hdiv_kernel_mod

  use argument_mod,            only : arg_type, func_type,       &
                                      GH_FIELD, GH_REAL,         &
                                      GH_WRITE, GH_READ,         &
                                      ANY_SPACE_9, GH_SCALAR,    &
                                      ANY_DISCONTINUOUS_SPACE_3, &
                                      GH_BASIS, GH_DIFF_BASIS,   &
                                      CELL_COLUMN, GH_EVALUATOR, &
                                      GH_INTEGER
  use constants_mod,           only : r_def, i_def
  use fs_continuity_mod,       only : W2
  use kernel_mod,              only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: convert_phys_to_hdiv_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                         &
         arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W2),                        &
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2),                        &
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2),                        &
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  W2),                        &
         arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_SPACE_9),               &
         arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)                              &
         /)
    type(func_type) :: meta_funcs(2) = (/                                       &
         func_type(W2, GH_BASIS),                                               &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                        &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: convert_phys_to_hdiv_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: convert_phys_to_hdiv_code

contains

!> @brief Computes a HDiv wind field from co-located physical components
!> @param[in]     nlayers        Number of layers in mesh
!> @param[in,out] u_hdiv         Computational wind field, in W2
!> @param[in]     u_lon          Zonal component of physical wind at W2 DoFs
!> @param[in]     u_lat          Meridional component of wind at W2 DoFs
!> @param[in]     u_up           Vertical component of physical wind at W2 DoFs
!> @param[in]     chi_1          1st coordinate field
!> @param[in]     chi_2          2nd coordinate field
!> @param[in]     chi_3          3rd coordinate field
!> @param[in]     panel_id       Field giving the ID for mesh panels
!> @param[in]     geometry       Integer indicating the domain geometry
!> @param[in]     ndf_w2         Number of DoFs per cell for W2
!> @param[in]     undf_w2        Number of DoFs for W2 for this partition
!> @param[in]     map_w2         Map of DoFs for lowest-layer cells for W2
!> @param[in]     basis_w2       Basis functions for W2 evaluated at W2 DoFs
!> @param[in]     ndf_chi        Number of DoFs per cell for Wchi
!> @param[in]     undf_chi       Number of DoFs for Wchi for this partition
!> @param[in]     map_chi        Map of DoFs for lowest-layer cells for Wchi
!> @param[in]     basis_chi      Basis functions for Wchi evaluated at W2 DoFs
!> @param[in]     diff_basis_chi Differential of the Wchi basis functions
!!                               evaluated at W2 DoFs
!> @param[in]     ndf_pid        Number of DoFs per cell for panel ID
!> @param[in]     undf_pid       Number of DoFs for panel ID for this partition
!> @param[in]     map_pid        Map of DoFs for lowest-layer cells for panel ID
subroutine convert_phys_to_hdiv_code( nlayers,        &
                                      u_hdiv,         &
                                      u_lon,          &
                                      u_lat,          &
                                      u_up,           &
                                      chi_1,          &
                                      chi_2,          &
                                      chi_3,          &
                                      panel_id,       &
                                      geometry,       &
                                      ndf_w2,         &
                                      undf_w2,        &
                                      map_w2,         &
                                      basis_w2,       &
                                      ndf_chi,        &
                                      undf_chi,       &
                                      map_chi,        &
                                      basis_chi,      &
                                      diff_basis_chi, &
                                      ndf_pid,        &
                                      undf_pid,       &
                                      map_pid )

  use sci_chi_transform_mod,      only : chi2llr
  use sci_coordinate_jacobian_mod, only : pointwise_coordinate_jacobian, &
                                          pointwise_coordinate_jacobian_inverse
  use coord_transform_mod,        only : sphere2cart_vector

  use base_mesh_config_mod,      only: topology, &
                                       geometry_spherical
  use finite_element_config_mod, only: coord_system
  use planet_config_mod,         only: scaled_radius

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, ndf_pid, ndf_chi
  integer(kind=i_def), intent(in) :: undf_w2, undf_pid, undf_chi
  integer(kind=i_def), intent(in) :: geometry

  integer(kind=i_def), intent(in) :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in) :: map_chi(ndf_chi)
  integer(kind=i_def), intent(in) :: map_pid(ndf_pid)

  real(kind=r_def),    intent(in) :: basis_w2(3,ndf_w2,ndf_w2)
  real(kind=r_def),    intent(in) :: basis_chi(1,ndf_chi,ndf_w2)
  real(kind=r_def),    intent(in) :: diff_basis_chi(3,ndf_chi,ndf_w2)

  real(kind=r_def),    intent(inout) :: u_hdiv(undf_w2)
  real(kind=r_def),    intent(in)    :: u_lon(undf_w2)
  real(kind=r_def),    intent(in)    :: u_lat(undf_w2)
  real(kind=r_def),    intent(in)    :: u_up(undf_w2)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  real(kind=r_def),    intent(in)    :: chi_1(undf_chi)
  real(kind=r_def),    intent(in)    :: chi_2(undf_chi)
  real(kind=r_def),    intent(in)    :: chi_3(undf_chi)

  ! Internal variables
  integer(kind=i_def) :: df_w2, df_chi, k, ipanel
  real(kind=r_def)    :: detj
  real(kind=r_def)    :: jac(3,3), jac_inv(3,3)
  real(kind=r_def)    :: u_spherical(3), u_cartesian(3), coords(3), llr(3)
  real(kind=r_def)    :: chi_1_cell(ndf_chi), chi_2_cell(ndf_chi), chi_3_cell(ndf_chi)

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Loop through faces
  do df_w2 = 1, ndf_w2

    ! Loop through layers
    do k = 0, nlayers - 1

      ! ---------------------------------------------------------------------- !
      ! Compute Jacobian for transformation
      ! ---------------------------------------------------------------------- !
      ! Create local array for chi field
      do df_chi = 1, ndf_chi
        chi_1_cell(df_chi) = chi_1(map_chi(df_chi) + k)
        chi_2_cell(df_chi) = chi_2(map_chi(df_chi) + k)
        chi_3_cell(df_chi) = chi_3(map_chi(df_chi) + k)
      end do

      ! Compute Jacobian at this W2 point
      call pointwise_coordinate_jacobian(coord_system, geometry, topology,     &
                                         scaled_radius, ndf_chi,               &
                                         chi_1_cell, chi_2_cell, chi_3_cell,   &
                                         ipanel,                               &
                                         basis_chi(:,:,df_w2),                 &
                                         diff_basis_chi(:,:,df_w2),            &
                                         jac, detj)
      ! Calculate inverse Jacobian
      jac_inv = pointwise_coordinate_jacobian_inverse(jac, detj)

      ! ---------------------------------------------------------------------- !
      ! Compute wind in Cartesian components
      ! ---------------------------------------------------------------------- !

      if ( geometry == geometry_spherical ) then
        ! Need to convert from spherical-polar components to Cartesian

        ! First need coordinates of this DoF
        coords(:) = 0.0_r_def
        do df_chi = 1, ndf_chi
          coords(1) = coords(1) + chi_1_cell(df_chi)*basis_chi(1,df_chi,df_w2)
          coords(2) = coords(2) + chi_2_cell(df_chi)*basis_chi(1,df_chi,df_w2)
          coords(3) = coords(3) + chi_3_cell(df_chi)*basis_chi(1,df_chi,df_w2)
        end do

        ! Convert coordinates from whatever coordinate system the model uses
        ! into spherical-polar coordinates
        call chi2llr(coords(1), coords(2), coords(3), &
                     ipanel, llr(1), llr(2), llr(3))

        u_spherical(1) = u_lon(map_w2(df_w2) + k)
        u_spherical(2) = u_lat(map_w2(df_w2) + k)
        u_spherical(3) = u_up(map_w2(df_w2) + k)

        ! Perform conversion of components
        u_cartesian = sphere2cart_vector(u_spherical, llr)

      else

        ! Planar geometry, so values at DoFs are already correcct
        u_cartesian(1) = u_lon(map_w2(df_w2) + k)
        u_cartesian(2) = u_lat(map_w2(df_w2) + k)
        u_cartesian(3) = u_up(map_w2(df_w2) + k)

      end if

      ! ---------------------------------------------------------------------- !
      ! Compute compuational wind
      ! ---------------------------------------------------------------------- !

      ! W2 transformation from computational (reference) values to physical is:
      ! u_physical = matmul(J, u_hdiv) / det(J)
      ! Therefore the inverse of this is
      ! u_hdiv = det(J)*matmul(J_inv, u_physical)
      ! Can pull off individual components of HDiv wind by dotting with W2 basis
      ! function (since these are orthogonal on the reference element)
      u_hdiv(map_w2(df_w2) + k) = detj*dot_product(basis_w2(:,df_w2,df_w2), &
                                                   matmul(jac_inv, u_cartesian))

    end do
  end do

end subroutine convert_phys_to_hdiv_code

end module sci_convert_phys_to_hdiv_kernel_mod
