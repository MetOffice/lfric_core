!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Compute an operator to project a scalar field (treated as a component of a vector)
!!        from  the scalar space Ws (= W3 or Wtheta) into a W2 field.
module sci_project_ws_to_w2_operator_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,       &
                                    GH_OPERATOR,               &
                                    GH_FIELD, GH_SCALAR,       &
                                    GH_REAL, GH_INTEGER,       &
                                    GH_READ, GH_WRITE,         &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    GH_BASIS, GH_DIFF_BASIS,   &
                                    CELL_COLUMN, GH_QUADRATURE_XYoZ
use constants_mod,           only : r_def, i_def
use fs_continuity_mod,       only : W2, Wchi

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: project_ws_to_w2_operator_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                               &
       arg_type(GH_OPERATOR, GH_REAL,    GH_WRITE,  W2, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD*3,  GH_REAL,    GH_READ,   Wchi),                          &
       arg_type(GH_FIELD,    GH_REAL,    GH_READ,   ANY_DISCONTINUOUS_SPACE_3),     &
       arg_type(GH_SCALAR,   GH_INTEGER, GH_READ)                                   &
       /)
  type(func_type) :: meta_funcs(3) = (/                              &
       func_type(W2,                        GH_BASIS),               &
       func_type(ANY_DISCONTINUOUS_SPACE_1, GH_BASIS),               &
       func_type(Wchi,                      GH_BASIS, GH_DIFF_BASIS) &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: project_ws_to_w2_operator_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: project_ws_to_w2_operator_code

contains

!> @brief Compute an operator to project a scalar field into a W2 field.
!> @details Compute an operator to project a scalar field
!!          from  the scalar space Ws (= W3 or Wtheta) into a W2 field.
!!          The component of the vector to be projected is given by the
!!          direction field.
!> @param[in]     cell                Horizontal cell index
!> @param[in]     nlayers             Number of layers
!> @param[in]     ncell               Total number of cells in the 3D mesh
!> @param[in,out] projection_operator Projection operator to compute
!> @param[in]     chi1                1st coordinate field in Wchi
!> @param[in]     chi2                2nd coordinate field in Wchi
!> @param[in]     chi3                3rd coordinate field in Wchi
!> @param[in]     panel_id            Field giving the ID for mesh panels.
!> @param[in]     direction           Index of the vector component (1,2 or 3) to project
!> @param[in]     ndf_w2              Number of degrees of freedom per cell for vector space
!> @param[in]     basis_w2            Basis functions for the vector space at quadrature points
!> @param[in]     ndf_ws              Number of degrees of freedom per cell for scalar space
!> @param[in]     basis_ws            Basis functions for the scalar space at quadrature points
!> @param[in]     ndf_wx              Number of degrees of freedom per cell for the coordinate space
!> @param[in]     undf_wx             Number of unique degrees of freedom for the coordinate space
!> @param[in]     map_wx              Dofmap for the cell at the base of the column for the coordinate space
!> @param[in]     basis_wx            Basis functions for the coordinate space at quadrature points
!> @param[in]     diff_basis_wx       Differential basis functions for the coordinate space at quadrature points
!> @param[in]     ndf_pid             Number of degrees of freedom per cell for panel_id
!> @param[in]     undf_pid            Number of unique degrees of freedom for panel_id
!> @param[in]     map_pid             Dofmap for the cell at the base of the column for panel_id
!> @param[in]     nqp_h               Number of horizontal quadrature points
!> @param[in]     nqp_v               Number of vertical quadrature points
!> @param[in]     wqp_h               Weights of horizontal quadrature points
!> @param[in]     wqp_v               Weights of vertical quadrature points
subroutine project_ws_to_w2_operator_code( cell, nlayers,              &
                                           ncell,                      &
                                           projection_operator,        &
                                           chi1, chi2, chi3,           &
                                           panel_id,                   &
                                           direction,                  &
                                           ndf_w2, basis_w2,           &
                                           ndf_ws, basis_ws,           &
                                           ndf_wx, undf_wx, map_wx,    &
                                           basis_wx, diff_basis_wx,    &
                                           ndf_pid, undf_pid, map_pid, &
                                           nqp_h, nqp_v, wqp_h, wqp_v )

  use sci_coordinate_jacobian_mod, only: pointwise_coordinate_jacobian

  use base_mesh_config_mod,      only: geometry, topology
  use finite_element_config_mod, only: coord_system
  use planet_config_mod,         only: scaled_radius

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in) :: cell, ncell, nlayers, nqp_h, nqp_v
  integer(kind=i_def),                     intent(in) :: ndf_w2, ndf_ws, ndf_wx, ndf_pid
  integer(kind=i_def),                     intent(in) :: undf_wx, undf_pid
  integer(kind=i_def), dimension(ndf_wx),  intent(in) :: map_wx
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in) :: basis_w2
  real(kind=r_def), dimension(1,ndf_ws,nqp_h,nqp_v), intent(in) :: basis_ws
  real(kind=r_def), dimension(1,ndf_wx,nqp_h,nqp_v), intent(in) :: basis_wx
  real(kind=r_def), dimension(3,ndf_wx,nqp_h,nqp_v), intent(in) :: diff_basis_wx

  real(kind=r_def), dimension(ncell, ndf_w2, ndf_ws),  intent(inout) :: projection_operator
  real(kind=r_def), dimension(undf_wx),                intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_pid),               intent(in)    :: panel_id

  integer(kind=i_def), intent(in) :: direction

  real(kind=r_def), intent(in) :: wqp_h(nqp_h)
  real(kind=r_def), intent(in) :: wqp_v(nqp_v)

  ! Internal variables
  integer(kind=i_def)                 :: df_x, df_s, df_2, ik, k, qp_h, qp_v
  real(kind=r_def), dimension(ndf_wx) :: chi1_e, chi2_e, chi3_e
  real(kind=r_def), dimension(3,3)    :: jac
  real(kind=r_def)                    :: detj
  real(kind=r_def), dimension(3)      :: a3d, v

  integer(kind=i_def) :: ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    do df_x = 1, ndf_wx
      chi1_e(df_x) = chi1(map_wx(df_x) + k)
      chi2_e(df_x) = chi2(map_wx(df_x) + k)
      chi3_e(df_x) = chi3(map_wx(df_x) + k)
    end do

    ik = 1 + k + (cell-1)*nlayers
    projection_operator(ik,:,:) = 0.0_r_def
    do qp_v = 1,nqp_v
      do qp_h = 1,nqp_h
        call pointwise_coordinate_jacobian(coord_system, geometry,          &
                                           topology, scaled_radius,         &
                                           ndf_wx, chi1_e, chi2_e, chi3_e,  &
                                           ipanel, basis_wx(:,:,qp_h,qp_v), &
                                           diff_basis_wx(:,:,qp_h,qp_v),    &
                                           jac, detj)
        a3d = 0.0_r_def
        do df_s = 1,ndf_ws
          a3d(direction) = basis_ws(1,df_s,qp_h,qp_v)

          do df_2 = 1,ndf_w2
            ! Projection from ws*3 to w2 = jac*v . a3d
            v = matmul(jac,basis_w2(:,df_2,qp_h,qp_v))

            projection_operator(ik, df_2, df_s) = projection_operator(ik, df_2, df_s) &
                                                + wqp_h(qp_h)*wqp_v(qp_v)*dot_product(v,a3d)
          end do
        end do
      end do
    end do
  end do

end subroutine project_ws_to_w2_operator_code

end module sci_project_ws_to_w2_operator_kernel_mod
