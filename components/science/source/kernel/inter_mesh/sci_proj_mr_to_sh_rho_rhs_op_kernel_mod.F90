!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the matrix for taking a mixing ratio from Wtheta to a mass
!> in shifted W3.
!> @details The transform matrix, is the matrix formed doing the integral
!> \f[ \int( \phi_i * \Phi_j \Psi_k dV \f]
!> where \phi_i are the basis functions of the W3 shifted space, \Phi_j are
!> the basis functions of Wtheta and \Psi_k are the basis functions for W3.
!> The end result is a tridiagonal matrix. To transform from mixing ratio in
!> Wtheta to density in shifted W3, we need also need to multiply by dry density.
!> This results in us having four matrix terms -- those on the upper diagonal,
!> lower diagonal and two contributing to the diagonal.
!> This kernel computes these matrix contributions. Each involves integrating
!> over a half-layer.
!> N.B. this will only work under the following conditions:
!> 1. We are at lowest order, so W3 basis functions are constant and Wtheta
!>    functions vary linearly with height.
!> 2. The shifted level is exactly halfway between two levels on the original
!>    mesh.
!> 3. We assume rehabilitation and do not divide by det J for W3 fields.
!>
module sci_proj_mr_to_sh_rho_rhs_op_kernel_mod

  use argument_mod,      only : arg_type, func_type,                    &
                                GH_FIELD, GH_REAL, GH_WRITE, GH_READ,   &
                                ANY_SPACE_9, ANY_DISCONTINUOUS_SPACE_3, &
                                GH_BASIS, GH_DIFF_BASIS,                &
                                CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  !>
  type, public, extends(kernel_type) :: proj_mr_to_sh_rho_rhs_op_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                      &
         arg_type(GH_FIELD*4, GH_REAL, GH_WRITE, W3),                        & ! I_lower/upper
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),               & ! chi_dl
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! panel_id
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta)                     & ! dummy_theta
         /)
    type(func_type) :: meta_funcs(2) = (/                                    &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS),                    &
         func_type(Wtheta,      GH_BASIS)                                    &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: proj_mr_to_sh_rho_rhs_op_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: proj_mr_to_sh_rho_rhs_op_code
contains

!> @brief Compute integrals for matrix transforming from mixing ratio to density.
!! @param[in] nlayers The number of layers of the original mesh.
!! @param[in,out] I_lower_i_ip1 The integral of the (i+1)-th Wtheta basis function on
!! the lower half of the original mesh. Is a W3 field.
!! @param[in,out] I_lower_i_i The integral of the i-th Wtheta basis function on
!! the lower half of the original mesh. Is a W3 field.
!! @param[in,out] I_upper_i_i The integral of the i-th Wtheta basis function on
!! the upper half of the original mesh. Is a W3 field.
!! @param[in,out] I_upper_i_im1 The integral of the (i-1)-th Wtheta basis function on
!! the upper half of the original mesh. Is a W3 field.
!! @param[in] chi_dl_1 The 1st coordinate field in Wchi for double level mesh.
!! @param[in] chi_dl_2 The 2nd coordinate field in Wchi for double level mesh.
!! @param[in] chi_dl_3 The 3rd coordinate field in Wchi for double level mesh.
!! @param[in] panel_id A field giving the ID for the mesh panels.
!! @param[in] dummy_theta An unused dummy variable in Wtheta.
!! @param[in] ndf_w3 The number of degrees of freedom per cell for w3
!! @param[in] undf_w3 The number of unique degrees of freedom for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] ndf_chi_dl The number of degrees of freedom per cell for chi
!! @param[in] undf_chi_dl The number of unique degrees of freedom for chi
!! @param[in] map_chi_dl Dofmap for the cell at the base of the column for chi
!! @param[in] chi_dl_basis 4-dim array for holding the WChi basis functions
!!                         evaluated at quadrature points.
!! @param[in] chi_diff_basis 4-dim array holding differential of the basis
!!                           functions evaluated at gaussian quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] ndf_wtheta The number of degrees of freedom per cell for wtheta
!! @param[in] undf_wtheta The number of unique degrees of freedom for wtheta
!! @param[in] map_wtheta Dofmap for the cell at the base of the column for wtheta
!! @param[in] wtheta_basis 4-dim array holding basis functions evaluated at
!!                         quadrature points
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine proj_mr_to_sh_rho_rhs_op_code(                                      &
                                          nlayers,                             &
                                          I_lower_i_ip1,                       &
                                          I_lower_i_i,                         &
                                          I_upper_i_i,                         &
                                          I_upper_i_im1,                       &
                                          chi_dl_1, chi_dl_2, chi_dl_3,        &
                                          panel_id,                            &
                                          dummy_theta,                         &
                                          ndf_w3, undf_w3, map_w3,             &
                                          ndf_chi_dl, undf_chi_dl, map_chi_dl, &
                                          chi_dl_basis,                        &
                                          chi_dl_diff_basis,                   &
                                          ndf_pid, undf_pid, map_pid,          &
                                          ndf_wtheta, undf_wtheta,             &
                                          map_wtheta, wtheta_basis,            &
                                          nqp_h, nqp_v, wqp_h, wqp_v           &
                                         )

  use sci_coordinate_jacobian_mod, only: coordinate_jacobian

  use finite_element_config_mod, only: coord_system
  use base_mesh_config_mod,      only: geometry, topology
  use planet_config_mod,         only: scaled_radius

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_w3, ndf_wtheta, ndf_chi_dl, ndf_pid
  integer(kind=i_def), intent(in) :: undf_w3, undf_wtheta, undf_chi_dl, undf_pid
  integer(kind=i_def), dimension(ndf_w3),     intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta
  integer(kind=i_def), dimension(ndf_chi_dl), intent(in) :: map_chi_dl
  integer(kind=i_def), dimension(ndf_pid),    intent(in) :: map_pid

  real(kind=r_def), dimension(1,ndf_wtheta,nqp_h,nqp_v), intent(in) :: wtheta_basis
  real(kind=r_def), dimension(3,ndf_chi_dl,nqp_h,nqp_v), intent(in) :: chi_dl_diff_basis
  real(kind=r_def), dimension(1,ndf_chi_dl,nqp_h,nqp_v), intent(in) :: chi_dl_basis

  real(kind=r_def), dimension(undf_w3),     intent(inout) :: I_lower_i_ip1
  real(kind=r_def), dimension(undf_w3),     intent(inout) :: I_lower_i_i
  real(kind=r_def), dimension(undf_w3),     intent(inout) :: I_upper_i_i
  real(kind=r_def), dimension(undf_w3),     intent(inout) :: I_upper_i_im1
  real(kind=r_def), dimension(undf_pid),    intent(in)   :: panel_id
  real(kind=r_def), dimension(undf_wtheta), intent(in)   :: dummy_theta
  real(kind=r_def), dimension(undf_chi_dl), intent(in)   :: chi_dl_1, chi_dl_2, chi_dl_3
  real(kind=r_def), dimension(nqp_h), intent(in)         :: wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)         :: wqp_v

  ! Internal variables
  integer(kind=i_def)                                    :: df, k
  integer(kind=i_def)                                    :: qp1, qp2
  real(kind=r_def), dimension(ndf_chi_dl)                :: lower_chi_1_e, &
                                                            lower_chi_2_e, &
                                                            lower_chi_3_e, &
                                                            upper_chi_1_e, &
                                                            upper_chi_2_e, &
                                                            upper_chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)               :: lower_dj, upper_dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v)           :: lower_jac, upper_jac
  real(kind=r_def), dimension(ndf_w3)                    :: I_lower_i_ip1_e,  &
                                                            I_lower_i_i_e, &
                                                            I_upper_i_i_e, &
                                                            I_upper_i_im1_e

  integer(kind=i_def) :: ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)


  do k = 0, nlayers-1
    ! Extract coordinates for lower and upper half cells
    ! Lower is lower half of orig layer and upper is upper half of orig layer
    do df = 1, ndf_chi_dl
      lower_chi_1_e(df) = chi_dl_1( map_chi_dl(df) + 2*k )
      lower_chi_2_e(df) = chi_dl_2( map_chi_dl(df) + 2*k )
      lower_chi_3_e(df) = chi_dl_3( map_chi_dl(df) + 2*k )
      upper_chi_1_e(df) = chi_dl_1( map_chi_dl(df) + 2*k + 1 )
      upper_chi_2_e(df) = chi_dl_2( map_chi_dl(df) + 2*k + 1 )
      upper_chi_3_e(df) = chi_dl_3( map_chi_dl(df) + 2*k + 1 )
    end do

    ! Get detj for lower and upper half cells
    call coordinate_jacobian(coord_system, geometry, topology, scaled_radius,        &
                             ndf_chi_dl, nqp_h, nqp_v, lower_chi_1_e, lower_chi_2_e, &
                             lower_chi_3_e, ipanel, chi_dl_basis, chi_dl_diff_basis, &
                             lower_jac, lower_dj)
    call coordinate_jacobian(coord_system, geometry, topology, scaled_radius,        &
                             ndf_chi_dl, nqp_h, nqp_v, upper_chi_1_e, upper_chi_2_e, &
                             upper_chi_3_e, ipanel, chi_dl_basis, chi_dl_diff_basis, &
                             upper_jac, upper_dj)

    ! Initialise values to zero
    do df = 1, ndf_w3
      I_lower_i_ip1_e(df) = 0.0_r_def
      I_lower_i_i_e(df) = 0.0_r_def
      I_upper_i_i_e(df) = 0.0_r_def
      I_upper_i_im1_e(df) = 0.0_r_def
    end do

    do qp1 = 1, nqp_h
      do qp2 = 1, nqp_v
        do df = 1, ndf_w3
          ! The quadrature points aren't correct for the wtheta basis, so we use
          ! the linearity of the functions to obtain the correct values
          ! N.B. Assume ndf_wtheta is only 2

          ! Lower half: wtheta_basis(2, 0.5*qp) -> 0.5*wtheta_basis(2, qp)
          I_lower_i_ip1_e(df) = I_lower_i_ip1_e(df) + wqp_h(qp1) * wqp_v(qp2) &
            * 0.5_r_def * wtheta_basis(1,2,qp1,qp2) * lower_dj(qp1,qp2)

          ! Lower half: wtheta_basis(1, 0.5*qp) -> 0.5 + 0.5*wtheta_basis(1, qp)
          I_lower_i_i_e(df) = I_lower_i_i_e(df) + wqp_h(qp1) * wqp_v(qp2) &
            * (0.5_r_def * wtheta_basis(1,1,qp1,qp2) + 0.5_r_def) * lower_dj(qp1,qp2)

          ! Upper half: wtheta_basis(2, 0.5+0.5*qp) -> 0.5 + 0.5*wtheta_basis(2, qp)
          I_upper_i_i_e(df) = I_upper_i_i_e(df) + wqp_h(qp1) * wqp_v(qp2) &
            * (0.5_r_def + 0.5_r_def * wtheta_basis(1,2,qp1,qp2)) * upper_dj(qp1,qp2)

          ! Upper half: wtheta_basis(1, 0.5+0.5*qp) -> 0.5*wtheta_basis(1, qp)
          I_upper_i_im1_e(df) = I_upper_i_im1_e(df) + wqp_h(qp1) * wqp_v(qp2) &
            * 0.5_r_def * wtheta_basis(1,1,qp1,qp2) * upper_dj(qp1,qp2)

        end do
      end do
    end do

    ! Return from reference element to physical element
    do df = 1, ndf_w3
      I_lower_i_ip1( map_w3(df) + k ) = I_lower_i_ip1_e(df)
      I_lower_i_i( map_w3(df) + k )   = I_lower_i_i_e(df)
      I_upper_i_i( map_w3(df) + k )   = I_upper_i_i_e(df)
      I_upper_i_im1( map_w3(df) + k ) = I_upper_i_im1_e(df)
    end do
  end do

end subroutine proj_mr_to_sh_rho_rhs_op_code

end module sci_proj_mr_to_sh_rho_rhs_op_kernel_mod
