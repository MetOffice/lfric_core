!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Compute the Coriolis operator to apply the rotation vector Omega to
!>        the wind fields.
!> @details The form of Coriolis operator is:
!> \f[ <v,2\Omega \times v> \f] where v is the test/trial function for
!> the velocity space and Omega is the rotation vector of the domain.
!> The Coriolis terms will then be applied by multiplying this operator
!> by a wind field.
!>
!
module compute_coriolis_matrix_kernel_mod

use constants_mod,           only: i_def, r_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, func_type,             &
                                   GH_OPERATOR, GH_FIELD,           &
                                   GH_READ, GH_WRITE,               &
                                   ANY_SPACE_9,                     &
                                   GH_BASIS, GH_DIFF_BASIS,         &
                                   CELLS, GH_QUADRATURE_XYoZ,       &
                                   ANY_DISCONTINUOUS_SPACE_3
use fs_continuity_mod,       only: W2

use coordinate_jacobian_mod, only: coordinate_jacobian
use base_mesh_config_mod,    only: geometry,           &
                                   geometry_spherical, &
                                   f_lat
use rotation_vector_mod,     only: rotation_vector_fplane,  &
                                   rotation_vector_sphere
use planet_config_mod,       only: scaled_omega
use cross_product_mod,       only: cross_product

implicit none
private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: compute_coriolis_matrix_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                &
      arg_type(GH_OPERATOR, GH_WRITE, W2, W2),                       &
      arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_9),                  &
      arg_type(GH_FIELD,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3)     &
       /)
  type(func_type) :: meta_funcs(2) = (/                              &
      func_type(W2, GH_BASIS),                                       &
      func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: compute_coriolis_matrix_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public compute_coriolis_matrix_code
contains

!> @brief Compute the Coriolis operator to apply the rotation vector Omega to
!>        the wind fields.
!!
!! @param[in] cell     Identifying number of cell.
!! @param[in] nlayers  Number of layers.
!! @param[in] ncell_3d ncell*ndf
!! @param[inout] mm    Local stencil or Coriolis operator.
!! @param[in] chi_sph_1 1st coordinate in spherical Wchi
!! @param[in] chi_sph_2 2nd coordinate in spherical Wchi
!! @param[in] chi_sph_3 3rd coordinate in spherical Wchi
!! @param[in] panel_id A field giving the ID for mesh panels.
!! @param[in] ndf      Degrees of freedom per cell.
!! @param[in] basis    Vector basis functions evaluated at quadrature points.
!! @param[in] ndf_chi_sph Number of degrees of freedom per cell for chi
!! @param[in] undf_chi_sph Number of unique degrees of freedom  for chi
!! @param[in] map_chi_sph Dofmap for the cell at the base of the column for chi
!! @param[in] chi_sph_basis Basis functions for spherical Wchi evaluated at
!!                          gaussian quadrature points
!! @param[in] chi_sph_diff_basis Differential of the spherical Wchi basis functions
!!                               evaluated at gaussian quadrature point
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h    Number of horizontal quadrature points.
!! @param[in] nqp_v    Number of vertical quadrature points.
!! @param[in] wqp_h    Horizontal quadrature weights.
!! @param[in] wqp_v    Vertical quadrature weights.
subroutine compute_coriolis_matrix_code(cell, nlayers, ncell_3d,           &
                                        mm,                                &
                                        chi_sph_1, chi_sph_2, chi_sph_3,   &
                                        panel_id,                          &
                                        ndf, basis,                        &
                                        ndf_chi_sph, undf_chi_sph,         &
                                        map_chi_sph,                       &
                                        basis_chi_sph, diff_basis_chi_sph, &
                                        ndf_pid, undf_pid, map_pid,        &
                                        nqp_h, nqp_v, wqp_h, wqp_v)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: cell, ndf, ndf_pid, ndf_chi_sph
  integer(kind=i_def), intent(in)    :: undf_pid, undf_chi_sph
  integer(kind=i_def), intent(in)    :: nqp_h, nqp_v
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ncell_3d

  real(kind=r_def), dimension(3,ndf,nqp_h,nqp_v), intent(in) :: basis

  integer(kind=i_def), intent(in)    :: map_chi_sph(ndf_chi_sph)
  integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)
  real(kind=r_def),    intent(inout) :: mm(ndf,ndf,ncell_3d)
  real(kind=r_def),    intent(in)    :: basis_chi_sph(1,ndf_chi_sph,nqp_h,nqp_v)
  real(kind=r_def),    intent(in)    :: diff_basis_chi_sph(3,ndf_chi_sph,nqp_h,nqp_v)
  real(kind=r_def),    intent(in)    :: chi_sph_1(undf_chi_sph)
  real(kind=r_def),    intent(in)    :: chi_sph_2(undf_chi_sph)
  real(kind=r_def),    intent(in)    :: chi_sph_3(undf_chi_sph)
  real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
  real(kind=r_def),    intent(in)    :: wqp_h(nqp_h)
  real(kind=r_def),    intent(in)    :: wqp_v(nqp_v)

  ! Internal variables
  integer(kind=i_def)                          :: df, df2, k, ik
  integer(kind=i_def)                          :: qp1, qp2

  real(kind=r_def), dimension(ndf_chi_sph)     :: chi_sph_1_e, chi_sph_2_e, chi_sph_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(3,nqp_h,nqp_v)   :: rotation_vector
  real(kind=r_def), dimension(3)               :: omega_cross_u
  real(kind=r_def), dimension(3)               :: jac_u, jac_v

  integer(kind=i_def) :: ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Loop over layers: Start from 1 as in this loop k is not an offset
  do k = 1, nlayers

     ! Indirect the chi coord field here
     do df = 1, ndf_chi_sph
        chi_sph_1_e(df) = chi_sph_1(map_chi_sph(df) + k - 1)
        chi_sph_2_e(df) = chi_sph_2(map_chi_sph(df) + k - 1)
        chi_sph_3_e(df) = chi_sph_3(map_chi_sph(df) + k - 1)
     end do

    ! Calculate rotation vector Omega = (0, 2*cos(lat), 2*sin(lat)) and Jacobian
    if ( geometry == geometry_spherical ) then
      call rotation_vector_sphere(ndf_chi_sph, nqp_h, nqp_v, chi_sph_1_e, chi_sph_2_e,       &
                                  chi_sph_3_e, ipanel, basis_chi_sph, rotation_vector)
    else
      call rotation_vector_fplane(nqp_h, nqp_v, scaled_omega, f_lat,           &
                                  rotation_vector)
    end if

    call coordinate_jacobian(ndf_chi_sph, nqp_h, nqp_v,                     &
                             chi_sph_1_e, chi_sph_2_e, chi_sph_3_e, ipanel, &
                             basis_chi_sph, diff_basis_chi_sph, jac, dj )

    ik = k + (cell-1)*nlayers
    mm(:,:,ik) = 0.0_r_def
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        do df2 = 1, ndf
          jac_u = matmul(jac(:,:,qp1,qp2),basis(:,df2,qp1,qp2))
          omega_cross_u = wqp_h(qp1)*wqp_v(qp2)                                &
                        *cross_product(rotation_vector(:,qp1,qp2),jac_u)       &
                        /dj(qp1,qp2)
          do df = 1, ndf
             jac_v = matmul(jac(:,:,qp1,qp2),basis(:,df,qp1,qp2))
             mm(df,df2,ik) = mm(df,df2,ik) - dot_product(jac_v,omega_cross_u)
          end do
        end do
      end do
    end do
  end do ! end of k loop

end subroutine compute_coriolis_matrix_code

end module compute_coriolis_matrix_kernel_mod
