!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Compute the values of dA at W2 locations.
!>
module sci_calc_dA_at_w2_kernel_mod

  use argument_mod,          only : arg_type, func_type,       &
                                    GH_FIELD, GH_REAL, GH_INC, &
                                    GH_READ, ANY_SPACE_1,      &
                                    GH_BASIS, GH_DIFF_BASIS,   &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    CELL_COLUMN, GH_EVALUATOR
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : W2
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : N, S, E, W, T, B

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: calc_dA_at_w2_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                    &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                       &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_1),              &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3) &
         /)
    type(func_type) :: meta_funcs(1) = (/                                  &
         func_type(ANY_SPACE_1, GH_BASIS, GH_DIFF_BASIS)                   &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: calc_dA_at_w2_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: calc_dA_at_w2_code

contains

!> @param[in]  nlayers        Integer the number of layers
!> @param[in,out] dA          The output field containing the dA values at W2 locations
!> @param[in]  chi1           The array of coordinates in the first direction
!> @param[in]  chi2           The array of coordinates in the second direction
!> @param[in]  chi3           The array of coordinates in the third direction
!> @param[in]  panel_id       A field giving the ID for mesh panels.
!> @param[in]  ndf_w2         The number of degrees of freedom per cell for the output field
!> @param[in]  undf_w2        The number of unique degrees of freedom for the output field
!> @param[in]  map_w2         Integer array holding the dofmap for the cell at the base of the column for the output field
!> @param[in]  ndf_chi        The number of degrees of freedom per cell for the coordinate field
!> @param[in]  undf_chi       The number of unique degrees of freedom for the coordinate field
!> @param[in]  map_chi        Integer array holding the dofmap for the cell at the base of the column for the coordinate field
!> @param[in]  basis_chi      The basis functions of the coordinate space evaluated at the nodal points
!> @param[in]  diff_basis_chi The diff basis functions of the coordinate space evaluated at the nodal points
!> @param[in]  ndf_pid        Number of degrees of freedom per cell for panel_id
!> @param[in]  undf_pid       Number of unique degrees of freedom for panel_id
!> @param[in]  map_pid        Dofmap for the cell at the base of the column for panel_id
subroutine calc_dA_at_w2_code( nlayers,                                  &
                               dA,                                       &
                               chi1, chi2, chi3,                         &
                               panel_id,                                 &
                               ndf_w2, undf_w2, map_w2,                  &
                               ndf_chi, undf_chi, map_chi,               &
                               basis_chi, diff_basis_chi,                &
                               ndf_pid, undf_pid, map_pid                &
                              )

  use sci_coordinate_jacobian_mod, only: coordinate_jacobian, coordinate_jacobian_inverse

  use finite_element_config_mod, only: coord_system
  use base_mesh_config_mod,      only: geometry, topology
  use planet_config_mod,         only: scaled_radius

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2
  integer(kind=i_def), intent(in)    :: undf_w2
  integer(kind=i_def), intent(in)    :: ndf_chi
  integer(kind=i_def), intent(in)    :: undf_chi
  integer(kind=i_def), intent(in)    :: ndf_pid
  integer(kind=i_def), intent(in)    :: undf_pid
  real(kind=r_def),    intent(inout) :: dA(undf_w2)
  real(kind=r_def),    intent(in)    :: chi1(undf_chi)
  real(kind=r_def),    intent(in)    :: chi2(undf_chi)
  real(kind=r_def),    intent(in)    :: chi3(undf_chi)
  real(kind=r_def), dimension(undf_pid), intent(in)         :: panel_id
  integer(kind=i_def), dimension(ndf_w2), intent(in)        :: map_w2
  integer(kind=i_def), dimension(ndf_chi), intent(in)       :: map_chi
  real(kind=r_def), dimension(3,ndf_chi,ndf_w2), intent(in) :: diff_basis_chi
  real(kind=r_def), dimension(1,ndf_chi,ndf_w2), intent(in) :: basis_chi
  integer(kind=i_def), dimension(ndf_pid), intent(in)       :: map_pid

  ! Internal variables
  integer(kind=i_def) :: df, k
  real(kind=r_def), dimension(ndf_chi)    :: chi1_e, chi2_e, chi3_e

  real(kind=r_def), dimension(ndf_w2)     :: dj
  real(kind=r_def), dimension(3,3,ndf_w2) :: jacobian
  real(kind=r_def), dimension(3,3,ndf_w2) :: jac_inv
  real(kind=r_def), dimension(3,3,ndf_w2) :: JTJinvT

  integer(kind=i_def) :: ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    do df = 1,ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do
    call coordinate_jacobian( coord_system, geometry, topology, scaled_radius, &
                              ndf_chi, ndf_w2, chi1_e, chi2_e, chi3_e,         &
                              ipanel, basis_chi, diff_basis_chi, jacobian, dj)

    call coordinate_jacobian_inverse(ndf_w2, jacobian, dj, jac_inv)

    do df = 1,ndf_w2
      JTJinvT(:,:,df) = matmul(jac_inv(:,:,df),transpose(jac_inv(:,:,df)))
    end do

    ! Face 1 = dx2 * dx3 = dy * dz (E,W)
    ! Face 2 = dx1 * dx3 = dx * dz (N,S)
    ! Face 3 = dx1 * dx2 = dx * dy (T,B)

    ! Lowest order only below, i.e. we only calculate on the face dofs
    dA(map_w2(W)+k) = dj(W)*sqrt(JTJinvT(1,1,W))
    dA(map_w2(E)+k) = dj(E)*sqrt(JTJinvT(1,1,E))
    dA(map_w2(N)+k) = dj(N)*sqrt(JTJinvT(2,2,N))
    dA(map_w2(S)+k) = dj(S)*sqrt(JTJinvT(2,2,S))
    dA(map_w2(T)+k) = dj(T)*sqrt(JTJinvT(3,3,T))
    dA(map_w2(B)+k) = dj(B)*sqrt(JTJinvT(3,3,B))

  end do

end subroutine calc_dA_at_w2_code

end module sci_calc_dA_at_w2_kernel_mod
