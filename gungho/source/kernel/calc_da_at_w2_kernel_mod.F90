!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Compute the values of dA at W2 locations.
!>
module calc_dA_at_w2_kernel_mod

  use argument_mod,          only : arg_type, func_type,        &
                                    GH_FIELD, GH_READ, GH_INC,  &
                                    ANY_SPACE_1, GH_DIFF_BASIS, &
                                    CELLS, GH_EVALUATOR
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : W2
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : N, S, E, W, T, B

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: calc_dA_at_w2_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/              &
        arg_type(GH_FIELD,    GH_INC,   W2),         &
        arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_1) &
        /)
    type(func_type) :: meta_funcs(1) = (/     &
        func_type(ANY_SPACE_1, GH_DIFF_BASIS) &
        /)
    integer(i_def) :: iterates_over = CELLS
    integer(i_def) :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass ::calc_dA_at_w2_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! Overload the default structure constructor for function space
  interface calc_dA_at_w2_kernel_type
    module procedure calc_dA_at_w2_kernel_constructor
  end interface

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public calc_dA_at_w2_code

contains

type(calc_dA_at_w2_kernel_type) function calc_dA_at_w2_kernel_constructor() result(self)
  implicit none
  return
end function calc_dA_at_w2_kernel_constructor

!> @param[in]  nlayers        Integer the number of layers
!> @param[out] dA             The output field containing the dA values at W2 locations
!> @param[in]  chi1           The array of coordinates in the first direction
!> @param[in]  chi2           The array of coordinates in the second direction
!> @param[in]  chi3           The array of coordinates in the third direction
!> @param[in]  ndf_w2         The number of degrees of freedom per cell for the output field
!> @param[in]  undf_w2        The number of unique degrees of freedom for the output field
!> @param[in]  map_w2         Integer array holding the dofmap for the cell at the base of the column for the output field
!> @param[in]  ndf_chi        The number of degrees of freedom per cell for the coordinate field
!> @param[in]  undf_chi       The number of unique degrees of freedom for the coordinate field
!> @param[in]  map_chi        Integer array holding the dofmap for the cell at the base of the column for the coordinate field
!> @param[in]  diff_basis_chi The diff basis functions of the coordinate space evaluated at the nodal points
subroutine calc_dA_at_w2_code( nlayers,                                  &
                               dA,                                       &
                               chi1, chi2, chi3,                         &
                               ndf_w2, undf_w2, map_w2,                  &
                               ndf_chi, undf_chi, map_chi,               &
                               diff_basis_chi                            )

  use coordinate_jacobian_mod, only: coordinate_jacobian, coordinate_jacobian_inverse

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2
  integer(kind=i_def), intent(in)    :: undf_w2
  integer(kind=i_def), intent(in)    :: ndf_chi
  integer(kind=i_def), intent(in)    :: undf_chi
  real(kind=r_def), dimension(undf_w2), intent(inout)       :: dA
  real(kind=r_def), dimension(undf_chi), intent(in)         :: chi1, chi2, chi3
  integer(kind=i_def), dimension(ndf_w2), intent(in)        :: map_w2
  integer(kind=i_def), dimension(ndf_chi), intent(in)       :: map_chi
  real(kind=r_def), dimension(3,ndf_chi,ndf_w2), intent(in) :: diff_basis_chi

  ! Internal variables
  integer(kind=i_def) :: df, k
  real(kind=r_def), dimension(ndf_chi)    :: chi1_e, chi2_e, chi3_e

  real(kind=r_def), dimension(ndf_w2)     :: dj
  real(kind=r_def), dimension(3,3,ndf_w2) :: jacobian
  real(kind=r_def), dimension(3,3,ndf_w2) :: jac_inv
  real(kind=r_def), dimension(3,3,ndf_w2) :: JTJinvT

  do k = 0, nlayers-1
    do df = 1,ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi, ndf_w2, 1, chi1_e, chi2_e, chi3_e,  &
                             diff_basis_chi, jacobian, dj)

    call coordinate_jacobian_inverse(ndf_w2, 1, jacobian, dj, jac_inv)

    do df = 1,ndf_w2
      JTJinvT(:,:,df) = matmul(jac_inv(:,:,df),transpose(jac_inv(:,:,df)))
    end do
      
    ! Lowest order only below, i.e. we only calculate on the face dofs
    dA(map_w2(N)+k) = dj(N)*sqrt(JTJinvT(1,1,N))
    dA(map_w2(S)+k) = dj(S)*sqrt(JTJinvT(1,1,S))
    dA(map_w2(W)+k) = dj(W)*sqrt(JTJinvT(2,2,W))
    dA(map_w2(E)+k) = dj(E)*sqrt(JTJinvT(2,2,E))
    dA(map_w2(T)+k) = dj(T)*sqrt(JTJinvT(3,3,T))
    dA(map_w2(B)+k) = dj(B)*sqrt(JTJinvT(3,3,B))

  end do

end subroutine calc_dA_at_w2_code

end module calc_dA_at_w2_kernel_mod
