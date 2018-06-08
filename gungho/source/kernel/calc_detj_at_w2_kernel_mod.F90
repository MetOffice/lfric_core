!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the values of detj at W2 locations.
!>
module calc_detj_at_w2_kernel_mod

  use argument_mod,      only : arg_type, func_type,       &
                                GH_FIELD, GH_READ, GH_INC, &
                                GH_DIFF_BASIS,             &
                                CELLS, GH_EVALUATOR
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W0, W2
  use kernel_mod,        only : kernel_type

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: calc_detj_at_w2_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/      &
        arg_type(GH_FIELD,    GH_INC,   W2), &
        arg_type(GH_FIELD*3,  GH_READ,  W0)  &
        /)
    type(func_type) :: meta_funcs(1) = (/ &
        func_type(W0, GH_DIFF_BASIS)      &
        /)
    integer :: iterates_over = CELLS
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass ::calc_detj_at_w2_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! Overload the default structure constructor for function space
  interface calc_detj_at_w2_kernel_type
    module procedure calc_detj_at_w2_kernel_constructor
  end interface

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public calc_detj_at_w2_code

contains

type(calc_detj_at_w2_kernel_type) function calc_detj_at_w2_kernel_constructor() result(self)
  implicit none
  return
end function calc_detj_at_w2_kernel_constructor

!> @param[in]  nlayers        Integer the number of layers
!> @param[out] detj_w2        The output field containing the detj values at W2 locations
!> @param[in]  chi1           The array of coordinates in the first direction
!> @param[in]  chi2           The array of coordinates in the second direction
!> @param[in]  chi3           The array of coordinates in the third direction
!> @param[in]  ndf_w2         The number of degrees of freedom per cell for the output field
!> @param[in]  undf_w2        The number of unique degrees of freedom for the output field
!> @param[in]  map_w2         Integer array holding the dofmap for the cell at the base of the column for the output field
!> @param[in]  nodal_basis_w2 The nodal basis functions evaluated at the nodal points for the W2 field
!> @param[in]  ndf_chi        The number of degrees of freedom per cell for the coordinate field
!> @param[in]  undf_chi       The number of unique degrees of freedom for the coordinate field
!> @param[in]  map_chi        Integer array holding the dofmap for the cell at the base of the column for the coordinate field
!> @param[in]  diff_basis_chi The diff basis functions of the coordinate space evaluated at the nodal points
subroutine calc_detj_at_w2_code( nlayers,                                  &
                                 detj_w2,                                  &
                                 chi1, chi2, chi3,                         &
                                 ndf_w2, undf_w2, map_w2,                  &
                                 ndf_chi, undf_chi, map_chi,               &
                                 diff_basis_chi                            )

  use coordinate_jacobian_mod, only: coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def),                            intent(in)    :: nlayers
  integer(kind=i_def),                            intent(in)    :: ndf_w2
  integer(kind=i_def),                            intent(in)    :: undf_w2
  integer(kind=i_def),                            intent(in)    :: ndf_chi
  integer(kind=i_def),                            intent(in)    :: undf_chi
  real(kind=r_def), dimension(undf_w2),           intent(inout) :: detj_w2
  real(kind=r_def), dimension(undf_chi),          intent(in)    :: chi1, chi2, chi3
  integer(kind=i_def), dimension(ndf_w2),         intent(in)    :: map_w2
  integer(kind=i_def), dimension(ndf_chi),        intent(in)    :: map_chi
  real(kind=r_def), dimension(3,ndf_chi,ndf_w2),  intent(in)    :: diff_basis_chi

  !Internal variables
  integer(kind=i_def)                  :: df, k
  real(kind=r_def), dimension(ndf_chi) :: chi1_e, chi2_e, chi3_e

  real(kind=r_def), dimension(ndf_w2)     :: dj
  real(kind=r_def), dimension(3,3,ndf_w2) :: jacobian

  do k = 0, nlayers-1
    do df = 1,ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi, ndf_w2, 1, chi1_e, chi2_e, chi3_e,  &
                             diff_basis_chi, jacobian, dj)

    do df = 1,ndf_w2
      detj_w2(map_w2(df)+k) = dj(df)
    end do
  end do

end subroutine calc_detj_at_w2_code

end module calc_detj_at_w2_kernel_mod
