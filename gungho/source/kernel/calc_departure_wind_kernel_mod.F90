!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel to compute the wind used for calculating the departure points
!! in the dimensionally split advection scheme. Input is the Piola wind and
!! output is the departure wind. The u_departure_wind variable is being used to
!! store the departure winds and we wish to have a positive u_departure_wind
!! representing a postive physical wind. Since u_departure_wind should not be
!! multiplied by the basis functions themselves in order to evaluate them, we
!! remove any dependency on the direction of the basis functions by multiplying
!! the u_piola wind values by the nodal basis functions.

module calc_departure_wind_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,                     &
                                    GH_FIELD, GH_READ, GH_WRITE,             &
                                    ANY_SPACE_9, W2,                         &
                                    GH_DIFF_BASIS, GH_BASIS,                 &
                                    CELLS, EVALUATOR, GH_EVALUATOR
use constants_mod,           only : r_def

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: calc_departure_wind_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD,    GH_WRITE, W2),                            &
       arg_type(GH_FIELD,    GH_READ,  W2),                            &
       arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_9)                    &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(W2, GH_BASIS),                                        &
       func_type(ANY_SPACE_9, GH_DIFF_BASIS)                           &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_EVALUATOR
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = EVALUATOR
contains
  procedure, nopass ::calc_departure_wind_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface calc_departure_wind_kernel_type
   module procedure calc_departure_wind_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public calc_departure_wind_code
contains

type(calc_departure_wind_kernel_type) function calc_departure_wind_kernel_constructor() result(self)
  return
end function calc_departure_wind_kernel_constructor

!> @param[in] nlayers Number of layers
!> @param[in] ndf Number of degrees of freedom per cell for the output field
!> @param[in] undf Number of unique degrees of freedom for the output field
!> @param[in] map Dofmap for the cell at the base of the column for the output field
!> @param[in] nodal_basis_u Basis functions evaluated at the nodal points for the W2 field
!> @param[inout] u_departure_wind Output field containing the departure wind used to calculate departure points
!> @param[in] u_piola Field for the Piola wind
!> @param[in] chi1 Coordinates in the first direction
!> @param[in] chi2 Coordinates in the second direction
!> @param[in] chi3 Coordinates in the third direction
!> @param[in] ndf_chi Number of degrees of freedom per cell for the coordinate field
!> @param[in] undf_chi Number of unique degrees of freedom for the coordinate field
!> @param[in] map_chi Dofmap for the cell at the base of the column for the coordinate field
!> @param[in] diff_basis_chi Differential basis functions of the coordinate space evaluated at the nodal points
subroutine calc_departure_wind_code(nlayers,                                  &
                                    u_departure_wind,                         &
                                    u_piola,                                  &
                                    chi1, chi2, chi3,                         &
                                    ndf, undf, map, nodal_basis_u,            &
                                    ndf_chi, undf_chi, map_chi,               &
                                    diff_basis_chi                            &
                                    )
  use coordinate_jacobian_mod, only: coordinate_jacobian
  implicit none
  !Arguments
  integer,                                    intent(in)    :: nlayers
  integer,                                    intent(in)    :: ndf, undf, &
                                                               ndf_chi, undf_chi
  integer,          dimension(ndf),           intent(in)    :: map
  real(kind=r_def), dimension(3,ndf,ndf),     intent(in)    :: nodal_basis_u
  integer,          dimension(ndf_chi),       intent(in)    :: map_chi
  real(kind=r_def), dimension(undf),          intent(in)    :: u_piola
  real(kind=r_def), dimension(undf_chi),      intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf),          intent(inout) :: u_departure_wind
  real(kind=r_def), dimension(3,ndf_chi,ndf), intent(in)    :: diff_basis_chi

  !Internal variables
  integer          :: df, k
  real(kind=r_def) :: jacobian(3,3,ndf,1), dj(ndf,1)
  real(kind=r_def), dimension(ndf_chi) :: chi1_e, chi2_e, chi3_e

  do k = 0, nlayers-1
    do df = 1,ndf_chi
      chi1_e(df) = chi1(map_chi(df) + k)
      chi2_e(df) = chi2(map_chi(df) + k)
      chi3_e(df) = chi3(map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi, ndf, 1, chi1_e, chi2_e, chi3_e,  &
                             diff_basis_chi, jacobian, dj)
    do df = 1,ndf
      u_departure_wind(map(df)+k) =                                           &
          dot_product(nodal_basis_u(:,df,df),abs(nodal_basis_u(:,df,df)))*    &
          u_piola(map(df)+k)/dj(df,1)
    end do
  end do

end subroutine calc_departure_wind_code

end module calc_departure_wind_kernel_mod
