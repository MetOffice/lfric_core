!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @brief Calculates the coefficients, a0,a1,a2, for 1D subgrid
!>        representation of rho, rho(x) = a0 + a1*x+a2*x**2 with 0<x<1,
!>        here x simply represents a local coordinate within a cell in either
!>        the chi1, chi2, or chi3 directions.

!> @details The kernel computes the coefficients a0,a1,a2 where rho is represented in 1D
!>         by the approximation rho(x) = a0+a1*x+a2*x**2
!>         Various cases for calculating a0,a1 and a2 are available, including 
!>         constant,linear and quadratic subgrid representations of rho.
!>         For linear representation there are several options. If no slope limiter is
!>         required then centered difference is used to estimate the slope.
!>         Slope limiters which are currently available are minmod and superbee.
!>         These slope limiters are extensively covered in the literature on slope limiters
!>         and have good performance.
!>         For quadratic representation of rho PPM is used and the options of
!>         positivity and monotonicity are available
!>
!>         Note that this kernel only works when rho is a W3 field at lowest order
!>         since it is assumed that ndf_w3 = 1 with stencil_map(1,:) containing
!>         the relevant dofmaps.
module subgrid_coeffs_kernel_mod

use argument_mod,       only : arg_type, func_type,        &
                               GH_FIELD, GH_INC, GH_WRITE, &
                               W3,                         &
                               GH_BASIS,                   &
                               CELLS
use constants_mod,      only : r_def, i_def
use subgrid_config_mod, only : subgrid_rho_approximation_constant_subgrid,     &
                               subgrid_rho_approximation_constant_positive,    &
                               subgrid_rho_approximation_ppm_no_limiter,       &
                               subgrid_rho_approximation_ppm_positive_only,    &
                               subgrid_rho_approximation_ppm_positive_monotone
use kernel_mod,         only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: subgrid_coeffs_kernel_type
  private
  type(arg_type) :: meta_args(1) = (/                                  &
       arg_type(GH_FIELD, GH_WRITE, W3)                                &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, public, nopass :: subgrid_coeffs_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface subgrid_coeffs_kernel_type
   module procedure subgrid_coeffs_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

type(subgrid_coeffs_kernel_type) function subgrid_coeffs_kernel_constructor() result(self)
  return
end function subgrid_coeffs_kernel_constructor

!> @brief Compute the subgrid reconstruction coeffiecients for a density field
!! @param[in] nlayers Number of layers
!! @param[in] subgridrho_option Option for which approximation to use
!! @param[in] ndf_w3 Number of degrees of freedom for W3 per cell
!! @param[in] undf_w3 Number unique of degrees of freedom for W3
!! @param[in] rho Density
!! @param[in] rho Orientation of cell
!! @param[in] stencil_length Local length of a stencil (5 for PPM)
!! @param[in] stencil_map Dofmap for the stencil
!! @param[in] direction Direction of cosmic update
!! @param[out] a0 Coefficient a0
!! @param[out] a1 Coefficient a1
!! @param[out] a2 Coefficient a2

subroutine subgrid_coeffs_code(                                               &
                                nlayers,                                      &
                                subgridrho_option,                            &
                                undf_w3,                                      &
                                rho,                                          &
                                cell_orientation,                             &
                                ndf_w3,                                       &
                                stencil_length,                               &
                                stencil_map,                                  &
                                direction,                                    &
                                a0,                                           &
                                a1,                                           &
                                a2                                            &
                                )

  use subgrid_rho_mod, only: minmod_function, maxmod_function, second_order_coeffs
  use cosmic_flux_mod,    only : stencil_ordering_and_orientation

  ! Arguments
  integer, intent(in)               :: nlayers
  integer, intent(in)               :: subgridrho_option
  integer, intent(in)               :: undf_w3
  real(kind=r_def), intent(in)      :: rho(undf_w3)
  real(kind=r_def), intent(in)      :: cell_orientation(undf_w3)
  integer, intent(in)               :: ndf_w3
  integer, intent(in)               :: stencil_length
  integer, intent(in)               :: stencil_map(1:ndf_w3,1:stencil_length)
  real(kind=r_def), intent(inout)   :: a0(undf_w3)
  real(kind=r_def), intent(inout)   :: a1(undf_w3)
  real(kind=r_def), intent(inout)   :: a2(undf_w3)
  integer, intent(in)               :: direction

  real(kind=r_def)               :: sigma1,sigma2
  real(kind=r_def)               :: coeffs(1:3)
  real(kind=r_def)               :: rho_local(1:stencil_length)

  integer :: k, ii
  integer :: stencil_ordering(1:stencil_length)

  logical :: positive,monotone

  ! If stencil is of length 5 for example, the stencil_ordering array returned is either
  ! (/ 4, 2, 1, 3, 5 /) or (/ 5, 3, 1, 2, 4 /) depending on the cell orientation and
  ! direction (x or y)
  call stencil_ordering_and_orientation(stencil_length,int(cell_orientation(stencil_map(1,1)),i_def),direction,stencil_ordering)


  do k=0,nlayers-1

    ! Rearrange the rho values such that array rho_local is an array of rho values
    ! from left to right with rho_local(1) the rho value for the left most cell,
    ! rho_local(2) the value one cell to the right and so on.
    do ii=1,stencil_length
      rho_local(ii) = rho( stencil_map(1,stencil_ordering(ii))+k )
    end do

    select case(subgridrho_option)
      case (subgrid_rho_approximation_constant_subgrid)
        a0(stencil_map(1,1)+k) = rho(stencil_map(1,1)+k)
        a1(stencil_map(1,1)+k) = 0.0_r_def
        a2(stencil_map(1,1)+k) = 0.0_r_def

      case (subgrid_rho_approximation_constant_positive)
        a0(stencil_map(1,1)+k) = max(rho(stencil_map(1,1)+k),0.0_r_def)
        a1(stencil_map(1,1)+k) = 0.0_r_def
        a2(stencil_map(1,1)+k) = 0.0_r_def

      case (subgrid_rho_approximation_ppm_no_limiter)
        positive=.false.
        monotone=.false.
        call second_order_coeffs(rho_local,coeffs,positive,monotone)
        a0(stencil_map(1,1)+k) = coeffs(1)
        a1(stencil_map(1,1)+k) = coeffs(2)
        a2(stencil_map(1,1)+k) = coeffs(3)

      case (subgrid_rho_approximation_ppm_positive_only, &
            subgrid_rho_approximation_ppm_positive_monotone)
        positive=.true.
        monotone=.false.
        if ( subgridrho_option == subgrid_rho_approximation_ppm_positive_monotone) monotone=.true.
        call second_order_coeffs(rho_local,coeffs,positive,monotone)
        a0(stencil_map(1,1)+k) = coeffs(1)
        a1(stencil_map(1,1)+k) = coeffs(2)
        a2(stencil_map(1,1)+k) = coeffs(3)

    end select

  end do

end subroutine subgrid_coeffs_code

end module subgrid_coeffs_kernel_mod
