!-----------------------------------------------------------------------------
! Copyright (c) 2021,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Add a random perturbation to initial potential temperature

module sci_random_perturb_kernel_mod

use, intrinsic :: iso_fortran_env, only: real32, real64

use argument_mod,               only : arg_type, func_type,            &
                                       GH_FIELD, GH_REAL, GH_SCALAR,   &
                                       GH_READ, GH_WRITE,              &
                                       GH_READWRITE,                   &
                                       ANY_SPACE_9, GH_BASIS,          &
                                       GH_DIFF_BASIS,                  &
                                       CELL_COLUMN, GH_EVALUATOR
use constants_mod,              only : i_def, rmdi
use kernel_mod,                 only : kernel_type
use fs_continuity_mod,          only : WTHETA, W3

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: random_perturb_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                      &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,       WTHETA), &
       arg_type(GH_SCALAR, GH_REAL, GH_READ ),             &
       arg_type(GH_SCALAR, GH_REAL, GH_READ ),             &
       arg_type(GH_SCALAR, GH_REAL, GH_READ )              &
       /)
  integer :: operates_on = CELL_COLUMN
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public random_perturb_code

! Generic interface for real32 and real64 types
interface random_perturb_code
  module procedure  &
    random_perturb_code_real32, &
    random_perturb_code_real64
end interface

contains

! ==================
! REAL32 PRECISION
! ==================
subroutine random_perturb_code_real32(  nlayers, theta, height_wtheta,     &
                        theta_pert_size, theta_pert_start, theta_pert_end, &
                        ndf_wtheta, undf_wtheta, map_wtheta)

implicit none

integer(kind=i_def), intent(in) :: nlayers
integer(kind=i_def), intent(in) :: ndf_wtheta
integer(kind=i_def), intent(in) :: undf_wtheta
integer(kind=i_def), intent(in),    dimension(ndf_wtheta)  :: map_wtheta
real(kind=real32),   intent(inout), dimension(undf_wtheta) :: theta
real(kind=real32),   intent(in),    dimension(undf_wtheta) :: height_wtheta
real(kind=real32),   intent(in) :: theta_pert_size
real(kind=real32),   intent(in) :: theta_pert_start
real(kind=real32),   intent(in) :: theta_pert_end

integer(kind=i_def) :: k
real(kind=real32)    :: pert(0:nlayers-1)

call random_number(pert)

pert(:) = 2.0_real32 * theta_pert_size * ( pert(:) - 0.5_real32 )

do k = 0, nlayers - 1
  if ( height_wtheta(map_wtheta(1) + k) <= theta_pert_end .and.    &
       height_wtheta(map_wtheta(1) + k) >= theta_pert_start ) then
    theta(map_wtheta(1) + k) = theta(map_wtheta(1) + k) + pert(k)
  end if
end do

end subroutine random_perturb_code_real32

! ==================
! REAL64 PRECISION
! ==================
subroutine random_perturb_code_real64(  nlayers, theta, height_wtheta,     &
                        theta_pert_size, theta_pert_start, theta_pert_end, &
                        ndf_wtheta, undf_wtheta, map_wtheta)

implicit none

integer(kind=i_def), intent(in) :: nlayers
integer(kind=i_def), intent(in) :: ndf_wtheta
integer(kind=i_def), intent(in) :: undf_wtheta
integer(kind=i_def), intent(in),    dimension(ndf_wtheta)  :: map_wtheta
real(kind=real64),   intent(inout), dimension(undf_wtheta) :: theta
real(kind=real64),   intent(in),    dimension(undf_wtheta) :: height_wtheta
real(kind=real64),   intent(in) :: theta_pert_size
real(kind=real64),   intent(in) :: theta_pert_start
real(kind=real64),   intent(in) :: theta_pert_end

integer(kind=i_def) :: k
real(kind=real64)    :: pert(0:nlayers-1)

call random_number(pert)

pert(:) = 2.0_real64 * theta_pert_size * ( pert(:) - 0.5_real64 )

do k = 0, nlayers - 1
  if ( height_wtheta(map_wtheta(1) + k) <= theta_pert_end .and.    &
       height_wtheta(map_wtheta(1) + k) >= theta_pert_start ) then
    theta(map_wtheta(1) + k) = theta(map_wtheta(1) + k) + pert(k)
  end if
end do

end subroutine random_perturb_code_real64

end module sci_random_perturb_kernel_mod
