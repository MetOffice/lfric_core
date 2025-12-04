!-----------------------------------------------------------------------------
! Copyright (c) 2021,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Add a random perturbation to initial potential temperature

module sci_random_perturb_kernel_mod

use argument_mod,               only : arg_type, func_type,            &
                                       GH_FIELD, GH_REAL,              &
                                       GH_READ, GH_WRITE,              &
                                       GH_READWRITE,                   &
                                       ANY_SPACE_9, GH_BASIS,          &
                                       GH_DIFF_BASIS,                  &
                                       CELL_COLUMN, GH_EVALUATOR
use constants_mod,              only : r_def, i_def, rmdi
use kernel_mod,                 only : kernel_type
use fs_continuity_mod,          only : WTHETA, W3

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: random_perturb_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                      &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,       WTHETA)  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass ::random_perturb_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public random_perturb_code
contains

subroutine random_perturb_code(nlayers, theta, height_wtheta, ndf_wtheta, undf_wtheta, map_wtheta)

implicit none

integer(kind=i_def), intent(in) :: nlayers
integer(kind=i_def), intent(in) :: ndf_wtheta
integer(kind=i_def), intent(in) :: undf_wtheta
integer(kind=i_def), intent(in), dimension(ndf_wtheta)  :: map_wtheta
real(kind=r_def), intent(inout), dimension(undf_wtheta) :: theta
real(kind=r_def), intent(in),    dimension(undf_wtheta) :: height_wtheta

integer(kind=i_def) :: k
real(kind=r_def)    :: pert(0:nlayers-1)

! These could become variables in the idealised namelist:
real(kind=r_def), parameter :: theta_pert_start = 5.0e3_r_def ! Height (m) at which perturbations start
real(kind=r_def), parameter :: theta_pert_end   = 7.0e3_r_def ! Height (m) at which perturbations end
real(kind=r_def), parameter :: theta_pert_size  = 0.5_r_def   ! Max size of perturbations (K)

call random_number(pert)

pert(:) = 2.0_r_def * theta_pert_size * ( pert(:) - 0.5_r_def )

do k = 0, nlayers - 1
  if ( height_wtheta(map_wtheta(1) + k) <= theta_pert_end .and.    &
       height_wtheta(map_wtheta(1) + k) >= theta_pert_start ) then
    theta(map_wtheta(1) + k) = theta(map_wtheta(1) + k) + pert(k)
  end if
end do

end subroutine random_perturb_code

end module sci_random_perturb_kernel_mod
