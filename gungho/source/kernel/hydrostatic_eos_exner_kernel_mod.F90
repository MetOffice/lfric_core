!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Computes exner from equation of state and hydrostatic balance

module hydrostatic_eos_exner_kernel_mod

use argument_mod,               only : arg_type, func_type,            &
                                       GH_FIELD, GH_READ, GH_WRITE,    &
                                       ANY_SPACE_9, ANY_SPACE_1,       &
                                       GH_BASIS, CELLS, GH_EVALUATOR
use constants_mod,              only : r_def, i_def
use planet_config_mod,          only : gravity, p_zero, kappa, rd, cp
use fs_continuity_mod,          only : WTHETA, W3
use kernel_mod,                 only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: hydrostatic_eos_exner_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                       &
       arg_type(GH_FIELD,   GH_WRITE, W3),                  &
       arg_type(GH_FIELD,   GH_READ,  W3),                  &
       arg_type(GH_FIELD,   GH_READ,  WTHETA),              &
       arg_type(GH_FIELD*3, GH_READ,  WTHETA),              &
       arg_type(GH_FIELD,   GH_READ,  W3)                   &
       /)
  type(func_type) :: meta_funcs(2) = (/                     &
       func_type(W3, GH_BASIS),                             &
       func_type(WTHETA, GH_BASIS)                          &
       /)
       integer :: iterates_over = CELLS
       integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass ::hydrostatic_eos_exner_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface hydrostatic_eos_exner_kernel_type
   module procedure hydrostatic_eos_exner_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public hydrostatic_eos_exner_code
contains

type(hydrostatic_eos_exner_kernel_type) &
   function hydrostatic_eos_exner_kernel_constructor() result(self)
implicit none
return
end function hydrostatic_eos_exner_kernel_constructor

!> @brief Computes density from equation of state
!! @param[in] nlayers Number of layers
!! @param[out] exner Exner pressure field
!! @param[in] rho Density field
!! @param[in] theta Potential temperature field
!! @param[in] height_w3 Height coordinate in w3
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number unique of degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] ndf_wt Number of degrees of freedom per cell for wtheta
!! @param[in] undf_wt Number unique of degrees of freedom  for wtheta
!! @param[in] map_wt Dofmap for the cell at the base of the column for wt
subroutine hydrostatic_eos_exner_code(nlayers, exner, rho, theta, &
                                 moist_dyn_gas, moist_dyn_tot, moist_dyn_fac, &
                                 height_w3,    &
                                 ndf_w3, undf_w3, map_w3, basis_w3, ndf_wt, &
                                 undf_wt, map_wt, basis_wt)
  
  use analytic_temperature_profiles_mod, only : analytic_temperature

  implicit none
  
  !Arguments
  integer, intent(in) :: nlayers, ndf_w3, undf_w3,  ndf_wt, undf_wt
  integer, dimension(ndf_w3), intent(in)  :: map_w3
  integer, dimension(ndf_wt), intent(in) :: map_wt

  real(kind=r_def), dimension(undf_w3),  intent(inout)       :: exner
  real(kind=r_def), dimension(undf_w3),  intent(in)          :: rho, height_w3
  real(kind=r_def), dimension(undf_wt),  intent(in)          :: moist_dyn_gas, moist_dyn_tot, moist_dyn_fac
  real(kind=r_def), dimension(undf_wt),  intent(in)          :: theta
  real(kind=r_def), dimension(1,ndf_w3,ndf_w3),  intent(in)  :: basis_w3
  real(kind=r_def), dimension(1,ndf_wt,ndf_w3),  intent(in)  :: basis_wt

  !Internal variables
  integer(kind=i_def)                  :: k, df, dft, df3
  real(kind=r_def), dimension(ndf_w3)  :: rho_e
  real(kind=r_def), dimension(ndf_wt)  :: theta_moist_e
  real(kind=r_def)                     :: rho_cell, theta_moist, dz

  !Compute exner from eqn of state in lowest level
  k = 0

  do df3 = 1, ndf_w3
    rho_e(df3) = rho( map_w3(df3) + k)
  end do

  do dft = 1, ndf_wt
    theta_moist_e(dft) = theta( map_wt(dft) + k) * moist_dyn_gas(map_wt(dft)+k)
  end do
    
  do df = 1, ndf_w3
    
    rho_cell = 0.0_r_def 
    do df3 = 1, ndf_w3
      rho_cell = rho_cell + rho_e(df3)*basis_w3(1,df3,df)
    end do
    
    theta_moist = 0.0_r_def
    do dft = 1, ndf_wt
      theta_moist = theta_moist + theta_moist_e(dft)*basis_wt(1,dft,df)
    end do

    exner(map_w3(df)+k) = (rd*rho_cell*theta_moist/p_zero)**(kappa/(1.0_r_def-kappa))
  end do
  
  ! Exner on other levels from hydrostatic balance
  do k = 1, nlayers-1
    dz = height_w3(map_w3(1)+k)-height_w3(map_w3(1)+k-1)
    theta_moist = moist_dyn_gas(map_wt(1)+k) * theta(map_wt(1)+k) /   &
                  moist_dyn_tot(map_wt(1)+k)
    exner(map_w3(1)+k) = exner(map_w3(1)+k-1) - gravity * dz / (cp * theta_moist)
  end do

end subroutine hydrostatic_eos_exner_code

end module hydrostatic_eos_exner_kernel_mod
