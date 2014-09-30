!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Module for computing a linear hydrostatially balanced reference state
module reference_profile_mod
use constants_mod, only: r_def, n_sq, gravity, cp, rd, kappa, p_zero

implicit none

contains
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> Subroutine Computes the reference profile for a single element
!! @param[in] ndf_w0     Integer. The size of the w0 arrays
!! @param[in] ndf_w3     Integer. The size of the w3 arrays
!! @param[in] exner_s    Real 1-dim array. Holds the exner reference profile
!! @param[in] rho_s      Real 1-dim array. Holds the rho reference profile
!! @param[in] theta_s    Real 1-dim array. Holds the theta reference profile
!! @param[in] z_w0       Real 1-dim array. Holds the z coordinate field for w0
!! @param[in] z_w3       Real 1-dim array. Holds the z coordinate field for w3
subroutine reference_profile(ndf_w0,ndf_w3,exner_s,rho_s,theta_s,z_w0,z_w3)

integer,       intent(in)     :: ndf_w0, ndf_w3
real(kind=r_def), intent(in)  :: z_w0(ndf_w0), z_w3(ndf_w3)
real(kind=r_def), intent(out) :: exner_s(ndf_w3), rho_s(ndf_w3), theta_s(ndf_w0)

real(kind=r_def), parameter :: theta_surf = 300.0_r_def
real(kind=r_def), parameter :: exner_surf = 1.0_r_def
real(kind=r_def), parameter :: rho_surf   = 1.0_r_def
real(kind=r_def)            :: theta_w3, nsq_over_g

integer :: df

nsq_over_g = n_sq/gravity

do df = 1, ndf_w0
  theta_s(df) = theta_surf * exp ( nsq_over_g * z_w0(df) )
end do
do df = 1, ndf_w3
  exner_s(df) = exner_surf - gravity**2/(cp*theta_surf*n_sq)   &
              * (1.0_r_def - exp ( - nsq_over_g * z_w3(df) ))
  
  theta_w3    = theta_surf * exp ( nsq_over_g * z_w3(df) )
  rho_s(df)   = p_zero/(rd*theta_w3) * exner_s(df) ** ((1.0_r_def - kappa)/kappa) 
end do

end subroutine reference_profile

end module reference_profile_mod
