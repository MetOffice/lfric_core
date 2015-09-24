!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Functions to compute a variety of analytic profiles
!> @details Collection of functions to return the value of a field at a given
!!          point based upon a specified analytic formula
module analytic_profiles_mod

use constants_mod, only: r_def, &
                         ZERO_WIND, &
                         SOLID_BODY_ROTATION_WIND, &
                         CONSTANT_UV_WIND, &
                         CONSTANT_SHEAR_UV_WIND
 
use log_mod,       only : log_event,         &
                   log_scratch_space, &
                   LOG_LEVEL_ERROR

implicit none

  
contains

!> @brief Compute an analytic wind field
!> @param[in] chi Position in physical coordinates
!> @param[in] choice Integer defining which specified formula to use
!> @param[in] num_option Number of sclaer options to supply
!> @param[in] option Array of real values used to generate the initial profile
!> @result u The result wind field vector (u,v,w)
function analytic_wind(chi, choice, num_options, option) result(u)
  implicit none
  real(kind=r_def), intent(in) :: chi(3)
  integer,          intent(in) :: choice, num_options
  real(kind=r_def), optional   :: option(num_options)
  real(kind=r_def)             :: u(3)

  if ( .not. present(option) ) option(:) = 0.0_r_def

  select case ( choice )

    case ( zero_wind )
      u(:) = 0.0_r_def
    case ( solid_body_rotation_wind )
      u(1) = option(1) * ( cos(chi(2))*cos(option(2)) &
                       + sin(chi(1))*sin(chi(2))*sin(option(2)) )
      u(2) = option(1) * cos(chi(1))*sin(option(2))
      u(3) = 0.0_r_def
    case ( constant_uv_wind )
      u(1) = option(1)
      u(2) = option(2)
      u(3) = 0.0_r_def
    case ( constant_shear_uv_wind )
      u(1) = option(1)*chi(3)/option(3)
      u(2) = option(2)*chi(3)/option(3)
      u(3) = 0.0_r_def 
    case default
      write( log_scratch_space, '(A)' )  'Invalid velocity profile choice, stopping'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select

end function analytic_wind  

end module analytic_profiles_mod
