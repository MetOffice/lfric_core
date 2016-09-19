!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Adds a Held-Suarez forcing

!> @details Kernel adds a Held-Suarez forcing based on Wedi and Smolarkiewicz 2009:
!> Wedi, N. P. and Smolarkiewicz, P. K. (2009), A framework for testing global 
!> non-hydrostatic models. Q.J.R. Meteorol. Soc., 135: 469-484. doi: 10.1002/qj.377

module held_suarez_kernel_mod
  
use kernel_mod,               only: kernel_type
use argument_mod,             only: arg_type, func_type,                 &
                                    GH_FIELD, GH_WRITE, GH_READ, GH_INC, &
                                    W0, W2, W3, ANY_SPACE_9,             &
                                    GH_BASIS, GH_DIFF_BASIS,             &
                                    CELLS
use constants_mod,            only: r_def
use coord_transform_mod,      only: xyz2llr
use calc_exner_pointwise_mod, only: calc_exner_pointwise
use planet_config_mod,        only: scaling_factor, kappa
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: held_suarez_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,   W2),                             &
       arg_type(GH_FIELD,   GH_INC,   W0),                             &
       arg_type(GH_FIELD,   GH_READ,  W2),                             &
       arg_type(GH_FIELD,   GH_READ,  W0),                             &
       arg_type(GH_FIELD,   GH_READ,  W3),                             &
       arg_type(GH_FIELD*3, GH_READ,  ANY_SPACE_9)                              &
       /)
  type(func_type) :: meta_funcs(4) = (/                                &
       func_type(W2, GH_BASIS),                        &
       func_type(W0, GH_BASIS),                         &
       func_type(W3, GH_BASIS),                                         &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                 &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: held_suarez_code
end type


!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface held_suarez_kernel_type
   module procedure held_suarez_kernel_constructor
end interface


!-------------------------------------------------------------------------------
! local parameters
!-------------------------------------------------------------------------------
! Held-Suarez parameters 
real(kind=r_def), parameter :: SIGMA_B = 0.7
real(kind=r_def), parameter :: KA = 1./86400.
real(kind=r_def), parameter :: KF = KA/40.0
real(kind=r_def), parameter :: KS = KA/4.0
real(kind=r_def), parameter :: T_MIN            = 200.0 ! Minimum/Stratospheric temperature 
real(kind=r_def), parameter :: T_SURF           = 315.0 ! surface temperature 
real(kind=r_def), parameter :: DT_EQ_POLE       = 60.0  ! Equator-Pole Temp diff (deltaT)_y
real(kind=r_def), parameter :: STATIC_STABILITY = 10.0  ! Static Stability temperature (delta \theta)_z


!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public held_suarez_code
contains

type(held_suarez_kernel_type) function held_suarez_kernel_constructor() result(self)
  return
end function held_suarez_kernel_constructor

!> @brief Adds a Held-Suarez forcing
!! @param[in] nlayers Number of layers
!! @param[inout] dtheta Theta increment
!! @param[in] theta Potential temperature
!! @param[inout] du Wind increment
!! @param[inout] u Wind
!! @param[in] rho Density
!! @param[in] chi_1 X component of the coordinate 
!! @param[in] chi_2 Y component of the coordinate 
!! @param[in] chi_3 Z component of the coordinate
!! @param[in] ndf_w0 Number of degrees of freedom per cell for w0
!! @param[in] undf_w0 Number of unique degrees of freedom for w0
!! @param[in] map_w0 Dofmap for the cell at the base of the column for w0
!! @param[in] basis_w0 Basis functions evaluated at gaussian quadrature points for w0
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2 Number of unique degrees of freedom for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] basis_w2 Basis functions evaluated at gaussian quadrature points for w2
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number of unique degrees of freedom for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] basis_w3 Basis functions evaluated at gaussian quadrature points for w3
!! @param[in] ndf_chi The number of degrees of freedom per cell for chi
!! @param[in] undf_chi The number of unique degrees of freedom for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] basis_chi Basis functions evaluated at gaussian quadrature points for chi
!! @param[in] diff_basis_chi Differential basis functions evaluated at gaussian quadrature points for chi
!! @param[in] nqp_h the number of horizontal quadrature points
!! @param[in] nqp_v the number of vertical quadrature points
!! @param[in] wqp_h the weights of the horizontal quadrature points
!! @param[in] wqp_v the weights of the vertical quadrature points
subroutine held_suarez_code(nlayers,                                           &
                            du, dtheta, u, theta, rho,                         &
                            chi_1, chi_2, chi_3,                               &
                            ndf_w2, undf_w2, map_w2, basis_w2,                 &
                            ndf_w0, undf_w0, map_w0, basis_w0,                 &
                            ndf_w3, undf_w3, map_w3, basis_w3,                 &
                            ndf_chi, undf_chi, map_chi, basis_chi, &
                            diff_basis_chi, nqp_h, nqp_v, wqp_h, wqp_v      &
                            )
  
  use coordinate_jacobian_mod,  only: coordinate_jacobian

  implicit none

  !Arguments
  integer, intent(in) :: nlayers

  integer, intent(in) :: ndf_w0, undf_w0  
  integer, intent(in) :: ndf_w2, undf_w2
  integer, intent(in) :: ndf_w3, undf_w3
  integer, intent(in) :: ndf_chi, undf_chi
  integer, intent(in) :: nqp_h, nqp_v
  real(kind=r_def), dimension(undf_w0), intent(inout) :: dtheta
  real(kind=r_def), dimension(undf_w2), intent(inout) :: du
  real(kind=r_def), dimension(undf_w0), intent(in)    :: theta
  real(kind=r_def), dimension(undf_w2), intent(in)    :: u
  real(kind=r_def), dimension(undf_w3), intent(in)    :: rho
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3

  integer, dimension(ndf_w0),  intent(in)                          :: map_w0
  real(kind=r_def), dimension(1, ndf_w0, nqp_h, nqp_v), intent(in) :: basis_w0

  integer, dimension(ndf_w2),  intent(in)                          :: map_w2
  real(kind=r_def), dimension(3, ndf_w2, nqp_h, nqp_v), intent(in) :: basis_w2 

  integer, dimension(ndf_w3),  intent(in)                          :: map_w3
  real(kind=r_def), dimension(1, ndf_w3, nqp_h, nqp_v), intent(in) :: basis_w3 

  integer, dimension(ndf_chi),  intent(in)                          :: map_chi
  real(kind=r_def), dimension(1, ndf_chi, nqp_h, nqp_v), intent(in) :: basis_chi
  real(kind=r_def), dimension(3, ndf_chi, nqp_h, nqp_v), intent(in) :: diff_basis_chi

  real(kind=r_def), dimension(nqp_h), intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) :: wqp_v

  !Internal variables
  integer               :: df, k, loc

  real(kind=r_def)            :: theta_eq, exner
  real(kind=r_def)            :: lon, r
  integer                     :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_at_dof, chi_2_at_dof, chi_3_at_dof
  real(kind=r_def), dimension(ndf_w3)          :: rho_at_dof
  real(kind=r_def), dimension(ndf_w2)          :: u_at_dof
  real(kind=r_def), dimension(ndf_w0)          :: theta_at_dof
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_w2)          :: du_at_dof
  real(kind=r_def), dimension(ndf_w0)          :: dtheta_at_dof
  real(kind=r_def)                             :: rho_at_quad, theta_at_quad
  real(kind=r_def)                             :: u_at_quad(3), jac_v(3), chi_at_quad(3), lat_at_quad

  real(kind=r_def) :: exner0 ! lowest level exner value
  real(kind=r_def) :: sigma  ! exner/exner0
             
  exner0 = 0.0
  dtheta(:) = 0.0
  du(:) = 0.0
  do k = 0, nlayers-1
    dtheta_at_dof(:) = 0.0
    du_at_dof(:) = 0.0

   ! Store the values for each degree of freedom
   do df = 1, ndf_chi
     loc = map_chi(df) + k
     chi_1_at_dof(df) = chi_1( loc )
     chi_2_at_dof(df) = chi_2( loc )
     chi_3_at_dof(df) = chi_3( loc )
   end do
   do df = 1, ndf_w0
     theta_at_dof(df) = theta( map_w0(df) + k )
   end do   
   do df = 1, ndf_w2
     u_at_dof(df) = u( map_w2(df) + k )
   end do  
   do df = 1, ndf_w3
     rho_at_dof(df) = rho( map_w3(df) + k )
   end do 
   ! Compute integrals over each cell       
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        rho_at_quad = 0.0_r_def 
        do df = 1, ndf_w3
          rho_at_quad  = rho_at_quad + rho_at_dof(df)*basis_w3(1,df,qp1,qp2) 
        end do
        theta_at_quad = 0.0_r_def
        chi_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w0
          theta_at_quad   = theta_at_quad                                      &
                          + theta_at_dof(df)*basis_w0(1,df,qp1,qp2)
        end do
        do df = 1, ndf_chi
          chi_at_quad(1)     = chi_at_quad(1) + chi_1_at_dof(df)*basis_chi(1,df,qp1,qp2)
          chi_at_quad(2)     = chi_at_quad(2) + chi_2_at_dof(df)*basis_chi(1,df,qp1,qp2)
          chi_at_quad(3)     = chi_at_quad(3) + chi_3_at_dof(df)*basis_chi(1,df,qp1,qp2)
        end do

        ! Now calculate rhs for theta 
        call xyz2llr(chi_at_quad(1), chi_at_quad(2), chi_at_quad(3), lon, lat_at_quad, r)
        u_at_quad(:) = 0.0_r_def 
        do df = 1, ndf_w2
          u_at_quad(:) = u_at_quad(:) &
             + u_at_dof(df)*basis_w2(:,df,qp1,qp2)
        end do
        exner = calc_exner_pointwise(rho_at_quad, theta_at_quad)

        if (k==0)then
          sigma = 1.0_r_def
          ! calculate mean exner value in bottom layer of quadrature points
          if (qp2==1)exner0 = exner0 + exner/real(nqp_v, r_def)
        else
          sigma = (exner/exner0)**(1.0/kappa)
        end if

        theta_eq = held_suarez_equilibrium_theta(exner, lat_at_quad) 
        do df = 1, ndf_w0
          dtheta_at_dof(df) = dtheta_at_dof(df) & 
             - held_suarez_newton_frequency(sigma, lat_at_quad) * &
             (theta_at_quad - theta_eq)
        end do

        ! Now calculate rhs for winds
        call coordinate_jacobian(ndf_chi, nqp_h, nqp_v,                      &
                                 chi_1_at_dof, chi_2_at_dof, chi_3_at_dof,  &
                                 diff_basis_chi, jac, dj)

        u_at_quad(:) = wqp_h(qp1)*wqp_v(qp2)&
           * matmul(jac(:,:,qp1,qp2), u_at_quad(:)*held_suarez_damping(sigma))/dj(qp1,qp2)

        do df = 1, ndf_w2
          jac_v = matmul(jac(:,:,qp1,qp2), basis_w2(:,df,qp1,qp2))
          du_at_dof(df) = du_at_dof(df) +  dot_product(jac_v,u_at_quad(:))
        end do
      end do
    end do
    do df = 1, ndf_w0
      dtheta(map_w0(df) + k) = dtheta(map_w0(df) + k) + dtheta_at_dof(df) 
    end do
    do df = 1, ndf_w2
      du(map_w2(df) + k) = du(map_w2(df) + k) + du_at_dof(df) 
    end do
  end do

end subroutine held_suarez_code

!> @brief Function to calculate equilibrium theta profile for held-suarez
!! @param[in] exner exner pressure
!! @param[in] lat latitude
function held_suarez_equilibrium_theta(exner, lat) result(theta_eq)

  implicit none

  real(kind=r_def), intent(in)  :: exner, lat
  real(kind=r_def) :: theta_eq         ! Equilibrium theta

  theta_eq = max(T_MIN, (T_SURF - DT_EQ_POLE*sin(lat)*sin(lat) &
     - STATIC_STABILITY*log(exner)*cos(lat)*cos(lat)/kappa))*exner

end function held_suarez_equilibrium_theta

!> @brief Function to calculate the newton relaxation frequency for held-suarez
!! @param[in] sigma nondimensional pressure p/p_surf
!! @param[in] lat latitude
function held_suarez_newton_frequency(sigma, lat) result(held_suarez_frequency)

  implicit none

  real(kind=r_def), intent(in) :: sigma 
  real(kind=r_def), intent(in) :: lat
  real(kind=r_def)             :: held_suarez_frequency
  real(kind=r_def)             :: sigma_func

  sigma_func = max((sigma - SIGMA_B)/(1.0 - SIGMA_B), 0.0_r_def)
  held_suarez_frequency = KA + (KS - KA)*sigma_func*(cos(lat)**4)

  ! If running on a scaled planet, then reduce the timescale...
  held_suarez_frequency = held_suarez_frequency*scaling_factor

end function held_suarez_newton_frequency
 
!> @brief Function to calculate the damping coefficent for held-suarez
!! @param[in] sigma nondimensional pressure p/p_surf
function held_suarez_damping(sigma) result(held_suarez_damping_rate)

  implicit none

  real(kind=r_def), intent(in) :: sigma 
  real(kind=r_def)             :: held_suarez_damping_rate
  real(kind=r_def) :: sigma_func

  sigma_func = max((sigma - SIGMA_B)/(1.0 - SIGMA_B), 0.0_r_def)
  held_suarez_damping_rate = -KF*sigma_func

  ! If running on a scaled planet, then reduce the timescale...
  held_suarez_damping_rate = held_suarez_damping_rate*scaling_factor

end function held_suarez_damping

end module held_suarez_kernel_mod
