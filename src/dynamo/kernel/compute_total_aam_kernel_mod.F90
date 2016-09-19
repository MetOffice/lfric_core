!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the cell integrated axial angular momentum

!> @details The kernel computes the  cell integrated axial angular momentum:
!> \f[ \int ( \rho * \hat{z} . [r \times { u + \Omega \times r } ] ) dV \f]
module compute_total_aam_kernel_mod
use argument_mod,      only : arg_type, func_type,         &
                              GH_FIELD, GH_READ, GH_WRITE, &
                              ANY_SPACE_9, W2, W3,                  &
                              GH_BASIS, GH_DIFF_BASIS,     &
                              CELLS
use constants_mod,     only : r_def
use kernel_mod,        only : kernel_type
use planet_config_mod, only : scaled_omega, scaled_radius

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: compute_total_aam_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, W3),                             &
       arg_type(GH_FIELD,   GH_READ,  W2),                             &
       arg_type(GH_FIELD,   GH_READ,  W3),                             &
       arg_type(GH_FIELD*3, GH_READ,  ANY_SPACE_9)                              &
       /)
  type(func_type) :: meta_funcs(3) = (/                                &
       func_type(W2, GH_BASIS),                                        &
       func_type(W3, GH_BASIS),                                        &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                          &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass ::compute_total_aam_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface compute_total_aam_kernel_type
   module procedure compute_total_aam_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public compute_total_aam_code
contains

type(compute_total_aam_kernel_type) function compute_total_aam_kernel_constructor() result(self)
  return
end function compute_total_aam_kernel_constructor

!> @brief The subroutine to compute the total axial angular momentum
!! @param[in] nlayers Number of layers
!! @param[out] aam Cell integrated axial angular momentum
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2 Number unique of degrees of freedom  for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Basis functions evaluated at quadrature points 
!! @param[in] u Velocity array
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number unique of degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Basis functions evaluated at gaussian quadrature points 
!! @param[in] rho density
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number unique of degrees of freedom  for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Basis functions evaluated at gaussian quadrature points 
!! @param[in] chi_diff_basis Differntial of the basis functions evaluated at gaussian quadrature point
!! @param[in] chi_1 Physical x coordinate in chi
!! @param[in] chi_2 Physical y coordinate in chi
!! @param[in] chi_3 Physical z coordinate in chi
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine compute_total_aam_code(                                                   &
                                  nlayers,                                           &
                                  aam, u, rho, chi_1, chi_2, chi_3,                  &
                                  ndf_w3, undf_w3, map_w3, w3_basis,                 &
                                  ndf_w2, undf_w2, map_w2, w2_basis,                 &
                                  ndf_chi, undf_chi, map_chi, chi_basis, chi_diff_basis,  &
                                  nqp_h, nqp_v, wqp_h, wqp_v                         &
                                 )

  use coordinate_jacobian_mod, only: coordinate_jacobian
  use coord_transform_mod,     only: xyz2llr, cart2sphere_vector
  use cross_product_mod,       only: cross_product

  !Arguments
  integer,                     intent(in)    :: nlayers, nqp_h, nqp_v
  integer,                     intent(in)    :: ndf_chi, ndf_w2, ndf_w3
  integer,                     intent(in)    :: undf_chi, undf_w2, undf_w3
  integer, dimension(ndf_chi), intent(in)   :: map_chi
  integer, dimension(ndf_w2),  intent(in)    :: map_w2
  integer, dimension(ndf_w3),  intent(in)    :: map_w3

  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in) :: w3_basis
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w3), intent(out)   :: aam
  real(kind=r_def), dimension(undf_w2), intent(in)    :: u
  real(kind=r_def), dimension(undf_w3), intent(in)    :: rho
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer               :: df, k, loc
  integer               :: qp1, qp2

  real(kind=r_def), dimension(ndf_chi)          :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_w3)          :: rho_e, aam_e
  real(kind=r_def), dimension(ndf_w2)          :: u_e

  real(kind=r_def) :: u_at_quad(3), scaled_omega_vec(3), r_vec(3), &
                      am(3), x_vec(3), u_vec(3), j_u(3)
  real(kind=r_def) :: rho_at_quad
  real(kind=r_def), parameter :: z_hat(3) = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)

  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_chi
      loc = map_chi(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             chi_diff_basis, jac, dj)
    do df = 1, ndf_w3
      rho_e(df) = rho( map_w3(df) + k )
      aam_e(df) = 0.0_r_def
    end do
    do df = 1, ndf_w2
      u_e(df) = u( map_w2(df) + k )
    end do
  ! compute the aam integrated over one cell
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        x_vec(:) = 0.0_r_def
        do df = 1, ndf_chi
          x_vec(1) = x_vec(1) + chi_1_e(df)*chi_basis(1,df,qp1,qp2)
          x_vec(2) = x_vec(2) + chi_2_e(df)*chi_basis(1,df,qp1,qp2)
          x_vec(3) = x_vec(3) + chi_3_e(df)*chi_basis(1,df,qp1,qp2)
        end do
        call xyz2llr(x_vec(1),x_vec(2),x_vec(3),r_vec(1),r_vec(2),r_vec(3))
        scaled_omega_vec(1) = 0.0_r_def
        scaled_omega_vec(2) = 2.0_r_def*scaled_omega*cos(r_vec(2))
        scaled_omega_vec(3) = 2.0_r_def*scaled_omega*sin(r_vec(2))

        rho_at_quad = 0.0_r_def 
        do df = 1, ndf_w3
          rho_at_quad  = rho_at_quad + rho_e(df)*w3_basis(1,df,qp1,qp2) 
        end do

        u_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w2
          u_at_quad(:) = u_at_quad(:) &
                       + u_e(df)*w2_basis(:,df,qp1,qp2)
        end do
   
        j_u = matmul(jac(:,:,qp1,qp2),u_at_quad)
        u_vec(:) = cart2sphere_vector(x_vec,j_u) &
                 + cross_product(scaled_omega_vec,r_vec)*dj(qp1,qp2)

        am(:) = cross_product(r_vec,u_vec)
        do df = 1, ndf_w3
          aam_e(df) = aam_e(df) &
                      + wqp_h(qp1)*wqp_v(qp2)*rho_at_quad*dot_product(z_hat,am)/scaled_radius**2
        end do
      end do
    end do
    do df = 1,ndf_w3
      aam(map_w3(df)+k) = aam_e(df)
    end do
  end do

end subroutine compute_total_aam_code

end module compute_total_aam_kernel_mod
