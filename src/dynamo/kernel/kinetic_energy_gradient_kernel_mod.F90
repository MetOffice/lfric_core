!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes the kinetic gradient component of the rhs of the momentum equation 
!>        for the nonlinear equations, written in the vector invariant form


!> @details The kernel computes the kinetic gradient component of the rhs of the
!>         momentum equation for the nonlinear equations,
!>         written in the vector invariant form
!>         This consists of four terms:
!>         This consists of four terms:
!>         Pressure gradient: \f[ cp*\theta*\nabla(\Pi)\f]
!>         geopotential gradient: \f[ \nabla(\Phi) ( \equiv g for some domains)\f]
!>         gradient of kinetic energy: \f[ \nabla(1/2*u.u) \f] 
!>         vorticity advection: \f[ \xi/\rho \times F (with vorticity \xi and mass flux F) \f]
!>         This results in:
!>         \f[ r_u = -\xi/\rho \times F - \nabla(\Phi + 1/2*u.u) - cp*\theta*\nabla(\Pi) \f]
module kinetic_energy_gradient_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,                 &
                                    GH_FIELD, GH_READ, GH_INC,           &
                                    ANY_SPACE_9, W2, GH_BASIS, GH_DIFF_BASIS,     &
                                    CELLS 
use constants_mod,           only : r_def

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: kinetic_energy_gradient_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W2),                              &
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9)                               &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(W2, GH_BASIS, GH_DIFF_BASIS),                         &
       func_type(ANY_SPACE_9,  GH_DIFF_BASIS)                          &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass ::kinetic_energy_gradient_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface kinetic_energy_gradient_kernel_type
   module procedure kinetic_energy_gradient_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public kinetic_energy_gradient_code
contains

type(kinetic_energy_gradient_kernel_type) function kinetic_energy_gradient_kernel_constructor() result(self)
  return
end function kinetic_energy_gradient_kernel_constructor

!> @brief Computes the kinetic gradient component of the rhs of the momentum equation 
!! @param[in] nlayers Number of layers
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2 Number unique of degrees of freedom  for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Basis functions evaluated at quadrature points 
!! @param[in] w2_diff_basis Differntial of the basis functions evaluated at  quadrature points
!! @param[inout] r_u Right hand side of momentum equation
!! @param[in] u Velocity
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number unique of degrees of freedom  for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_diff_basis Differntial of the basis functions evaluated at gaussian quadrature point
!! @param[in] chi_1 Physical x coordinate in chi
!! @param[in] chi_2 Physical y coordinate in chi
!! @param[in] chi_3 Physical z coordinate in chi
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine kinetic_energy_gradient_code(nlayers,                                          &
                                        r_u, u, chi_1, chi_2, chi_3,                      &
                                        ndf_w2, undf_w2, map_w2, w2_basis, w2_diff_basis, &
                                        ndf_chi, undf_chi, map_chi, chi_diff_basis,       &
                                        nqp_h, nqp_v, wqp_h, wqp_v                        &
                                        )
                           
  use coordinate_jacobian_mod,  only: coordinate_jacobian
  
  !Arguments
  integer, intent(in) :: nlayers,nqp_h, nqp_v
  integer, intent(in) :: ndf_chi, ndf_w2
  integer, intent(in) :: undf_chi, undf_w2
  integer, dimension(ndf_chi), intent(in) :: map_chi
  integer, dimension(ndf_w2), intent(in) :: map_w2

  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_basis 
  real(kind=r_def), dimension(1,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_diff_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w2), intent(inout) :: r_u
  real(kind=r_def), dimension(undf_w2), intent(in)    :: u
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer               :: df, k, loc 
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_w2)          :: ru_e, u_e

  real(kind=r_def) :: u_at_quad(3)
  real(kind=r_def) :: ke_at_quad, dv
  
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

    do df = 1, ndf_w2
      u_e(df) = u( map_w2(df) + k )
      ru_e(df) = 0.0_r_def
    end do    
  ! compute the RHS integrated over one cell    
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
! k.e term
        u_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w2
          u_at_quad(:) = u_at_quad(:) &
                       + u_e(df)*w2_basis(:,df,qp1,qp2)
        end do
        ke_at_quad = 0.5_r_def*dot_product(matmul(jac(:,:,qp1,qp2),u_at_quad), &
                                           matmul(jac(:,:,qp1,qp2),u_at_quad))/(dj(qp1,qp2)**2)

        do df = 1, ndf_w2
          dv = w2_diff_basis(1,df,qp1,qp2)
          ru_e(df) = ru_e(df) +  wqp_h(qp1)*wqp_v(qp2)*dv*ke_at_quad                                                       

        end do
      end do
    end do
    do df = 1, ndf_w2
      r_u( map_w2(df) + k ) =  r_u( map_w2(df) + k ) + ru_e(df)
    end do 
  end do
  
end subroutine kinetic_energy_gradient_code

end module kinetic_energy_gradient_kernel_mod
