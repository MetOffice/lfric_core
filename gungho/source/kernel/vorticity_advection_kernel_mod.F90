!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the vorticity component of the rhs of the momentum 
!>        equation for the nonlinear equations, written in the vector invariant form


!> @details The kernel computes the  vorticity component of the rhs of the momentum equation 
!>         for the nonlinear equations, written in the vector invariant form
!>         This consists of four terms:
!>         Pressure gradient: \f[ cp*\theta*\nabla(\Pi) \f]
!>         geopotential gradient: \f[ \nabla(\Phi) ( \equiv g for some domains) \f]
!>         gradient of kinetic energy: \f[ \nabla(1/2*u.u) \f] 
!>         vorticity advection: \f[ \xi/\rho \times F (with vorticity \xi and mass flux F) \f]
!>         This results in:
!>         \f[ r_u = -\xi/\rho \times F - \nabla(\Phi + 1/2*u.u) - cp*\theta*\nabla(\Pi) \f]
module vorticity_advection_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,                 &
                                    GH_FIELD, GH_READ, GH_INC,           &
                                    ANY_SPACE_9, W1, W2, W3,             &
                                    GH_BASIS, GH_DIFF_BASIS,             &
                                    CELLS, GH_QUADRATURE_XYoZ,           &
                                    QUADRATURE_XYoZ
use constants_mod,           only : r_def
use cross_product_mod,       only : cross_product

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: vorticity_advection_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                  &
       arg_type(GH_FIELD,   GH_INC,  W2),                              &
       arg_type(GH_FIELD,   GH_READ, W2),                              &
       arg_type(GH_FIELD,   GH_READ, W3),                              &
       arg_type(GH_FIELD,   GH_READ, W1),                              &
       arg_type(GH_FIELD*3, GH_READ, any_space_9)                      &
       /)
  type(func_type) :: meta_funcs(4) = (/                                &
       func_type(W2, GH_BASIS),                                        &
       func_type(W3, GH_BASIS),                                        &
       func_type(W1, GH_BASIS),                                        &
       func_type(ANY_SPACE_9, GH_DIFF_BASIS)                           &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = QUADRATURE_XYoZ
contains
  procedure, nopass ::vorticity_advection_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface vorticity_advection_kernel_type
   module procedure vorticity_advection_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public vorticity_advection_code
contains

type(vorticity_advection_kernel_type) function vorticity_advection_kernel_constructor() result(self)
  return
end function vorticity_advection_kernel_constructor

!> @brief Compute the advection of the wind field by the vorticity
!! @param[in] nlayers Number of layers
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2 Number unique of degrees of freedom  for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Basis functions evaluated at quadrature points 
!! @param[inout] r_u Right hand side of the momentum equation
!! @param[in] mass_flux Mass flux
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number unique of degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Basis functions evaluated at gaussian quadrature points 
!! @param[in] rho Density
!! @param[in] ndf_w1 Number of degrees of freedom per cell for w1
!! @param[in] undf_w1 Number unique of degrees of freedom  for w1
!! @param[in] map_w1 Dofmap for the cell at the base of the column for w1
!! @param[in] w1_basis Basis functions evaluated at gaussian quadrature points 
!! @param[in] xi Vorticity
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi
!! @param[in] undf_chi Number unique of degrees of freedom  for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_diff_basis Differential of the basis functions evaluated at gaussian quadrature point
!! @param[in] chi_1 Physical x coordinate in chi
!! @param[in] chi_2 Physical y coordinate in chi
!! @param[in] chi_3 Physical z coordinate in chi
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine vorticity_advection_code(nlayers,                                           &
                                    r_u, mass_flux, rho, xi,                           &
                                    chi_1, chi_2, chi_3,                               &
                                    ndf_w2, undf_w2, map_w2, w2_basis,                 &
                                    ndf_w3, undf_w3, map_w3, w3_basis,                 &
                                    ndf_w1, undf_w1, map_w1, w1_basis,                 &
                                    ndf_chi, undf_chi, map_chi, chi_diff_basis,        &
                                    nqp_h, nqp_v, wqp_h, wqp_v                         &
                                    )
                           
  use coordinate_jacobian_mod,  only: coordinate_jacobian, &
                                      coordinate_jacobian_inverse
  
  !Arguments
  integer, intent(in) :: nlayers,nqp_h, nqp_v
  integer, intent(in) :: ndf_chi, ndf_w1, ndf_w2, ndf_w3
  integer, intent(in) :: undf_chi, undf_w1, undf_w2, undf_w3
  integer, dimension(ndf_chi), intent(in) :: map_chi
  integer, dimension(ndf_w1), intent(in) :: map_w1
  integer, dimension(ndf_w2), intent(in) :: map_w2
  integer, dimension(ndf_w3), intent(in) :: map_w3

  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in) :: w3_basis  
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_basis 
  real(kind=r_def), dimension(3,ndf_w1,nqp_h,nqp_v), intent(in) :: w1_basis 
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w2), intent(inout) :: r_u
  real(kind=r_def), dimension(undf_w2), intent(in)    :: mass_flux
  real(kind=r_def), dimension(undf_w3), intent(in)    :: rho
  real(kind=r_def), dimension(undf_w1), intent(in)    :: xi
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer               :: df, k, loc 
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac, jac_inv
  real(kind=r_def), dimension(ndf_w3)          :: rho_e
  real(kind=r_def), dimension(ndf_w2)          :: ru_e, f_e
  real(kind=r_def), dimension(ndf_w1)          :: xi_e

  real(kind=r_def), dimension(3) :: xi_at_quad, f_at_quad, jac_v, &
                                    vorticity, v, j_xi, j_f
  real(kind=r_def) :: rho_at_quad, vorticity_term
  
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
    call coordinate_jacobian_inverse(nqp_h, nqp_v, jac, dj, jac_inv)   

    do df = 1, ndf_w3
      rho_e(df) = rho( map_w3(df) + k )
    end do    
    do df = 1, ndf_w2
      f_e(df) = mass_flux( map_w2(df) + k )
      ru_e(df) = 0.0_r_def
    end do    
    do df = 1, ndf_w1
      xi_e(df) = xi( map_w1(df) + k )
    end do    
  ! compute the RHS integrated over one cell    
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        rho_at_quad = 0.0_r_def 
        do df = 1, ndf_w3
          rho_at_quad  = rho_at_quad + rho_e(df)*w3_basis(1,df,qp1,qp2) 
        end do

        f_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w2
          f_at_quad(:) = f_at_quad(:) &
                       + f_e(df)*w2_basis(:,df,qp1,qp2)
        end do

! Vorticity advection term
        xi_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w1
          xi_at_quad(:) = xi_at_quad(:) &
                        + xi_e(df)*w1_basis(:,df,qp1,qp2)
        end do
        j_xi = matmul(transpose(jac_inv(:,:,qp1,qp2)),xi_at_quad)
        j_f  = matmul(jac(:,:,qp1,qp2),f_at_quad)
        vorticity = cross_product(j_xi, j_f)
        vorticity = vorticity/(dj(qp1,qp2)*rho_at_quad)

        do df = 1, ndf_w2
          v  = w2_basis(:,df,qp1,qp2)
          jac_v = matmul(jac(:,:,qp1,qp2),v)

          vorticity_term = - dot_product(jac_v,vorticity)

          ru_e(df) = ru_e(df) +  wqp_h(qp1)*wqp_v(qp2)*vorticity_term

        end do
      end do
    end do
    do df = 1, ndf_w2
      r_u( map_w2(df) + k ) =  r_u( map_w2(df) + k ) + ru_e(df)
    end do 
  end do
  
end subroutine vorticity_advection_code

end module vorticity_advection_kernel_mod
