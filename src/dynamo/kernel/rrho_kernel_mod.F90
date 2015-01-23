!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes rhs of the continuity equation

!> @detail The kernel computes thr rhs of the continuity equation for the linear equations with
!>         no advection
module rrho_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, &          ! the type
                                    gh_read, gh_write, w0, w2, w3, fe, cells 
use reference_profile_mod,   only : reference_profile
use constants_mod,           only : n_sq, gravity, r_def
use mesh_generator_mod,      only : xyz2llr, sphere2cart_vector
use mesh_mod,                only : l_spherical

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: rrho_kernel_type
  private
  type(arg_type) :: meta_args(5) = [  &
       arg_type(gh_write,w3,fe,.true., .false.,.false.,.true.),       &
       arg_type(gh_read ,w2,fe,.true., .true. ,.false., .false.),     &
       arg_type(gh_read ,w0,fe,.false.,.true. ,.false., .false.),     &
       arg_type(gh_read ,w0,fe,.false.,.false.,.false.,.false.),      &
       arg_type(gh_read ,w0,fe,.false.,.false.,.false.,.false.)       &
       ]
  integer :: iterates_over = cells
contains
  procedure, nopass ::rrho_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface rrho_kernel_type
   module procedure rrho_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public rrho_code
contains

type(rrho_kernel_type) function rrho_kernel_constructor() result(self)
  return
end function rrho_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_w3 The number of degrees of freedom per cell for w3
!! @param[in] undf_w3 The number of (local) unique degrees of freedom
!! @param[in] map_w3 Integer array holding the dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Real 4-dim array holding basis functions evaluated at quadrature points 
!! @param[inout] r_rho Real array the data 
!! @param[in] ndf_w2 The number of degrees of freedom per cell for w2
!! @param[in] undf_w2 The number of (local) unique degrees of freedom
!! @param[in] map_w2 Integer array holding the dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis Real 4-dim array holding basis functions evaluated at quadrature points 
!! @param[in] w2_diff_basis Real 4-dim array holding differential basis functions evaluated at quadrature points 
!! @param[in] orientation Integer array holding the orientation of the fs
!! @param[in] u Real array. The velocity data
!! @param[in] ndf_w0 The number of degrees of freedom per cell for w0
!! @param[in] undf_w0 The number of (local) unique degrees of freedom
!! @param[in] map_w0 Integer array holding the dofmap for the cell at the base of the column for w0
!! @param[in] w0_basis Real 4-dim array holding the basis functions for w0 evaluated at  quadrature points 
!! @param[in] w0_diff_basis Real 4-dim array holding differential of the basis functions for w0 evaluated at quadrature points 
!! @param[in] chi_1 Real array. the physical x coordinate in w0
!! @param[in] chi_2 Real array. the physical y coordinate in w0
!! @param[in] chi_3 Real array. the physical z coordinate in w0
!! @param[in] nqp_h Integer, number of quadrature points in the horizontal
!! @param[in] nqp_v Integer, number of quadrature points in the vertical
!! @param[in] wqp_h Real array. Quadrature weights horizontal
!! @param[in] wqp_v Real array. Quadrature weights vertical
subroutine rrho_code(nlayers,                                                  &
                     ndf_w3, undf_w3, map_w3, w3_basis, r_rho,                 &
                     ndf_w2, undf_w2, map_w2, w2_basis, w2_diff_basis,         &
                     orientation, u,                                           &
                     ndf_w0, undf_w0, map_w0, w0_basis, w0_diff_basis,         &
                     chi_1, chi_2, chi_3,                                      &
                     nqp_h, nqp_v, wqp_h, wqp_v         )
                             
  use coordinate_jacobian_mod, only: coordinate_jacobian
  use reference_profile_mod,   only: reference_profile                          
  
  !Arguments
  integer, intent(in) :: nlayers, nqp_h, nqp_v
  integer, intent(in) :: ndf_w0, ndf_w2, ndf_w3
  integer, intent(in) :: undf_w0, undf_w2, undf_w3
  integer, dimension(ndf_w3), intent(in) :: map_w3
  integer, dimension(ndf_w2), intent(in) :: map_w2
  integer, dimension(ndf_w0), intent(in) :: map_w0
  integer, dimension(ndf_w2), intent(in) :: orientation

  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in) :: w3_basis  
  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_basis 
  real(kind=r_def), dimension(1,ndf_w0,nqp_h,nqp_v), intent(in) :: w0_basis 
  real(kind=r_def), dimension(1,ndf_w2,nqp_h,nqp_v), intent(in) :: w2_diff_basis
  real(kind=r_def), dimension(3,ndf_w0,nqp_h,nqp_v), intent(in) :: w0_diff_basis 

  real(kind=r_def), dimension(undf_w3), intent(inout) :: r_rho
  real(kind=r_def), dimension(undf_w0), intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_w2), intent(in)    :: u

  real(kind=r_def), dimension(nqp_h), intent(in)      ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

  !Internal variables
  integer               :: df, k, loc 
  integer               :: qp1, qp2
  
  real(kind=r_def), dimension(ndf_w0)           :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(ndf_w2)           :: u_e
  real(kind=r_def), dimension(nqp_h,nqp_v)        :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v)    :: jac
  real(kind=r_def), dimension(ndf_w3) :: rrho_e
  real(kind=r_def) :: rho_s_at_quad, exner_s_at_quad, &
                      theta_s_at_quad, div_u_at_quad, &
                      div_term, buoy_term, vec_term
  real(kind=r_def) :: u_at_quad(3), k_sphere(3), k_cart(3), &
                      x_at_quad(3), llr(3)

  k_sphere = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def /)

  do k = 0, nlayers-1
  ! Extract element arrays of chi
    do df = 1, ndf_w0
      loc = map_w0(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do
    do df = 1, ndf_w3
      rrho_e(df) = 0.0_r_def
    end do
    call coordinate_jacobian(ndf_w0, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             w0_diff_basis, jac, dj)
    do df = 1, ndf_w2
      u_e(df) = u( map_w2(df) + k )*real(orientation(df))
    end do
  ! compute the RHS integrated over one cell
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        x_at_quad(:) = 0.0_r_def
        do df = 1, ndf_w0
          x_at_quad(1) = x_at_quad(1) + chi_1_e(df)*w0_basis(1,df,qp1,qp2) 
          x_at_quad(2) = x_at_quad(2) + chi_2_e(df)*w0_basis(1,df,qp1,qp2)
          x_at_quad(3) = x_at_quad(3) + chi_3_e(df)*w0_basis(1,df,qp1,qp2)
        end do
        call reference_profile(exner_s_at_quad, rho_s_at_quad, & 
                               theta_s_at_quad, x_at_quad)
        u_at_quad(:) = 0.0_r_def
        div_u_at_quad = 0.0_r_def
        do df = 1, ndf_w2
          u_at_quad(:)  = u_at_quad(:)  + u_e(df)*w2_basis(:,df,qp1,qp2)
          div_u_at_quad = div_u_at_quad + u_e(df)*w2_diff_basis(1,df,qp1,qp2)
        end do
       
        if ( l_spherical ) then
          call xyz2llr(x_at_quad(1), x_at_quad(2), x_at_quad(3), &
                       llr(1), llr(2), llr(3))
          k_cart(:) = sphere2cart_vector(k_sphere, llr)   
        else
          k_cart(:) = k_sphere(:)
        end if     

        div_term  =  - rho_s_at_quad*div_u_at_quad 
        vec_term  = dot_product(matmul(jac(:,:,qp1,qp2),u_at_quad),k_cart)
        buoy_term = n_sq/gravity*rho_s_at_quad*vec_term

        do df = 1, ndf_w3
          rrho_e(df) = rrho_e(df) + wqp_h(qp1)*wqp_v(qp2)*w3_basis(1,df,qp1,qp2)*( buoy_term + div_term )
        end do
      end do
    end do
    do df = 1, ndf_w3
      r_rho( map_w3(df) + k ) =  rrho_e(df)
    end do 
  end do
  
end subroutine rrho_code

end module rrho_kernel_mod
