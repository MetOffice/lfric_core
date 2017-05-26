!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Compute the mass matrix for W3 weighted by 1/rho

module weighted_m3_rho_kernel_mod
use constants_mod,           only: r_def, i_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, func_type,                      &
                                   GH_OPERATOR, GH_FIELD, GH_READ, GH_WRITE, &
                                   ANY_SPACE_1, W3, GH_BASIS, GH_DIFF_BASIS, &
                                   CELLS, GH_QUADRATURE_XYoZ, QUADRATURE_XYoZ
use coordinate_jacobian_mod, only: coordinate_jacobian
implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public, extends(kernel_type) :: weighted_m3_rho_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_OPERATOR, GH_WRITE, W3, W3),                        &
       arg_type(GH_FIELD,    GH_READ,  W3),                            &
       arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_1)                    &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(W3, GH_BASIS),                                        &
       func_type(ANY_SPACE_1, GH_DIFF_BASIS)                           &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_QUADRATURE_XYoZ
  ! gh_shape replaces evaluator_shape and will be removed by #1066
  integer :: evaluator_shape = QUADRATURE_XYoZ
contains
  procedure, nopass :: weighted_m3_rho_code
end type weighted_m3_rho_kernel_type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface weighted_m3_rho_kernel
   module procedure weighted_m3_rho_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public weighted_m3_rho_code
contains

type(weighted_m3_rho_kernel_type) function weighted_m3_rho_constructor() result(self)
  return
end function weighted_m3_rho_constructor
  
!> @brief Computes the mass matrix for the w3 space weighted by the reference
!!        density
!! @param[in] cell Cell number
!! @param[in] nlayers Number of layers.
!! @param[in] ncell_3d ncell*ndf
!! @param[inout] mm Mass matrix data array
!! @param[in] rho Density
!! @param[in] chi1 Chi in the first dir
!! @param[in] chi2 Chi in the 2nd dir
!! @param[in] chi3 Chi in the 3rd dir
!! @param[in] ndf_w3 Number of degrees of freedom per cell for the operator space.
!! @param[in] undf_w3 Total number of degrees of freedom for the W3 space
!! @param[in] map_w3 Dofmap for the bottom layer in the W3 space
!! @param[in] basis_w3 Basis functions evaluated at quadrature points.
!! @param[in] ndf_chi Number of degrees of freedom per cell for the coordinate field.
!! @param[in] undf_chi Number of unique degrees of freedum  for chi field
!! @param[in] map_chi Dofmap for the cell at the base of the column.
!! @param[in] diff_basis_chi Differential basis functions evaluated at quadrature points.
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine weighted_m3_rho_code(cell, nlayers, ncell_3d,            &
                                mm,                                 &
                                rho,                                &
                                chi1, chi2, chi3,                   &
                                ndf_w3, undf_w3, map_w3, basis_w3,  &
                                ndf_chi, undf_chi,                  &
                                map_chi, diff_basis_chi,            & 
                                nqp_h, nqp_v, wqp_h, wqp_v )

  implicit none

  !Arguments
  integer(kind=i_def), intent(in)     :: cell, nqp_h, nqp_v
  integer(kind=i_def), intent(in)     :: nlayers, ndf_w3, ndf_chi, undf_chi, undf_w3
  integer(kind=i_def), intent(in)     :: ncell_3d

  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3

  real(kind=r_def), dimension(ndf_w3,ndf_w3,ncell_3d),  intent(inout)  :: mm

  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: diff_basis_chi
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v),  intent(in) :: basis_w3

  real(kind=r_def), dimension(undf_w3),  intent(in)           :: rho
  real(kind=r_def), dimension(undf_chi), intent(in)           :: chi1
  real(kind=r_def), dimension(undf_chi), intent(in)           :: chi2
  real(kind=r_def), dimension(undf_chi), intent(in)           :: chi3

  real(kind=r_def), dimension(nqp_h), intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) :: wqp_v

  !Internal variables
  integer(kind=i_def)                          :: df, df1, df2, k, ik, loc
  integer(kind=i_def)                          :: qp1, qp2
  real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
  real(kind=r_def), dimension(ndf_w3)          :: rho_e   
  real(kind=r_def)                             :: rho_quad
  real(kind=r_def)                             :: integrand
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac

  do k = 0, nlayers-1
    ik = 1 + k + (cell-1)*nlayers
     
    do df = 1, ndf_chi
      loc = map_chi(df) + k
      chi1_e(df) = chi1(loc)
      chi2_e(df) = chi2(loc)
      chi3_e(df) = chi3(loc)
    end do

    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                             diff_basis_chi, jac, dj)
    do df = 1, ndf_w3
      rho_e(df) = rho(map_w3(df) + k)
    end do

    do df2 = 1, ndf_w3
      do df1 = 1, ndf_w3 
        mm(df1,df2,ik) = 0.0_r_def
        do qp2 = 1, nqp_v
          do qp1 = 1, nqp_h 
            rho_quad = 0.0_r_def
            do df = 1,ndf_w3
              rho_quad = rho_quad + rho_e(df)*basis_w3(1,df,qp1,qp2)
            end do
            integrand = wqp_h(qp1)*wqp_v(qp2) & 
                       *basis_w3(1,df1,qp1,qp2)*basis_w3(1,df2,qp1,qp2)  &
                       /rho_quad*dj(qp1,qp2) 
            mm(df1,df2,ik) = mm(df1,df2,ik) + integrand
          end do
        end do
      end do
    end do  
  end do ! end of k loop

end subroutine weighted_m3_rho_code

end module weighted_m3_rho_kernel_mod
