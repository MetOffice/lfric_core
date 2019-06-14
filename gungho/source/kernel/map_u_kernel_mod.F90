!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the rhs for the mapping of the wind field.
!>
!> The kernel computes a very crude approximation to the rhs of the equation u = u0 
!> where u0 is the physical wind field. The computational wind field is 
!> projected onto using Galerkin projection.
!>
module map_u_kernel_mod

  use argument_mod,            only : arg_type, func_type,       &
                                      GH_FIELD, GH_INC, GH_READ, &
                                      ANY_SPACE_9,               &
                                      GH_BASIS, GH_DIFF_BASIS,   &
                                      CELLS, GH_QUADRATURE_XYoZ, &
                                      GH_REAL
  use constants_mod,           only : r_def, i_def
  use fs_continuity_mod,       only : W2, W3, WTHETA 
  use kernel_mod,              only : kernel_type

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: map_u_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/             &
        arg_type(GH_FIELD,   GH_INC,  W2),          &
        arg_type(GH_FIELD,   GH_READ,  W3),         &
        arg_type(GH_FIELD,   GH_READ,  W3),         &
        arg_type(GH_FIELD,   GH_READ,  WTHETA),     &
        arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9)  &
        /)
    type(func_type) :: meta_funcs(4) = (/               &
        func_type(W2, GH_BASIS),                        &
        func_type(W3, GH_BASIS),                        &
        func_type(WTHETA, GH_BASIS),                    &
        func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS) &
        /)
    integer :: iterates_over = CELLS
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, public, nopass :: map_u_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! overload the default structure constructor for function space
  interface map_u_kernel_type
    module procedure map_u_kernel_constructor
  end interface

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
contains

type(map_u_kernel_type) function map_u_kernel_constructor() result(self)
  implicit none
  return
end function map_u_kernel_constructor

!> @brief Compute the right hand side to mapise the wind field.
!! @param[in] nlayers Number of layers
!! @param[inout] rhs Right hand side field to compute
!! @param[in] chi_1 X component of the coordinate field
!! @param[in] chi_2 Y component of the coordinate field
!! @param[in] chi_3 Z component of the coordinate field
!> @param[in] time Time (timestep multiplied by dt)
!! @param[in] ndf Number of degrees of freedom per cell
!! @param[in] undf Total number of degrees of freedom
!! @param[in] map Dofmap for the cell at the base of the column
!! @param[in] basis Basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_chi Number of dofs per cell for the coordinate field
!! @param[in] undf_chi Total number of degrees of freedom
!! @param[in] map_chi Dofmap for the coordinate field
!! @param[in] chi_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] chi_diff_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine map_u_code(nlayers, &
                      rhs, & 
                      u_lon, u_lat, u_up, & 
                      chi_1, chi_2, chi_3,  &
                      ndf_w2, undf_w2, &
                      map_w2, basis_w2, &
                      ndf_w3, undf_w3, map_w3, basis_w3, &
                      ndf_wth, undf_wth, map_wth, basis_wt, &
                      ndf_chi, undf_chi, &
                      map_chi, chi_basis, chi_diff_basis, &
                      nqp_h, nqp_v, wqp_h, wqp_v &
                      )

  use base_mesh_config_mod,       only : geometry, &
                                         geometry_spherical
  use coordinate_jacobian_mod,    only : coordinate_jacobian
  use coord_transform_mod,        only : sphere2cart_vector, xyz2llr

  implicit none

  !Arguments
  integer, intent(in) :: nlayers, ndf_w2, ndf_chi, ndf_wth, ndf_w3
  integer, intent(in) :: undf_w2, undf_chi, undf_wth, undf_w3
  integer, intent(in) :: nqp_h, nqp_v

  integer, dimension(ndf_w2),  intent(in) :: map_w2
  integer, dimension(ndf_chi), intent(in) :: map_chi
  integer, dimension(ndf_w3),  intent(in) :: map_w3
  integer, dimension(ndf_wth), intent(in) :: map_wth

  real(kind=r_def), intent(in), dimension(3,ndf_w2, nqp_h,nqp_v)  :: basis_w2 
  real(kind=r_def), intent(in), dimension(1,ndf_w3, nqp_h,nqp_v)  :: basis_w3 
  real(kind=r_def), intent(in), dimension(1,ndf_wth, nqp_h,nqp_v) :: basis_wt 
  real(kind=r_def), intent(in), dimension(3,ndf_chi,nqp_h,nqp_v)  :: chi_diff_basis
  real(kind=r_def), intent(in), dimension(1,ndf_chi,nqp_h,nqp_v)  :: chi_basis

  real(kind=r_def), dimension(undf_w2),  intent(inout) :: rhs
  real(kind=r_def), dimension(undf_w3),  intent(inout) :: u_lat, u_lon
  real(kind=r_def), dimension(undf_wth), intent(inout) :: u_up
  real(kind=r_def), dimension(undf_chi), intent(in)    :: chi_1, chi_2, chi_3

  real(kind=r_def), dimension(nqp_h), intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) ::  wqp_v

  !Internal variables
  integer(kind=i_def)                          :: df, k, qp1, qp2
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jacobian
  real(kind=r_def), dimension(ndf_chi)         :: chi_1_cell, chi_2_cell, chi_3_cell
  real(kind=r_def), dimension(3)               :: u_physical, u_spherical, xyz, llr
  real(kind=r_def)                             :: integrand

  do k = 0, nlayers-1
    do df = 1, ndf_chi
      chi_1_cell(df) = chi_1( map_chi(df) + k)
      chi_2_cell(df) = chi_2( map_chi(df) + k)
      chi_3_cell(df) = chi_3( map_chi(df) + k)
    end do
    call coordinate_jacobian(ndf_chi, &
                             nqp_h, &
                             nqp_v, &
                             chi_1_cell, &
                             chi_2_cell, &
                             chi_3_cell, &
                             chi_diff_basis, &
                             jacobian, &
                             dj)
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        ! Compute analytical vector wind in physical space
        xyz(:) = 0.0_r_def
        do df = 1, ndf_chi
          xyz(1) = xyz(1) + chi_1_cell(df)*chi_basis(1,df,qp1,qp2)
          xyz(2) = xyz(2) + chi_2_cell(df)*chi_basis(1,df,qp1,qp2)
          xyz(3) = xyz(3) + chi_3_cell(df)*chi_basis(1,df,qp1,qp2)
        end do

        if ( geometry == geometry_spherical ) then
          call xyz2llr(xyz(1), xyz(2), xyz(3), llr(1), llr(2), llr(3))
          u_spherical(1) = u_lon(map_w3(1)+k)*basis_w3(1,1,qp1,qp2)
          u_spherical(2) = u_lat(map_w3(1)+k)*basis_w3(1,1,qp1,qp2)
          u_spherical(3) = u_up(map_wth(1)+k)*basis_wt(1,1,qp1,qp2) &
                         + u_up(map_wth(2)+k)*basis_wt(1,2,qp1,qp2)
          u_physical     = sphere2cart_vector(u_spherical,llr) 
        else
          u_physical(1) = u_lon(map_w3(1)+k)*basis_w3(1,1,qp1,qp2)
          u_physical(2) = u_lat(map_w3(1)+k)*basis_w3(1,1,qp1,qp2)
          u_physical(3) = u_up(map_wth(1)+k)*basis_wt(1,1,qp1,qp2) &
                        + u_up(map_wth(2)+k)*basis_wt(1,2,qp1,qp2)
         end if

        do df = 1, ndf_w2 
          integrand = dot_product(matmul(jacobian(:,:,qp1,qp2),&
                                         basis_w2(:,df,qp1,qp2)),u_physical)
          rhs(map_w2(df) + k) = rhs(map_w2(df) + k) &
                               + wqp_h(qp1)*wqp_v(qp2)*integrand
        end do
      end do
    end do
  end do

end subroutine map_u_code

end module map_u_kernel_mod
