!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Module for computing the Jacobian matrix, its determinant and
!> inverse for a coordinate field
module coordinate_jacobian_mod

use constants_mod, only: r_def

implicit none

contains
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> @brief Subroutine Computes the element Jacobian of the coordinate transform from
!! reference space \f[ \hat{\chi} \f] to physical space chi assuming the coordinate
!! field is in W0
!> @details Compute the Jacobian of the coordinate transform from
!> reference space \f[ \hat{\chi} \f] to physical space \f[ \chi \f]
!> \f[ J^{i,j} = \frac{\partial \chi_i} / {\partial \hat{\chi_j}} \f] 
!> and the derterminant det(J)
!! @param[in] ndf        Size of the chi arrays
!! @param[in] ngp_h      Number of quadrature points in horizontal direction
!! @param[in] ngp_v      Number of quadrature points in vertical direction
!! @param[in] chi_1      Coordinate field
!! @param[in] chi_2      Coordinate field
!! @param[in] chi_3      Coordinate field
!! @param[in] diff_basis Grad of W0 basis functions
!! @param[out] jac       Jacobian on quadrature points
!! @param[out] dj        Determinant of the Jacobian on quadrature points
subroutine coordinate_jacobian(ndf, ngp_h, ngp_v, chi_1, chi_2, chi_3, diff_basis, jac, dj)
!-------------------------------------------------------------------------------
! Compute the Jacobian J^{i,j} = d chi_i / d \hat{chi_j} and the 
! derterminant det(J)
!-------------------------------------------------------------------------------

integer,          intent(in)  :: ndf, ngp_h, ngp_v
real(kind=r_def), intent(in)  :: chi_1(ndf), chi_2(ndf), chi_3(ndf)
real(kind=r_def), intent(in)  :: diff_basis(3,ndf,ngp_h,ngp_v)
real(kind=r_def), intent(out) :: jac(3,3,ngp_h,ngp_v)
real(kind=r_def), intent(out) :: dj(ngp_h,ngp_v)

! Hardwired values for cartesian domain
!real(kind=r_def) :: dx = 6000.0_r_def, &
!                    dy = 2000.0_r_def, &
!                    dz = 2000.0_r_def

integer :: i, j, df, dir

do j = 1,ngp_v
  do i = 1,ngp_h
    jac(:,:,i,j) = 0.0_r_def
    do df = 1,ndf
      do dir = 1,3
        jac(1,dir,i,j) = jac(1,dir,i,j) + chi_1(df)*diff_basis(dir,df,i,j)
        jac(2,dir,i,j) = jac(2,dir,i,j) + chi_2(df)*diff_basis(dir,df,i,j)
        jac(3,dir,i,j) = jac(3,dir,i,j) + chi_3(df)*diff_basis(dir,df,i,j)
      end do
    end do

! Hard wired values for cartesian biperiodic domain this needs correcting
!    jac(:,:,i,j) = 0.0_r_def
!    jac(1,1,i,j) = dx
!    jac(2,2,i,j) = dy
!    jac(3,3,i,j) = dz

    dj(i,j) = jac(1,1,i,j)*(jac(2,2,i,j)*jac(3,3,i,j)        &
                          - jac(2,3,i,j)*jac(3,2,i,j))       &
            - jac(1,2,i,j)*(jac(2,1,i,j)*jac(3,3,i,j)        &
                          - jac(2,3,i,j)*jac(3,1,i,j))       &
            + jac(1,3,i,j)*(jac(2,1,i,j)*jac(3,2,i,j)        &
                          - jac(2,2,i,j)*jac(3,1,i,j)) 
  end do
end do

end subroutine coordinate_jacobian

!> @brief Subroutine Computes the inverse of the Jacobian of the coordinate transform from
!! reference space \f[\hat{\chi}\f] to physical space \f[ \chi \f] 
!> @details Compute the inverse of the Jacobian 
!> \f[ J^{i,j} = \frac{\partial \chi_i} / {\partial \hat{\chi_j}} \f] 
!> and the derterminant det(J)
!! @param[in] ngp_h      Number of quadrature points in horizontal direction
!! @param[in] ngp_v      Number of quadrature points in vertical direction
!! @param[in] jac        Jacobian on quadrature points
!! @param[in] dj         Determinant of the Jacobian
!! @param[out] jac_inv   Inverse of the Jacobian on quadrature points
subroutine coordinate_jacobian_inverse(ngp_h, ngp_v, jac, dj, jac_inv)

use matrix_invert_mod, only: matrix_invert_3x3

implicit none

integer,          intent(in)  :: ngp_h, ngp_v
real(kind=r_def), intent(in)  :: jac(3,3,ngp_h,ngp_v)
real(kind=r_def), intent(in)  :: dj(ngp_h,ngp_v)
real(kind=r_def), intent(out) :: jac_inv(3,3,ngp_h,ngp_v)

real(r_def) :: dummy
integer :: i, k

!> @todo This is here to maintain the API. If it turns out we don't want this it should
!> be removed.
dummy = dj(1,1)

! Calculates the inverse Jacobian from the analytic inversion formula
do k = 1,ngp_v
   do i = 1,ngp_h
     jac_inv(:,:,i,k) = matrix_invert_3x3(jac(:,:,i,k))
   end do
end do

end subroutine coordinate_jacobian_inverse

end module coordinate_jacobian_mod
