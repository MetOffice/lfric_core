!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Module for computing the Jacobian matrix, its determinant and
!> inverse for a coordinate field
module coordinate_jacobian_mod

use constants_mod, only: r_def, i_def

implicit none

private

public :: coordinate_jacobian
public :: coordinate_jacobian_inverse
public :: pointwise_coordinate_jacobian
public :: pointwise_coordinate_jacobian_inverse

contains

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  !> @brief Subroutine Computes the element Jacobian of the coordinate transform from
  !! reference space \f$ \hat{\chi} \f$ to physical space chi
  !> @details Compute the Jacobian of the coordinate transform from
  !> reference space \f[ \hat{\chi} \f] to physical space \f[ \chi \f]
  !> \f[ J^{i,j} = \frac{\partial \chi_i} / {\partial \hat{\chi_j}} \f]
  !> and the determinant det(J)
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
  ! determinant det(J)
  !-------------------------------------------------------------------------------
    implicit none

    integer(kind=i_def), intent(in)  :: ndf, ngp_h, ngp_v

    real(kind=r_def), intent(in)  :: chi_1(ndf), chi_2(ndf), chi_3(ndf)
    real(kind=r_def), intent(in)  :: diff_basis(3,ndf,ngp_h,ngp_v)
    real(kind=r_def), intent(out) :: jac(3,3,ngp_h,ngp_v)
    real(kind=r_def), intent(out) :: dj(ngp_h,ngp_v)

    real(kind=r_def) :: dx, dy, dz

    integer(kind=i_def) :: i, j, df, dir

    jac(:,:,:,:) = 0.0_r_def
    do j = 1,ngp_v
      do i = 1,ngp_h
        do df = 1,ndf
          do dir = 1,3
            jac(1,dir,i,j) = jac(1,dir,i,j) + chi_1(df)*diff_basis(dir,df,i,j)
            jac(2,dir,i,j) = jac(2,dir,i,j) + chi_2(df)*diff_basis(dir,df,i,j)
            jac(3,dir,i,j) = jac(3,dir,i,j) + chi_3(df)*diff_basis(dir,df,i,j)
          end do
        end do

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
  !! reference space \f$\hat{\chi}\f$ to physical space \f$ \chi \f$
  !> @details Compute the inverse of the Jacobian
  !> \f[ J^{i,j} = \frac{\partial \chi_i} / {\partial \hat{\chi_j}} \f]
  !> and the determinant det(J)
  !! @param[in] ngp_h      Number of quadrature points in horizontal direction
  !! @param[in] ngp_v      Number of quadrature points in vertical direction
  !! @param[in] jac        Jacobian on quadrature points
  !! @param[in] dj         Determinant of the Jacobian
  !! @param[out] jac_inv   Inverse of the Jacobian on quadrature points
  subroutine coordinate_jacobian_inverse(ngp_h, ngp_v, jac, dj, jac_inv)

    use matrix_invert_mod, only: matrix_invert_3x3

    implicit none

    integer(kind=i_def), intent(in)  :: ngp_h, ngp_v

    real(kind=r_def), intent(in)  :: jac(3,3,ngp_h,ngp_v)
    real(kind=r_def), intent(in)  :: dj(ngp_h,ngp_v)
    real(kind=r_def), intent(out) :: jac_inv(3,3,ngp_h,ngp_v)

    real(kind=r_def) :: dummy
    integer(kind=i_def) :: i, k

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

  !> @brief Subroutine Computes the element Jacobian of the coordinate transform from
  !! reference space \f$ \hat{\chi} \f$ to physical space chi for a single point
  !> @details Compute the Jacobian of the coordinate transform from
  !> reference space \f[ \hat{\chi} \f] to physical space \f[ \chi \f]
  !> \f[ J^{i,j} = \frac{\partial \chi_i} / {\partial \hat{\chi_j}} \f]
  !> and the determinant det(J) for a single point
  !! @param[in] ndf        Size of the chi arrays
  !! @param[in] chi_1      Coordinate field
  !! @param[in] chi_2      Coordinate field
  !! @param[in] chi_3      Coordinate field
  !! @param[in] diff_basis Grad of W0 basis functions
  !! @param[out] jac       Jacobian on quadrature points
  !! @param[out] dj        Determinant of the Jacobian on quadrature points
  subroutine pointwise_coordinate_jacobian(ndf, chi_1, chi_2, chi_3, diff_basis, jac, dj)

    implicit none

    integer(kind=i_def), intent(in)  :: ndf

    real(kind=r_def), intent(in)  :: chi_1(ndf), chi_2(ndf), chi_3(ndf)
    real(kind=r_def), intent(in)  :: diff_basis(3,ndf)
    real(kind=r_def), intent(out) :: jac(3,3)
    real(kind=r_def), intent(out) :: dj

    integer(kind=i_def) :: df, dir

    real(kind=r_def) :: dx, dy, dz

    jac(:,:) = 0.0_r_def
    do df = 1,ndf
      do dir = 1,3
        jac(1,dir) = jac(1,dir) + chi_1(df)*diff_basis(dir,df)
        jac(2,dir) = jac(2,dir) + chi_2(df)*diff_basis(dir,df)
        jac(3,dir) = jac(3,dir) + chi_3(df)*diff_basis(dir,df)
      end do
    end do

    dj = jac(1,1)*(jac(2,2)*jac(3,3)        &
                 - jac(2,3)*jac(3,2))       &
       - jac(1,2)*(jac(2,1)*jac(3,3)        &
                 - jac(2,3)*jac(3,1))       &
       + jac(1,3)*(jac(2,1)*jac(3,2)        &
                 - jac(2,2)*jac(3,1))

  end subroutine pointwise_coordinate_jacobian

  !> @brief Subroutine Computes the inverse of the Jacobian of the coordinate transform from
  !! reference space \f$\hat{\chi}\f$ to physical space \f$ \chi \f$ for a
  !! single point
  !> @details Compute the inverse of the Jacobian
  !> \f[ J^{i,j} = \frac{\partial \chi_i} / {\partial \hat{\chi_j}} \f]
  !> and the determinant det(J)
  !! @param[in] jac        Jacobian on quadrature points
  !! @param[in] dj         Determinant of the Jacobian
  !! @result    jac_inv    Inverse of the Jacobian on quadrature points
  function pointwise_coordinate_jacobian_inverse(jac, dj) result(jac_inv)
    implicit none

    real(kind=r_def)              :: jac_inv(3,3)
    real(kind=r_def), intent(in)  :: jac(3,3)
    real(kind=r_def), intent(in)  :: dj

    real(kind=r_def) :: idj

    idj = 1.0_r_def/dj

    jac_inv(1,1) =  (jac(2,2)*jac(3,3) - jac(2,3)*jac(3,2))*idj
    jac_inv(1,2) = -(jac(1,2)*jac(3,3) - jac(1,3)*jac(3,2))*idj
    jac_inv(1,3) =  (jac(1,2)*jac(2,3) - jac(1,3)*jac(2,2))*idj
    jac_inv(2,1) = -(jac(2,1)*jac(3,3) - jac(2,3)*jac(3,1))*idj
    jac_inv(2,2) =  (jac(1,1)*jac(3,3) - jac(1,3)*jac(3,1))*idj
    jac_inv(2,3) = -(jac(1,1)*jac(2,3) - jac(1,3)*jac(2,1))*idj
    jac_inv(3,1) =  (jac(2,1)*jac(3,2) - jac(2,2)*jac(3,1))*idj
    jac_inv(3,2) = -(jac(1,1)*jac(3,2) - jac(1,2)*jac(3,1))*idj
    jac_inv(3,3) =  (jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1))*idj

  end function pointwise_coordinate_jacobian_inverse

end module coordinate_jacobian_mod
