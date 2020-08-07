!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module get_unit_test_w3nodal_basis_mod
! A module containing a collection of helper routines that provide canned
! basis functions (and differential basis functions) evaluated on w3 nodal
! points for use when writing unit tests of kernels.
!
! Not all the combinations of (diff) basis functions have been added. There are
! many of them and only a few will ever be used, so they are being added as
! necessary.
!
! If you require a (diff) basis function that you think should be in this file
! but isn't, then feel free to add it.

  use constants_mod, only : r_def

  implicit none

  private

  public :: get_w0_w3nodal_diff_basis,     &
            get_wtheta_w3nodal_basis,      &
            get_w0_w3nodal_basis

  contains

!---------------------------------------------------------------------

  subroutine get_w0_w3nodal_diff_basis(diff_basis_w0)
    ! Return the diff basis function for a field on a w0 function space
    ! evaluated on w3 nodal points
    implicit none
    real(r_def), allocatable, intent(out) :: diff_basis_w0(:,:,:)

    allocate(diff_basis_w0(3,8,1))
    diff_basis_w0 = reshape( [-0.25_r_def, -0.25_r_def, -0.25_r_def, &
                               0.25_r_def, -0.25_r_def, -0.25_r_def, &
                               0.25_r_def,  0.25_r_def, -0.25_r_def, &
                              -0.25_r_def,  0.25_r_def, -0.25_r_def, &
                              -0.25_r_def, -0.25_r_def,  0.25_r_def, &
                               0.25_r_def, -0.25_r_def,  0.25_r_def, &
                               0.25_r_def,  0.25_r_def,  0.25_r_def, &
                              -0.25_r_def,  0.25_r_def,  0.25_r_def], [3,8,1] )
  end subroutine get_w0_w3nodal_diff_basis

!---------------------------------------------------------------------

  subroutine get_wtheta_w3nodal_basis(basis_wtheta)
    ! Return the basis function for a field on a wtheta function space
    ! evaluated on w3 nodal points

    implicit none

    real(r_def), allocatable :: basis_wtheta(:,:,:)

    ! Lowest order scalar. 2 Wtheta basis functions. 1 W3 nodal points.
    allocate(basis_wtheta(1,2,1))
    basis_wtheta(:, 1, 1) =  0.50_r_def
    basis_wtheta(:, 2, 1) =  0.50_r_def

  end subroutine get_wtheta_w3nodal_basis

!---------------------------------------------------------------------

  subroutine get_w0_w3nodal_basis(basis_w0)
    ! Return the basis function for a field on w0 function space
    ! evaluated on w3 nodal points

    implicit none

    real(r_def), allocatable :: basis_w0(:,:,:)
    allocate(basis_w0(1,8,1))
    basis_w0(:, :, :) =  0.125_r_def

  end subroutine get_w0_w3nodal_basis

!---------------------------------------------------------------------

end module get_unit_test_w3nodal_basis_mod
