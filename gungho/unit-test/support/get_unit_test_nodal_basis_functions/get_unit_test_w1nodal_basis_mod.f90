!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module get_unit_test_w1nodal_basis_mod
! A module containing a collection of helper routines that provide canned 
! basis functions (and differential basis functions) evaluated on w1 nodal
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

  public :: get_w1_w1nodal_basis, &
            get_w0_w1nodal_diff_basis
            
  contains
  
!---------------------------------------------------------------------

  subroutine get_w1_w1nodal_basis(basis_w1) 
    ! Return the basis function for a field on a w1 function space
    ! evaluated on w1 nodal points 
    implicit none
    real(r_def), allocatable, intent(out) :: basis_w1(:,:,:)

    allocate(basis_w1(3,12,12))
    basis_w1 = reshape( [ &
        0.0_r_def,  1.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def,  1.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.5_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  1.0_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def, &
        0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.5_r_def,  0.0_r_def,  1.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.5_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  1.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  1.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def, &
        0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  1.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def,  0.5_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        1.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def,  1.0_r_def,  0.0_r_def, &
        0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def, &
        0.0_r_def,  1.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def, &
        0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  1.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def,  1.0_r_def, &
        0.0_r_def,  0.0_r_def],  [3,12,12] )
  end subroutine get_w1_w1nodal_basis

!---------------------------------------------------------------------

  subroutine get_w0_w1nodal_diff_basis(diff_basis_w0) 
    ! Return the diff basis function for a field on a w0 function space 
    ! evaluated on w1 nodal points 
    implicit none
    real(r_def), allocatable, intent(out) :: diff_basis_w0(:,:,:)

    allocate(diff_basis_w0(3,8,12))
    diff_basis_w0 = reshape( [ &
       -0.5_r_def, -1.0_r_def, -0.5_r_def,  0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def, -0.5_r_def, &
        1.0_r_def, -0.5_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def, -1.0_r_def, &
       -0.5_r_def, -0.5_r_def,  1.0_r_def, -0.5_r_def, -0.5_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def, -0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.5_r_def, -1.0_r_def, -0.5_r_def,  0.5_r_def, &
        1.0_r_def, -0.5_r_def, -0.5_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.5_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def, -0.5_r_def,  0.0_r_def, &
        0.0_r_def, -0.5_r_def,  0.0_r_def,  1.0_r_def,  0.5_r_def, &
       -0.5_r_def, -1.0_r_def,  0.5_r_def, -0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def, &
        0.5_r_def, -0.5_r_def, -0.5_r_def, -1.0_r_def,  0.5_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def, -0.5_r_def, -0.5_r_def, &
        1.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def, &
       -0.5_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def, -0.5_r_def, &
       -1.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def, -0.5_r_def,  0.0_r_def,  0.0_r_def, &
        0.5_r_def, -0.5_r_def,  1.0_r_def,  0.0_r_def,  0.5_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def, -0.5_r_def,  0.0_r_def, &
        0.5_r_def,  0.5_r_def, -1.0_r_def, -0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
       -0.5_r_def,  0.0_r_def,  0.5_r_def,  0.5_r_def,  1.0_r_def, &
       -0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, -0.5_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def, &
        0.0_r_def,  0.0_r_def, -0.5_r_def,  0.5_r_def, -1.0_r_def, &
        0.0_r_def, -0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def, -0.5_r_def, &
        0.5_r_def,  1.0_r_def,  0.0_r_def,  0.0_r_def, -0.5_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def, -0.5_r_def, -0.5_r_def, &
       -1.0_r_def,  0.5_r_def,  0.5_r_def,  0.0_r_def,  0.0_r_def, &
        0.5_r_def,  0.0_r_def,  0.0_r_def, -0.5_r_def,  1.0_r_def, &
        0.5_r_def,  0.0_r_def,  0.0_r_def, -0.5_r_def,  0.0_r_def, &
        0.0_r_def, -0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def, -1.0_r_def, -0.5_r_def, &
        0.5_r_def,  1.0_r_def, -0.5_r_def,  0.5_r_def,  0.0_r_def, &
        0.5_r_def,  0.0_r_def,  0.0_r_def,  0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
       -0.5_r_def,  0.0_r_def,  0.0_r_def, -0.5_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def, -0.5_r_def,  0.0_r_def,  0.0_r_def, &
        0.5_r_def, -1.0_r_def,  0.5_r_def,  0.5_r_def,  1.0_r_def, &
        0.5_r_def, -0.5_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def,  0.0_r_def, &
        0.0_r_def,  0.0_r_def, -0.5_r_def,  0.0_r_def,  0.0_r_def, &
       -0.5_r_def,  0.0_r_def, -0.5_r_def,  0.0_r_def,  0.0_r_def, &
       -0.5_r_def,  0.0_r_def,  1.0_r_def,  0.5_r_def,  0.5_r_def, &
       -1.0_r_def,  0.5_r_def,  0.5_r_def], [3,8,12] )
  end subroutine get_w0_w1nodal_diff_basis

!---------------------------------------------------------------------

end module get_unit_test_w1nodal_basis_mod
