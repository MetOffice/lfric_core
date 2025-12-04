!-----------------------------------------------------------------------------
! Copyright (c) 2019,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

module get_unit_test_wthetanodal_basis_mod
! A module containing helper routines that provide canned basis functions
! for a range of function spaces computed at the location of the w2 nodal
! points.

  use constants_mod,                 only : i_def, r_def

  implicit none

  private

  public :: get_w2_wthetanodal_basis,     &
            get_w0_wthetanodal_basis,     &
            get_w0_wthetanodal_diff_basis

contains

  subroutine get_w2_wthetanodal_basis(basis_w2)
    ! Provide w2 basis functions computed at wtheta nodal points

    implicit none

    real(r_def), allocatable :: basis_w2(:,:,:)

    ! Lowest order vector. 6 W2 basis functions. 2 Wtheta nodal points.
    allocate(basis_w2(3,6,2))
    basis_w2(:, 1, 1) = (/ 0.5_r_def,  0.0_r_def,  0.0_r_def /)
    basis_w2(:, 2, 1) = (/ 0.0_r_def, -0.5_r_def,  0.0_r_def /)
    basis_w2(:, 3, 1) = (/ 0.5_r_def,  0.0_r_def,  0.0_r_def /)
    basis_w2(:, 4, 1) = (/ 0.0_r_def, -0.5_r_def,  0.0_r_def /)
    basis_w2(:, 5, 1) = (/ 0.0_r_def,  0.0_r_def,  1.0_r_def /)
    basis_w2(:, 6, 1) = (/ 0.0_r_def,  0.0_r_def,  0.0_r_def /)
    basis_w2(:, 1, 2) = (/ 0.5_r_def,  0.0_r_def,  0.0_r_def /)
    basis_w2(:, 2, 2) = (/ 0.0_r_def, -0.5_r_def,  0.0_r_def /)
    basis_w2(:, 3, 2) = (/ 0.5_r_def,  0.0_r_def,  0.0_r_def /)
    basis_w2(:, 4, 2) = (/ 0.0_r_def, -0.5_r_def,  0.0_r_def /)
    basis_w2(:, 5, 2) = (/-0.0_r_def, -0.0_r_def, -0.0_r_def /)
    basis_w2(:, 6, 2) = (/ 0.0_r_def,  0.0_r_def,  1.0_r_def /)

  end subroutine get_w2_wthetanodal_basis

  subroutine get_w0_wthetanodal_basis(basis_w0)
    ! Provide w0 basis functions computed at wtheta nodal points

    implicit none

    real(r_def), allocatable :: basis_w0(:,:,:)

    ! Six W2 basis functions. Two Wtheta nodal points.
    allocate(basis_w0(1,8,2))
    basis_w0(:, 1, 1) =  0.25_r_def
    basis_w0(:, 2, 1) =  0.25_r_def
    basis_w0(:, 3, 1) =  0.25_r_def
    basis_w0(:, 4, 1) =  0.25_r_def
    basis_w0(:, 5, 1) =  0.00_r_def
    basis_w0(:, 6, 1) =  0.00_r_def
    basis_w0(:, 7, 1) =  0.00_r_def
    basis_w0(:, 8, 1) =  0.00_r_def
    basis_w0(:, 1, 2) = -0.00_r_def
    basis_w0(:, 2, 2) = -0.00_r_def
    basis_w0(:, 3, 2) = -0.00_r_def
    basis_w0(:, 4, 2) = -0.00_r_def
    basis_w0(:, 5, 2) =  0.25_r_def
    basis_w0(:, 6, 2) =  0.25_r_def
    basis_w0(:, 7, 2) =  0.25_r_def
    basis_w0(:, 8, 2) =  0.25_r_def
  end subroutine get_w0_wthetanodal_basis

  subroutine get_w0_wthetanodal_diff_basis(diff_basis_w0)
    ! Provide w0 diff basis functions computed at wtheta nodal points

    implicit none

    real(r_def), allocatable :: diff_basis_w0(:,:,:)

    ! Lowest order vector. 8 W0 basis functions. 2 Wtheta points.
    allocate(diff_basis_w0(3,8,2))
    diff_basis_w0(:, 1, 1) = (/ -0.5_r_def, -0.50_r_def, -0.25_r_def /)
    diff_basis_w0(:, 2, 1) = (/  0.5_r_def, -0.50_r_def, -0.25_r_def /)
    diff_basis_w0(:, 3, 1) = (/  0.5_r_def,  0.50_r_def, -0.25_r_def /)
    diff_basis_w0(:, 4, 1) = (/ -0.5_r_def,  0.50_r_def, -0.25_r_def /)
    diff_basis_w0(:, 5, 1) = (/ -0.0_r_def, -0.00_r_def,  0.25_r_def /)
    diff_basis_w0(:, 6, 1) = (/  0.0_r_def, -0.00_r_def,  0.25_r_def /)
    diff_basis_w0(:, 7, 1) = (/  0.0_r_def,  0.00_r_def,  0.25_r_def /)
    diff_basis_w0(:, 8, 1) = (/ -0.0_r_def,  0.00_r_def,  0.25_r_def /)
    diff_basis_w0(:, 1, 2) = (/  0.0_r_def,  0.00_r_def, -0.25_r_def /)
    diff_basis_w0(:, 2, 2) = (/ -0.0_r_def,  0.00_r_def, -0.25_r_def /)
    diff_basis_w0(:, 3, 2) = (/ -0.0_r_def, -0.00_r_def, -0.25_r_def /)
    diff_basis_w0(:, 4, 2) = (/  0.0_r_def, -0.00_r_def, -0.25_r_def /)
    diff_basis_w0(:, 5, 2) = (/ -0.5_r_def, -0.50_r_def,  0.25_r_def /)
    diff_basis_w0(:, 6, 2) = (/  0.5_r_def, -0.50_r_def,  0.25_r_def /)
    diff_basis_w0(:, 7, 2) = (/  0.5_r_def,  0.50_r_def,  0.25_r_def /)
    diff_basis_w0(:, 8, 2) = (/ -0.5_r_def,  0.50_r_def,  0.25_r_def /)
  end subroutine get_w0_wthetanodal_diff_basis


end module get_unit_test_wthetanodal_basis_mod
