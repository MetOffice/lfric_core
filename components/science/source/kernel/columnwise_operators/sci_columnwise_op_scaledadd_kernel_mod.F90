!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which adds a columnwise operator to another one.
!> @details Calculates C = alpha * A + beta * B.

module sci_columnwise_op_scaledadd_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type,                          &
                                    GH_COLUMNWISE_OPERATOR, GH_SCALAR, &
                                    GH_REAL, GH_READ, GH_WRITE,        &
                                    ANY_SPACE_1, ANY_SPACE_2,          &
                                    CELL_COLUMN

use constants_mod,           only : r_solver, i_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: columnwise_op_scaledadd_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                                 &
       arg_type(GH_COLUMNWISE_OPERATOR, GH_REAL, GH_READ,  ANY_SPACE_1, ANY_SPACE_2), &
       arg_type(GH_COLUMNWISE_OPERATOR, GH_REAL, GH_READ,  ANY_SPACE_1, ANY_SPACE_2), &
       arg_type(GH_COLUMNWISE_OPERATOR, GH_REAL, GH_WRITE, ANY_SPACE_1, ANY_SPACE_2), &
       arg_type(GH_SCALAR,              GH_REAL, GH_READ),                            &
       arg_type(GH_SCALAR,              GH_REAL, GH_READ)                             &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: columnwise_op_scaledadd_kernel_code
end type columnwise_op_scaledadd_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: columnwise_op_scaledadd_kernel_code

contains

  !> @brief The subroutine which is called directly from the PSy layer and
  !!        calculates operation \f$C = \alpha A + \beta B\f$.
  !>
  !> @param[in] cell the horizontal cell index
  !> @param[in] ncell_2d total number of cells in 2d grid
  !> @param[in] columnwise_matrix_A banded matrix \f$A\f$
  !> @param[in] nrow_A number of rows in the banded matrix A
  !> @param[in] ncol_A number of columns in the banded matrix A
  !> @param[in] bandwidth_A bandwidth of the banded matrix
  !> @param[in] alpha_A banded matrix parameter \f$\alpha\f$
  !> @param[in] beta_A banded matrix parameter \f$\beta\f$
  !> @param[in] gamma_m_A banded matrix parameter \f$\gamma_-\f$
  !> @param[in] gamma_p_A banded matrix parameter \f$\gamma_+\f$
  !> @param[in] columnwise_matrix_B banded matrix \f$B\f$
  !> @param[in] nrow_B number of rows in the banded matrix B
  !> @param[in] ncol_B number of columns in the banded matrix B
  !> @param[in] bandwidth_B bandwidth of the banded matrix
  !> @param[in] alpha_B banded matrix parameter \f$\alpha\f$
  !> @param[in] beta_B banded matrix parameter \f$\beta\f$
  !> @param[in] gamma_m_B banded matrix parameter \f$\gamma_-\f$
  !> @param[in] gamma_p_B banded matrix parameter \f$\gamma_+\f$
  !> @param[in,out] columnwise_matrix_C banded matrix \f$C\f$
  !> @param[in] nrow_C number of rows in the banded matrix C
  !> @param[in] ncol_C number of columns in the banded matrix C
  !> @param[in] bandwidth_C bandwidth of the banded matrix
  !> @param[in] alpha_C banded matrix parameter \f$\alpha\f$
  !> @param[in] beta_C banded matrix parameter \f$\beta\f$
  !> @param[in] gamma_m_C banded matrix parameter \f$\gamma_-\f$
  !> @param[in] gamma_p_C banded matrix parameter \f$\gamma_+\f$
  !> @param[in] alpha Scaling parameter \f$\alpha\f$
  !> @param[in] beta Scaling parameter \f$\beta\f$
  subroutine columnwise_op_scaledadd_kernel_code(cell,                      &
                                                 ncell_2d,                  &
                                                 columnwise_matrix_A,       &
                                                 nrow_A,                    &
                                                 ncol_A,                    &
                                                 bandwidth_A,               &
                                                 alpha_A,                   &
                                                 beta_A,                    &
                                                 gamma_m_A,                 &
                                                 gamma_p_A,                 &
                                                 columnwise_matrix_B,       &
                                                 nrow_B,                    &
                                                 ncol_B,                    &
                                                 bandwidth_B,               &
                                                 alpha_B,                   &
                                                 beta_B,                    &
                                                 gamma_m_B,                 &
                                                 gamma_p_B,                 &
                                                 columnwise_matrix_C,       &
                                                 nrow_C,                    &
                                                 ncol_C,                    &
                                                 bandwidth_C,               &
                                                 alpha_C,                   &
                                                 beta_C,                    &
                                                 gamma_m_C,                 &
                                                 gamma_p_C,                 &
                                                 alpha,                     &
                                                 beta)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: cell, ncell_2d
    integer(kind=i_def), intent(in) :: nrow_A, ncol_A
    integer(kind=i_def), intent(in) :: nrow_B, ncol_B
    integer(kind=i_def), intent(in) :: nrow_C, ncol_C
    integer(kind=i_def), intent(in) :: bandwidth_A, bandwidth_B, bandwidth_C
    real(kind=r_solver), dimension(bandwidth_A,nrow_A,ncell_2d), intent(in)    :: columnwise_matrix_A
    real(kind=r_solver), dimension(bandwidth_B,nrow_B,ncell_2d), intent(in)    :: columnwise_matrix_B
    real(kind=r_solver), dimension(bandwidth_C,nrow_C,ncell_2d), intent(inout) :: columnwise_matrix_C

    integer(kind=i_def), intent(in) :: alpha_A, beta_A, gamma_m_A, gamma_p_A
    integer(kind=i_def), intent(in) :: alpha_B, beta_B, gamma_m_B, gamma_p_B
    integer(kind=i_def), intent(in) :: alpha_C, beta_C, gamma_m_C, gamma_p_C
    real(kind=r_solver), intent(in) :: alpha, beta

    ! Internal parameters
    integer(kind=i_def) :: i, j ! Row and column index index
    ! Smallest index in a particular row
    integer(kind=i_def) :: j_minus_A, j_minus_B, j_minus_C, j_plus_A, j_plus_B

    ! Number of rows (same for all matrices)
    integer(kind=i_def) :: nrow

    ! The variables bandwidth_B, gamma_m_B and ncol_B are not used in the code
    ! below and will trigger a compiler warning, which will abort the
    ! compilation. To avoid this, the following line uses those variables, but
    ! the result of the computation is irrelevant.

    i = bandwidth_B + gamma_m_B + ncol_B
    nrow = nrow_A
    columnwise_matrix_C(:,:,cell) = 0.0_r_solver
    ! Add matrix A
    do i = 1, nrow
       j_minus_A = ceiling((alpha_A*i-gamma_p_A)/(1.0_r_solver*beta_A), i_def)
       j_plus_A = floor((alpha_A*i+gamma_m_A)/(1.0_r_solver*beta_A), i_def)
       j_minus_B = ceiling((alpha_B*i-gamma_p_B)/(1.0_r_solver*beta_B), i_def)
       j_plus_B = floor((alpha_B*i+gamma_m_B)/(1.0_r_solver*beta_B), i_def)
       j_minus_C = ceiling((alpha_C*i-gamma_p_C)/(1.0_r_solver*beta_C), i_def)
       do j = MAX(1,j_minus_A), MIN(ncol_A,j_plus_A)
          columnwise_matrix_C(j-j_minus_C+1,i,cell)             &
            = columnwise_matrix_C(j-j_minus_C+1,i,cell)         &
            + alpha * columnwise_matrix_A(j-j_minus_A+1,i,cell)
       end do
       do j = MAX(1,j_minus_B), MIN(ncol_A,j_plus_B)
          columnwise_matrix_C(j-j_minus_C+1,i,cell)             &
            = columnwise_matrix_C(j-j_minus_C+1,i,cell)         &
            + beta * columnwise_matrix_B(j-j_minus_B+1,i,cell)
       end do
    end do

  end subroutine columnwise_op_scaledadd_kernel_code

end module sci_columnwise_op_scaledadd_kernel_mod
