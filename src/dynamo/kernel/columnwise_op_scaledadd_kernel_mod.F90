!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which adds a columnwise operator to another one
!> @detail calculates op_B -> op_B + omega * op_A

module columnwise_op_scaledadd_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type,                               &
                                    GH_COLUMNWISE_OPERATOR,                 &
                                    GH_READ, GH_INC,  GH_REAL,              &
                                    ANY_SPACE_1, ANY_SPACE_2,               &
                                    CELLS 

use constants_mod,           only : r_def, i_def

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: columnwise_op_scaledadd_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                       &
       arg_type(GH_COLUMNWISE_OPERATOR, GH_READ, ANY_SPACE_1, ANY_SPACE_2), &
       arg_type(GH_COLUMNWISE_OPERATOR, GH_INC, ANY_SPACE_1, ANY_SPACE_2),  &
       arg_type(GH_REAL, GH_READ)                                           &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: columnwise_op_scaledadd_kernel_code
end type columnwise_op_scaledadd_kernel_type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface columnwise_op_scaledadd_kernel_type
   module procedure columnwise_op_scaledadd_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public columnwise_op_scaledadd_kernel_code
contains
  
  type(columnwise_op_scaledadd_kernel_type) function columnwise_op_scaledadd_kernel_constructor() result(self)
    implicit none
    return
  end function columnwise_op_scaledadd_kernel_constructor

  !> @brief The subroutine which is called directly from the PSY layer and
  !> applies the operator as lhs += A.x
  !>
  !> @param [in] cell the horizontal cell index
  !> @param [in] ncell_2d total number of cells in 2d grid
  !> @param [in] columnwise_matrix_A banded matrix op_A
  !> @param [inout] columnwise_matrix_B banded matrix op_B
  !> @param [in] nrow_A number of rows in the banded matrix A
  !> @param [in] ncol_A number of columns in the banded matrix A
  !> @param [in] bandwidth_A bandwidth of the banded matrix
  !> @param [in] alpha_A banded matrix parameter \f$\alpha\f$
  !> @param [in] beta_A banded matrix parameter \f$\beta\f$
  !> @param [in] gamma_m_A banded matrix parameter \f$\gamma_-\f$
  !> @param [in] gamma_p_A banded matrix parameter \f$\gamma_+\f$
  !> @param [in] nrow_B number of rows in the banded matrix B
  !> @param [in] ncol_B number of columns in the banded matrix B
  !> @param [in] bandwidth_B bandwidth of the banded matrix
  !> @param [in] alpha_B banded matrix parameter \f$\alpha\f$
  !> @param [in] beta_B banded matrix parameter \f$\beta\f$
  !> @param [in] gamma_m_B banded matrix parameter \f$\gamma_-\f$
  !> @param [in] gamma_p_B banded matrix parameter \f$\gamma_+\f$
  !> @param [in] omega Scaling parameter
  subroutine columnwise_op_scaledadd_kernel_code(cell,                      &
                                                 ncell_2d,                  &
                                                 columnwise_matrix_A,       &
                                                 columnwise_matrix_B,       &
                                                 nrow_A,                    &
                                                 ncol_A,                    &
                                                 bandwidth_A,               &
                                                 alpha_A,                   &
                                                 beta_A,                    &
                                                 gamma_m_A,                 &
                                                 gamma_p_A,                 &
                                                 nrow_B,                    &
                                                 ncol_B,                    &
                                                 bandwidth_B,               &
                                                 alpha_B,                   &
                                                 beta_B,                    &
                                                 gamma_m_B,                 &
                                                 gamma_p_B,                 &
                                                 omega)
    implicit none
    
    ! Arguments
    integer(kind=i_def), intent(in) :: cell,  ncell_2d
    integer(kind=i_def), intent(in) :: nrow_A, ncol_A
    integer(kind=i_def), intent(in) :: nrow_B, ncol_B
    integer(kind=i_def), intent(in) :: bandwidth_A, bandwidth_B
    real(kind=r_def), dimension(bandwidth_A,nrow_A,ncell_2d), intent(in) :: columnwise_matrix_A
    real(kind=r_def), dimension(bandwidth_B,nrow_B,ncell_2d), intent(inout) :: columnwise_matrix_B

    integer(kind=i_def), intent(in) :: alpha_A, beta_A, gamma_m_A, gamma_p_A
    integer(kind=i_def), intent(in) :: alpha_B, beta_B, gamma_m_B, gamma_p_B
    real(kind=r_def), intent(in) :: omega

    ! Internal parameters
    integer(kind=i_def) :: i,j ! Row and column index index
    ! Smallest index in a particular row
    integer(kind=i_def) :: j_minus_A, j_minus_B, j_plus_A

    ! The variables bandwidth_B, gamma_m_B and ncol_B are not used in the code
    ! below and will trigger a compiler warning, which will abort the
    ! compilation. To avoid this, the following line uses those variables, but
    ! the result of the computation is irrelevant.
    i = bandwidth_B + gamma_m_B + ncol_B

    do i=1, nrow_A
       j_minus_A = ceiling((alpha_A*i-gamma_p_A)/(1.0_r_def*beta_A),i_def)
       j_plus_A = floor((alpha_A*i+gamma_m_A)/(1.0_r_def*beta_A),i_def)
       j_minus_B = ceiling((alpha_B*i-gamma_p_B)/(1.0_r_def*beta_B),i_def)
       do j=MAX(1,j_minus_A), MIN(ncol_A,j_plus_A)
          columnwise_matrix_B(j-j_minus_B+1,i,cell)             &
            = columnwise_matrix_B(j-j_minus_B+1,i,cell)         &
            + omega * columnwise_matrix_A(j-j_minus_A+1,i,cell)
       end do
    end do

  end subroutine columnwise_op_scaledadd_kernel_code

end module columnwise_op_scaledadd_kernel_mod
