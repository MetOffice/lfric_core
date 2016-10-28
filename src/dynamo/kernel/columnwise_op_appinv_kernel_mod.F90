!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which applies the inverse of a columnwise assembled operator
!>
!> @detail Given a field \f$x\f$, this calculates \f$y+=A^{-1}x\f$, i.e.
!> solves \f$A\delta y=x\f$ for \f$\delta y\f$ and adds this to \f$y\f$.
!> Only works if the assembled matrix is square and tridiagonal, and the 
!> Thomas algorithm can be used. The PSY layer should check whether the
!> matrix has the correct properties.
!>
!> Reference for Thomas algorithm:
!>   Numerical Recipes, Second Edition, CAMBRIDGE UNIVERSITY PRESS (1992)

module columnwise_op_appinv_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,                    &
                                    GH_FIELD, GH_COLUMNWISE_OPERATOR,       &
                                    GH_READ, GH_INC,                        &
                                    ANY_SPACE_1, ANY_SPACE_2,               &
                                    GH_COLUMN_INDIRECTION_DOFMAP,           &
                                    CELLS 

use constants_mod,           only : r_def, i_def

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: columnwise_op_appinv_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                      &
       arg_type(GH_FIELD,    GH_INC,  ANY_SPACE_1),                        &  
       arg_type(GH_FIELD,    GH_READ, ANY_SPACE_2),                        &
       arg_type(GH_COLUMNWISE_OPERATOR, GH_READ, ANY_SPACE_1, ANY_SPACE_2) &
       /)
  type(func_type) :: meta_funcs(2) = (/                                    &
       func_type(ANY_SPACE_1, GH_COLUMN_INDIRECTION_DOFMAP),               &
       func_type(ANY_SPACE_2, GH_COLUMN_INDIRECTION_DOFMAP)                &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: columnwise_op_appinv_kernel_code
end type columnwise_op_appinv_kernel_type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface columnwise_op_appinv_kernel_type
   module procedure columnwise_op_appinv_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public columnwise_op_appinv_kernel_code
contains
  
  type(columnwise_op_appinv_kernel_type) function columnwise_op_appinv_kernel_constructor() result(self)
    implicit none
    return
  end function columnwise_op_appinv_kernel_constructor

  !> @brief The subroutine which is called directly from the PSY layer and
  !> applies the operator as lhs += A^{-1}.x
  !>
  !> @param [in] cell the horizontal cell index
  !> @param [in] ncell_2d number of cells in 2d grid
  !> @param [inout] lhs Resulting field lhs += A^{-1}.x
  !> @param [in] x input field
  !> @param [in] columnwise_matrix banded matrix to assemble into
  !> @param [in] ndf1 number of degrees of freedom per cell for the to-space
  !> @param [in] undf1 unique number of degrees of freedom  for the to-space
  !> @param [in] map1 dofmap for the to-space
  !> @param [in] ndf2 number of degrees of freedom per cell for the from-space
  !> @param [in] undf2 unique number of degrees of freedom for the from-space 
  !> @param [in] map2 dofmap for the from-space
  !> @param [in] nrow number of rows in the banded matrix
  !> @param [in] ncol number of columns in the banded matrix
  !> @param [in] bandwidth bandwidth of the banded matrix
  !> @param [in] alpha banded matrix parameter \f$\alpha\f$
  !> @param [in] beta banded matrix parameter \f$\beta\f$
  !> @param [in] gamma_m banded matrix parameter \f$\gamma_-\f$
  !> @param [in] gamma_p banded matrix parameter \f$\gamma_+\f$
  !> @param [in] indirection_dofmap_to indirection map for to-space
  !> @param [in] indirection_dofmap_from indirection map for from-space
  subroutine columnwise_op_appinv_kernel_code(cell,              &
                                              ncell_2d,          &
                                              lhs, x,            & 
                                              columnwise_matrix, &
                                              ndf1, undf1, map1, &
                                              ndf2, undf2, map2, &
                                              nrow,              &
                                              ncol,              &
                                              bandwidth,         &
                                              alpha,             &
                                              beta,              &
                                              gamma_m,           &
                                              gamma_p,           &
                                              indirection_dofmap_to, &
                                              indirection_dofmap_from)
    implicit none
    
    ! Arguments
    integer(kind=i_def), intent(in) :: cell, ncell_2d
    integer(kind=i_def), intent(in) :: nrow, ncol, bandwidth
    integer(kind=i_def), intent(in) :: undf1, ndf1
    integer(kind=i_def), intent(in) :: undf2, ndf2
    real(kind=r_def), dimension(undf1), intent(inout) :: lhs
    real(kind=r_def), dimension(undf2), intent(in) :: x
    real(kind=r_def), dimension(bandwidth,nrow,ncell_2d), intent(in) :: columnwise_matrix
    integer(kind=i_def), dimension(ndf1), intent(in) :: map1
    integer(kind=i_def), dimension(ndf2), intent(in) :: map2

    integer(kind=i_def), intent(in) :: alpha, beta, gamma_m, gamma_p
    integer(kind=i_def), dimension(nrow), intent(in) :: indirection_dofmap_to
    integer(kind=i_def), dimension(ncol), intent(in) :: indirection_dofmap_from

    ! Internal parameters
    integer(kind=i_def) :: i, mu_i ! Row and column index

    ! Arrays c' and d' used in Thomas algorithm
    real(kind=r_def), dimension(nrow) :: c_prime, d_prime
    ! inverse denominator
    real(kind=r_def) :: inv_denom
    
    ! Spurious instructions to avoid 'unused variable' warnings
    i = alpha + beta + gamma_m + gamma_p

    ! Step 1: Forward sweep, loop over all rows
    do i=1, nrow
       mu_i = map2(1) + indirection_dofmap_from(i) - 1
       if (i == 1) then 
          ! First row
          inv_denom = 1.0_r_def/columnwise_matrix(2,i,cell)
          c_prime(i) = inv_denom * columnwise_matrix(3,i,cell)
          d_prime(i) = inv_denom * x(mu_i)
       else 
          ! Subsequent rows 2,...,nrow-1
          inv_denom = 1.0_r_def / ( columnwise_matrix(2,i,cell) &
                    - columnwise_matrix(1,i,cell) * c_prime(i-1) )
          if (i < nrow) then
             ! We don't need c' in the last row
             c_prime(i) = inv_denom * columnwise_matrix(3,i,cell)
          end if
          d_prime(i) = inv_denom * ( x(mu_i) &
                     - columnwise_matrix(1,i,cell) * d_prime(i-1) )
       end if
    end do
    ! Step 2: Backward sweep (substitution), loop over all rows backwards
    do i=nrow,1,-1
       ! Overwrite d' with solution and then copy to correct position in vector
       mu_i = map1(1) + indirection_dofmap_to(i) - 1
       if (i<nrow) then 
          d_prime(i) = d_prime(i) - c_prime(i) * d_prime(i+1)
       end if
       lhs(mu_i) = d_prime(i)
    end do

  end subroutine columnwise_op_appinv_kernel_code

end module columnwise_op_appinv_kernel_mod
