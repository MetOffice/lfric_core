!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which constructs a CMA representation of the diagonally
!!        lumped inverse of the vertical velocity mass matrix.
!> @details Extract the diagonal of a locally assembled matrix (LMA) for the
!!          vertical velocity mass matrix. The inverse of this diagonal is then
!!          assembled into a CMA. Note that this CMA is diagonal, i.e. it has
!!          parameters \f$\alpha=\beta=1\f$ and \f$\gamma_-=\gamma_+=0\f$.
!>

module sci_columnwise_op_asm_m2v_lumped_inv_krnl_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type,                            &
                                    GH_OPERATOR, GH_COLUMNWISE_OPERATOR, &
                                    GH_REAL, GH_READ, GH_WRITE,          &
                                    ANY_SPACE_1,                         &
                                    CELL_COLUMN

use constants_mod,           only : r_solver, i_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: columnwise_op_asm_m2v_lumped_inv_krnl_type
  private
  type(arg_type) :: meta_args(2) = (/                                                 &
       arg_type(GH_OPERATOR,            GH_REAL, GH_READ,  ANY_SPACE_1, ANY_SPACE_1), &
       arg_type(GH_COLUMNWISE_OPERATOR, GH_REAL, GH_WRITE, ANY_SPACE_1, ANY_SPACE_1)  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: columnwise_op_asm_m2v_lumped_inv_krnl_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: columnwise_op_asm_m2v_lumped_inv_krnl_code

contains

!> @brief The subroutine which is called directly from the PSy layer and
!!        assembles the LMA into a CMA operator.
!> @details Given an LMA representation of local operator for the vertical
!!          velocity mass matrix, assemble the columnwise matrix which
!!          represents the inverse lumped mass matrix.
!>
!> @param[in]  cell Horizontal cell index
!> @param[in]  nlayers Number of vertical layers
!> @param[in]  ncell_2d Number of cells in 2D grid
!> @param[in]  ncell_3d Total number of cells
!> @param[in]  local_stencil Locally assembled matrix
!> @param[in,out] columnwise_matrix Banded matrix to assemble into
!> @param[in]  nrow Number of rows (and columns) in the banded matrix
!> @param[in]  bandwidth Bandwidth of the banded matrix
!> @param[in]  alpha Banded matrix parameter \f$\alpha=1\f$
!> @param[in]  beta Banded matrix parameter \f$\beta=1\f$
!> @param[in]  gamma_m Banded matrix parameter \f$\gamma_-=0\f$
!> @param[in]  gamma_p Banded matrix parameter \f$\gamma_+=0\f$
!> @param[in]  ndf Number of degrees of freedom per cell for the function space
!> @param[in]  column_banded_dofmap List of offsets for function space
subroutine columnwise_op_asm_m2v_lumped_inv_krnl_code(cell,                    &
                                                      nlayers,                 &
                                                      ncell_2d,                &
                                                      ncell_3d,                &
                                                      local_stencil,           &
                                                      columnwise_matrix,       &
                                                      nrow,                    &
                                                      bandwidth,               &
                                                      alpha,                   &
                                                      beta,                    &
                                                      gamma_m,                 &
                                                      gamma_p,                 &
                                                      ndf,                     &
                                                      column_banded_dofmap)

  implicit none

  ! Arguments
  integer(kind=i_def),                                        intent(in)  :: cell
  integer(kind=i_def),                                        intent(in)  :: nlayers
  integer(kind=i_def),                                        intent(in)  :: ncell_3d
  integer(kind=i_def),                                        intent(in)  :: ncell_2d
  integer(kind=i_def),                                        intent(in)  :: nrow
  integer(kind=i_def),                                        intent(in)  :: bandwidth
  integer(kind=i_def),                                        intent(in)  :: ndf
  integer(kind=i_def),                                        intent(in)  :: alpha
  integer(kind=i_def),                                        intent(in)  :: beta
  integer(kind=i_def),                                        intent(in)  :: gamma_m
  integer(kind=i_def),                                        intent(in)  :: gamma_p
  integer(kind=i_def),    dimension(ndf,nlayers),             intent(in)  :: column_banded_dofmap
  real   (kind=r_solver), dimension(ncell_3d,ndf,ndf),        intent(in)  :: local_stencil
  real   (kind=r_solver), dimension(bandwidth,nrow,ncell_2d), intent(inout) :: columnwise_matrix

  ! Internal parameters
  integer(kind=i_def) :: df1    ! Loop index for dofs
  integer(kind=i_def) :: i      ! Row and column index
  integer(kind=i_def) :: ik     ! ncell3d counter
  integer(kind=i_def) :: k      ! nlayers counter

  k = alpha + beta + gamma_m + gamma_p

  ! Initialise matrix to zero
  columnwise_matrix( :, :, cell ) = 0.0_r_solver
  ! Loop over all vertical layers add add up diagonal entries
  do k = 1, nlayers
    ik = (cell-1)*nlayers + k ! Cell index in 3D
    do df1 = 1, ndf
      i = column_banded_dofmap( df1, k )
      columnwise_matrix( 1, i, cell ) = columnwise_matrix( 1, i, cell ) &
                                      + real(local_stencil( ik, df1, df1 ), r_solver)
    end do
  end do

  ! Calculate inverse of diagonal entries
  where (columnwise_matrix(:,:,cell) /= 0.0_r_solver)
     columnwise_matrix(:,:,cell) = 1.0_r_solver/columnwise_matrix(:,:,cell)
  end where

end subroutine columnwise_op_asm_m2v_lumped_inv_krnl_code

end module sci_columnwise_op_asm_m2v_lumped_inv_krnl_mod
