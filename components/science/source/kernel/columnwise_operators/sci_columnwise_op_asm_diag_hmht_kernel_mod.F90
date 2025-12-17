!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which calculates the diagonal contribution to the term
!!        D_h*M_{2v,lumped,inv}*D_h.
!> @details Takes the operator D_h in LMA representation and the field
!!          representation of the diagonally lumped horizontal velocity mass
!!          matrix. Based on this, the kernel assembles a CMA which contains
!!          the diagonal couplings of the term D_h*M_{2v,lumped,inv}*D_h^T.

module sci_columnwise_op_asm_diag_hmht_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type,                 &
                                    GH_OPERATOR, GH_REAL,     &
                                    GH_COLUMNWISE_OPERATOR,   &
                                    GH_READ, GH_WRITE,        &
                                    ANY_SPACE_1, ANY_SPACE_2, &
                                    CELL_COLUMN

use constants_mod,           only : r_solver, i_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: columnwise_op_asm_diag_hmht_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                                 &
       arg_type(GH_OPERATOR,            GH_REAL, GH_READ,  ANY_SPACE_1, ANY_SPACE_2), &
       arg_type(GH_OPERATOR,            GH_REAL, GH_READ,  ANY_SPACE_2, ANY_SPACE_2), &
       arg_type(GH_COLUMNWISE_OPERATOR, GH_REAL, GH_WRITE, ANY_SPACE_1, ANY_SPACE_1)  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: columnwise_op_asm_diag_hmht_kernel_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: columnwise_op_asm_diag_hmht_kernel_code

contains

!> @brief The subroutine which is called directly from the PSy layer and
!!        assembles the LMA into a CMA operator.
!> @details Given an LMA representation of the operator mapping between two
!!          horizontally discontinuous spaces, assemble the columnwise matrix
!!          representation of the operator.
!>
!> @param[in]  cell Horizontal cell index
!> @param[in]  nlayers Number of vertical layers
!> @param[in]  ncell_2d Number of cells in 2D grid
!> @param[in]  ncell_3d Total number of cells
!> @param[in]  local_stencil_Dh Locally assembled matrix for \f$D_h\f$
!> @param[in]  ncell_3d_tmp Total number of cells (unused duplicate)
!> @param[in]  local_stencil_M2h Locally assembled matrix for \f$M_{2h}\f$
!> @param[in,out] columnwise_matrix Banded matrix to assemble into
!> @param[in]  nrow Number of rows (and columns) in the banded matrix
!> @param[in]  bandwidth Bandwidth of the banded matrix
!> @param[in]  alpha banded Matrix parameter \f$\alpha\f$
!> @param[in]  beta banded Matrix parameter \f$\beta\f$
!> @param[in]  gamma_m Banded matrix parameter \f$\gamma_-\f$
!> @param[in]  gamma_p Banded matrix parameter \f$\gamma_+\f$
!> @param[in]  ndf_w3 Number of dofs per cell for the W_3 space
!> @param[in]  column_banded_dofmap List of offsets for W3-space
!> @param[in]  ndf_w2h Number of dofs per cell for the W_{2h} space
subroutine columnwise_op_asm_diag_hmht_kernel_code(cell,                 &
                                                   nlayers,              &
                                                   ncell_2d,             &
                                                   ncell_3d,             &
                                                   local_stencil_Dh,     &
                                                   ncell_3d_tmp,         &
                                                   local_stencil_M2h,    &
                                                   columnwise_matrix,    &
                                                   nrow,                 &
                                                   bandwidth,            &
                                                   alpha,                &
                                                   beta,                 &
                                                   gamma_m,              &
                                                   gamma_p,              &
                                                   ndf_w3,               &
                                                   column_banded_dofmap, &
                                                   ndf_w2h               )

  implicit none

  ! Arguments
  integer(kind=i_def),                                        intent(in)  :: cell
  integer(kind=i_def),                                        intent(in)  :: nlayers
  integer(kind=i_def),                                        intent(in)  :: ncell_3d
  integer(kind=i_def),                                        intent(in)  :: ncell_3d_tmp
  integer(kind=i_def),                                        intent(in)  :: ncell_2d
  integer(kind=i_def),                                        intent(in)  :: alpha
  integer(kind=i_def),                                        intent(in)  :: beta
  integer(kind=i_def),                                        intent(in)  :: gamma_m
  integer(kind=i_def),                                        intent(in)  :: gamma_p
  integer(kind=i_def),                                        intent(in)  :: nrow
  integer(kind=i_def),                                        intent(in)  :: bandwidth
  integer(kind=i_def),                                        intent(in)  :: ndf_w3
  integer(kind=i_def),                                        intent(in)  :: ndf_w2h
  integer(kind=i_def), dimension(ndf_w3,nlayers),             intent(in)  :: column_banded_dofmap
  real   (kind=r_solver), dimension(ncell_3d,ndf_w3,ndf_w2h), intent(in)  :: local_stencil_Dh
  real   (kind=r_solver), dimension(ncell_3d,ndf_w2h,ndf_w2h),intent(in)  :: local_stencil_M2h
  real   (kind=r_solver), dimension(bandwidth,nrow,ncell_2d), intent(inout) :: columnwise_matrix


  ! Internal parameters
  integer(kind=i_def) :: df1, df2, df3  ! Loop indices for dofs
  integer(kind=i_def) :: i,j            ! Row and column index index
  integer(kind=i_def) :: j_minus        ! First column in a row
  integer(kind=i_def) :: ik             ! ncell3d counter
  integer(kind=i_def) :: k              ! nlayers counter
  real   (kind=r_solver) :: tmp            ! Local contribution

  k = gamma_m

  ! Initialise matrix to zero
  columnwise_matrix( :, :, cell ) = 0.0_r_solver
  ! Loop over all vertical layers
  do k = 1, nlayers
    ik = (cell-1)*nlayers + k ! Cell index in 3D
    do df1 = 1, ndf_w3
      i = column_banded_dofmap( df1, k )
      j_minus = ceiling((alpha*i-gamma_p)/(1.0_r_solver*beta), i_def)
      do df2 = 1, ndf_w3
        tmp = 0.0_r_solver
        do df3 = 1, ndf_w2h
          tmp = tmp                               &
              + local_stencil_Dh ( ik, df1, df3 ) &
              * local_stencil_Dh ( ik, df2, df3 ) &
              / local_stencil_M2h( ik, df3, df3 )
        end do
        j = column_banded_dofmap( df2, k )
        columnwise_matrix( j-j_minus+1, i, cell ) &
           = columnwise_matrix( j-j_minus+1, i, cell ) + real(tmp, r_solver)
      end do
    end do
  end do

end subroutine columnwise_op_asm_diag_hmht_kernel_code

end module sci_columnwise_op_asm_diag_hmht_kernel_mod
