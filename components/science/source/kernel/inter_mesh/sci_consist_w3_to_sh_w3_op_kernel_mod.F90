!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Forms the operator for mapping from W3 to shifted W3.
!> @details Calculates the operator for mapping from W3 to shifted W3.
!!          The resulting (non-square) matrix representing this operation is
!!          bi-diagonal. We store the diagonal components in W3 fields.
!!          This only works for the lowest-order finite element spaces.
module sci_consist_w3_to_sh_w3_op_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    GH_READ, GH_WRITE,         &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    GH_OPERATOR, CELL_COLUMN
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : W3
  use kernel_mod,            only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !! Psy layer.
  !!
  type, public, extends(kernel_type) :: consist_w3_to_sh_w3_op_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                       &
         arg_type(GH_FIELD*2,  GH_REAL, GH_WRITE, W3),                        &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3,  &
                                                  ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ,  W3, W3)                     &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: consist_w3_to_sh_w3_op_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: consist_w3_to_sh_w3_op_code

contains

!> @brief Forms the operator for mapping from W3 to shifted W3.
!> @param[in] cell         Column index
!> @param[in] nlayers      Number of layers in the original mesh
!> @param[in,out] T_ip1    Below-diagonal components of the transform matrix.
!>                         Is a field in W3.
!> @param[in,out] T_i      Above-diagonal components of the transform matrix.
!>                         Is a field in W3.
!> @param[in] ncells_sh    Total number of cells on the shifted mesh
!> @param[in] mm_w3_sh_inv Inverse of the shifted W3 mass matrix
!> @param[in] ncells       Total number of cells on the original mesh
!> @param[in] mm_w3        W3 mass matrix
!> @param[in] ndf_w3       Number of DoFs per cell for W3
!> @param[in] undf_w3      Number of universal DoFs for W3
!> @param[in] map_w3       DoFmap for W3
!> @param[in] ndf_w3_sh    Number of DoFs for W3 shifted
subroutine consist_w3_to_sh_w3_op_code( cell,           &
                                        nlayers,        &
                                        T_ip1,          &
                                        T_i,            &
                                        ncells_sh,      &
                                        mm_w3_sh_inv,   &
                                        ncells,         &
                                        mm_w3,          &
                                        ndf_w3,         &
                                        undf_w3,        &
                                        map_w3,         &
                                        ndf_w3_sh       &
                                      )

!  use sci_coordinate_jacobian_mod, only: coordinate_jacobian


  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in) :: cell, nlayers
  integer(kind=i_def),                           intent(in) :: ncells_sh, ncells
  integer(kind=i_def),                           intent(in) :: ndf_w3_sh, ndf_w3
  integer(kind=i_def),                           intent(in) :: undf_w3
  integer(kind=i_def), dimension(ndf_w3),        intent(in) :: map_w3

  real(kind=r_def),    dimension(undf_w3),    intent(inout) :: T_ip1, T_i

  real(kind=r_def), dimension(ncells_sh,ndf_w3_sh,ndf_w3_sh), intent(in) :: mm_w3_sh_inv
  real(kind=r_def), dimension(ncells,ndf_w3,ndf_w3),          intent(in) :: mm_w3

  ! Internal variables
  integer(kind=i_def) :: df, k, ik, ik_sh, nlayers_sh

  ! Lowest order so only consider 1 DoF
  df = 1
  nlayers_sh = nlayers + 1

  do k = 0, nlayers-1
    ik_sh = (cell-1)*nlayers_sh + k + 1  ! Index of cell on shifted mesh
    ik = (cell-1)*nlayers + k + 1        ! Index of cell on original mesh

    ! Formula is matmul(mm_w3_sh_inv, matmul(halves, mm_w3))
    ! halves is a bidiagonal matrix mapping from W3 to shifted W3,
    ! whose non-zero entries are 1/2

    ! T_i corresponds to the lower-half layers of the original mesh
    T_i(map_w3(df)+k) = 0.5_r_def * mm_w3(ik,df,df) * mm_w3_sh_inv(ik_sh,df,df)

    ! T_ip1 corresponds to the upper-half layers of the original mesh
    T_ip1(map_w3(df)+k) = 0.5_r_def * mm_w3(ik,df,df) * mm_w3_sh_inv(ik_sh+1,df,df)

  end do

end subroutine consist_w3_to_sh_w3_op_code

end module sci_consist_w3_to_sh_w3_op_kernel_mod
