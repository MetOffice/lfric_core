!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Stores the diagonal elements of a mass matrix

!> @details Stores the diagonal elements of a mass matrix M into a field D
!>          i.e D(df) = M(df,df)


module sci_mm_diagonal_kernel_mod
  use argument_mod,            only : arg_type,              &
                                      GH_FIELD, GH_OPERATOR, &
                                      GH_READ, GH_INC,       &
                                      GH_REAL, ANY_SPACE_1,  &
                                      CELL_COLUMN
  use constants_mod,           only : i_def, r_single, r_double
  use kernel_mod,              only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  type, public, extends(kernel_type) :: mm_diagonal_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                    &
         arg_type(GH_FIELD,    GH_REAL, GH_INC,  ANY_SPACE_1),             &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ, ANY_SPACE_1, ANY_SPACE_1) &
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: mm_diagonal_kernel_code

  ! Generic interface for real32 and real64 types
  interface mm_diagonal_kernel_code
    module procedure  &
      mm_diagonal_kernel_code_r_single, &
      mm_diagonal_kernel_code_r_double
  end interface

contains

  !> @brief Stores the diagonal of a mass_matrix
  !> @brief real32 and real64 variants
  !> @param[in] cell Horizontal cell index
  !! @param[in] nlayers Number of layers
  !! @param[inout] mm_diag Field array to store the diagonal entries
  !!               of the mass matrix (reals)
  !! @param[in] ncell_3d Total number of cells
  !! @param[in] mass_matrix Array holding mass matrix values (reals)
  !! @param[in] ndf Number of degrees of freedom per cell
  !! @param[in] undf Unique number of degrees of freedom
  !! @param[in] map Dofmap for the cell at the base of the column

  ! R_SINGLE PRECISION
  ! ==================
  subroutine mm_diagonal_kernel_code_r_single(cell,        &
                                              nlayers,     &
                                              mm_diag,     &
                                              ncell_3d,    &
                                              mass_matrix, &
                                              ndf, undf, map)

    implicit none

    ! Arguments
    integer(kind=i_def),                    intent(in)  :: cell, nlayers, ndf
    integer(kind=i_def),                    intent(in)  :: undf, ncell_3d
    integer(kind=i_def), dimension(ndf),    intent(in)  :: map
    real(kind=r_single), dimension(undf), intent(inout) :: mm_diag
    real(kind=r_single), dimension(ncell_3d,ndf,ndf), intent(in) :: mass_matrix

    ! Internal variables
    integer(kind=i_def) :: df, k, ik

    do k = 0, nlayers-1
      ik = (cell-1)*nlayers + k + 1
      do df = 1,ndf
        mm_diag(map(df)+k) = mm_diag(map(df)+k) + mass_matrix(ik,df,df)
      end do
    end do

  end subroutine mm_diagonal_kernel_code_r_single


  ! R_DOUBLE PRECISION
  ! ==================
  subroutine mm_diagonal_kernel_code_r_double(cell,        &
                                              nlayers,     &
                                              mm_diag,     &
                                              ncell_3d,    &
                                              mass_matrix, &
                                              ndf, undf, map)

   implicit none

   ! Arguments
   integer(kind=i_def),                    intent(in)  :: cell, nlayers, ndf
   integer(kind=i_def),                    intent(in)  :: undf, ncell_3d
   integer(kind=i_def), dimension(ndf),    intent(in)  :: map
   real(kind=r_double), dimension(undf), intent(inout) :: mm_diag
   real(kind=r_double), dimension(ncell_3d,ndf,ndf), intent(in) :: mass_matrix

   ! Internal variables
   integer(kind=i_def) :: df, k, ik

   do k = 0, nlayers-1
     ik = (cell-1)*nlayers + k + 1
     do df = 1,ndf
        mm_diag(map(df)+k) = mm_diag(map(df)+k) + mass_matrix(ik,df,df)
     end do
   end do

  end subroutine mm_diagonal_kernel_code_r_double

end module sci_mm_diagonal_kernel_mod
