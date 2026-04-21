!-------------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Copy a field to specified halo depth

module copy_field_halo_kernel_mod

  use argument_mod,         only: arg_type,                  &
                                  GH_FIELD, GH_REAL,         &
                                  GH_READ, GH_WRITE,         &
                                  OWNED_AND_HALO_CELL_COLUMN,&
                                  ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,        only: r_def, i_def
  use kernel_mod,           only: kernel_type

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: copy_field_halo_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1)  &
         /)
    integer :: operates_on = OWNED_AND_HALO_CELL_COLUMN
  contains
    procedure, nopass :: copy_field_halo_code
  end type copy_field_halo_kernel_type

  public :: copy_field_halo_code

contains

  !> @brief Copy field to specified halo depth
  !> @param[in]     nlayers       The number of layers
  !> @param[in,out] field_out     Output field
  !> @param[in]     field_in      Input field
  !> @param[in]     ndf           Number of degrees of freedom per cell
  !> @param[in]     undf          Number of total degrees of freedom
  !> @param[in]     map           Dofmap for the cell at the base of the column
  subroutine copy_field_halo_code(nlayers,          &
                                  field_out,        &
                                  field_in,         &
                                  ndf, undf, map)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf, undf
    integer(kind=i_def), intent(in), dimension(ndf) :: map

    ! Arguments passed explicitly from algorithm
    real(kind=r_def), intent(in), dimension(undf) :: field_in
    real(kind=r_def), intent(inout), dimension(undf) :: field_out

    ! Local arguments
    integer(kind=i_def) :: k, dof

    do dof = 1, ndf
      do k = 0, nlayers-1
        field_out(map(dof)+k) = field_in(map(dof)+k)
      end do
    end do

  end subroutine copy_field_halo_code

end module copy_field_halo_kernel_mod
