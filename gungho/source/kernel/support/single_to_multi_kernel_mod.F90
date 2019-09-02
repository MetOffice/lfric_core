!-------------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Place a single-level field into a higher-order field
!> @details Temporary infrastructure required to place a single-level field
!>          from IO data into a higher-order field currently being used
!>          to handle multi-dimensional fields. It will be retired when
!>          multi-dimensional fields are properly implemented.

module single_to_multi_kernel_mod

  use argument_mod,  only: arg_type, GH_FIELD, GH_INTEGER, GH_INC, GH_READ, &
                           ANY_SPACE_1, ANY_SPACE_2, CELLS
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone 
  type, public, extends(kernel_type) :: single_to_multi_kernel_type
      private
      type(arg_type) :: meta_args(3) = (/               &
          arg_type(GH_FIELD,   GH_INC,   ANY_SPACE_1),  & ! multi dim field
          arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_2),  & ! single lev field
          arg_type(GH_INTEGER, GH_READ               )  & ! dim to update
          /)
      integer :: iterates_over = CELLS
  contains
      procedure, nopass :: single_to_multi_code
  end type

  ! Overload the default structure constructor for function space
  interface single_to_multi_kernel_type
      module procedure single_to_multi_kernel_constructor
  end interface

  public single_to_multi_code

contains

  type (single_to_multi_kernel_type) function &
  single_to_multi_kernel_constructor() result(self)
    implicit none
    return
  end function single_to_multi_kernel_constructor

  !> @param[in]     nlayers       The number of layers
  !> @param[in,out] multi_field   Multi-dimensional field to update
  !> @param[in]     single_field  Single-level field to update it with
  !> @param[in]     update_dim    Dimension to update
  !> @param[in]     ndf_multi     Number of DOFs per cell for multi-dim field
  !> @param[in]     undf_multi    Number of total DOFs for multi-dim field
  !> @param[in]     map_multi     Dofmap for cell for multi-dim fields
  !> @param[in]     ndf_2d        Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d       Number of total DOFs for 2D fields
  !> @param[in]     map_2d        Dofmap for cell for 2D fields
  subroutine single_to_multi_code(nlayers,                          &
                                  multi_field,                      &
                                  single_field,                     &
                                  update_dim,                       &
                                  ndf_multi, undf_multi, map_multi, &
                                  ndf_2d, undf_2d, map_2d)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, update_dim
    integer(kind=i_def), intent(in) :: ndf_multi, undf_multi
    integer(kind=i_def), intent(in) :: map_multi(ndf_multi)
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)

    real(kind=r_def), intent(inout) :: multi_field(undf_multi)

    real(kind=r_def), intent(in) :: single_field(undf_2d)

    ! Copy 2D field to correct location in multi-dim field
    multi_field(map_multi(update_dim)) = single_field(map_2d(1))

  end subroutine single_to_multi_code

end module single_to_multi_kernel_mod
