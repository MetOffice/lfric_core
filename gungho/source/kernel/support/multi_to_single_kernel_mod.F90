!-------------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Extract a single-level field from a higher-order field
!> @details Temporary infrastructure required to extract single-level fields
!>          for IO purposes from higher-order fields currently being used
!>          to handle multi-dimensional fields. It will be retired when
!>          multi-dimensional fields are properly implemented.

module multi_to_single_kernel_mod

  use argument_mod,  only: arg_type, GH_FIELD, GH_INTEGER, GH_READ, &
                           ANY_SPACE_1, ANY_SPACE_2, CELLS, GH_WRITE
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone 
  type, public, extends(kernel_type) :: multi_to_single_kernel_type
      private
      type(arg_type) :: meta_args(3) = (/               &
          arg_type(GH_FIELD,   GH_WRITE, ANY_SPACE_1),  & ! single lev field
          arg_type(GH_FIELD,   GH_READ,  ANY_SPACE_2),  & ! multi dim field
          arg_type(GH_INTEGER, GH_READ               )  & ! dim to update
          /)
      integer :: iterates_over = CELLS
  contains
      procedure, nopass :: multi_to_single_code
  end type

  ! Overload the default structure constructor for function space
  interface multi_to_single_kernel_type
      module procedure multi_to_single_kernel_constructor
  end interface

  public multi_to_single_code

contains

  type (multi_to_single_kernel_type) function &
  multi_to_single_kernel_constructor() result(self)
    implicit none
    return
  end function multi_to_single_kernel_constructor

  !> @param[in]     nlayers       The number of layers
  !> @param[out]    single_field  Single-level field to update
  !> @param[in]     multi_field   Multi-dimensional field to update it with
  !> @param[in]     update_dim    Dimension to update
  !> @param[in]     ndf_2d        Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d       Number of total DOFs for 2D fields
  !> @param[in]     map_2d        Dofmap for cell for 2D fields
  !> @param[in]     ndf_multi     Number of DOFs per cell for multi-dim field
  !> @param[in]     undf_multi    Number of total DOFs for multi-dim field
  !> @param[in]     map_multi     Dofmap for cell for multi-dim fields
  subroutine multi_to_single_code(nlayers,                          &
                                  single_field,                     &
                                  multi_field,                      &
                                  update_dim,                       &
                                  ndf_2d, undf_2d, map_2d,          &    
                                  ndf_multi, undf_multi, map_multi)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, update_dim
    integer(kind=i_def), intent(in) :: ndf_multi, undf_multi
    integer(kind=i_def), intent(in) :: map_multi(ndf_multi)
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d)

    real(kind=r_def), intent(in) :: multi_field(undf_multi)

    real(kind=r_def), intent(out) :: single_field(undf_2d)

    ! Extract 2D field from correct location in multi-dim field
    single_field(map_2d(1)) = multi_field(map_multi(update_dim))

  end subroutine multi_to_single_code

end module multi_to_single_kernel_mod
