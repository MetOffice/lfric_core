!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Convert a multi_data field to a higher-order field
!> @details Temporary infrastructure required to convert a multi-data field
!>          into a higher-order field for IO purposes. It will be retired
!>          when multi-dimensional fields are properly implemented.
!>          Also used without propoer psyclone support of multi-data fields
!>          see https://github.com/stfc/PSyclone/issues/868

module multi_to_high_kernel_mod

  use argument_mod,  only: arg_type, CELLS,           &
                           GH_FIELD, GH_INTEGER,      &
                           GH_READ, GH_WRITE,         &
                           ANY_DISCONTINUOUS_SPACE_1, &
                           ANY_DISCONTINUOUS_SPACE_2
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: multi_to_high_kernel_type
      private
      type(arg_type) :: meta_args(3) = (/                             &
          arg_type(GH_FIELD,   GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  & ! high-order field
          arg_type(GH_FIELD,   GH_READ,  ANY_DISCONTINUOUS_SPACE_2),  & ! multi-data field
          arg_type(GH_INTEGER, GH_READ               )                & ! ndata
          /)
      integer :: iterates_over = CELLS
  contains
      procedure, nopass :: multi_to_high_code
  end type

  public multi_to_high_code

contains

  !> @param[in]     nlayers       The number of layers
  !> @param[out]    high_field    Higher-order field to write to
  !> @param[in]     multi_field   Multi-data field to write from
  !> @param[in]     ndata         Dimension of multi-data field
  !> @param[in]     ndf_high      Number of DOFs per cell for high-order field
  !> @param[in]     undf_high     Number of total DOFs for high-order field
  !> @param[in]     map_high      Dofmap for cell for high-order fields
  !> @param[in]     ndf_multi     Number of DOFs per cell for multi-data field
  !> @param[in]     undf_multi    Number of total DOFs for multi-data field
  !> @param[in]     map_multi     Dofmap for cell for multi-data fields
  subroutine multi_to_high_code(nlayers,                          &
                                high_field,                       &
                                multi_field,                      &
                                ndata,                            &
                                ndf_high, undf_high, map_high,    &
                                ndf_multi, undf_multi, map_multi)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, ndata
    integer(kind=i_def), intent(in) :: ndf_multi, undf_multi
    integer(kind=i_def), intent(in) :: map_multi(ndf_multi)
    integer(kind=i_def), intent(in) :: ndf_high, undf_high
    integer(kind=i_def), intent(in) :: map_high(ndf_high)

    real(kind=r_def), intent(in) :: multi_field(undf_multi)
    real(kind=r_def), intent(out)  :: high_field(undf_high)

    integer(kind=i_def) :: i

    ! Convert multi-data field to high-order field
    do i = 1, ndata
      high_field(map_high(i)) = multi_field(map_multi(1)+i-1)
    end do

  end subroutine multi_to_high_code

end module multi_to_high_kernel_mod
