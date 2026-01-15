!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the average of nearest WTHETA space to W0 space.
!> @details Kernel to average a WTHETA lower-level field to the horizontal parts
!!          of a W0 field. The method is valid for the bottom DoFs, horizontal
!!          and lower level of the lowest-order finite elements on a cubed-sphere mesh.

module sci_wth_to_w0_average_kernel_mod

  use argument_mod,      only : arg_type,             &
                                GH_FIELD, GH_REAL,    &
                                GH_INC, GH_READ,      &
                                CELL_COLUMN
  use constants_mod,     only : i_def, r_def
  use fs_continuity_mod, only : W0, WTHETA
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: wth_to_w0_average_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/               &
         arg_type(GH_FIELD, GH_REAL, GH_INC,  W0),    & ! field_w0
         arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),& ! field_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W0)    & ! rmultiplicity_w0
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: wth_to_w0_average_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: wth_to_w0_average_code

  contains

  !> @brief Computes a 1-2-1 Filter from WTHETA to W0 space.
  !> @param[in]     nlayers          Number of layers
  !> @param[in,out] field_w0         Output field from Filter on W0 space
  !> @param[in]     field_wth        Input field for filter on WTHETA space
  !> @param[in]     rmultiplicity_w0 Reciprocal of how many times the dof has been visited in total on w0 space
  !> @param[in]     ndf_w0           Number of degrees of freedom per cell for W0
  !> @param[in]     undf_w0          Number of unique degrees of freedom for W0
  !> @param[in]     map_w0           Dofmap for the cell at the base of the column for W0
  !> @param[in]     ndf_wtheta       Number of degrees of freedom per cell for WTHETA
  !> @param[in]     undf_wtheta      Number of unique degrees of freedom for WTHETA
  !> @param[in]     map_wtheta       Dofmap for the cell at the base of the column for WTHETA
  subroutine wth_to_w0_average_code(nlayers,                 &
                                    field_w0, field_wth,     &
                                    rmultiplicity_w0,        &
                                    ndf_w0, undf_w0, map_w0, &
                                    ndf_wtheta, undf_wtheta, &
                                    map_wtheta)


  implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w0, undf_w0
    integer(kind=i_def), intent(in) :: ndf_wtheta, undf_wtheta

    integer(kind=i_def), intent(in), dimension(ndf_w0)     :: map_w0
    integer(kind=i_def), intent(in), dimension(ndf_wtheta) :: map_wtheta

    real(kind=r_def), intent(inout), dimension(undf_w0)     :: field_w0
    real(kind=r_def), intent(in),    dimension(undf_wtheta) :: field_wth
    real(kind=r_def), intent(in),    dimension(undf_w0)     :: rmultiplicity_w0

    ! Internal variables
    integer(kind=i_def) :: df, k

    do df = 1, 4 ! Use bottom four W0 DoFs in cell
      do k = 0, nlayers
        ! Use rmultiplicity from bottom layer to ensure appropriate average
        field_w0(map_w0(df) + k) = field_w0(map_w0(df) + k)                    &
          + field_wth(map_wtheta(1) + k)*rmultiplicity_w0(map_w0(df))
      end do
    end do

  end subroutine wth_to_w0_average_code

end module sci_wth_to_w0_average_kernel_mod
