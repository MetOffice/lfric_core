!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Divides the Piola wind values by detJ at W2 dofs.
!>
module cosmic_departure_wind_kernel_mod

  use argument_mod,       only : arg_type, func_type,       &
                                 GH_FIELD, GH_READ, GH_INC, &
                                 ANY_SPACE_9,               &
                                 GH_DIFF_BASIS, GH_BASIS,   &
                                 CELLS, GH_EVALUATOR
  use constants_mod,      only : r_def, i_def
  use flux_direction_mod, only : x_direction, y_direction
  use fs_continuity_mod,  only : W2
  use kernel_mod,         only : kernel_type

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: cosmic_departure_wind_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                     &
        arg_type(GH_FIELD,    GH_INC,   W2),                &
        arg_type(GH_FIELD,    GH_READ,  W2),                &
        arg_type(GH_FIELD*3,  GH_READ,  ANY_SPACE_9)        &
        /)
    type(func_type) :: meta_funcs(2) = (/                   &
        func_type(W2,          GH_BASIS),                   &
        func_type(ANY_SPACE_9, GH_DIFF_BASIS)               &
        /)
    integer :: iterates_over = CELLS
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass ::cosmic_departure_wind_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! Overload the default structure constructor for function space
  interface cosmic_departure_wind_kernel_type
    module procedure cosmic_departure_wind_kernel_constructor
  end interface

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public cosmic_departure_wind_code

contains

  type(cosmic_departure_wind_kernel_type) function cosmic_departure_wind_kernel_constructor() result(self)
    implicit none
    return
  end function cosmic_departure_wind_kernel_constructor

  !> @brief Divides the Piola wind values by detJ at W2 dofs.
  !> @details The Piola wind values are divided by detJ at the W2 dofs which rescales
  !>          the wind values to the computational grid.
  !> @param[in] nlayers Number of layers
  !> @param[inout] u_departure_wind Output field containing the departure wind used to calculate departure points
  !> @param[in] u_piola Field for the Piola wind
  !> @param[in] detj_at_w2 Input field containing the detj values at W2 locations
  !> @param[in] ndf Number of degrees of freedom per cell for the output field
  !> @param[in] undf Number of unique degrees of freedom for the output field
  !> @param[in] map Dofmap for the cell at the base of the column for the output field
  !> @param[in] direction The direction in which the winds are calculated
  subroutine cosmic_departure_wind_code(nlayers,                                  &
                                        u_departure_wind,                         &
                                        u_piola,                                  &
                                        detj_at_w2,                               &
                                        ndf, undf, map,                           &
                                        direction                                 &
                                        )

    use coordinate_jacobian_mod, only: coordinate_jacobian

    implicit none

    ! Arguments
    integer(kind=i_def),                        intent(in)    :: nlayers
    integer(kind=i_def),                        intent(in)    :: ndf, undf
    integer(kind=i_def), dimension(ndf),        intent(in)    :: map
    real(kind=r_def), dimension(undf),          intent(in)    :: u_piola
    real(kind=r_def), dimension(undf),          intent(in)    :: detj_at_w2
    real(kind=r_def), dimension(undf),          intent(inout) :: u_departure_wind
    integer(kind=i_def),                        intent(in)    :: direction

    ! Internal variables
    integer(kind=i_def)                  :: df, k
    real(kind=r_def)                     :: mult_factor

    ! Change the sign of the output winds depending on whether it is an x or y
    ! direction update
    if (direction == x_direction ) then
      mult_factor = 1.0_r_def
    elseif (direction == y_direction ) then
      mult_factor = -1.0_r_def
    end if

    ! Divide the input Piola wind values by the corresponding detJ value.
    do k = 0, nlayers-1
      do df = 1,ndf
        u_departure_wind(map(df)+k) = mult_factor*u_piola(map(df)+k)/detj_at_w2(map(df)+k)
      end do
    end do

  end subroutine cosmic_departure_wind_code

end module cosmic_departure_wind_kernel_mod
