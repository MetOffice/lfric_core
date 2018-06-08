!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Corrects the sign of the winds used to calculate departure points in
!>        Cosmic.

module correct_cosmic_wind_kernel_mod

  use argument_mod,      only : arg_type, func_type,       &
                                GH_FIELD, GH_READ, GH_INC, &
                                CELLS
  use constants_mod,     only : i_def, r_def
  use cosmic_flux_mod,   only : dof_to_update
  use fs_continuity_mod, only : W3, W2
  use kernel_mod,        only : kernel_type

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: correct_cosmic_wind_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                    &
        arg_type(GH_FIELD,  GH_INC,   W2),                 &
        arg_type(GH_FIELD,  GH_READ,  W2),                 &
        arg_type(GH_FIELD,  GH_READ,  W3)                  &
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass ::correct_cosmic_wind_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! Overload the default structure constructor for function space
  interface correct_cosmic_wind_kernel_type
    module procedure correct_cosmic_wind_kernel_constructor
  end interface

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public correct_cosmic_wind_code

contains

type(correct_cosmic_wind_kernel_type) function correct_cosmic_wind_kernel_constructor() result(self)
  implicit none
  return
end function correct_cosmic_wind_kernel_constructor

  !> @brief Corrects the sign of the winds used to calculate departure points in
  !>        Cosmic.
  !> @details The sign of the winds need correcting for use in Cosmic since the 
  !>          Cosmic scheme is a finite-volume scheme and we require the physical
  !>          values of the wind, rather than the Piola wind values.
  !> @param[in] nlayers Integer the number of layers
  !> @param[in] wind_out The output field for the wind
  !> @param[in] wind_in The input field for the wind
  !> @param[in] orientation The orientation of the cells, in particular in the halo
  !> @param[in] undf_w2 The number of unique degrees of freedom for the wind fields
  !> @param[in] ndf_w2 The number of degrees of freedom per cell for the wind fields
  !> @param[in] map_w2 Integer array holding the dofmap for the cell at the base of the column for the wind fields
  !> @param[in] undf_w3 The number of unique degrees of freedom for the cell orientation field
  !> @param[in] ndf_w3 The number of degrees of freedom per cell for the cell orientation field
  !> @param[in] map_w3 Integer array holding the dofmap for the cell at the base of the column for the cell orientation field
  !> @param[in] direction The direction in which the winds are corrected
  subroutine correct_cosmic_wind_code(nlayers,                                 &
                                      wind_out,                                &
                                      wind_in,                                 &
                                      orientation,                             &
                                      undf_w2,ndf_w2, map_w2,                  &
                                      undf_w3,ndf_w3, map_w3,                  &
                                      direction                                &
                                      )


    use coordinate_jacobian_mod, only: coordinate_jacobian
    use log_mod,                 only: log_event, log_scratch_space, LOG_LEVEL_INFO
    use flux_direction_mod,      only: x_direction, y_direction

    implicit none

    ! Arguments
    integer(kind=i_def),                        intent(in)    :: nlayers
    integer(kind=i_def),                        intent(in)    :: ndf_w3, undf_w3
    integer(kind=i_def), dimension(ndf_w3),     intent(in)    :: map_w3
    integer(kind=i_def),                        intent(in)    :: ndf_w2, undf_w2, direction
    integer(kind=i_def), dimension(ndf_w2),     intent(in)    :: map_w2
    real(kind=r_def), dimension(undf_w3),       intent(in)    :: orientation
    real(kind=r_def), dimension(undf_w2),       intent(in)    :: wind_in
    real(kind=r_def), dimension(undf_w2),       intent(inout) :: wind_out

    ! Internal variables
    integer(kind=i_def)   :: df, k, local_dofs_x(1:2), local_dofs_y(1:2)
    integer(kind=i_def)   :: int_orientation

    int_orientation = int(orientation(map_w3(1)),i_def)

    if (int_orientation > 0_i_def .and. int_orientation < 5_i_def) then

      local_dofs_x = dof_to_update(int_orientation,x_direction)
      local_dofs_y = dof_to_update(int_orientation,y_direction)

      do k = 0, nlayers-1

          if (direction == x_direction ) then
            if (int_orientation == 3 .or. int_orientation == 2) then
              wind_out(map_w2(local_dofs_x(1)) + k) = -1.0_r_def*wind_in(map_w2(local_dofs_x(1)) + k)
              wind_out(map_w2(local_dofs_x(2)) + k) = -1.0_r_def*wind_in(map_w2(local_dofs_x(2)) + k)
            else
              wind_out(map_w2(local_dofs_x(1)) + k) = wind_in(map_w2(local_dofs_x(1)) + k)
              wind_out(map_w2(local_dofs_x(2)) + k) = wind_in(map_w2(local_dofs_x(2)) + k)
            end if
          elseif (direction == y_direction) then
            if (int_orientation == 3 .or. int_orientation == 4) then
              wind_out(map_w2(local_dofs_y(1)) + k) = -1.0_r_def*wind_in(map_w2(local_dofs_y(1)) + k)
              wind_out(map_w2(local_dofs_y(2)) + k) = -1.0_r_def*wind_in(map_w2(local_dofs_y(2)) + k)
            else
              wind_out(map_w2(local_dofs_y(1)) + k) = wind_in(map_w2(local_dofs_y(1)) + k)
              wind_out(map_w2(local_dofs_y(2)) + k) = wind_in(map_w2(local_dofs_y(2)) + k)
            end if
          endif

      end do

    end if

  end subroutine correct_cosmic_wind_code

end module correct_cosmic_wind_kernel_mod
