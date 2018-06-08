!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief  Extracts the y direction component from a W2 field.
!>
!> @details Extracts the y direction component from a W2 field and uses the
!>          cell orientation field, currently held in W3, to do this.
!>          The extraction of the y component of the W2 field is done for the
!>          COSMIC transport scheme which acts in the two horizontal directions
!>          separately.
!>
module extract_y_kernel_mod

  use argument_mod,      only : arg_type, func_type,       &
                                GH_FIELD, GH_READ, GH_INC, &
                                GH_BASIS, CELLS, GH_INTEGER
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W2, W3
  use kernel_mod,        only : kernel_type

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: extract_y_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/     &
        arg_type(GH_FIELD,   GH_READ,  W3), &
        arg_type(GH_FIELD,   GH_READ,  W2), &
        arg_type(GH_FIELD,   GH_INC,   W2), &
        arg_type(GH_INTEGER, GH_READ     )  &
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass :: extract_y_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! Overload the default structure constructor for function space
  interface extract_y_kernel_type
    module procedure extract_y_kernel_constructor
  end interface

  public extract_y_code

contains

type(extract_y_kernel_type) function extract_y_kernel_constructor() result(self)
  implicit none
  return
end function extract_y_kernel_constructor

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> @brief  Extracts the y direction component from a W2 field.
!> @details Extracts the y direction component from a W2 field and uses the
!>          cell orientation field, currently held in W3, to do this.
!>          The extraction of the y component of the W2 field is done for the
!>          COSMIC transport scheme which acts in the two horizontal directions
!>          separately.
!! @param[in]  nlayers           The number of model levels
!! @param[in]  cell_orientation  Cell orientation values held in W3 field
!! @param[in]  w2_field_in       Input W2 field
!! @param[out] y_field           Extracted y-direction component of input W2 field
!! @param[in]  undf_w3           Number of unique degrees of freedom for W3
!! @param[in]  ndf_w3            Number of degrees of freedom per cell in W3
!! @param[in]  map_w3            Dofmap for the cell at the base of the column
!! @param[in]  undf_w2           Number of unique degrees of freedom for W2
!! @param[in]  ndf_w2            Number of degrees of freedom per cell in W2
!! @param[in]  map_w2            Dofmap for the cell at the base of the column
subroutine extract_y_code( nlayers,                      &
                           cell_orientation,             &
                           w2_field_in,                  &
                           y_field,                      &
                           undf_w3,                      &
                           ndf_w3,                       &
                           map_w3,                       &
                           undf_w2,                      &
                           ndf_w2,                       &
                           map_w2 )


  use cosmic_flux_mod, only : w2_dof

  implicit none

  integer(kind=i_def), intent(in)                     :: nlayers
  integer(kind=i_def), intent(in)                     :: undf_w3
  integer(kind=i_def), intent(in)                     :: undf_w2
  real(kind=r_def), intent(in)                        :: cell_orientation(1:undf_w3)
  real(kind=r_def), intent(in)                        :: w2_field_in(1:undf_w2)
  real(kind=r_def), intent(out)                       :: y_field(1:undf_w2)
  integer(kind=i_def), intent(in)                     :: ndf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in)  :: map_w3
  integer(kind=i_def), intent(in)                     :: ndf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in)  :: map_w2

  integer(kind=i_def) :: k
  integer(kind=i_def) :: orientation_of_cell
  integer(kind=i_def) :: dof_vector(1:4)

  orientation_of_cell = int(cell_orientation(map_w3(1)),i_def)

  if (orientation_of_cell > 0_i_def .and. orientation_of_cell < 5_i_def ) then
    dof_vector = (/ w2_dof(orientation_of_cell,2),     &
                    w2_dof(orientation_of_cell,3),     &
                    w2_dof(orientation_of_cell,4),     &
                    w2_dof(orientation_of_cell,1)      /)

    do k=0, nlayers-1

      y_field(map_w2(dof_vector(1))+k) = w2_field_in(map_w2(dof_vector(1))+k)
      y_field(map_w2(dof_vector(2))+k) = -9999.0_r_def
      y_field(map_w2(dof_vector(3))+k) = w2_field_in(map_w2(dof_vector(3))+k)
      y_field(map_w2(dof_vector(4))+k) = -9999.0_r_def

    end do
  end if

end subroutine extract_y_code

end module extract_y_kernel_mod
