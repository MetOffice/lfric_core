!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office GitHub Copilot Enterprise.
!-----------------------------------------------------------------------------

!> @brief Defines an object to pair 2D field arrays with a unique identifier.
module id_r64_field_array_2d_pair_mod

  use constants_mod,         only: i_def
  use field_real64_mod,      only: field_real64_type
  use function_space_mod,    only: function_space_type
  use id_abstract_pair_mod,  only: id_abstract_pair_type

  implicit none

  private

  ! ========================================================================== !
  ! ID-Field 2D Array Pair
  ! ========================================================================== !

  !> @brief An object pairing a 2D field array with a unique identifier
  !>
  type, public, extends(id_abstract_pair_type) :: id_r64_field_array_2d_pair_type

    private

    type(field_real64_type), allocatable :: field_array_(:,:)

  contains

    procedure, public :: initialise
    procedure, public :: copy_initialise
    procedure, public :: get_field_array

    final :: destructor

  end type id_r64_field_array_2d_pair_type

contains

  !> @brief Initialises the id_r64_field_array_2d_pair object with new fields
  !> @param[in] fs           The function space of the new fields
  !> @param[in] id           The integer ID to pair with the field array
  !> @param[in] dim1         The size of the first dimension of the field array
  !> @param[in] dim2         The size of the second dimension of the field array
  !> @param[in] halo_depth   Optional halo depth for field (to overwrite the
  !!                         default halo depth)
  subroutine initialise(self, fs, id, dim1, dim2, halo_depth)

    implicit none

    class(id_r64_field_array_2d_pair_type),   intent(inout) :: self
    type(function_space_type),       pointer, intent(in)    :: fs
    integer(kind=i_def),                      intent(in)    :: id
    integer(kind=i_def),                      intent(in)    :: dim1
    integer(kind=i_def),                      intent(in)    :: dim2
    integer(kind=i_def),            optional, intent(in)    :: halo_depth

    integer(kind=i_def) :: i, j

    allocate(self%field_array_(dim1, dim2))

    do j = 1, dim2
      do i = 1, dim1
        call self%field_array_(i,j)%initialise(fs, halo_depth=halo_depth)
      end do
    end do

    call self%set_id(id)

  end subroutine initialise

  !> @brief Initialises the id_r64_field_array_2d_pair object by copying in fields
  !> @param[in] field_array   The fields to be stored in the paired object
  !> @param[in] id            The integer ID to pair with the field_array
  subroutine copy_initialise(self, field_array, id)

    implicit none

    class(id_r64_field_array_2d_pair_type), intent(inout) :: self
    type(field_real64_type),                intent(in)    :: field_array(:,:)
    integer(kind=i_def),                    intent(in)    :: id

    integer(kind=i_def) :: i, j

    allocate(self%field_array_(size(field_array,1), size(field_array,2)))

    do j = 1, size(field_array,2)
      do i = 1, size(field_array,1)
        call self%field_array_(i,j)%initialise(field_array(i,j)%get_function_space(), &
                                               halo_depth=field_array(i,j)%get_field_halo_depth())
        call field_array(i,j)%copy_field_serial(self%field_array_(i,j))
      end do
    end do

    call self%set_id(id)

  end subroutine copy_initialise

  !> @brief Get the field_array corresponding to the paired object
  !> @param[in] self     The paired object
  !> @return             The field_array
  function get_field_array(self) result(field_array)

    implicit none

    class(id_r64_field_array_2d_pair_type), target, intent(in) :: self
    type(field_real64_type),                        pointer    :: field_array(:,:)

    field_array => self%field_array_

  end function get_field_array

  !> @brief Loop over fields in the 2D array and call finaliser on each
  !> Finaliser for the 2D field array.
  subroutine destructor(self)
    implicit none
    type(id_r64_field_array_2d_pair_type), intent(inout) :: self
    integer(kind=i_def) :: i, j
    if ( allocated(self%field_array_) ) then
      do j = 1, size(self%field_array_, 2)
        do i = 1, size(self%field_array_, 1)
          call self%field_array_(i,j)%field_final()
        end do
      end do
      deallocate(self%field_array_)
    end if
  end subroutine destructor

end module id_r64_field_array_2d_pair_mod
