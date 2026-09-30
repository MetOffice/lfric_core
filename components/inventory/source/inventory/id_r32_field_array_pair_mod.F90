!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Defines an object to pair field arrays with a unique identifier.
module id_r32_field_array_pair_mod

  use constants_mod,         only: i_def
  use field_real32_mod,      only: field_real32_type
  use function_space_mod,    only: function_space_type
  use id_abstract_pair_mod,  only: id_abstract_pair_type

  implicit none

  private

  ! ========================================================================== !
  ! ID-Field Array Pair
  ! ========================================================================== !

  !> @brief An object pairing a field array with a unique identifier
  !>
  type, public, extends(id_abstract_pair_type) :: id_r32_field_array_pair_type

    private

    type(field_real32_type), allocatable :: field_array_(:)

  contains

    procedure, public :: initialise
    procedure, public :: copy_initialise
    procedure, public :: get_field_array

    final :: destructor

  end type id_r32_field_array_pair_type

contains

  !> @brief Initialises the id_r32_field_pair object with a new field
  !> @param[in] fs           The function space of the new field
  !> @param[in] id           The integer ID to pair with the field
  !> @param[in] array_size   The size of the field array
  !> @param[in] halo_depth   Optional halo depth for field (to overwrite the
  !!                         default halo depth)
  subroutine initialise(self, fs, id, array_size, halo_depth)

    implicit none

    class(id_r32_field_array_pair_type),      intent(inout) :: self
    type(function_space_type),       pointer, intent(in)    :: fs
    integer(kind=i_def),                      intent(in)    :: id
    integer(kind=i_def),                      intent(in)    :: array_size
    integer(kind=i_def),            optional, intent(in)    :: halo_depth

    integer(kind=i_def) :: i

    allocate(self%field_array_(array_size))

    do i = 1, array_size
      call self%field_array_(i)%initialise(fs, halo_depth = halo_depth)
    end do

    call self%set_id(id)

  end subroutine initialise

  !> @brief Initialises the id_r32_field_array_pair object by copying in fields
  !> @param[in] field_array   The fields to be stored in the paired object
  !> @param[in] id            The integer ID to pair with the field_array
  subroutine copy_initialise(self, field_array, id)

    implicit none

    class(id_r32_field_array_pair_type), intent(inout) :: self
    type(field_real32_type),             intent(in)    :: field_array(:)
    integer(kind=i_def),                 intent(in)    :: id

    integer(kind=i_def) :: i

    allocate(self%field_array_(size(field_array)))

    do i = 1, size(field_array)
      call self%field_array_(i)%initialise(field_array(i)%get_function_space(), &
                                           halo_depth=field_array(i)%get_field_halo_depth())
      call field_array(i)%copy_field_serial(self%field_array_(i))
    end do

    call self%set_id(id)

  end subroutine copy_initialise

  !> @brief Get the field_array corresponding to the paired object
  !> @param[in] self     The paired object
  !> @return             The field_array
  function get_field_array(self) result(field_array)

    implicit none

    class(id_r32_field_array_pair_type), target, intent(in) :: self
    type(field_real32_type),                     pointer    :: field_array(:)

    field_array => self%field_array_

  end function get_field_array

  !> @brief Loop over fields in array and call finaliser on each
  !> Finaliser for the field array.
  subroutine destructor(self)
    implicit none
    type(id_r32_field_array_pair_type), intent(inout) :: self
    integer(kind=i_def) :: i
    if( allocated(self%field_array_) ) then
      do i = 1, size(self%field_array_)
        call self%field_array_(i)%field_final()
      end do
      deallocate(self%field_array_)
    end if
  end subroutine destructor

end module id_r32_field_array_pair_mod
