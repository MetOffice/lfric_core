!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Holds and manages fields in a collection
!>
!> @details A container that holds a collection of fields. Fields that are
!>          presented to the field_collection through the add_field() method are
!>          copied, so when the original goes out of scope, the copy in the
!>          field_collection will continue to be maintained.
!
module field_collection_mod

  use constants_mod,      only: i_def, l_def
  use field_mod,          only: field_type
  use log_mod,            only: log_event, log_scratch_space, &
                                LOG_LEVEL_ERROR, LOG_LEVEL_INFO
  use linked_list_mod,    only: linked_list_type, &
                                linked_list_item_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Type that holds a collection of fields in a linked list
  !-----------------------------------------------------------------------------
  type, public :: field_collection_type
    private
    type(linked_list_type) :: field_list

  contains
    procedure, public :: add_field
    procedure, public :: get_field

    procedure, public :: clear
    final             :: field_collection_destructor
  end type field_collection_type
  !-----------------------------------------------------------------------------

  interface field_collection_type
    module procedure field_collection_constructor
  end interface

contains

!> Constructor for a field collection
function field_collection_constructor() result(self)

  implicit none

  type(field_collection_type) :: self

  self%field_list = linked_list_type()

end function field_collection_constructor


!> Adds a field to the collection. The field maintained in the collection will
!> be a copy of the original.
!> @param [in] field The field that is to be copied into the collection
subroutine add_field(self, field)

  implicit none

  class(field_collection_type), intent(inout) :: self
  type(field_type), intent(in) :: field

  call self%field_list%insert_item( field )

end subroutine add_field


!> Access a field from the collection
!> @param [in] field_name The name of the field to be accessed
!> @return field Pointer to the field that is extracted
function get_field(self, field_name) result(field)

  implicit none

  class(field_collection_type), intent(inout) :: self

  character(*), intent(in) :: field_name
  type(field_type), pointer :: field

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type),pointer :: loop => null()

  ! start at the head of the mesh collection linked list
  loop => self%field_list%get_head()

  do
    ! If list is empty or we're at the end of list and we didn't find the
    ! field, return a null pointer
    if ( .not. associated(loop) ) then
      nullify(field)
      exit
    end if
    ! otherwise search list for the name of field we want

    ! 'cast' to the field_type 
    select type(f => loop%payload)
      type is (field_type)
      if ( trim(field_name) == trim(f%get_name()) ) then
          field => f
          exit
      end if
    end select

    loop => loop%next
  end do

end function get_field


!> Clears all items from the field collection linked list
subroutine clear(self)

  implicit none

  class(field_collection_type), intent(inout) :: self

  call self%field_list%clear()

  return
end subroutine clear

!> Destructor for the field collection
subroutine field_collection_destructor(self)

  implicit none

  type (field_collection_type), intent(inout) :: self

  call self%clear()

  return
end subroutine field_collection_destructor

end module field_collection_mod
