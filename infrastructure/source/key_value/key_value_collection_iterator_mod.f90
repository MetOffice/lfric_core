!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Provides functionality for iterating over all members of a key_value
!>        collection
!>
!> @details Provides functionality for iteratively returning every member
!>          of a key_value collection. The order of the key_values returned is
!>          not defined and can change if the implementation of the key_value
!>          collection is changed
!
module key_value_collection_iterator_mod

  use constants_mod,            only: l_def
  use key_value_mod,            only: key_value_type
  use linked_list_mod,          only: linked_list_item_type
  use key_value_collection_mod, only: key_value_collection_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Type that iterates through a key_value collection
  !-----------------------------------------------------------------------------
  type, public :: key_value_collection_iterator_type
    private
    !> A pointer to the key_value collection being iterated over
    type(key_value_collection_type), pointer :: collection
    !> A pointer to the linked list item within the collection that will
    !> contain the next key_value to be returned
    type(linked_list_item_type), pointer :: current
  contains
    procedure, public :: initialise
    procedure, public :: next
    procedure, public :: has_next
  end type key_value_collection_iterator_type

contains

!> Initialise a key_value collection iterator
!> @param [in] collection The collection to iterate over
subroutine initialise(self, collection)

  implicit none

  class(key_value_collection_iterator_type) :: self
  type(key_value_collection_type), target :: collection

  ! Store a pointer to the collection being iterated over
  self%collection => collection

  ! Start the iterator at the beginning of the key_value list.
  nullify(self%current)
  self%current => self%collection%get_next_item(self%current)

end subroutine initialise

!> Returns the next key_value from the collection
!> @return A polymorphic key_value pointer to the key_value that is next in
!>         the collection
function next(self) result (key_value)

  implicit none

  class(key_value_collection_iterator_type), intent(inout), target :: self
  class(key_value_type), pointer :: key_value

  ! Empty lists are valid
  !
  if (.not. associated(self%current)) then
    key_value => null()
    return
  end if

  ! Extract a pointer to the current key_value in the collection
  select type(listkey_value => self%current%payload)
    class is (key_value_type)
      key_value => listkey_value
  end select

  ! Move the current item pointer onto the next key_value in the collection
  self%current => self%collection%get_next_item(self%current)

end function next

!> Checks if there are any further key_values in the collection being iterated over
!> @return next true if there is another key_value in the collection, and false if
!> there isn't.
function has_next(self) result(next)
  implicit none
  class(key_value_collection_iterator_type), intent(in) :: self
  logical(l_def) :: next
  next = .true.
  if(.not.associated(self%current)) next = .false.
end function has_next

end module key_value_collection_iterator_mod
