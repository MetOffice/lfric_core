!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Provides functionality for iterating over all members of a defined
!>        namelist (bar) collection.
!>
!> @details Provides functionality for iteratively returning every member
!>          of the defined namelist (bar) collection. The order of
!>          the namelists returned is not defined and can change if the
!>          implementation of the namelist collection is changed.
!
module bar_nml_iterator_mod

  use constants_mod,   only: l_def
  use linked_list_mod, only: linked_list_type, &
                             linked_list_item_type

  use bar_nml_mod, only: bar_nml_type

  implicit none

  private
  public :: bar_nml_iterator_type

  !-----------------------------------------------------------------------------
  ! Type that iterates through a linked list of only bar_nml_type
  !-----------------------------------------------------------------------------
  type :: bar_nml_iterator_type

    private

    !> A pointer to the namelist list being iterated over
    type(linked_list_type), pointer :: bar_list

    !> A pointer to the linked list item within the
    !> linked list that will contain the next namelist
    !> to be returned
    type(linked_list_item_type), pointer :: current

  contains

    procedure, public :: initialise
    procedure, public :: next
    procedure, public :: has_next

  end type bar_nml_iterator_type

contains

!> @brief Initialise a bar namelist collection iterator
!> @param [in] nml_list Linked list containing only
!>                      bar_nml_types to iterate over.
subroutine initialise(self, nml_list)

  implicit none

  class(bar_nml_iterator_type), intent(inout) :: self
  type(linked_list_type), intent(in), target :: nml_list

  ! Store a pointer to the collection being iterated over
  self%bar_list => nml_list

  ! Start the iterator at the beginning of the nml_list.
  nullify(self%current)
  self%current => self%bar_list%get_head()

end subroutine initialise

!> @brief  Returns the next bar namelist from the collection
!> @return A pointer to the next bar namelist in the collection
function next(self) result (nml_obj)

  implicit none

  class(bar_nml_iterator_type), intent(inout), target :: self
  type(bar_nml_type), pointer :: nml_obj

  nml_obj => null()

  ! Empty lists are valid
  !
  if (.not. associated(self%current)) return

  ! Extract a pointer to the current namelist
  select type(list_nml => self%current%payload)
    type is (bar_nml_type)
      nml_obj => list_nml
  end select

  ! Move the current item pointer onto the next item
  self%current => self%current%next

end function next

!> @brief Checks if there are any further namelists in the collection
!>        being iterated over.
!> @return next  .true. if there is another namelist in the collection.
function has_next(self) result(next)

  implicit none

  class(bar_nml_iterator_type), intent(in) :: self
  logical(l_def) :: next

  next = .true.
  if (.not.associated(self%current)) next = .false.

end function has_next

end module bar_nml_iterator_mod
