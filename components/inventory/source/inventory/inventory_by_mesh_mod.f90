!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! GitHub Copilot (Claude Sonnet 5).
!-----------------------------------------------------------------------------
!
!> @brief Holds and manages objects that are paired with meshes in an inventory
!>
!> @details A container that holds a collection of objects that are paired with
!!          meshes.
!
module inventory_by_mesh_mod

  use constants_mod,                    only: i_def, l_def, str_def, r_single, &
                                              r_double
  use field_mod,                        only: field_type
  use field_real32_mod,                 only: field_real32_type
  use field_real64_mod,                 only: field_real64_type
  use function_space_mod,               only: function_space_type
  use integer_field_mod,                only: integer_field_type
  use log_mod,                          only: log_event, log_scratch_space, &
                                              LOG_LEVEL_ERROR,              &
                                              LOG_LEVEL_DEBUG
  use linked_list_data_mod,             only: linked_list_data_type
  use linked_list_mod,                  only: linked_list_type, &
                                              linked_list_item_type
  use mesh_mod,                         only: mesh_type
  use operator_mod,                     only: operator_type
  use operator_real32_mod,              only: operator_real32_type
  use operator_real64_mod,              only: operator_real64_type
  use id_abstract_pair_mod,             only: id_abstract_pair_type
  use id_r32_field_pair_mod,            only: id_r32_field_pair_type
  use id_r64_field_pair_mod,            only: id_r64_field_pair_type
  use id_integer_field_pair_mod,        only: id_integer_field_pair_type
  use id_r32_field_array_pair_mod,      only: id_r32_field_array_pair_type
  use id_r64_field_array_pair_mod,      only: id_r64_field_array_pair_type
  use id_integer_field_array_pair_mod,  only: id_integer_field_array_pair_type
  use id_r32_field_array_2d_pair_mod,   only: id_r32_field_array_2d_pair_type
  use id_r64_field_array_2d_pair_mod,   only: id_r64_field_array_2d_pair_type
  use id_integer_pair_mod,              only: id_integer_pair_type
  use id_integer_array_pair_mod,        only: id_integer_array_pair_type
  use id_real32_pair_mod,               only: id_real32_pair_type
  use id_real64_pair_mod,               only: id_real64_pair_type
  use id_logical_pair_mod,              only: id_logical_pair_type
  use id_r32_operator_pair_mod,         only: id_r32_operator_pair_type
  use id_r64_operator_pair_mod,         only: id_r64_operator_pair_type

  implicit none

  private

  ! Set the default table length -- this number should be higher than the
  ! expected number of meshes. If this limit is actually reached, then
  ! it can be increased in the future
  integer(kind=i_def), parameter :: default_table_len = 20

  !-----------------------------------------------------------------------------
  ! Type that holds an inventory of paired objects in a linked list
  !-----------------------------------------------------------------------------
  type, extends(linked_list_data_type), public :: inventory_by_mesh_type
    private
    !> The name of the inventory if provided.
    character(str_def)     :: name = 'unnamed_inventory'
    !> A hash table of linked lists of objects contained within the databse
    type(linked_list_type), allocatable :: paired_object_list(:)
    !> The size of the hash table to use.
    ! (Default to a value that represents an uninitialised hash table)
    integer(kind=i_def) :: table_len = 0
  contains
    ! Generic routines -- the same between different types of inventory
    procedure, public :: initialise
    procedure, public :: add_paired_object
    procedure, public :: get_length
    procedure, public :: get_name
    procedure, public :: get_table_len
    procedure, public :: is_initialised
    procedure, public :: clear
    ! Routines that need to discriminate between different types of object
    ! To support new objects, these need editing
    procedure, public :: get_paired_object
    procedure, public :: paired_object_exists
    procedure, public :: remove_paired_object
    ! Specific routines for adding different objects to the inventory
    ! To support new objects, add more routines here
    procedure, public :: add_r32_field
    procedure, public :: add_r64_field
    procedure, public :: add_integer_field
    procedure, public :: add_r32_intermesh_field
    procedure, public :: add_r64_intermesh_field
    generic           :: add_field => add_r32_field,           &
                                      add_r64_field,           &
                                      add_integer_field,       &
                                      add_r32_intermesh_field, &
                                      add_r64_intermesh_field
    procedure, public :: add_r32_field_array
    procedure, public :: add_r64_field_array
    procedure, public :: add_integer_field_array
    procedure, public :: add_r32_field_array_2d
    procedure, public :: add_r64_field_array_2d
    generic           :: add_field_array => add_r32_field_array,     &
                                            add_r64_field_array,     &
                                            add_integer_field_array, &
                                            add_r32_field_array_2d,  &
                                            add_r64_field_array_2d
    procedure, public :: add_integer
    procedure, public :: add_integer_array
    procedure, public :: add_logical
    procedure, public :: add_real32
    procedure, public :: add_real64
    generic           :: add_real => add_real32, &
                                     add_real64
    procedure, public :: add_r32_operator
    procedure, public :: add_r64_operator
    generic           :: add_operator => add_r32_operator, &
                                         add_r64_operator
    ! Specific routines for copying different objects into the inventory
    ! To support new objects, add more routines here
    procedure, public :: copy_r32_field
    procedure, public :: copy_r64_field
    procedure, public :: copy_integer_field
    procedure, public :: copy_r32_intermesh_field
    procedure, public :: copy_r64_intermesh_field
    generic           :: copy_field => copy_r32_field,           &
                                       copy_r64_field,           &
                                       copy_integer_field,       &
                                       copy_r32_intermesh_field, &
                                       copy_r64_intermesh_field
    procedure, public :: copy_r32_field_array
    procedure, public :: copy_r64_field_array
    procedure, public :: copy_integer_field_array
    procedure, public :: copy_r32_field_array_2d
    procedure, public :: copy_r64_field_array_2d
    generic           :: copy_field_array => copy_r32_field_array,     &
                                             copy_r64_field_array,     &
                                             copy_integer_field_array, &
                                             copy_r32_field_array_2d,  &
                                             copy_r64_field_array_2d
    procedure, public :: copy_r32_operator
    procedure, public :: copy_r64_operator
    generic           :: copy_operator => copy_r32_operator, &
                                          copy_r64_operator
    ! Specific routines for getting different objects from the inventory
    ! To support new objects, add more routines here
    procedure, public :: get_r32_field
    procedure, public :: get_r64_field
    procedure, public :: get_integer_field
    procedure, public :: get_r32_intermesh_field
    procedure, public :: get_r64_intermesh_field
    generic           :: get_field => get_r32_field,           &
                                      get_r64_field,           &
                                      get_integer_field,       &
                                      get_r32_intermesh_field, &
                                      get_r64_intermesh_field
    procedure, public :: get_r32_field_array
    procedure, public :: get_r64_field_array
    procedure, public :: get_integer_field_array
    procedure, public :: get_r32_field_array_2d
    procedure, public :: get_r64_field_array_2d
    generic           :: get_field_array => get_r32_field_array,     &
                                            get_r64_field_array,     &
                                            get_integer_field_array, &
                                            get_r32_field_array_2d,  &
                                            get_r64_field_array_2d
    procedure, public :: get_integer
    procedure, public :: get_integer_array
    procedure, public :: get_logical
    procedure, public :: get_real32
    procedure, public :: get_real64
    generic           :: get_real => get_real32, &
                                     get_real64
    procedure, public :: get_r32_operator
    procedure, public :: get_r64_operator
    generic           :: get_operator => get_r32_operator, &
                                         get_r64_operator

    ! Overloaded routine (which will trigger a failure)
    procedure, private :: inventory_copy_constructor

    procedure, public :: compute_intermesh_id

    generic, public   :: assignment(=) => inventory_copy_constructor
    final             :: inventory_by_mesh_destructor
  end type inventory_by_mesh_type

contains

! ============================================================================ !
! GENERIC INVENTORY ROUTINES -- these appear in each type of inventory
! ============================================================================ !
!> @brief Initialises an inventory
!> @param[in] name The name given to the inventory
subroutine initialise(self, name, table_len)

  implicit none

  class(inventory_by_mesh_type), intent(inout) :: self
  character(*),        optional, intent(in)    :: name
  integer(kind=i_def), optional, intent(in)    :: table_len

  integer(kind=i_def) :: i

  if( present(table_len)) then
    self%table_len = table_len
  else
    self%table_len = default_table_len
  end if

  ! Create the hash table for the inventory
  allocate(self%paired_object_list(0:self%table_len-1))
  do i = 0, self%table_len-1
    self%paired_object_list(i) = linked_list_type()
  end do

  if (present(name)) self%name = trim(name)

end subroutine initialise

!> @brief Adds a paired object to the inventory.
!> @param[in] paired_object  The object that is to be copied into the inventory
subroutine add_paired_object(self, paired_object)

  implicit none

  class(inventory_by_mesh_type), intent(inout) :: self
  class(id_abstract_pair_type),  intent(in)    :: paired_object
  integer(kind=i_def) :: hash, id

  ! Extract ID
  id = paired_object%get_id()

  ! Check if object exists in inventory already
  ! If it does, throw an error
  if ( self%paired_object_exists( id ) ) then
    write(log_scratch_space, '(A,I8,5A)')                                      &
        'Paired object on mesh [', id, '] corresponds to a hash that',         &
        ' already exists in inventory_by_mesh: ', trim(self%name),             &
        '. If this object corresponds to a new mesh, you may need to ',        &
        'increase the table length of the inventory'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    call self%remove_paired_object(id)
  end if

  ! Finished checking - so the object must be good to add - so add it
  hash = mod(id, self%get_table_len())
  call self%paired_object_list(hash)%insert_item( paired_object )
  write(log_scratch_space, '(A,I8,2A)') &
    'Adding object on mesh [', id, '] to inventory_by_mesh: ', trim(self%name)
  call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

end subroutine add_paired_object

!> @brief Returns the number of entries in the inventory
!> @return The number of entries in the inventory
function get_length(self) result(length)

  implicit none

  class(inventory_by_mesh_type), intent(in) :: self
  integer(kind=i_def) :: length
  integer(kind=i_def) :: i

  length = 0
  do i = 0, self%get_table_len()-1
    length = length + self%paired_object_list(i)%get_length()
  end do

end function get_length

!> @brief Returns the name of the inventory
!> @return The name of the inventory
function get_name(self) result(name)

  implicit none

  class(inventory_by_mesh_type), intent(in) :: self
  character(str_def) :: name

  name = self%name

end function get_name

!> @brief Returns the length of hash table for this inventory
!> @return The length of the inventory's hash table
function get_table_len(self) result(table_len)

  implicit none

  class(inventory_by_mesh_type), intent(in) :: self
  integer(kind=i_def) :: table_len

  if ( self%table_len == 0 ) then
    call log_event("inventory_by_mesh: Attempt to use uninitialised collection", &
                    LOG_LEVEL_ERROR)
  end if

  table_len = self%table_len

end function get_table_len

!> @brief Returns whether the inventory has been initialised
!> @return Logical for whether the inventory is initialised
function is_initialised(self) result(initialised)

  implicit none

  class(inventory_by_mesh_type), intent(in) :: self
  logical(kind=l_def) :: initialised

  initialised = allocated(self%paired_object_list)

end function is_initialised

!> DEPRECATED: Assignment operator between inventory_by_mesh_type pairs.
!> Currently, this routine generates a (hopefully) useful message, then
!> forces an error.
!>
!> @param[out] self
!> @param[in]  source
subroutine inventory_copy_constructor(self, source)

  implicit none
  class(inventory_by_mesh_type), intent(inout) :: self
  type(inventory_by_mesh_type),  intent(in)    :: source

  write(log_scratch_space,'(A,A)')&
     '"inventory_by_mesh2 = inventory_by_mesh1" syntax not supported. '  // &
     'Use "call inventory_by_mesh1%copy_inventory(inventory_by_mesh2)". '// &
     'Inventory: ', source%get_name()
  call log_event( log_scratch_space, LOG_LEVEL_ERROR )

end subroutine inventory_copy_constructor

!> @brief Clears all items from the inventory
subroutine clear(self)

  implicit none

  class(inventory_by_mesh_type), intent(inout) :: self
  integer(i_def) :: i

  if (allocated(self%paired_object_list)) then
    do i = 0, self%get_table_len()-1
      call self%paired_object_list(i)%clear()
    end do
    deallocate(self%paired_object_list)
  end if

  return

end subroutine clear

!> @brief Destructor for the inventory
subroutine inventory_by_mesh_destructor(self)

  implicit none

  type (inventory_by_mesh_type), intent(inout) :: self

  call self%clear()

  return

end subroutine inventory_by_mesh_destructor

! ============================================================================ !
! ROUTINES TO CYCLE THROUGH LINKED LIST
! ============================================================================ !

!> @brief Access a paired object from the inventory
!> @param[in] id             The ID of the object to be accessed
!> @returns   paired_object  Pointer to the paired_object to be set
function get_paired_object(self, id) result(paired_object)

  implicit none

  class(inventory_by_mesh_type), intent(in) :: self
  integer(kind=i_def),           intent(in) :: id

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type),  pointer :: loop => null()
  class(id_abstract_pair_type), pointer :: paired_object
  integer(kind=i_def) :: hash, loop_id

  hash = mod(id, self%get_table_len())

  ! start at the head of the inventory linked list
  loop => self%paired_object_list(hash)%get_head()

  do
    ! If inventory is empty or we're at the end of list and we didn't find the
    ! object, fail with an error
    if ( .not. associated(loop) ) then
      write(log_scratch_space, '(A,I8,2A)') 'get_paired_object: No object on mesh [', &
         id, '] in inventory_by_mesh: ', trim(self%name)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if

    ! otherwise search list for the ID of object we want
    select type(list_paired_object => loop%payload)
      type is (id_r32_field_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      type is (id_r64_field_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      type is (id_integer_field_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      type is (id_r32_operator_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      type is (id_r64_operator_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      type is (id_r32_field_array_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      type is (id_r64_field_array_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      type is (id_r32_field_array_2d_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      type is (id_r64_field_array_2d_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      type is (id_integer_field_array_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      type is (id_integer_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      type is (id_integer_array_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      type is (id_logical_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      type is (id_real32_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      type is (id_real64_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          paired_object => list_paired_object
          exit
        end if
      class default
        call log_event('Type of ID paired object not supported', LOG_LEVEL_ERROR)
    end select

    loop => loop%next
  end do

end function get_paired_object

!> @brief Check if a paired object is present in the inventory
!> @param[in] id The ID for the paired object
!> @return exists Flag stating if paired object is present or not
function paired_object_exists(self, id) result(exists)

  implicit none

  class(inventory_by_mesh_type), intent(in) :: self

  integer(kind=i_def), intent(in) :: id
  integer(kind=i_def)             :: loop_id
  logical(kind=l_def)             :: exists
  integer(kind=i_def)             :: hash

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type),  pointer :: loop => null()

  ! Calculate the hash of the object being searched for
  hash = mod(id, self%get_table_len())

  ! start at the head of the inventory linked list
  loop => self%paired_object_list(hash)%get_head()

  do
    ! If list is empty or we're at the end of list and we didn't find the
    ! object, set 'exists' to be false
    if ( .not. associated(loop) ) then
      exists=.false.
      exit
    end if

    ! otherwise search list for the object we want
    select type(list_paired_object => loop%payload)
      type is (id_r32_field_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      type is (id_r64_field_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      type is (id_integer_field_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      type is (id_r32_operator_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      type is (id_r64_operator_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      type is (id_r32_field_array_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      type is (id_r64_field_array_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      type is (id_r32_field_array_2d_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      type is (id_r64_field_array_2d_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      type is (id_integer_field_array_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      type is (id_integer_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      type is (id_integer_array_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      type is (id_logical_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      type is (id_real32_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      type is (id_real64_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          exists=.true.
          exit
        end if
      class default
        call log_event('Type of ID paired object not supported', LOG_LEVEL_ERROR)
    end select

    loop => loop%next
  end do

end function paired_object_exists

!> @brief Remove a paired_id type from the inventory
!> @param[in] id The ID of the paired object to be removed
subroutine remove_paired_object(self, id)

  implicit none

  class(inventory_by_mesh_type), intent(inout) :: self
  integer(kind=i_def),           intent(in)    :: id

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type),  pointer :: loop => null()
  integer(kind=i_def) :: hash, loop_id

  ! Calculate the hash of the field being removed
  hash = mod(id, self%get_table_len())

  ! start at the head of the inventory linked list
  loop => self%paired_object_list(hash)%get_head()

  do
    ! If list is empty or we're at the end of list and we didn't find the
    ! object, fail with an error
    if ( .not. associated(loop) ) then
      write(log_scratch_space, '(A,I8,2A)') 'remove_paired_object: No object paired to mesh [', &
         id, '] in the inventory: ', trim(self%name)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR)
    end if

    ! otherwise search list for the object we want
    select type(list_paired_object => loop%payload)
      type is (id_r32_field_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      type is (id_r64_field_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      type is (id_integer_field_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      type is (id_r32_operator_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      type is (id_r64_operator_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      type is (id_r32_field_array_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      type is (id_r64_field_array_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      type is (id_r32_field_array_2d_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      type is (id_r64_field_array_2d_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      type is (id_integer_field_array_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      type is (id_integer_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      type is (id_integer_array_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      type is (id_logical_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      type is (id_real32_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      type is (id_real64_pair_type)
        loop_id = list_paired_object%get_id()
        if ( id == loop_id ) then
          call self%paired_object_list(hash)%remove_item(loop)
          exit
        end if
      class default
        call log_event('Type of ID paired object not supported', LOG_LEVEL_ERROR)
    end select
    loop => loop%next
  end do

end subroutine remove_paired_object

! ============================================================================ !
! ADD OBJECT ROUTINES -- these are specific to an inventory_by_mesh
! ============================================================================ !

!> @brief Adds an r32 field to the inventory and returns a pointer to it
!> @param[out] field       Pointer to the field that is to be added
!> @param[in]  fs          Function space for the field to be created
!> @param[in]  mesh        The mesh to pair the field with
!> @param[in]  name        Optional name for the field
!> @param[in]  halo_depth  Optional halo depth for field (to overwrite the
!!                         default halo depth)
subroutine add_r32_field(self, field, fs, mesh, name, halo_depth)

  implicit none

  class(inventory_by_mesh_type),      intent(inout) :: self
  type(field_real32_type),   pointer, intent(out)   :: field
  type(function_space_type), pointer, intent(in)    :: fs
  type(mesh_type),                    intent(in)    :: mesh
  character(*),             optional, intent(in)    :: name
  integer(i_def),           optional, intent(in)    :: halo_depth
  type(id_r32_field_pair_type)                      :: paired_object
  character(str_def)                                :: local_name

  if ( present(name) ) then
    local_name = name
  else
    local_name = 'none'
  end if

  ! Set up the paired_object
  call paired_object%initialise(fs, mesh%get_id(), name=local_name, &
                                halo_depth=halo_depth)
  call self%add_paired_object(paired_object)
  call self%get_r32_field(mesh, field)

end subroutine add_r32_field

!> @brief Adds an r64 field to the inventory and returns a pointer to it
!> @param[out] field       Pointer to the field that is to be added
!> @param[in]  fs          Function space for the field to be created
!> @param[in]  mesh        The mesh to pair the field with
!> @param[in]  name        Optional name for the field
!> @param[in]  halo_depth  Optional halo depth for field (to overwrite the
!!                         default halo depth)
subroutine add_r64_field(self, field, fs, mesh, name, halo_depth)

  implicit none

  class(inventory_by_mesh_type),      intent(inout) :: self
  type(field_real64_type),   pointer, intent(out)   :: field
  type(function_space_type), pointer, intent(in)    :: fs
  type(mesh_type),                    intent(in)    :: mesh
  character(*),             optional, intent(in)    :: name
  integer(i_def),           optional, intent(in)    :: halo_depth
  type(id_r64_field_pair_type)                      :: paired_object
  character(str_def)                                :: local_name

  if ( present(name) ) then
    local_name = name
  else
    local_name = 'none'
  end if

  ! Set up the paired_object
  call paired_object%initialise(fs, mesh%get_id(), name=local_name, &
                                halo_depth=halo_depth)
  call self%add_paired_object(paired_object)
  call self%get_r64_field(mesh, field)

end subroutine add_r64_field

!> @brief Adds an r32 field to the inventory, that is used for intermesh mapping
!!        and returns a pointer to it
!> @param[out] field        Pointer to the field that is to be added
!> @param[in]  fs           Function space for the field to be created
!> @param[in]  source_mesh  The source mesh for the intermesh transform
!> @param[in]  target_mesh  The target mesh for the intermesh transform
subroutine add_r32_intermesh_field(self, field, fs, source_mesh, target_mesh)

  implicit none

  class(inventory_by_mesh_type),      intent(inout) :: self
  type(field_real32_type),   pointer, intent(out)   :: field
  type(function_space_type), pointer, intent(in)    :: fs
  type(mesh_type),                    intent(in)    :: source_mesh
  type(mesh_type),                    intent(in)    :: target_mesh
  type(id_r32_field_pair_type)                      :: paired_object
  integer(kind=i_def)                               :: intermesh_id

  ! Make a unique ID for transform between source and target mesh
  intermesh_id = self%compute_intermesh_id(source_mesh, target_mesh)
  ! Set up the paired_object
  call paired_object%initialise(fs, intermesh_id)
  call self%add_paired_object(paired_object)
  call self%get_r32_intermesh_field(source_mesh, target_mesh, field)

end subroutine add_r32_intermesh_field

!> @brief Adds an r64 field to the inventory, that is used for intermesh mapping
!!        and returns a pointer to it
!> @param[out] field        Pointer to the field that is to be added
!> @param[in]  fs           Function space for the field to be created
!> @param[in]  source_mesh  The source mesh for the intermesh transform
!> @param[in]  target_mesh  The target mesh for the intermesh transform
subroutine add_r64_intermesh_field(self, field, fs, source_mesh, target_mesh)

  implicit none

  class(inventory_by_mesh_type),      intent(inout) :: self
  type(field_real64_type),   pointer, intent(out)   :: field
  type(function_space_type), pointer, intent(in)    :: fs
  type(mesh_type),                    intent(in)    :: source_mesh
  type(mesh_type),                    intent(in)    :: target_mesh
  type(id_r64_field_pair_type)                      :: paired_object
  integer(kind=i_def)                               :: intermesh_id

  ! Make a unique ID for transform between source and target mesh
  intermesh_id = self%compute_intermesh_id(source_mesh, target_mesh)
  ! Set up the paired_object
  call paired_object%initialise(fs, intermesh_id)
  call self%add_paired_object(paired_object)
  call self%get_r64_intermesh_field(source_mesh, target_mesh, field)

end subroutine add_r64_intermesh_field

!> @brief Adds an integer field to the inventory and returns a pointer to it
!> @param[out] field        Pointer to the field that is to be added
!> @param[in]  fs           Function space for the field to be created
!> @param[in]  mesh         The mesh to pair the field with
!> @param[in]  name         Optional name for the field
!> @param[in]  halo_depth   Optional halo depth for field (to overwrite the
!!                          default halo depth)
subroutine add_integer_field(self, field, fs, mesh, name, halo_depth)

  implicit none

  class(inventory_by_mesh_type),      intent(inout) :: self
  type(integer_field_type),  pointer, intent(out)   :: field
  type(function_space_type), pointer, intent(in)    :: fs
  type(mesh_type),                    intent(in)    :: mesh
  character(*),             optional, intent(in)    :: name
  integer(i_def),           optional, intent(in)    :: halo_depth
  type(id_integer_field_pair_type)                  :: paired_object
  character(str_def)                                :: local_name

  if ( present(name) ) then
    local_name = name
  else
    local_name = 'none'
  end if

  ! Set up the paired_object
  call paired_object%initialise(fs, mesh%get_id(), name=local_name, &
                                halo_depth=halo_depth)
  call self%add_paired_object(paired_object)
  call self%get_integer_field(mesh, field)

end subroutine add_integer_field

!> @brief Adds an r32 field_array to the inventory and returns a pointer to it
!> @param[out] field_array  Pointer to the field_array that is to be added
!> @param[in]  fs           Function space for the fields to be created
!> @param[in]  array_size   Size of array of fields to be created
!> @param[in]  mesh         The mesh to pair the field_array with
!> @param[in]  halo_depth   Optional halo depth for field (to overwrite the
!!                          default halo depth)
subroutine add_r32_field_array(self, field_array, fs, array_size, mesh, halo_depth)

  implicit none

  class(inventory_by_mesh_type),      intent(inout) :: self
  type(field_real32_type),   pointer, intent(out)   :: field_array(:)
  type(function_space_type), pointer, intent(in)    :: fs
  integer(kind=i_def),                intent(in)    :: array_size
  type(mesh_type),                    intent(in)    :: mesh
  integer(i_def),           optional, intent(in)    :: halo_depth
  type(id_r32_field_array_pair_type)                :: paired_object

  ! Set up the paired_object
  call paired_object%initialise(fs, mesh%get_id(), array_size, &
                                halo_depth=halo_depth)
  call self%add_paired_object(paired_object)
  call self%get_r32_field_array(mesh, field_array)

end subroutine add_r32_field_array

!> @brief Adds an r64 field_array to the inventory and returns a pointer to it
!> @param[out] field_array  Pointer to the field_array that is to be added
!> @param[in]  fs           Function space for the fields to be created
!> @param[in]  array_size   Size of array of fields to be created
!> @param[in]  mesh         The mesh to pair the field_array with
!> @param[in]  halo_depth   Optional halo depth for field (to overwrite the
!!                          default halo depth)
subroutine add_r64_field_array(self, field_array, fs, array_size, mesh, halo_depth)

  implicit none

  class(inventory_by_mesh_type),      intent(inout) :: self
  type(field_real64_type),   pointer, intent(out)   :: field_array(:)
  type(function_space_type), pointer, intent(in)    :: fs
  integer(kind=i_def),                intent(in)    :: array_size
  type(mesh_type),                    intent(in)    :: mesh
  integer(i_def),           optional, intent(in)    :: halo_depth
  type(id_r64_field_array_pair_type)                :: paired_object

  ! Set up the paired_object
  call paired_object%initialise(fs, mesh%get_id(), array_size, &
                                  halo_depth=halo_depth)
  call self%add_paired_object(paired_object)
  call self%get_r64_field_array(mesh, field_array)

end subroutine add_r64_field_array

!> @brief Adds a 2D r32 field_array to the inventory and returns a pointer to it
!> @param[out] field_array  Pointer to the 2D field_array that is to be added
!> @param[in]  fs           Function space for the fields to be created
!> @param[in]  dim1         Size of the first dimension of the field array
!> @param[in]  dim2         Size of the second dimension of the field array
!> @param[in]  mesh         The mesh to pair the field_array with
!> @param[in]  halo_depth   Optional halo depth for field (to overwrite the
!!                          default halo depth)
subroutine add_r32_field_array_2d(self, field_array, fs, dim1, dim2, mesh, halo_depth)

  implicit none

  class(inventory_by_mesh_type),      intent(inout) :: self
  type(field_real32_type),   pointer, intent(out)   :: field_array(:,:)
  type(function_space_type), pointer, intent(in)    :: fs
  integer(kind=i_def),                intent(in)    :: dim1
  integer(kind=i_def),                intent(in)    :: dim2
  type(mesh_type),                    intent(in)    :: mesh
  integer(i_def),           optional, intent(in)    :: halo_depth
  type(id_r32_field_array_2d_pair_type)             :: paired_object

  ! Set up the paired_object
  call paired_object%initialise(fs, mesh%get_id(), dim1, dim2, &
                                halo_depth=halo_depth)
  call self%add_paired_object(paired_object)
  call self%get_r32_field_array_2d(mesh, field_array)

end subroutine add_r32_field_array_2d

!> @brief Adds a 2D r64 field_array to the inventory and returns a pointer to it
!> @param[out] field_array  Pointer to the 2D field_array that is to be added
!> @param[in]  fs           Function space for the fields to be created
!> @param[in]  dim1         Size of the first dimension of the field array
!> @param[in]  dim2         Size of the second dimension of the field array
!> @param[in]  mesh         The mesh to pair the field_array with
!> @param[in]  halo_depth   Optional halo depth for field (to overwrite the
!!                          default halo depth)
subroutine add_r64_field_array_2d(self, field_array, fs, dim1, dim2, mesh, halo_depth)

  implicit none

  class(inventory_by_mesh_type),      intent(inout) :: self
  type(field_real64_type),   pointer, intent(out)   :: field_array(:,:)
  type(function_space_type), pointer, intent(in)    :: fs
  integer(kind=i_def),                intent(in)    :: dim1
  integer(kind=i_def),                intent(in)    :: dim2
  type(mesh_type),                    intent(in)    :: mesh
  integer(i_def),           optional, intent(in)    :: halo_depth
  type(id_r64_field_array_2d_pair_type)             :: paired_object

  ! Set up the paired_object
  call paired_object%initialise(fs, mesh%get_id(), dim1, dim2, &
                                halo_depth=halo_depth)
  call self%add_paired_object(paired_object)
  call self%get_r64_field_array_2d(mesh, field_array)

end subroutine add_r64_field_array_2d

!> @brief Adds an integer field_array to the inventory and returns a pointer to it
!> @param[out] field_array  Pointer to the field_array that is to be added
!> @param[in]  fs           Function space for the fields to be created
!> @param[in]  array_size   Size of array of fields to be created
!> @param[in]  mesh         The mesh to pair the field_array with
!> @param[in]  halo_depth   Optional halo depth for field (to overwrite the
!!                          default halo depth)
subroutine add_integer_field_array(self, field_array, fs, array_size, mesh, halo_depth)

  implicit none

  class(inventory_by_mesh_type),      intent(inout) :: self
  type(integer_field_type),  pointer, intent(out)   :: field_array(:)
  type(function_space_type), pointer, intent(in)    :: fs
  integer(kind=i_def),                intent(in)    :: array_size
  type(mesh_type),                    intent(in)    :: mesh
  integer(i_def),           optional, intent(in)    :: halo_depth
  type(id_integer_field_array_pair_type)            :: paired_object

  ! Set up the paired_object
  call paired_object%initialise(fs, mesh%get_id(), array_size, &
                                halo_depth=halo_depth)
  call self%add_paired_object(paired_object)
  call self%get_integer_field_array(mesh, field_array)

end subroutine add_integer_field_array

!> @brief Adds an integer to the inventory
!> @param[in] number      The integer that is to be added to the inventory
!> @param[in] mesh        The mesh to pair the integer with
subroutine add_integer(self, number, mesh)

  implicit none

  class(inventory_by_mesh_type), intent(inout) :: self
  integer(kind=i_def),           intent(in)    :: number
  type(mesh_type),               intent(in)    :: mesh
  type(id_integer_pair_type)                   :: paired_object

  ! Set up the paired_object
  call paired_object%initialise(number, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine add_integer

!> @brief Adds an integer_array to the inventory
!> @param[in] numbers  The integer_array that is to be added to the inventory
!> @param[in] mesh     The mesh to pair the integer_array with
subroutine add_integer_array(self, numbers, mesh)

  implicit none

  class(inventory_by_mesh_type), intent(inout) :: self
  integer(kind=i_def),           intent(in)    :: numbers(:)
  type(mesh_type),               intent(in)    :: mesh
  type(id_integer_array_pair_type)             :: paired_object

  ! Set up the paired_object
  call paired_object%initialise(numbers, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine add_integer_array

!> @brief Adds a real32 to the inventory
!> @param[in] number      The real that is to be added to the inventory
!> @param[in] mesh        The mesh to pair the real with
subroutine add_real32(self, number, mesh)

  implicit none

  class(inventory_by_mesh_type), intent(inout) :: self
  real(kind=r_single),           intent(in)    :: number
  type(mesh_type),               intent(in)    :: mesh
  type(id_real32_pair_type)                    :: paired_object

  ! Set up the paired_object
  call paired_object%initialise(number, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine add_real32

!> @brief Adds a real64 to the inventory
!> @param[in] number      The real that is to be added to the inventory
!> @param[in] mesh        The mesh to pair the real with
subroutine add_real64(self, number, mesh)

  implicit none

  class(inventory_by_mesh_type), intent(inout) :: self
  real(kind=r_double),           intent(in)    :: number
  type(mesh_type),               intent(in)    :: mesh
  type(id_real64_pair_type)                    :: paired_object

  ! Set up the paired_object
  call paired_object%initialise(number, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine add_real64

!> @brief Adds a  logical to the inventory
!> @param[in] bool_flag   The logical that is to be added to the inventory
!> @param[in] mesh        The mesh to pair the logical with
subroutine add_logical(self, bool_flag, mesh)

  implicit none

  class(inventory_by_mesh_type), intent(inout) :: self
  logical(kind=i_def),           intent(in)    :: bool_flag
  type(mesh_type),               intent(in)    :: mesh
  type(id_logical_pair_type)                   :: paired_object

  ! Set up the paired_object
  call paired_object%initialise(bool_flag, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine add_logical

!> @brief Adds an r32 operator to the inventory and returns a pointer to it
!> @param[out] operator_out Pointer to the operator that is to be added
!> @param[in]  fs_target    Function space for the target field of the operator
!> @param[in]  fs_source    Function space for the source field of the operator
!> @param[in]  mesh         The mesh to pair the operator with
subroutine add_r32_operator(self, operator_out, fs_target, fs_source, mesh)

  implicit none

  class(inventory_by_mesh_type),      intent(inout) :: self
  type(operator_real32_type),pointer, intent(out)   :: operator_out
  type(function_space_type), pointer, intent(in)    :: fs_target
  type(function_space_type), pointer, intent(in)    :: fs_source
  type(mesh_type),                    intent(in)    :: mesh
  type(id_r32_operator_pair_type)                   :: paired_object

  ! Set up the paired_object
  call paired_object%initialise(fs_target, fs_source, mesh%get_id())
  call self%add_paired_object(paired_object)
  call self%get_r32_operator(mesh, operator_out)

end subroutine add_r32_operator

!> @brief Adds an r64 operator to the inventory and returns a pointer to it
!> @param[out] operator_out Pointer to the operator that is to be added
!> @param[in]  fs_target    Function space for the target field of the operator
!> @param[in]  fs_source    Function space for the source field of the operator
!> @param[in]  mesh         The mesh to pair the operator with
subroutine add_r64_operator(self, operator_out, fs_target, fs_source, mesh)

  implicit none

  class(inventory_by_mesh_type),      intent(inout) :: self
  type(operator_real64_type),pointer, intent(out)   :: operator_out
  type(function_space_type), pointer, intent(in)    :: fs_target
  type(function_space_type), pointer, intent(in)    :: fs_source
  type(mesh_type),                    intent(in)    :: mesh
  type(id_r64_operator_pair_type)                   :: paired_object

  ! Set up the paired_object
  call paired_object%initialise(fs_target, fs_source, mesh%get_id())
  call self%add_paired_object(paired_object)
  call self%get_r64_operator(mesh, operator_out)

end subroutine add_r64_operator

! ============================================================================ !
! COPY OBJECT ROUTINES -- these are specific to an inventory_by_mesh
! ============================================================================ !

!> @brief Copies an r32 field into the inventory
!> @param[in] field       The field that is to be copied into the inventory
!> @param[in] mesh        The mesh to pair the field with
subroutine copy_r32_field(self, field, mesh)

  implicit none

  class(inventory_by_mesh_type), intent(inout) :: self
  type(field_real32_type),          intent(in)    :: field
  type(mesh_type),               intent(in)    :: mesh
  type(id_r32_field_pair_type)                 :: paired_object

  ! Set up the paired_object
  call paired_object%copy_initialise(field, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine copy_r32_field

!> @brief Copies an r64 field into the inventory
!> @param[in] field       The field that is to be copied into the inventory
!> @param[in] mesh        The mesh to pair the field with
subroutine copy_r64_field(self, field, mesh)

  implicit none

  class(inventory_by_mesh_type), intent(inout) :: self
  type(field_real64_type),       intent(in)    :: field
  type(mesh_type),               intent(in)    :: mesh
  type(id_r64_field_pair_type)                 :: paired_object

  ! Set up the paired_object
  call paired_object%copy_initialise(field, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine copy_r64_field

!> @brief Copies an r32 field into the inventory, that is used for intermesh mapping
!> @param[in] field        The field that is to be copied into the inventory
!> @param[in] source_mesh  The source mesh for the intermesh transform
!> @param[in] target_mesh  The target mesh for the intermesh transform
subroutine copy_r32_intermesh_field(self, field, source_mesh, target_mesh)

  implicit none

  class(inventory_by_mesh_type), intent(inout) :: self
  type(field_real32_type),       intent(in)    :: field
  type(mesh_type),               intent(in)    :: source_mesh
  type(mesh_type),               intent(in)    :: target_mesh
  type(id_r32_field_pair_type)                 :: paired_object
  integer(kind=i_def)                          :: intermesh_id

  ! Make a unique ID for transform between source and target mesh
  intermesh_id = self%compute_intermesh_id(source_mesh, target_mesh)
  ! Set up the paired_object
  call paired_object%copy_initialise(field, intermesh_id)
  call self%add_paired_object(paired_object)

end subroutine copy_r32_intermesh_field

!> @brief Copies an r64 field into the inventory, that is used for intermesh mapping
!> @param[in] field        The field that is to be copied into the inventory
!> @param[in] source_mesh  The source mesh for the intermesh transform
!> @param[in] target_mesh  The target mesh for the intermesh transform
subroutine copy_r64_intermesh_field(self, field, source_mesh, target_mesh)

  implicit none

  class(inventory_by_mesh_type), intent(inout) :: self
  type(field_real64_type),       intent(in)    :: field
  type(mesh_type),               intent(in)    :: source_mesh
  type(mesh_type),               intent(in)    :: target_mesh
  type(id_r64_field_pair_type)                 :: paired_object
  integer(kind=i_def)                          :: intermesh_id

  ! Make a unique ID for transform between source and target mesh
  intermesh_id = self%compute_intermesh_id(source_mesh, target_mesh)
  ! Set up the paired_object
  call paired_object%copy_initialise(field, intermesh_id)
  call self%add_paired_object(paired_object)

end subroutine copy_r64_intermesh_field

!> @brief Copies an integer field into the inventory
!> @param[in] field       The field that is to be copied into the inventory
!> @param[in] mesh        The mesh to pair the field with
subroutine copy_integer_field(self, field, mesh)

  implicit none

  class(inventory_by_mesh_type), intent(inout) :: self
  type(integer_field_type),      intent(in)    :: field
  type(mesh_type),               intent(in)    :: mesh
  type(id_integer_field_pair_type)             :: paired_object

  ! Set up the paired_object
  call paired_object%copy_initialise(field, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine copy_integer_field

!> @brief Copies an r32 field_array into the inventory
!> @param[in] field_array The field_array that is to be copied into the inventory
!> @param[in] mesh        The mesh to pair the field_array with
subroutine copy_r32_field_array(self, field_array, mesh)

  implicit none

  class(inventory_by_mesh_type),     intent(inout) :: self
  type(field_real32_type),           intent(in)    :: field_array(:)
  type(mesh_type),                   intent(in)    :: mesh
  type(id_r32_field_array_pair_type)               :: paired_object

  ! Set up the paired_object
  call paired_object%copy_initialise(field_array, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine copy_r32_field_array

!> @brief Copies an r64 field_array into the inventory
!> @param[in] field_array The field_array that is to be copied into the inventory
!> @param[in] mesh        The mesh to pair the field_array with
subroutine copy_r64_field_array(self, field_array, mesh)

  implicit none

  class(inventory_by_mesh_type),     intent(inout) :: self
  type(field_real64_type),           intent(in)    :: field_array(:)
  type(mesh_type),                   intent(in)    :: mesh
  type(id_r64_field_array_pair_type)               :: paired_object

  ! Set up the paired_object
  call paired_object%copy_initialise(field_array, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine copy_r64_field_array

!> @brief Copies a 2D r32 field_array into the inventory
!> @param[in] field_array The 2D field_array that is to be copied into the inventory
!> @param[in] mesh        The mesh to pair the field_array with
subroutine copy_r32_field_array_2d(self, field_array, mesh)

  implicit none

  class(inventory_by_mesh_type),        intent(inout) :: self
  type(field_real32_type),              intent(in)    :: field_array(:,:)
  type(mesh_type),                      intent(in)    :: mesh
  type(id_r32_field_array_2d_pair_type)               :: paired_object

  ! Set up the paired_object
  call paired_object%copy_initialise(field_array, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine copy_r32_field_array_2d

!> @brief Copies a 2D r64 field_array into the inventory
!> @param[in] field_array The 2D field_array that is to be copied into the inventory
!> @param[in] mesh        The mesh to pair the field_array with
subroutine copy_r64_field_array_2d(self, field_array, mesh)

  implicit none

  class(inventory_by_mesh_type),        intent(inout) :: self
  type(field_real64_type),              intent(in)    :: field_array(:,:)
  type(mesh_type),                      intent(in)    :: mesh
  type(id_r64_field_array_2d_pair_type)               :: paired_object

  ! Set up the paired_object
  call paired_object%copy_initialise(field_array, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine copy_r64_field_array_2d

!> @brief Copies an integer field_array into the inventory
!> @param[in] field_array The field_array that is to be copied into the inventory
!> @param[in] mesh        The mesh to pair the field_array with
subroutine copy_integer_field_array(self, field_array, mesh)

  implicit none

  class(inventory_by_mesh_type),         intent(inout) :: self
  type(integer_field_type),              intent(in)    :: field_array(:)
  type(mesh_type),                       intent(in)    :: mesh
  type(id_integer_field_array_pair_type)               :: paired_object

  ! Set up the paired_object
  call paired_object%copy_initialise(field_array, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine copy_integer_field_array

!> @brief Copies an r32 operator into the inventory.
!> @param[in] operator_in The operator that is to be copied into the inventory
!> @param[in] mesh        The mesh to pair the operator with
subroutine copy_r32_operator(self, operator_in, mesh)

  implicit none

  class(inventory_by_mesh_type),     intent(inout) :: self
  type(operator_real32_type),        intent(in)    :: operator_in
  type(mesh_type),                   intent(in)    :: mesh
  type(id_r32_operator_pair_type)                  :: paired_object

  ! Set up the paired_object
  call paired_object%copy_initialise(operator_in, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine copy_r32_operator

!> @brief Copies an r64 operator into the inventory.
!> @param[in] operator_in The operator that is to be copied into the inventory
!> @param[in] mesh        The mesh to pair the operator with
subroutine copy_r64_operator(self, operator_in, mesh)

  implicit none

  class(inventory_by_mesh_type),     intent(inout) :: self
  type(operator_real64_type),        intent(in)    :: operator_in
  type(mesh_type),                   intent(in)    :: mesh
  type(id_r64_operator_pair_type)                  :: paired_object

  ! Set up the paired_object
  call paired_object%copy_initialise(operator_in, mesh%get_id())
  call self%add_paired_object(paired_object)

end subroutine copy_r64_operator

! ============================================================================ !
! GET OBJECT ROUTINES -- these are specific to an inventory_by_mesh
! ============================================================================ !

!> @brief Sets a pointer to an r32 field from the inventory
!> @param[in]  mesh    The mesh of the r32 field to be accessed
!> @param[out] field   Field pointer to the r32 field to be accessed
subroutine get_r32_field(self, mesh, field)

  implicit none

  class(inventory_by_mesh_type),    intent(in)  :: self
  type(mesh_type),                  intent(in)  :: mesh
  type(field_real32_type), pointer, intent(out) :: field
  class(id_abstract_pair_type), pointer         :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_r32_field_pair_type)
      field => this%get_field()
    class default
      call log_event('Paired ID object must be of r32 field type', LOG_LEVEL_ERROR)
  end select

end subroutine get_r32_field

!> @brief Sets a pointer to an r64 field from the inventory
!> @param[in]  mesh    The mesh of the r64 field to be accessed
!> @param[out] field   Field pointer to the r64 field to be accessed
subroutine get_r64_field(self, mesh, field)

  implicit none

  class(inventory_by_mesh_type),    intent(in)  :: self
  type(mesh_type),                  intent(in)  :: mesh
  type(field_real64_type), pointer, intent(out) :: field
  class(id_abstract_pair_type), pointer         :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_r64_field_pair_type)
      field => this%get_field()
    class default
      call log_event('Paired ID object must be of r64 field type', LOG_LEVEL_ERROR)
  end select

end subroutine get_r64_field

!> @brief Sets a pointer to an r32 field used for intermesh transforms
!> @param[in]  source_mesh   The source mesh for the intermesh transform
!> @param[in]  target_mesh   The target mesh for the intermesh transform
!> @param[out] field         Field pointer to the r32 field to be accessed
subroutine get_r32_intermesh_field(self, source_mesh, target_mesh, field)

  implicit none

  class(inventory_by_mesh_type),    intent(in)  :: self
  type(mesh_type),                  intent(in)  :: source_mesh
  type(mesh_type),                  intent(in)  :: target_mesh
  type(field_real32_type), pointer, intent(out) :: field
  class(id_abstract_pair_type), pointer         :: paired_object
  integer(kind=i_def)                           :: intermesh_id

  ! Make a unique ID for transform between source and target mesh
  intermesh_id = self%compute_intermesh_id(source_mesh, target_mesh)

  paired_object => self%get_paired_object(intermesh_id)

  select type(this => paired_object)
    type is (id_r32_field_pair_type)
      field => this%get_field()
    class default
      call log_event('Paired ID object must be of r32 field type', LOG_LEVEL_ERROR)
  end select

end subroutine get_r32_intermesh_field

!> @brief Sets a pointer to an r64 field used for intermesh transforms
!> @param[in]  source_mesh   The source mesh for the intermesh transform
!> @param[in]  target_mesh   The target mesh for the intermesh transform
!> @param[out] field         Field pointer to the r64 field to be accessed
subroutine get_r64_intermesh_field(self, source_mesh, target_mesh, field)

  implicit none

  class(inventory_by_mesh_type),    intent(in)  :: self
  type(mesh_type),                  intent(in)  :: source_mesh
  type(mesh_type),                  intent(in)  :: target_mesh
  type(field_real64_type), pointer, intent(out) :: field
  class(id_abstract_pair_type), pointer         :: paired_object
  integer(kind=i_def)                           :: intermesh_id

  ! Make a unique ID for transform between source and target mesh
  intermesh_id = self%compute_intermesh_id(source_mesh, target_mesh)

  paired_object => self%get_paired_object(intermesh_id)

  select type(this => paired_object)
    type is (id_r64_field_pair_type)
      field => this%get_field()
    class default
      call log_event('Paired ID object must be of r64 field type', LOG_LEVEL_ERROR)
  end select

end subroutine get_r64_intermesh_field

!> @brief Sets a pointer to an integer field from the inventory
!> @param[in]     mesh    The mesh of the integer field to be accessed
!> @param[out] field   Field pointer to the integer field to be accessed
subroutine get_integer_field(self, mesh, field)

  implicit none

  class(inventory_by_mesh_type),     intent(in)  :: self
  type(mesh_type),                   intent(in)  :: mesh
  type(integer_field_type), pointer, intent(out) :: field
  class(id_abstract_pair_type),      pointer    :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_integer_field_pair_type)
      field => this%get_field()
    class default
      call log_event('Paired ID object must be of integer field type', LOG_LEVEL_ERROR)
  end select

end subroutine get_integer_field

!> @brief Sets a pointer to an array of r32 fields from the inventory
!> @param[in]  mesh          The mesh of the r32 field array to be accessed
!> @param[out] field_array   Pointer to the array of r32 fields to be accessed
subroutine get_r32_field_array(self, mesh, field_array)

  implicit none

  class(inventory_by_mesh_type),    intent(in)  :: self
  type(mesh_type),                  intent(in)  :: mesh
  type(field_real32_type), pointer, intent(out) :: field_array(:)
  class(id_abstract_pair_type), pointer         :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_r32_field_array_pair_type)
      field_array => this%get_field_array()
    class default
      call log_event('Paired ID object must be of r32 field_array type', LOG_LEVEL_ERROR)
  end select

end subroutine get_r32_field_array

!> @brief Sets a pointer to an array of r64 fields from the inventory
!> @param[in]  mesh          The mesh of the r64 field array to be accessed
!> @param[out] field_array   Pointer to the array of r64 fields to be accessed
subroutine get_r64_field_array(self, mesh, field_array)

  implicit none

  class(inventory_by_mesh_type),    intent(in)  :: self
  type(mesh_type),                  intent(in)  :: mesh
  type(field_real64_type), pointer, intent(out) :: field_array(:)
  class(id_abstract_pair_type), pointer         :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_r64_field_array_pair_type)
      field_array => this%get_field_array()
    class default
      call log_event('Paired ID object must be of r64 field_array type', LOG_LEVEL_ERROR)
  end select

end subroutine get_r64_field_array

!> @brief Sets a pointer to a 2D array of r32 fields from the inventory
!> @param[in]  mesh          The mesh of the r32 field array to be accessed
!> @param[out] field_array   Pointer to the 2D array of r32 fields to be accessed
subroutine get_r32_field_array_2d(self, mesh, field_array)

  implicit none

  class(inventory_by_mesh_type),    intent(in)  :: self
  type(mesh_type),                  intent(in)  :: mesh
  type(field_real32_type), pointer, intent(out) :: field_array(:,:)
  class(id_abstract_pair_type), pointer         :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_r32_field_array_2d_pair_type)
      field_array => this%get_field_array()
    class default
      call log_event('Paired ID object must be of r32 2D field_array type', LOG_LEVEL_ERROR)
  end select

end subroutine get_r32_field_array_2d

!> @brief Sets a pointer to a 2D array of r64 fields from the inventory
!> @param[in]  mesh          The mesh of the r64 field array to be accessed
!> @param[out] field_array   Pointer to the 2D array of r64 fields to be accessed
subroutine get_r64_field_array_2d(self, mesh, field_array)

  implicit none

  class(inventory_by_mesh_type),    intent(in)  :: self
  type(mesh_type),                  intent(in)  :: mesh
  type(field_real64_type), pointer, intent(out) :: field_array(:,:)
  class(id_abstract_pair_type), pointer         :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_r64_field_array_2d_pair_type)
      field_array => this%get_field_array()
    class default
      call log_event('Paired ID object must be of r64 2D field_array type', LOG_LEVEL_ERROR)
  end select

end subroutine get_r64_field_array_2d

!> @brief Sets a pointer to an array of integer fields from the inventory
!> @param[in]  mesh          The mesh of the integer field array to be accessed
!> @param[out] field_array   Pointer to the array of integer fields to be accessed
subroutine get_integer_field_array(self, mesh, field_array)

  implicit none

  class(inventory_by_mesh_type),     intent(in)  :: self
  type(mesh_type),                   intent(in)  :: mesh
  type(integer_field_type), pointer, intent(out) :: field_array(:)
  class(id_abstract_pair_type),      pointer     :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_integer_field_array_pair_type)
      field_array => this%get_field_array()
    class default
      call log_event('Paired ID object must be of integer field_array type', LOG_LEVEL_ERROR)
  end select

end subroutine get_integer_field_array

!> @brief Sets a pointer to an integer from the inventory
!> @param[in]  mesh    The mesh of the integer to be accessed
!> @param[out] number  Pointer to the integer to be accessed
subroutine get_integer(self, mesh, number)

  implicit none

  class(inventory_by_mesh_type), intent(in)  :: self
  type(mesh_type),               intent(in)  :: mesh
  integer(kind=i_def),  pointer, intent(out) :: number
  class(id_abstract_pair_type),  pointer     :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_integer_pair_type)
      number => this%get_integer()
    class default
      call log_event('Paired ID object must be of integer type', LOG_LEVEL_ERROR)
  end select

end subroutine get_integer

!> @brief Sets a pointer to an integer_array from the inventory
!> @param[in]  mesh     The mesh of the integer array to be accessed
!> @param[out] numbers  Pointer to the integer_array to be accessed
subroutine get_integer_array(self, mesh, numbers)

  implicit none

  class(inventory_by_mesh_type), intent(in)  :: self
  type(mesh_type),               intent(in)  :: mesh
  integer(kind=i_def),  pointer, intent(out) :: numbers(:)
  class(id_abstract_pair_type),  pointer     :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_integer_array_pair_type)
      numbers => this%get_integer_array()
    class default
      call log_event('Paired ID object must be of integer type', LOG_LEVEL_ERROR)
  end select

end subroutine get_integer_array

!> @brief Sets a pointer to a real32 from the inventory
!> @param[in]  mesh    The mesh of the real to be accessed
!> @param[out] number  Pointer to the real to be accessed
subroutine get_real32(self, mesh, number)

  implicit none

  class(inventory_by_mesh_type), intent(in)  :: self
  type(mesh_type),               intent(in)  :: mesh
  real(kind=r_single),  pointer, intent(out) :: number
  class(id_abstract_pair_type),  pointer     :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_real32_pair_type)
      number => this%get_real32()
    class default
      call log_event('Paired ID object must be of real32 type', LOG_LEVEL_ERROR)
  end select

end subroutine get_real32

!> @brief Sets a pointer to a real64 from the inventory
!> @param[in]  mesh    The mesh of the real to be accessed
!> @param[out] number  Pointer to the real to be accessed
subroutine get_real64(self, mesh, number)

  implicit none

  class(inventory_by_mesh_type), intent(in)  :: self
  type(mesh_type),               intent(in)  :: mesh
  real(kind=r_double),  pointer, intent(out) :: number
  class(id_abstract_pair_type),  pointer     :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_real64_pair_type)
      number => this%get_real64()
    class default
      call log_event('Paired ID object must be of real64 type', LOG_LEVEL_ERROR)
  end select

end subroutine get_real64

!> @brief Sets a pointer to a logical from the inventory
!> @param[in]  mesh       The mesh of the logical to be accessed
!> @param[out] bool_flag  Pointer to the logical to be accessed
subroutine get_logical(self, mesh, bool_flag)

  implicit none

  class(inventory_by_mesh_type), intent(in)  :: self
  type(mesh_type),               intent(in)  :: mesh
  logical(kind=l_def),  pointer, intent(out) :: bool_flag
  class(id_abstract_pair_type),  pointer     :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_logical_pair_type)
      bool_flag => this%get_logical()
    class default
      call log_event('Paired ID object must be of logical type', LOG_LEVEL_ERROR)
  end select

end subroutine get_logical

!> @brief Sets a pointer to an r32 operator from the inventory
!> @param[in]  mesh          The mesh of the operator to be accessed
!> @param[out] operator_out  Pointer to be set to the desired operator
subroutine get_r32_operator(self, mesh, operator_out)

  implicit none

  class(inventory_by_mesh_type),       intent(in)  :: self
  type(mesh_type),                     intent(in)  :: mesh
  type(operator_real32_type), pointer, intent(out) :: operator_out
  class(id_abstract_pair_type), pointer            :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_r32_operator_pair_type)
      operator_out => this%get_operator()
    class default
      call log_event('Paired ID object must be of r32 operator type', LOG_LEVEL_ERROR)
  end select

end subroutine get_r32_operator

!> @brief Sets a pointer to an r64 operator from the inventory
!> @param[in]  mesh          The mesh of the operator to be accessed
!> @param[out] operator_out  Pointer to be set to the desired operator
subroutine get_r64_operator(self, mesh, operator_out)

  implicit none

  class(inventory_by_mesh_type),       intent(in)  :: self
  type(mesh_type),                     intent(in)  :: mesh
  type(operator_real64_type), pointer, intent(out) :: operator_out
  class(id_abstract_pair_type), pointer            :: paired_object

  paired_object => self%get_paired_object(mesh%get_id())

  select type(this => paired_object)
    type is (id_r64_operator_pair_type)
      operator_out => this%get_operator()
    class default
      call log_event('Paired ID object must be of r64 operator type', LOG_LEVEL_ERROR)
  end select

end subroutine get_r64_operator

!> @brief Computes a hash of two integers to create an intermesh ID
!> @param[in] source_mesh_id  ID of source mesh
!> @param[in] target_mesh_id  ID of target mesh
!> @returns   intermesh_id    ID for the combination of source and target meshes
function compute_intermesh_id(self, source_mesh, target_mesh) result(intermesh_id)

  implicit none

  class(inventory_by_mesh_type),    intent(in)  :: self
  type(mesh_type),                  intent(in)  :: source_mesh
  type(mesh_type),                  intent(in)  :: target_mesh

  integer(kind=i_def) :: source_mesh_id
  integer(kind=i_def) :: target_mesh_id
  integer(kind=i_def) :: intermesh_id

  source_mesh_id = source_mesh%get_id()
  target_mesh_id = target_mesh%get_id()

  if (source_mesh_id < 1) then
    write(log_scratch_space, '(A,I8,A)') &
      'source_mesh_id', source_mesh_id, 'is not a positive integer'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if
  if (target_mesh_id < 1) then
    write(log_scratch_space, '(A,I8,A)') &
      'target_mesh_id', target_mesh_id, 'is not a positive integer'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  ! Hash the source and target mesh IDs using the Hopcroft and Ullman pairing
  intermesh_id = (source_mesh_id + target_mesh_id - 2) * &
                 (source_mesh_id + target_mesh_id - 1) / 2 + source_mesh_id

end function compute_intermesh_id

end module inventory_by_mesh_mod
