!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!>
!> @brief Holds and manages function spaces created during a model run.
!>
!> @details A container which holds type definition of a collection of
!>          function spaces. The collection holds function spaces as 
!>          singletons. It will handle the creation and storing of 
!>          requested function spaces.
!
module function_space_collection_mod

  use constants_mod,      only: i_def, l_def
  use function_space_mod, only: function_space_type
  use fs_continuity_mod,  only: W0, W1, W2, W3, Wtheta, W2V, W2H, Wchi, fs_name
  use log_mod,            only: log_event, log_scratch_space                  &
                              , LOG_LEVEL_ERROR, LOG_LEVEL_TRACE
  use linked_list_mod,    only : linked_list_type, &
                                 linked_list_item_type
  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Type which is is a collection of function spaces held in a linked list
  !-----------------------------------------------------------------------------
  type, public :: function_space_collection_type
    private
    type(linked_list_type)        :: fs_list
  contains
    procedure, public :: get_fs
    procedure, public :: get_fs_collection_size
    procedure, public :: clear
    final             :: function_space_collection_destructor
  end type function_space_collection_type
  !-----------------------------------------------------------------------------

  interface function_space_collection_type
    module procedure function_space_collection_constructor
  end interface

  ! Module level variable to make the function space collection
  ! globally available
  type(function_space_collection_type), public    ::    function_space_collection

contains
!-----------------------------------------------------------------------------
! Construct the function space collection
!-----------------------------------------------------------------------------
!> Function to construct a function space collection

type(function_space_collection_type) function &
                                     function_space_collection_constructor()

  function_space_collection_constructor%fs_list = linked_list_type()

end function function_space_collection_constructor


!-----------------------------------------------------------------------------
! Get or create a function space
!-----------------------------------------------------------------------------
!> Function to get an instance of a function space from the linked list
!> or create it if it doesn't exist
function get_fs(self, mesh_id, element_order, dynamo_fs) result(fs)

  implicit none

  class(function_space_collection_type), intent(inout) :: self
  integer(i_def), intent(in)                           :: mesh_id
  integer(i_def), intent(in)                           :: element_order
  integer(i_def), intent(in)                           :: dynamo_fs

  type(function_space_type), pointer      :: fs
  type(linked_list_item_type),pointer     :: loop ! temp pointer for looping

  integer(i_def) :: fs_id

  select case (dynamo_fs)

  case (W0,W1,W2,W3,WTHETA,W2V,W2H, WCHI)
  case default
    write(log_scratch_space, '(2(A,I0),A)')                                    &
      'Function space type not defined for Dynamo. Available types are '     //&
      '[W0 | W1 | W2 | W3 | WTHETA | W2V | W2H | WCHI]'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    return

  end select

  if (element_order < 0) then
    write(log_scratch_space, '(A,I0)')                                         &
      'Function space element order must be >= 0   ',element_order
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    return
  end if

  ! Generate id for requested function space
  ! can use the passed mesh_id
  fs_id = 1000000*mesh_id + (1000*element_order) + dynamo_fs

  ! point at the head of the fs linked list
  loop => self%fs_list%get_head()


  ! loop through list
  do
    if ( .not. associated(loop) ) then
      ! List is empty or at the end of list and we didn't find it
      ! so create new function space and add it 
      call self%fs_list%insert_item(                                           &
                            function_space_type(mesh_id, element_order, dynamo_fs))

      write(log_scratch_space, '(A,2(I0,A))')                                  &
        'Generated order-',element_order,' '//trim(fs_name(dynamo_fs))//  &
        '-function space singleton (', fs_id,')'
      call log_event(log_scratch_space, LOG_LEVEL_TRACE)


      ! At this point the desired function space is the tail of the list
      ! so just retrieve it and exit loop

      loop => self%fs_list%get_tail()

      ! 'cast' to the function_space_type
      select type(v => loop%payload)
        type is (function_space_type)
            fs => v
      end select
      exit
    
    end if
    ! otherwise search list for the id we want
    if ( fs_id == loop%payload%get_id() ) then
      ! 'cast' to the function_space_type 
      select type(v => loop%payload)
        type is (function_space_type)
          fs => v
      end select
      exit
    end if
    loop => loop%next
  end do

  return
end function get_fs

!----------------------------------------------------------------------------
! Get the size of the function space collection
! (only really used in unit tests)
!-----------------------------------------------------------------------------
!> Function to return the number of function spaces currently
!> held in the collection

function get_fs_collection_size(self) result(fs_list_length)

  class(function_space_collection_type), intent(in)   :: self

  integer(i_def) :: fs_list_length

  fs_list_length = self%fs_list%get_length()
  
  return

end function get_fs_collection_size


!-----------------------------------------------------------------------------
! Clear the function space collection
!-----------------------------------------------------------------------------
!> Function to clear all items from the function space collection
!> linked list
subroutine clear(self)

  implicit none

  class(function_space_collection_type), intent(inout) :: self

  call self%fs_list%clear()

end subroutine clear

!-----------------------------------------------------------------------------
! Function space collection destructor
!-----------------------------------------------------------------------------

subroutine function_space_collection_destructor(self)

  implicit none

  type (function_space_collection_type), intent(inout) :: self

  call self%clear()

end subroutine function_space_collection_destructor


end module function_space_collection_mod
