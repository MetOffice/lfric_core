!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   Holds and manages a collection of fieldspec objects
!>
!> @details A singleton container which holds a collection of fieldspec
!>          objects that are used to manage information about fields that
!>          will be created during a run

module fieldspec_collection_mod

  use constants_mod,             only: i_def, str_def
  use fs_continuity_mod,         only: Wtheta
  use fieldspec_mod,             only: fieldspec_type
  use linked_list_mod,           only: linked_list_type, linked_list_item_type
  use log_mod,                   only: log_event, log_scratch_space, &
                                       LOG_LEVEL_ERROR, LOG_LEVEL_INFO

  implicit none

  private

  type, public :: fieldspec_collection_type

    private

    ! List of the global_mesh_type objects in this collection.
    type(linked_list_type) :: fieldspec_list

  contains
    procedure, public  :: add_fieldspec
    procedure, public  :: get_fieldspec
!    procedure, public  :: populate
    procedure, public  :: clear
    final              :: fieldspec_collection_destructor

  end type fieldspec_collection_type

  interface fieldspec_collection_type
    module procedure fieldspec_collection_constructor
  end interface

  !> @brief Singleton instance of a fieldspec_collection_type object.
  !>
  type(fieldspec_collection_type), public, allocatable :: &
      fieldspec_collection

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Constructs an empty fieldspec collection object.
  !> @return The constructed fieldspec collection object.
  !>
  function fieldspec_collection_constructor() result(self)

    implicit none

    type(fieldspec_collection_type) :: self

    self%fieldspec_list = linked_list_type()


  end function fieldspec_collection_constructor

  !===========================================================================
  !> @brief Destructor tears down object prior to being freed.
  !>
  subroutine fieldspec_collection_destructor( self )

    ! Object finalizer
    implicit none

    type (fieldspec_collection_type), intent(inout) :: self

    call self%clear()

    return
  end subroutine fieldspec_collection_destructor

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Creates a new fieldspec object and adds to the collection
  !>        Checks that the fieldspec unique id is not already in use
  !>
  subroutine add_fieldspec( self, unique_id, mesh_id, domain, &
                            order, field_kind, field_type )

    implicit none

    class(fieldspec_collection_type), intent(inout) :: self
    character(*),                     intent(in)    :: unique_id
    integer(i_def),                   intent(in)    :: mesh_id
    integer(i_def),                   intent(in)    :: domain
    integer(i_def),                   intent(in)    :: order
    integer(i_def),                   intent(in)    :: field_kind
    integer(i_def),                   intent(in)    :: field_type

    type(fieldspec_type) :: fieldspec_to_add

    ! Check this unique id is not already used and if it is, abort
    ! with error

    if ( associated( self%get_fieldspec( trim(unique_id) ) ) ) then

      write(log_scratch_space, '(3A)') &
            'The fieldspec unique id "', trim(unique_id), &
            '" is already in use'

      call log_event( log_scratch_space, LOG_LEVEL_ERROR)

    end if

    ! Create new fieldspec object

    fieldspec_to_add = fieldspec_type( unique_id,       &
                                       mesh_id,         &
                                       domain,          &
                                       order,           &
                                       field_kind,      &
                                       field_type)
    ! Add it to the list
    call self%fieldspec_list%insert_item( fieldspec_to_add )


    return
  end subroutine add_fieldspec

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Populates the fieldspec collection
  !> this will be enabled and populated on a subsequent task (#2014)
  !>
!  subroutine populate( self )
!
!    implicit none
!
!    class(fieldspec_collection_type), intent(inout) :: self
!
!
!    call log_event( 'Populating the fieldspec collection', LOG_LEVEL_INFO )
!
!    ! Hardcode all the prognostic fields we want to create
!
!    ! Eventually this will be read in from the iodef.xml
!
!    ! #ToDo: load from xml...
!
!    call self%add_fieldspec( "theta", 1_i_def, wtheta, 0_i_def, 0_i_def, 0_i_def )
!
!    return
!  end subroutine populate


  !===========================================================================
  !> @brief Requests a fieldspec object by unique id
  !>
  !> @param[in] unique_id of fieldspec object
  !>
  !> @return A pointer to the fieldspec object
  !>
  function get_fieldspec( self, unique_id ) result( fieldspec )

    implicit none

    class(fieldspec_collection_type)    :: self
    character(*),   intent(in)          :: unique_id

    type(fieldspec_type), pointer :: fieldspec

    ! Pointer to linked list - used for looping through the list
    type(linked_list_item_type),pointer :: loop => null()

    ! nullify this here as you can't do it on the declaration because it's also on the function def
    fieldspec => null()

    ! start at the head of the mesh collection linked list
    loop => self%fieldspec_list%get_head()

    do
      ! If list is empty or we're at the end of list and we didn't find the
      ! unique id, return a null pointer
      if ( .not. associated(loop) ) then
        nullify(fieldspec)
        exit
      end if

      ! 'cast' to the fieldspec_type and then we can get at the
      ! information in the payload
      select type( m => loop%payload )
        type is (fieldspec_type)
          fieldspec => m
      end select

      ! Now check for the id we want and keep looping until we
      ! find it or get to the end of the list
      if ( trim(unique_id) /= trim(fieldspec%get_unique_id()) ) then

        nullify(fieldspec)
        loop => loop%next

      else

        exit

      end if

    end do


    nullify(loop)

  end function get_fieldspec

  !===========================================================================
  !> @brief Forced clear of all the fieldspec objects in the collection.
  !>
  !> This routine should not need to be called manually except (possibly) in
  !> pfunit tests
  !>
  subroutine clear( self )

    ! Clear all items from the linked list in the collection
    implicit none

    class(fieldspec_collection_type), intent(inout) :: self

    call self%fieldspec_list%clear()

    return
  end subroutine clear

end module fieldspec_collection_mod
