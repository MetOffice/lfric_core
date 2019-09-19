!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!>
!> @brief Holds and manages halo_routing objects created during a model run.
!>
!> @details A container which a collection of halo_routing objects.
!>          The collection holds halo_routing objects as 
!>          singletons. It will handle the creation and storing of 
!>          requested halo_routing objects.
!
module halo_routing_collection_mod

  use constants_mod,      only: i_def, l_def
  use halo_routing_mod,   only: halo_routing_type
  use linked_list_mod,    only: linked_list_type, &
                                linked_list_item_type
  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Type which is is a collection of function spaces held in a linked list
  !-----------------------------------------------------------------------------
  type, public :: halo_routing_collection_type
    private
    !> Linked list which will hold the halo_routing objects
    type(linked_list_type) :: halo_routing_list
  contains
    !> Extracts a specific halo_routing object from the list
    procedure, public :: get_halo_routing
    !> Clears the link list of any content
    procedure, public :: clear
    ! Destructor, calls clear.
    final             :: halo_routing_collection_destructor
  end type halo_routing_collection_type
  !-----------------------------------------------------------------------------

  interface halo_routing_collection_type
    module procedure halo_routing_collection_constructor
  end interface

  ! Module level variable to make the function space collection
  ! globally available
  type(halo_routing_collection_type), public, allocatable :: &
      halo_routing_collection

contains

!> Function to construct a function space collection
!> @return self The new halo_routing_collection object
function halo_routing_collection_constructor() result(self)

  implicit none

  type(halo_routing_collection_type) :: self

  self%halo_routing_list = linked_list_type()

end function halo_routing_collection_constructor

!> Function to get an instance of a halo_routing object from the linked list
!> or create it if it doesn't exist (and return the newly created one)
!> @param [in] mesh_id Id of the mesh for which this information will
!>                     be valid for
!> @param [in] element_order The element order for which this information
!>                           will be valid
!> @param [in] lfric_fs The function space continuity type for which this
!>                      information will be valid
!> @return The halo_routing object that matches the input parameters
function get_halo_routing(self, mesh_id, element_order, lfric_fs) &
    result(halo_routing)

  implicit none

  class(halo_routing_collection_type), intent(inout) :: self

  type(halo_routing_type), pointer :: halo_routing

  integer(i_def), intent(in) :: mesh_id
  integer(i_def), intent(in) :: element_order
  integer(i_def), intent(in) :: lfric_fs

  halo_routing => get_halo_routing_from_list( self, &
                                              mesh_id, &
                                              element_order, &
                                              lfric_fs )

  if (.not. associated(halo_routing)) then

    call self%halo_routing_list%insert_item( halo_routing_type( mesh_id, &
                                                                element_order, &
                                                                lfric_fs ) )

    halo_routing => get_halo_routing_from_list( self, &
                                                mesh_id, &
                                                element_order, &
                                                lfric_fs )

  end if

  return
end function get_halo_routing

!------------------------------------------------------------------------------
! Note: This routine is not part of the API - but is used by get_halo_routing
! Function to scan the halo_routing collection for a
! halo_routing object with the given signature and return a pointer to it.
! A null pointer is returned if the requested halo_routing object does not exist.
!
function get_halo_routing_from_list(self, mesh_id, element_order, lfric_fs) &
    result(instance)

  implicit none

  class(halo_routing_collection_type), intent(inout) :: self
  integer(i_def), intent(in) :: mesh_id
  integer(i_def), intent(in) :: element_order
  integer(i_def), intent(in) :: lfric_fs

  type(halo_routing_type),   pointer  :: instance

  type(linked_list_item_type), pointer  :: loop

  ! Point to head of the function space linked list
  loop => self%halo_routing_list%get_head()

  ! Loop through the linked list
  do
    if ( .not. associated(loop) ) then
      ! Have reached the end of the list so either
      ! the list is empty or at the end of list.
      instance => null()

      loop => self%halo_routing_list%get_tail()
      exit
    end if

    ! 'cast' to the halo_routing_type 
    select type(listhalo_routing => loop%payload)
      type is (halo_routing_type)
      if ( mesh_id == listhalo_routing%get_mesh_id() .and. &
           element_order == listhalo_routing%get_element_order() .and. &
           lfric_fs == listhalo_routing%get_lfric_fs() ) then
        instance => listhalo_routing
        exit
      end if
    end select

    loop => loop%next
  end do

  nullify(loop)
  return
end function get_halo_routing_from_list

!> Function to clear all items from the field info collection linked list
subroutine clear(self)

  implicit none

  class(halo_routing_collection_type), intent(inout) :: self

  call self%halo_routing_list%clear()

  return
end subroutine clear

!> Destructor that will be called by the Fortran runtime environment when the 
!> field info collection object goes out of scope
subroutine halo_routing_collection_destructor(self)

  implicit none

  type (halo_routing_collection_type), intent(inout) :: self

  call self%clear()

  return
end subroutine halo_routing_collection_destructor

end module halo_routing_collection_mod
