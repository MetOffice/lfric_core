!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
!
!>
!> @brief   Holds and manages the multiple global meshes used to setup a model
!>          run.
!>
!> @details A container which holds a collection of global meshes
!>          It will handle the creation and storing of requested global meshes.
!
module global_mesh_collection_mod

  use constants_mod,   only: r_def, i_def, imdi, str_max_filename
  use log_mod,         only: log_event, log_scratch_space, &
                             LOG_LEVEL_ERROR, LOG_LEVEL_TRACE
  use linked_list_mod, only: linked_list_type, linked_list_item_type
  use global_mesh_mod, only: global_mesh_type

  implicit none

  private

  type, public :: global_mesh_collection_type

    private

    !> Linked list of the global_mesh_type objects in this collection
    type(linked_list_type) :: global_mesh_list
    !>
    !> @brief Number of panels in the mesh layout. npanels is set to
    !>        be the same as the 1st global mesh loaded into the collection.
    !>        All subsequent global meshes added should have been specified
    !>        with the same nume of panles
    !> @deprecated  Once multiple global meshes and associated mappings
    !>              are available in ugrid files.

    ! At present, global_mesh_type objects which are described from ugrid
    ! files which contain details for only one global mesh and thus contain no
    ! information about mapping between global meshes at different resolutions.
    ! As a consequence of this, mappings are calculated by the
    ! global_mesh_collection as between subsequent global meshes as they are
    ! added to the global mesh_collection. This calculation requires that all
    ! meshes have the same number of panels in the mesh. npanels is set to
    ! be the same as the 1st global mesh loaded into the collection.
    integer(i_def)         :: npanels = imdi
    !>
    !> An unused allocatable integer that prevents an intenal compiler error
    !> with the Gnu Fortran compiler. Adding an allocatable forces the compiler
    !> to accept that the object has a finaliser. It gets confused without it.
    !> This is a workaround for GCC bug id 61767 - when this bug is fixed, the
    !> integer can be removed.
    integer(i_def), allocatable :: dummy_for_gnu

  contains
    !> @brief Adds a global mesh object to the collection from a
    !>        ugrid file which contains a single global mesh.
    !> @param [in] filename  File containing details of a single global
    !>                       mesh object.
    !> @return ID of the global mesh added to collection.

    ! Adds a global mesh object to the collection from ugrid files
    ! which contain a single global mesh per file. Maps between
    ! each global mesh are created on the fly between subsequent
    ! global meshes in the collection in the order that they are
    ! added.
    !
    ! At present, there is no way to uniquely identify each global
    ! mesh that is read in. Global meshes are currently identified
    ! by an integer id counter which is assigned on creation of the
    ! global mesh object. When multiple global meshes per file are
    ! possible there will be a need to use a hash function to
    ! identify each global mesh object read in.
    procedure, public  :: add_new_global_mesh
    !>
    !> @brief Adds hardwired global mesh object to the collection
    !>        (for unit testing).
    !> @return ID of the global mesh added to collection
    procedure, public  :: add_unit_test_global_mesh
    !>
    !> @brief Requests global mesh object with specified global_mesh_id from
    !>        the collection.
    !> @param [in] global_mesh_id  Integer id of global mesh object requested.
    !> @return Global mesh object with requested global_mesh_id if present in
    !>         collection, a null pointer is returned if there is no global
    !>         mesh with the requested id.
    procedure, public  :: get_global_mesh
    !>
    !> @brief Forced clear of all the global mesh objects in the collection.
    !>        This routine should not need to be called manually except
    !>        (possibly) in pfunit tests
    procedure, public  :: clear
    !>
    !> @brief Internal routine to generate cell mappings (if possible) between
    !>        subsequent global meshes added to the collection. This
    !>        routine is present as mapping data is not currently available
    !>        in the global mesh ugrid files. Mappings are created on a
    !>        "per mesh panel" basis and is checked against the variable
    !>        npanels.
    !> @deprecated When ugrid files with multiple meshes/mappings are
    !>             available
    procedure, private :: map_global_meshes
    !>
    !> @brief Finalizer routine, should be called automatically by
    !>        code when the object is out of scope
    final              :: global_mesh_collection_destructor

  end type global_mesh_collection_type

  interface global_mesh_collection_type
    module procedure global_mesh_collection_constructor
  end interface

  ! Module variable allows access to the single mesh collection
  type(global_mesh_collection_type), public, allocatable :: &
      global_mesh_collection

contains

!> @brief Constructs the mesh collection object
!> @return The constructed mesh collection object
function global_mesh_collection_constructor() result(self)

  implicit none

  type(global_mesh_collection_type) :: self

  self%global_mesh_list = linked_list_type()

end function global_mesh_collection_constructor


!==============================================================================
function add_new_global_mesh( self, filename, npanels ) result (global_mesh_id)

  implicit none

  class(global_mesh_collection_type), intent(inout) :: self
  character(len=str_max_filename), intent(in)       :: filename
  integer(i_def),                  intent(in)       :: npanels

  integer(i_def)          :: n_global_meshes
  integer(i_def)          :: global_mesh_id
  type (global_mesh_type) :: global_mesh_to_add
  type (global_mesh_type), pointer :: global_mesh_at_tail => null()

  integer(i_def) :: source_global_mesh_id
  integer(i_def) :: target_global_mesh_id

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type), pointer :: list_item => null()

  global_mesh_to_add = global_mesh_type(trim(filename))
  global_mesh_id     = global_mesh_to_add%get_id()

  n_global_meshes = self%global_mesh_list%get_length()

  ! If there is at least one other mesh in the collection, attempt to
  ! create a map from it to the global mesh being added
  if (n_global_meshes >= 1) then

    ! In order to map the global meshes the two global
    ! meshes should have the same number of panels
    if (self%npanels == npanels) then

      list_item => self%global_mesh_list%get_tail()
      select type(m => list_item%payload)
      type is (global_mesh_type)
        global_mesh_at_tail => m
      end select

      source_global_mesh_id = global_mesh_at_tail%get_id()
      target_global_mesh_id = global_mesh_id

      call self%global_mesh_list%insert_item( global_mesh_to_add )
      call self%map_global_meshes( source_global_mesh_id, &
                                   target_global_mesh_id )

    else
       write(log_scratch_space,'(A,I0)')                                     &
           "This global mesh collection is for global meshes of npanels = ", &
           self%npanels
       call log_event(log_scratch_space,LOG_LEVEL_ERROR)
       return
     end if

  else

    self%npanels = npanels
    call self%global_mesh_list%insert_item( global_mesh_to_add )

  end if

  return
end function add_new_global_mesh


!==============================================================================
subroutine map_global_meshes( self,                   &
                              source_global_mesh_id,  &
                              target_global_mesh_id )

  implicit none

  class(global_mesh_collection_type), intent(inout) :: self
  integer(i_def), intent(in) :: source_global_mesh_id
  integer(i_def), intent(in) :: target_global_mesh_id

  type (global_mesh_type), pointer :: source_mesh => null()
  type (global_mesh_type), pointer :: target_mesh => null()
  type (global_mesh_type), pointer :: coarse_mesh => null()
  type (global_mesh_type), pointer :: fine_mesh   => null()

  integer(i_def) :: i,j,n,count       ! counters
  integer(i_def) :: coarse_panel_start_id
  integer(i_def) :: fine_panel_start_id
  integer(i_def) :: coarse_panel_end_id
  integer(i_def) :: fine_panel_end_id
  integer(i_def) :: coarse_id
  integer(i_def) :: fine_id

  integer(i_def) :: coarse_ndiv
  integer(i_def) :: fine_ndiv
  integer(i_def) :: coarse_ncells
  integer(i_def) :: fine_ncells
  integer(i_def) :: source_ncells
  integer(i_def) :: target_ncells

  integer(i_def) :: factor ! ratio of ndivs

  integer(i_def) :: start_x
  integer(i_def) :: start_y
  integer(i_def) :: end_x
  integer(i_def) :: end_y

  integer(i_def), allocatable :: coarse_panel_ids(:,:)
  integer(i_def), allocatable :: fine_panel_ids(:,:)

  integer(i_def), allocatable :: coarse_to_fine_gid_map(:,:)
  integer(i_def), allocatable :: fine_to_coarse_gid_map(:,:)
  integer(i_def), allocatable :: tmp_panel_ids(:)


  ! Check for identical source and target global meshes
  if (source_global_mesh_id == target_global_mesh_id) then
    write(log_scratch_space,'(A)')                                           &
         "Cannot have identical consective global mesh ids"
    call log_event(log_scratch_space,LOG_LEVEL_TRACE)
    return
  end if

  source_mesh => self%get_global_mesh(source_global_mesh_id)
  target_mesh => self%get_global_mesh(target_global_mesh_id)

  ! Figure out which of the meshes has fewer
  ! cells i.e. the coarse one
  source_ncells = source_mesh%get_ncells()
  target_ncells = target_mesh%get_ncells()

  if ( mod(source_ncells, target_ncells) == 0 .or. &
       mod(target_ncells, source_ncells) == 0 ) then

    if (source_ncells == target_ncells) then
      write(log_scratch_space, '(A)')                                          &
          'Meshes have equal number of cells, direct mapping, '//              &
          'no mapping created'
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      return

    else if (target_ncells > source_ncells) then
      coarse_mesh   => source_mesh
      fine_mesh     => target_mesh
      coarse_ncells = source_ncells
      fine_ncells   = target_ncells

    else
      coarse_mesh   => target_mesh
      fine_mesh     => source_mesh
      coarse_ncells = target_ncells
      fine_ncells   = source_ncells

    end if

  else
    ! The number of cells in one of these meshes
    ! needs to be a factor of the number cells in the other.
    write(log_scratch_space, '(2(A,I0))')                                      &
        'Unable to generate intermesh connectivity for global meshes: id:',    &
        source_global_mesh_id, ';id:', target_global_mesh_id
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    return

  end if

  coarse_id = coarse_mesh%get_id()
  fine_id   = fine_mesh%get_id()

  ! Now we know which of the source and target are the
  ! coarser and finer meshes, generate a coarse_to_fine
  ! and fine_to_coarse map per panel
  coarse_ndiv = INT(SQRT( real(coarse_ncells/self%npanels, r_def) ), i_def)
  fine_ndiv   = INT(SQRT( real(fine_ncells/self%npanels,   r_def) ), i_def)
  factor      = fine_ndiv/coarse_ndiv

  allocate( coarse_panel_ids ( coarse_ndiv, coarse_ndiv ))
  allocate( fine_panel_ids   ( fine_ndiv,   fine_ndiv   ))
  allocate( coarse_to_fine_gid_map( factor**2, coarse_ncells ) )
  allocate( fine_to_coarse_gid_map( 1, fine_ncells ))

  do n=1, self%npanels

    ! Get the id of the initial cell for the panel
    coarse_panel_start_id = ((n-1)*coarse_ncells / self%npanels) + 1
    fine_panel_start_id   = ((n-1)*fine_ncells   / self%npanels) + 1

    ! =======================================================
    ! Populate arrays with ids and reshape to panel layout
    ! =======================================================

    ! For coarse mesh
    allocate(tmp_panel_ids(coarse_ncells))
    count=1
    coarse_panel_end_id = coarse_panel_start_id + coarse_ncells-1
    do i=coarse_panel_start_id, coarse_panel_end_id
      tmp_panel_ids(count) = i
      count=count+1
    end do

    coarse_panel_ids = reshape(tmp_panel_ids,(/coarse_ndiv,coarse_ndiv/))
    deallocate(tmp_panel_ids)

    ! For fine mesh
    allocate(tmp_panel_ids(fine_ncells))
    count=1
    fine_panel_end_id =  fine_panel_start_id + fine_ncells-1
    do i=fine_panel_start_id, fine_panel_end_id
      tmp_panel_ids(count) = i
      count=count+1
    end do
    fine_panel_ids = reshape(tmp_panel_ids,(/fine_ndiv,fine_ndiv/))
    deallocate(tmp_panel_ids)

    ! ===============================
    ! Populate intermesh gid maps
    ! ===============================

    ! Coarse to Fine
    do j=1, coarse_ndiv
      start_y = ((j-1)*factor) + 1
      end_y   = start_y + factor - 1
      do i=1, coarse_ndiv
        start_x = ((i-1)*factor) + 1
        end_x   = start_x + factor - 1

        coarse_to_fine_gid_map(:,coarse_panel_ids(i,j)) =                      &
            reshape( fine_panel_ids( start_x:end_x, start_y:end_y ),           &
                     (/factor**2/) )
      end do
    end do

    ! Fine to Coarse
    coarse_panel_end_id = coarse_panel_start_id                                &
                        + (coarse_ncells/self%npanels) - 1

    do j=coarse_panel_start_id, coarse_panel_end_id
      do i=1, factor**2
        fine_to_coarse_gid_map(1, coarse_to_fine_gid_map(i,j)) = j
      end do
    end do

  end do

  ! Now we have the intermesh gid maps from coarse to fine
  ! and fine to coarse meshes. Create a global_mesh_map for each
  ! direction.
  call coarse_mesh % add_global_mesh_map( fine_mesh, &
                                          coarse_to_fine_gid_map )
  call fine_mesh   % add_global_mesh_map( coarse_mesh, &
                                          fine_to_coarse_gid_map )

  deallocate( coarse_panel_ids )
  deallocate( fine_panel_ids   )
  deallocate( coarse_to_fine_gid_map )
  deallocate( fine_to_coarse_gid_map )

  return
end subroutine map_global_meshes


!==============================================================================
function add_unit_test_global_mesh(self) result(global_mesh_id)

  implicit none

  class(global_mesh_collection_type), intent(inout) :: self
  integer(i_def) :: global_mesh_id

  type (global_mesh_type) :: global_mesh

  integer(i_def) :: n_global_meshes

  n_global_meshes = self%global_mesh_list%get_length()

  if (n_global_meshes == 0) then
    self%npanels   = 1
    global_mesh    = global_mesh_type()
    global_mesh_id = global_mesh%get_id()
    call self%global_mesh_list%insert_item( global_mesh )
  else
    write(log_scratch_space,'(A)')                                             &
         "Test global mesh must only be added to an "                        //&
          "empty global mesh collection."
    call log_event(log_scratch_space,LOG_LEVEL_ERROR)
    return
  end if

  return
end function add_unit_test_global_mesh


!==============================================================================
function get_global_mesh( self, global_mesh_id ) result( global_mesh )

  implicit none

  class(global_mesh_collection_type) :: self
  integer(i_def), intent(in) :: global_mesh_id

  type(global_mesh_type), pointer :: global_mesh

  ! Pointer to linked list - used for looping through the list
  type(linked_list_item_type),pointer :: loop => null()

  ! start at the head of the mesh collection linked list
  loop => self%global_mesh_list%get_head()

  do
    ! If list is empty or we're at the end of list and we didn't find the
    ! mesh_id, return a null pointer
    if ( .not. associated(loop) ) then
      nullify(global_mesh)
      exit
    end if

    ! Otherwise search list for the id we want
    if ( global_mesh_id == loop%payload%get_id() ) then
      ! 'cast' to the global_mesh_type
      select type(m => loop%payload)
        type is (global_mesh_type)
          global_mesh => m
      end select
      exit
    end if
    loop => loop%next
  end do

end function get_global_mesh

!==============================================================================
subroutine clear(self)

  ! Clear all items from the linked list in the collection
  implicit none

  class(global_mesh_collection_type), intent(inout) :: self

  call self%global_mesh_list%clear()

  return
end subroutine clear


!==============================================================================
subroutine global_mesh_collection_destructor(self)

  ! Object finalizer
  implicit none

  type (global_mesh_collection_type), intent(inout) :: self

  return
end subroutine global_mesh_collection_destructor



end module global_mesh_collection_mod
