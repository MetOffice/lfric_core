module exchange_map_collection_mod

  use constants_mod,          only: i_def, i_halo_index, l_def
  use function_space_mod,     only: function_space_type
  use function_space_collection_mod, &
                              only: function_space_collection_type, &
                                    function_space_collection
  use halo_comms_mod,         only: exchange_map_type
  use mesh_mod,               only: mesh_type
  use log_mod,                only: log_event, LOG_LEVEL_ERROR
  use lfric_mpi_mod,          only: global_mpi, get_lfric_datatype
  use linked_list_mod,        only: linked_list_type, &
                                    linked_list_item_type
  use linked_list_data_mod,   only: linked_list_data_type

  implicit none

  private

  !collection type
  type, public :: exchange_map_collection_type
    private
    type(linked_list_type) :: exchange_map_list
  contains
    procedure, public :: get_exchange_map
  end type exchange_map_collection_type
  interface exchange_map_collection_type
    module procedure exchange_map_collection_constructor
  end interface
contains

!> \brief Constructor for an exchange_map_collection object.
!>
!> Initializes a new exchange_map_collection_type instance and
!> allocates an empty linked list inside the object.
!>
!> \return A fully initialised exchange_map_collection_type object.
function exchange_map_collection_constructor() result(self)

  implicit none
  !> Constructed exchange_map_collection instance.
  type(exchange_map_collection_type) :: self

  self%exchange_map_list = linked_list_type()

end function exchange_map_collection_constructor

!> \brief Retrieve or construct an exchange map for the given mesh
!!        and finite-element parameters.
!>
!! This routine searches the collection for an existing
!! exchange_map_type object that matches the supplied parameters.
!! If no matching exchange map is found, a new one is constructed,
!! initialised, inserted into the collection, and then returned.
!>
!> \param[in]     mesh            The mesh object associated with the exchange map.
!> \param[in]     element_order_h Horizontal element order.
!> \param[in]     element_order_v Vertical element order.
!> \param[in]     lfric_fs        Identifier for the LFRic function space.
!> \param[in]     ndata           Number of data items per degree of freedom.
!> \param[in]     halo_depth      Number of halo layers to include.
!>
!> \return A pointer to the corresponding exchange_map_type object.
function get_exchange_map( self, &
                   mesh, &
                   element_order_h, &
                   element_order_v, &
                   lfric_fs, &
                   ndata, &
                   halo_depth )  result(exchange_maps)
  implicit none

  class(exchange_map_collection_type), intent(inout) :: self

  type(exchange_map_type), pointer :: exchange_maps

  type(mesh_type), intent(in), pointer :: mesh

  integer(i_def), intent(in) :: ndata
  integer(i_def), intent(in) :: halo_depth
  integer(i_def), intent(in) :: element_order_v, element_order_h
  integer(i_def), intent(in) :: lfric_fs

  type(function_space_type), pointer :: function_space
  integer(i_halo_index), allocatable :: global_dof_id(:)
  integer(i_def), allocatable :: halo_start(:)
  integer(i_def), allocatable :: halo_finish(:)
  integer(i_def) :: idepth
  integer(i_def) :: last_owned_dof
  integer(i_def) :: mesh_id

  nullify( function_space )

  exchange_maps => get_exchange_maps_from_list( self, &
                                              mesh, &
                                              element_order_h, &
                                              element_order_v, &
                                              lfric_fs, &
                                              ndata, &
                                              halo_depth )

  if (.not. associated(exchange_maps)) then

    !Get indices of owned and halo cells
    function_space => function_space_collection%get_fs( mesh, element_order_h, &
                                                        element_order_v, &
                                                        lfric_fs, &
                                                        ndata)

    last_owned_dof = function_space%get_last_dof_owned()

    ! Set up the global dof index array
    call function_space%get_global_dof_id(global_dof_id)

    ! Set up the boundaries of the different depths of halo
    allocate( halo_start(halo_depth) )
    allocate( halo_finish(halo_depth) )

    do idepth = 1, halo_depth

      halo_start(idepth)  = function_space%get_last_dof_owned()+1
      halo_finish(idepth) = function_space%get_last_dof_halo(idepth)
      ! The above assumes there is a halo cell following the last owned cell.
      ! This might not be true (e.g. in a serial run), so fix the start/finish
      ! points when that happens
      if ( halo_start(idepth) > function_space%get_last_dof_halo(idepth) ) then
        halo_start(idepth)  = function_space%get_last_dof_halo(idepth)
        halo_finish(idepth) = halo_start(idepth) - 1
      end if

    end do

    mesh_id = mesh%get_id()
    call self%exchange_map_list%insert_item( exchange_map_type( global_dof_id,&
                                                                last_owned_dof,&
                                                                halo_start, &
                                                                halo_finish, &
                                                                mesh_id, &
                                                                element_order_h, &
                                                                element_order_v, &
                                                                lfric_fs, &
                                                                ndata, &
                                                                halo_depth ))
    deallocate( halo_start, halo_finish, global_dof_id )

    exchange_maps => get_exchange_maps_from_list( self, &
                                                mesh, &
                                                element_order_h, &
                                                element_order_v, &
                                                lfric_fs, &
                                                ndata, &
                                                halo_depth )

  end if

  return
end function get_exchange_map

!> \brief Search the exchange-map collection for a matching entry.
!!
!! This routine walks through the linked list stored inside the
!! exchange-map collection and checks whether an existing
!! exchange map matches the supplied data.
!! If a matching map is found, a pointer to it is returned.
!! If no match is found, the function returns a null pointer.
!>
!> \param[in,out] self            The exchange-map collection to search.
!> \param[in]     mesh            The mesh used to identify the map.
!> \param[in]     element_order_h Horizontal element order.
!> \param[in]     element_order_v Vertical element order.
!> \param[in]     lfric_fs        Identifier for the LFRic function space.
!> \param[in]     ndata           Number of data items per degree of freedom.
!> \param[in]     halo_depth      Depth of the halo region.
!>
!> \return Pointer to an existing exchange map, or a null pointer if none match.
function get_exchange_maps_from_list(self, &
                                    mesh, &
                                    element_order_h, &
                                    element_order_v, &
                                    lfric_fs, &
                                    ndata, &
                                    halo_depth) &
    result(instance)

  implicit none

  class(exchange_map_collection_type), intent(inout) :: self

  type(mesh_type), intent(in), pointer :: mesh

  integer(i_def),  intent(in) :: ndata
  integer(i_def),  intent(in) :: halo_depth
  integer(i_def),  intent(in) :: element_order_h, element_order_v
  integer(i_def),  intent(in) :: lfric_fs

  type(exchange_map_type),   pointer  :: instance

  type(linked_list_item_type), pointer  :: loop

  integer(i_def) :: mesh_id

  mesh_id = mesh%get_id()
  ! Point to head of the exchange map linked list
  loop => self%exchange_map_list%get_head()

  ! Loop through the linked list
  do
    if ( .not. associated(loop) ) then
      ! Have reached the end of the list so either
      ! the list is empty or at the end of list.
      instance => null()

      loop => self%exchange_map_list%get_tail()
      exit
    end if

    ! 'cast' to the halo_routing_type
    select type(listhalo_routing => loop%payload)
      type is (exchange_map_type)
      if ( mesh_id == listhalo_routing%get_exchange_map_mesh_id() .and. &
           element_order_h == listhalo_routing%get_exchange_map_element_order_h() .and. &
           element_order_v == listhalo_routing%get_exchange_map_element_order_v() .and. &
           lfric_fs == listhalo_routing%get_exchange_map_lfric_fs() .and. &
           ndata == listhalo_routing%get_exchange_map_ndata() .and. &
           halo_depth == listhalo_routing%get_exchange_map_halo_depth() ) then
        instance => listhalo_routing
        exit
      end if
    end select

    loop => loop%next
  end do

  nullify(loop)
  return
end function get_exchange_maps_from_list

end module exchange_map_collection_mod
