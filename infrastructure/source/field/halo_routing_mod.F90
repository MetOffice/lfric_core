!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Holds and manages information about how to route halo exchange
!>        communicatoins of a particular type of field
!>
!> @details Holds information that can be attached to a number of different
!>          fields that describes how to perform halo exchanges. This
!>          information is in the form of a YAXT redistribution map.
!
module halo_routing_mod

  use constants_mod,         only: i_def, i_native, i_halo_index, &
                                   r_def, real_type
  use function_space_mod,    only: function_space_type
  use function_space_collection_mod, &
                             only: function_space_collection_type, &
                                   function_space_collection
  use linked_list_data_mod,  only: linked_list_data_type
  use mesh_mod,              only: mesh_type
  use mpi_mod,               only: generate_redistribution_map, get_mpi_datatype
  use yaxt,                  only: xt_redist

  implicit none

  private

  type, extends(linked_list_data_type), public :: halo_routing_type
    private
    !> Id of the mesh used in the function space that this information
    !> is valid for
    integer(i_def) :: mesh_id
    !> Order of the function space that this information is valid for
    integer(i_def) :: element_order
    !> Enumerated value representing the continutity of the function space
    !> that this information is valid for
    integer(i_native) :: lfric_fs   
    !> YAXT redistribution map 
    type(xt_redist), allocatable :: redist(:)
  contains
    !> Gets the mesh_id for which the halo_routing object is valid
    procedure, public :: get_mesh_id
    !> Gets the element_order for which the halo_routing object is valid
    procedure, public :: get_element_order
    !> Gets the function space continuity type for which the halo_routing
    !> object is valid
    procedure, public :: get_lfric_fs
    !> Gets a YAXT redistribution map for halo swapping
    procedure, public :: get_redist
  end type halo_routing_type
  !-----------------------------------------------------------------------------

  interface halo_routing_type
    module procedure halo_routing_constructor
  end interface

contains

!> Constructor for a halo_routing object
!> @param [in] mesh_id Id of the mesh for which this information will
!>                     be valid
!> @param [in] element_order The element order for which this information
!>                           will be valid
!> @param [in] lfric_fs The function space continuity type for which this
!>                      information will be valid
!> @return The new halo_routing object
function halo_routing_constructor( mesh_id, element_order, lfric_fs )  result(self)

  implicit none

  integer(i_def), intent(in) :: mesh_id
  integer(i_def), intent(in) :: element_order
  integer(i_def), intent(in) :: lfric_fs

  type(halo_routing_type) :: self

  type(function_space_type), pointer :: function_space => null()
  type(mesh_type), pointer :: mesh => null()
  integer(i_halo_index), allocatable :: global_dof_id(:)
  integer(i_def) :: idepth
  integer(i_def) :: halo_start, halo_finish

  self%mesh_id = mesh_id
  self%element_order = element_order
  self%lfric_fs = lfric_fs

  function_space => function_space_collection%get_fs( mesh_id, &
                                                      element_order, &
                                                      lfric_fs )

  ! set up the global dof index array
  allocate( global_dof_id( function_space%get_ndof_glob() ) )
  call function_space%get_global_dof_id(global_dof_id)

  mesh => function_space%get_mesh()

  ! Set up YAXT redistribution map for halo exchanges
  allocate( self%redist(mesh%get_halo_depth()) )

  do idepth = 1, mesh%get_halo_depth()

    halo_start  = function_space%get_last_dof_owned()+1
    halo_finish = function_space%get_last_dof_halo(idepth)
    ! If this is a serial run (no halos), halo_start will be out of bounds,
    ! so re-initialise halo_start and halo_finish specifically for a serial run
    if(halo_start > function_space%get_last_dof_halo(idepth))then
      halo_start  = function_space%get_last_dof_halo(idepth)
      halo_finish = function_space%get_last_dof_halo(idepth)-1
    end if

    ! Get the redistribution map objects for doing halo exchanges later
    self%redist(idepth) = generate_redistribution_map( &
                     global_dof_id(1:function_space%get_last_dof_owned()), &
                     global_dof_id( halo_start:halo_finish ), &
                     get_mpi_datatype( real_type, r_def ) )
  end do
end function halo_routing_constructor

!> Gets the mesh_id for which this object is valid
!> @return Id of the mesh that this information is valid for
function get_mesh_id(self) result (mesh_id)
  implicit none
  class(halo_routing_type), intent(in) :: self
  integer(i_def) :: mesh_id
  mesh_id = self%mesh_id
  return
end function get_mesh_id

!> Gets the element_order for which this object is valid
!> @return The element order that this information is valid for
function get_element_order(self) result (element_order)
  implicit none
  class(halo_routing_type), intent(in) :: self
  integer(i_def) :: element_order
  element_order = self%element_order
  return
end function get_element_order

!> Gets the function space continuity type for which this object is valid
!> @return The function space continuity type that this information is valid for
function get_lfric_fs(self) result (lfric_fs)
  implicit none
  class(halo_routing_type), intent(in) :: self
  integer(i_def) :: lfric_fs
  lfric_fs = self%lfric_fs
  return
end function get_lfric_fs

!> Gets a YAXT redistribution map for halo swapping
!> @param [in] depth The depth of halo exchange that the redistribution map
!>                   will be used for
!> @return The YAXT redistribution map for a halo exchange of a particular
!>         depth of halo, valid for a particualar field type
function get_redist(self, depth) result (redist)

  implicit none
  class(halo_routing_type) :: self
  integer(i_def) , intent(in) :: depth

  type(xt_redist) :: redist

  redist = self%redist(depth)

  return
end function get_redist

end module halo_routing_mod
