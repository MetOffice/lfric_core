!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module for the fieldspec class
!>
!> @details A fieldspec contains all relevant information to create a field


module fieldspec_mod

  use constants_mod,        only: i_def, str_def
  use linked_list_data_mod, only: linked_list_data_type

  implicit none

  private

  type, extends(linked_list_data_type), public :: fieldspec_type

    private

    !> Unique id used by the diagnostic system to identify the field
    character(str_def) :: unique_id
    !> Other information used to create a field
    integer(i_def)     :: mesh_id
    integer(i_def)     :: domain
    integer(i_def)     :: order
    integer(i_def)     :: field_kind
    integer(i_def)     :: field_type

  contains

    !> Getter to return the unique_id
    procedure, public :: get_unique_id

    !> Getter to return the mesh_id
    procedure, public :: get_mesh_id

    !> Getter to return the domain
    procedure, public :: get_domain

    !> Getter to return the order
    procedure, public :: get_order

    !> Getter to return the kind
    procedure, public :: get_kind

    !> Getter to return the type
    procedure, public :: get_type


  end type fieldspec_type

  interface fieldspec_type
    module procedure fieldspec_constructor
  end interface


contains


  !> Construct a <code>fieldspec_type</code> object.
  !>
  !> @param [in] unique_id A unique identifer for the field
  !> @param [in] mesh_id The mesh to create field with
  !> @param [in] domain The domain (function space) to create field with
  !> @param [in] order The order of the field
  !> @param [in] field_kind The kind of the field
  !> @param [in] field_type The type of the field
  !> @return self the fieldspec object
  !>
  function fieldspec_constructor( unique_id, mesh_id,     &
                                  domain, order,          &
                                  field_kind, field_type ) &
                                 result(self)

    use log_mod,         only : log_event, &
                                LOG_LEVEL_ERROR
    implicit none

    character(*),               intent(in)    :: unique_id
    integer(i_def),             intent(in)    :: mesh_id
    integer(i_def),             intent(in)    :: domain
    integer(i_def),             intent(in)    :: order
    integer(i_def),             intent(in)    :: field_kind
    integer(i_def),             intent(in)    :: field_type

    type(fieldspec_type), target :: self

    self%unique_id  = trim(unique_id)
    self%mesh_id    = mesh_id
    self%domain     = domain
    self%order      = order
    self%field_kind = field_kind
    self%field_type = field_type

  end function fieldspec_constructor


  !> Getter for unique_id
  !> @param[in]  self  fieldspec_type
  !> @return unique_id

  function get_unique_id(self) result(unique_id)

    implicit none

    class(fieldspec_type), intent(in) :: self
    character(str_def) :: unique_id

    unique_id = trim(self%unique_id)

  end function get_unique_id

  !> Getter for mesh_id
  !> @param[in]  self  fieldspec_type
  !> @return mesh_id
  function get_mesh_id(self) result(mesh_id)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: mesh_id

    mesh_id = self%mesh_id

  end function get_mesh_id

  !> Getter for domain
  !> @param[in]  self  fieldspec_type
  !> @return domain
  function get_domain(self) result(domain)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: domain

    domain = self%domain

  end function get_domain

  !> Getter for order
  !> @param[in]  self  fieldspec_type
  !> @return order
  function get_order(self) result(order)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: order

    order = self%order

  end function get_order

  !> Getter for field_kind
  !> @param[in]  self  fieldspec_type
  !> @return field_kind
  function get_kind(self) result(field_kind)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: field_kind

    field_kind = self%field_kind

  end function get_kind

  !> Getter for field_type
  !> @param[in]  self  fieldspec_type
  !> @return field_type
  function get_type(self) result(field_type)

    implicit none

    class(fieldspec_type), intent(in) :: self
    integer(i_def) :: field_type

    field_type = self%field_type

  end function get_type


end module fieldspec_mod
