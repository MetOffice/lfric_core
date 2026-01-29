!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>  @brief Module containing an XIOS interface for fields
!>
module lfric_xios_field_mod

  use constants_mod,        only: str_def
  use field_collection_iterator_mod, &
                            only: field_collection_iterator_type
  use field_mod,            only: field_type
  use field_real32_mod,     only: field_real32_type
  use field_real64_mod,     only: field_real64_type
  use field_parent_mod,     only: field_parent_type
  use function_space_mod,   only: function_space_type
  use fs_continuity_mod,    only: W0, W1, W2, W3, Wtheta, W2H, W2V, &
                                  name_from_functionspace
  use integer_field_mod,    only: integer_field_type
  use log_mod,              only: log_event, log_level_trace, log_level_error
  use mesh_mod,             only: mesh_type
  use xios,                 only: xios_domain, xios_axis, xios_grid,          &
                                  xios_gridgroup, xios_is_valid_grid,         &
                                  xios_field,xios_fieldgroup, xios_add_child, &
                                  xios_get_handle, xios_set_attr, xios_get_attr

  implicit none

private

!> @brief Container for field properties need by XIOS
!>
type, public :: lfric_xios_field_type
  private

  !> XIOS handle for field object
  type(xios_field)                  :: handle
  !> Unique identifier for XIOS field handle
  character(str_def)                :: xios_id
   !> Name of the XIOS field as it appears in the output file
  character(str_def)                :: xios_name
  !> Unique identifier for XIOS fieldgroup handle
  character(str_def)                :: fieldgroup_id
  !> Pointer to model field associated with XIOS field
  class(field_parent_type), pointer :: model_field

contains

  procedure, public :: register
  procedure, public :: recv
  procedure, public :: send
  procedure, public :: get_model_field
  procedure, public :: set_model_field
  procedure, public :: get_xios_name
  procedure, public :: set_xios_name
  final             :: lfric_xios_field_final

end type lfric_xios_field_type

interface lfric_xios_field_type
  module procedure lfric_xios_field_constructor
end interface

contains

!> Constructor for the lfric_xios_field_type object.
!>
!> @param[in] model_field   A pointer to the model field associated with this
!>                          XIOS field
!> @param[in] xios_id       Optional: The XIOS ID for this field
!> @param[in] fieldgroup_id Optional: The XIOS ID for the field group that this
!>                          field belongs to
function lfric_xios_field_constructor(model_field, xios_id, fieldgroup_id) result(self)

  implicit none

  type(lfric_xios_field_type) :: self

  class(field_parent_type), pointer, intent(in) :: model_field
  character(len=*), optional,        intent(in) :: xios_id
  character(len=*), optional,        intent(in) :: fieldgroup_id

  ! Point object to model data field
  self%model_field => model_field

  ! Create unique XIOS ID for XIOS field
  if (present(xios_id)) then
    if (present(fieldgroup_id)) then
      self%xios_id = trim(adjustl(xios_id))//"__"//trim(adjustl(fieldgroup_id))
    else
      self%xios_id = trim(adjustl(xios_id))
    end if
  else
    if (present(fieldgroup_id)) then
      self%xios_id = trim(adjustl(model_field%get_name()))//"__"//trim(adjustl(fieldgroup_id))
    else
      self%xios_id = trim(adjustl(model_field%get_name()))
    end if
  end if

  if (present(fieldgroup_id)) then
    self%fieldgroup_id = trim(adjustl(fieldgroup_id))
  else
    self%fieldgroup_id = "field_definition"
  end if

  return

end function lfric_xios_field_constructor

!> Registers a representation of the model field with the associated XIOS field
!> group
subroutine register(self, field_read_access)

  implicit none

  class(lfric_xios_field_type), intent(inout) :: self
  logical,                      intent(in)    :: field_read_access

  type(xios_fieldgroup) :: fieldgroup_hdl
  type(xios_domain)     :: domain
  type(xios_axis)       :: vert_axis, ndata_axis
  type(xios_gridgroup)  :: grid_definition
  type(xios_grid)       :: new_grid
  character(str_def)    :: domain_id, axis_id, grid_id

  type(mesh_type),           pointer :: mesh => null()
  type(function_space_type), pointer :: vspace => null()

  call log_event( "Registering XIOS field ["//trim(self%xios_id)//      &
                  "] with field group ["//trim(self%fieldgroup_id)//"]", &
                  log_level_trace )

  ! Get field group handle and add field
  call xios_get_handle(trim(adjustl(self%fieldgroup_id)), fieldgroup_hdl)
  call xios_add_child(fieldgroup_hdl, self%handle, trim(self%xios_id))
  call xios_set_attr( self%handle,                                     &
                      name=trim(adjustl(self%model_field%get_name())), &
                      read_access=field_read_access )

  ! Set up dimensions of output field
  vspace => self%model_field%get_function_space()
  select case(self%model_field%which_function_space())
  case (W0)
    domain_id="node"
    axis_id  ="vert_axis_full_levels"
  case (W2H)
    domain_id="edge"
    axis_id  ="vert_axis_half_levels"
  case (W2V)
    domain_id="face"
    axis_id  ="vert_axis_full_levels"
  case (Wtheta)
    domain_id="face"
    axis_id  ="vert_axis_full_levels"
  case (W3)
    domain_id="face"
    axis_id  ="vert_axis_half_levels"
  end select

  mesh => self%model_field%get_mesh()
  ! If field had only a single layer of DoFs then define this using a domain,
  ! otherwise use a grid.
  if ((mesh%get_nlayers() == 1) .and. (vspace%get_ndata() == 1)) then
    call xios_set_attr(self%handle, domain_ref=trim(domain_id))

  else
    call xios_get_handle("grid_definition", grid_definition)

    ! Create grid ID using function space information.
    if (vspace%get_ndata() == 1) then
      grid_id = trim(name_from_functionspace(self%model_field%which_function_space()))//"_grid"
    else if (mesh%get_nlayers() == 1) then
      grid_id = trim(domain_id)//"_"//char(vspace%get_ndata())//"_grid"
    else
      grid_id = trim(name_from_functionspace(self%model_field%which_function_space()))// &
                "_ndata_"//char(vspace%get_ndata())//"_grid"
    end if

    ! If this grid does not already exist, create it.
    if (.not. xios_is_valid_grid(trim(grid_id))) then
      call xios_add_child(grid_definition, new_grid, trim(grid_id))

      call xios_add_child(new_grid, domain)
      call xios_set_attr(domain, domain_ref=trim(domain_id))

      if (mesh%get_nlayers() > 1) then
        call xios_add_child(new_grid, vert_axis)
        call xios_set_attr(vert_axis, axis_ref=trim(axis_id))
      end if

      if (vspace%get_ndata() > 1) then
        call xios_add_child(new_grid, ndata_axis)
        call xios_set_attr(ndata_axis, n_glo=vspace%get_ndata())
      end if

    end if

    call xios_set_attr( self%handle, grid_ref=trim(grid_id) )

  end if

  ! Set output field precision
  select type(fld => self%model_field)
  type is (field_real32_type)
    call xios_set_attr(self%handle, prec=4)

  type is (field_real64_type)
    call xios_set_attr(self%handle, prec=8)

  type is (integer_field_type)
    call xios_set_attr(self%handle, prec=2)

  class default
    call log_event( "Field ["//trim(self%model_field%get_name())//"] does " // &
                    "not have a valid type for XIOS I/O", log_level_error )

  end select

end subroutine register

!> Recieves data from the associated XIOS field and sends it to the
!> corresponding model field
subroutine recv(self)

  implicit none

  class(lfric_xios_field_type), intent(inout) :: self

  call log_event( "Recieving data from XIOS field ["//trim(self%xios_id)//     &
                  "] to model field ["//trim(self%model_field%get_name())//"]",&
                  log_level_trace )

  select type(fld => self%model_field)
  type is (field_real32_type)
    call fld%read_field(self%xios_id)

  type is (field_real64_type)
    call fld%read_field(self%xios_id)

  type is (integer_field_type)
    call fld%read_field(self%xios_id)

  end select

end subroutine recv

!> Gets the data from the associated model field and sends it to the
!> corresponding XIOS field
subroutine send(self)

  implicit none

  class(lfric_xios_field_type), intent(inout) :: self

  call log_event( "Sending data from model field ["//           &
                  trim(self%model_field%get_name())//           &
                  "] to XIOS field ["//trim(self%xios_id)//"]", &
                  log_level_trace )

  select type(fld => self%model_field)
  type is (field_real32_type)
    call fld%write_field(self%xios_id)

  type is (field_real64_type)
    call fld%write_field(self%xios_id)

  type is (integer_field_type)
    call fld%write_field(self%xios_id)

  end select

end subroutine send

!> Returns a pointer to the associated model field
function get_model_field(self) result(model_field_ptr)

  implicit none

  class(lfric_xios_field_type), intent(in) :: self
  class(field_parent_type), pointer :: model_field_ptr

  model_field_ptr => self%model_field

end function get_model_field

!> Sets the associated model field
!> @param[in] new_model_field  A pointer to the new model field to associate with
!>                             this XIOS field
subroutine set_model_field(self, new_model_field)

  implicit none

  class(lfric_xios_field_type),         intent(inout) :: self
  class(field_parent_type),     target, intent(in)    :: new_model_field

  self%model_field => new_model_field

end subroutine set_model_field

!> Returns the associated XIOS field name
function get_xios_name(self) result(xios_name_out)

  implicit none

  class(lfric_xios_field_type), intent(in) :: self
  character(len=str_def), allocatable :: xios_name_out

  xios_name_out = trim(adjustl(self%xios_name))

end function get_xios_name

!> Sets the associated XIOS field name
!> @param[in] xios_name_in  The new XIOS field name to associate with this field
subroutine set_xios_name(self, xios_name_in)

  implicit none

  class(lfric_xios_field_type), intent(inout) :: self
  character(len=*),             intent(in)    :: xios_name_in

  self%xios_name = trim(adjustl(xios_name_in))

end subroutine set_xios_name

!> Finaliser for the lfric_xios_field_type object
subroutine lfric_xios_field_final(self)

  implicit none

  type(lfric_xios_field_type), intent(inout) :: self

  nullify(self%model_field)

end subroutine lfric_xios_field_final

end module lfric_xios_field_mod
