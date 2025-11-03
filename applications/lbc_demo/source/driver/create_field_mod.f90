!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>
!>
module create_field_mod

  ! Infrastructure
  use constants_mod,                 only: i_def
  use field_parent_mod,              only: read_interface, write_interface
  use field_mod,                     only: field_type
  use integer_field_mod,             only: integer_field_type
  use field_collection_mod,          only: field_collection_type
  use function_space_mod,            only: function_space_type
  use function_space_collection_mod, only: function_space_collection
  use fs_continuity_mod,             only: W0, W2H, W2V, Wtheta, W3
  use mesh_mod,                      only: mesh_type

  ! I/O methods
  use lfric_xios_read_mod,  only: read_field_generic
  use lfric_xios_write_mod, only: write_field_generic

  implicit none

  private

  public :: create_field

  interface create_field
    procedure :: create_real_field
    procedure :: create_integer_field
  end interface create_field

contains

!> @brief Creates real fields, assigns their IO behaviours and adds them
!>        to a field collection.
!> @param[in,out] field_collection  Field collection to add to.
!> @param[in,out] fld               Real field object.
!> @param[in]     field_name        Field name to be created
!> @param[in]     mesh              Mesh for function space
!> @param[in]     fs_id             Function space identifier
!> @param[in]     element_order_h   Function-space order in horizontal.
!> @param[in]     element_order_v   Function-space order in vertical.
!> @param[in]     ndata_per_dof     Size of the field's multi-data axis.
subroutine create_real_field( field_collection, &
                              fld,              &
                              field_name,       &
                              mesh,             &
                              fs_id,            &
                              element_order_h,  &
                              element_order_v,  &
                              ndata_per_dof )

  implicit none

  ! Arguments
  type(field_collection_type), intent(inout) :: field_collection

  type(field_type), intent(inout) :: fld
  character(len=*), intent(in)    :: field_name
  integer(i_def),   intent(in)    :: fs_id
  integer(i_def),   intent(in)    :: element_order_h
  integer(i_def),   intent(in)    :: element_order_v
  integer(i_def),   intent(in)    :: ndata_per_dof

  type(mesh_type), pointer, intent(in) :: mesh

  ! Local variables

  ! Pointers
  procedure(read_interface),  pointer :: tmp_read_ptr
  procedure(write_interface), pointer :: tmp_write_ptr
  type(function_space_type),  pointer :: vector_space

  nullify(tmp_read_ptr)
  nullify(tmp_write_ptr)
  nullify(vector_space)

  ! Set up function space
  vector_space =>                                            &
      function_space_collection%get_fs( mesh,                &
                                        element_order_h,     &
                                        element_order_v,     &
                                        fs_id,               &
                                        ndata=ndata_per_dof, &
                                        ndata_first=.true. )

  ! Set up I/O methods
  tmp_read_ptr  => read_field_generic
  tmp_write_ptr => write_field_generic

  ! Initialise field object from specifications
  if (mesh%get_halo_depth() == 0) then
    call fld%initialise( vector_space, name=field_name, halo_depth=0 )
  else
    call fld%initialise( vector_space, name=field_name )
  end if

  call fld%set_write_behaviour( tmp_write_ptr )
  call fld%set_read_behaviour( tmp_read_ptr )

  ! Add field to group
  call field_collection%add_field( fld )

end subroutine create_real_field


!> @brief Creates integer fields, assigns their IO behaviours and adds them
!>        to a field collection.
!> @param[in,out] field_collection  Field collection to add to.
!> @param[in,out] fld               Integer field object.
!> @param[in]     field_name        Field name to be created.
!> @param[in]     mesh              Mesh for function space.
!> @param[in]     fs_id             Fucntion space identifier.
!> @param[in]     element_order_h   Function-space order in horizontal.
!> @param[in]     element_order_v   Function-space order in vertical.
!> @param[in]     ndata_per_dof     Size of the field's multi-data axis.
subroutine create_integer_field( field_collection, &
                                 fld,              &
                                 field_name,       &
                                 mesh,             &
                                 fs_id,            &
                                 element_order_h,  &
                                 element_order_v,  &
                                 ndata_per_dof )

  implicit none

  ! Arguments
  type(field_collection_type), intent(inout) :: field_collection
  type(integer_field_type),    intent(inout) :: fld

  character(len=*), intent(in) :: field_name
  integer(i_def),   intent(in) :: fs_id
  integer(i_def),   intent(in) :: element_order_h
  integer(i_def),   intent(in) :: element_order_v
  integer(i_def),   intent(in) :: ndata_per_dof

  type(mesh_type),  intent(in), pointer :: mesh

  ! Local variables

  ! Pointers
  procedure(read_interface),  pointer :: tmp_read_ptr
  procedure(write_interface), pointer :: tmp_write_ptr
  type(function_space_type),  pointer :: vector_space

  nullify(tmp_read_ptr)
  nullify(tmp_write_ptr)
  nullify(vector_space)

  vector_space =>                                            &
      function_space_collection%get_fs( mesh,                &
                                        element_order_h,     &
                                        element_order_v,     &
                                        fs_id,               &
                                        ndata=ndata_per_dof, &
                                        ndata_first=.true. )
  ! Set up I/O methods
  tmp_read_ptr  => read_field_generic
  tmp_write_ptr => write_field_generic

  ! Initialise field object from specifications
  if (mesh%get_halo_depth() == 0) then
    call fld%initialise( vector_space, name=field_name, halo_depth=0 )
  else
    call fld%initialise( vector_space, name=field_name )
  end if

  call fld%set_write_behaviour( tmp_write_ptr )
  call fld%set_read_behaviour( tmp_read_ptr )

  ! Add field to group
  call field_collection%add_field( fld )

end subroutine create_integer_field

end module create_field_mod
