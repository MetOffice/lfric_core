!-------------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Setup of fields to be read in for the multifile IO demo
!> @details Handles the setup of all the fields that will be passed to XIOS
!>          to be read in, creating the fields and adding them to the collection
!>          and setting the read behaviour.
module multifile_field_setup_mod

  use constants_mod,                 only: i_def, str_def
  use driver_modeldb_mod,            only: modeldb_type
  use field_collection_mod,          only: field_collection_type
  use field_mod,                     only: field_type
  use field_parent_mod,              only: read_interface
  use fs_continuity_mod,             only: Wtheta
  use function_space_collection_mod, only: function_space_collection
  use function_space_mod,            only: function_space_type
  use lfric_xios_read_mod,           only: read_field_generic
  use mesh_mod,                      only: mesh_type
  use mesh_collection_mod,           only: mesh_collection

  implicit none

  public create_multifile_io_fields

contains

  !> @details Creates the fields needed for the multifile IO
  !> @param[in,out] modeldb The model database in which to store model data.
  subroutine create_multifile_io_fields(modeldb)

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    type(field_collection_type), pointer :: multifile_io_fields
    type(field_type)                     :: multifile_io_field
    procedure(read_interface), pointer   :: tmp_ptr

    type(function_space_type), pointer :: fs
    type(mesh_type),           pointer :: mesh

    character(str_def) :: prime_mesh_name

    integer(i_def) :: order_h
    integer(i_def) :: order_v

    prime_mesh_name = modeldb%config%base_mesh%prime_mesh_name()

    order_h = modeldb%config%finite_element%element_order_h()
    order_v = modeldb%config%finite_element%element_order_v()

    mesh => mesh_collection%get_mesh(prime_mesh_name)
    fs   => function_space_collection%get_fs(mesh, order_h, order_v, Wtheta)

    call modeldb%fields%add_empty_field_collection("multifile_io_fields")
    multifile_io_fields => modeldb%fields%get_field_collection("multifile_io_fields")

    call multifile_io_field%initialise(fs, name="multifile_field")

    tmp_ptr => read_field_generic
    call multifile_io_field%set_read_behaviour(tmp_ptr)
    call multifile_io_fields%add_field(multifile_io_field)

  end subroutine create_multifile_io_fields

end module multifile_field_setup_mod
